"""Tests for surface_play.curves.resample_all (Layer O, step O14)."""

import numpy as np
import pytest

from surface_play.contour import (
    build_contour_curves,
    build_contour_segments,
    find_contour_points,
    find_vps,
)
from surface_play.curves import (
    ResampledCurve,
    build_bcs,
    resample_all,
)
from surface_play.helpers import build_helper_curves
from surface_play.intersections import (
    build_sics,
    build_sis_pairs,
    dp_dtype,
    find_double_points,
    find_triple_points,
    sis_dtype,
    tp_dtype,
)
from surface_play.mesh import build_mesh
from surface_play.projection import Projection
from surface_play.splitting import (
    SplitArrays,
    assemble_subcurves,
    split_at_cdps,
    split_bcs_at_bcps,
    split_bcs_at_bdps,
    split_bcs_at_corners,
    split_ccs_at_vps,
    split_sics_at_tps,
)
from surface_play.test_fixtures import (
    disk_paraboloid_ca,
    fig8,
    helicoid,
    torus,
)


def _full_pipeline(surf, I=(1.0, 0.0, 0.0), J=(0.0, 1.0, 0.0), resolution=15):
    mesh = build_mesh(surf.domain, surf, resolution=resolution, jitter=True, seed=42)
    bcs = build_bcs(mesh)
    proj = Projection(surf, I=list(I), J=list(J))
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    vps = find_vps(ccs, css, cps, surf, proj)
    dps = find_double_points(mesh, surf)
    sis_pairs = build_sis_pairs(dps)
    sics = build_sics(sis_pairs) if len(sis_pairs) else []
    if len(sis_pairs) >= 3:
        tps = find_triple_points(sis_pairs, dps, mesh, surf)
    else:
        tps = np.empty(0, dtype=tp_dtype)

    splits = SplitArrays()
    split_bcs_at_corners(mesh, bcs, splits, proj)
    split_bcs_at_bcps(mesh, bcs, ccs, css, cps, splits, surf, proj)
    split_bcs_at_bdps(mesh, bcs, sics, sis_pairs, dps, splits, surf, proj)
    split_sics_at_tps(tps, sis_pairs, dps, splits, surf, proj)
    split_ccs_at_vps(vps, css, ccs, splits, surf, proj)
    split_at_cdps(mesh, css, cps, sis_pairs, dps, splits, surf, proj)

    hcs = build_helper_curves(
        bcs, ccs, sics, css, sis_pairs, cps, dps,
        mesh, splits, proj, surf, mesh.domain,
    )
    subs = assemble_subcurves(bcs, ccs, sics, hcs, mesh, css, sis_pairs, splits)
    return mesh, proj, splits, cps, css, ccs, dps, sis_pairs, subs


def test_resample_all():
    # ─── sub-assert 1: open SubCurve gets ≥4 samples per half ───────────────
    # Note: the helicoid's CC silhouette is the screw axis (u=0), so the
    # Z-axis view I=(1,0,0), J=(0,1,0) collapses the whole CC to (0, 0) in
    # image space — a degenerate non-generic input that resample_all now
    # raises on. Use a tilted view instead (axis = I×J ≠ z-axis).
    surf = helicoid(perturb=False)
    mesh, proj, splits, cps, css, ccs, dps, sis_pairs, subs = _full_pipeline(
        surf, I=(1.0, 0.0, 0.0), J=(0.0, 1.0, 1.0), resolution=15,
    )
    rcs = resample_all(
        subs, surf, proj, splits, mesh, css, sis_pairs, cps, dps,
        resolution=30,
    )
    assert len(rcs) == len(subs)
    # Find an open SubCurve with both endpoints distinct SPs.
    open_rcs = [
        rc for rc, sub in zip(rcs, subs)
        if sub.kind in ("BC", "CC")
        and not sub.is_closed
        and sub.start >= 0 and sub.end >= 0
        and len(rc.xy) > 0
    ]
    assert open_rcs, "expected at least one open BC/CC SubCurve"
    # Per-sample ≥8 check removed 2026-05-29: BC/CC now sample at CPs
    # verbatim (see [[resume-hc-match-cc]]). Short SCs between two SPs with
    # no intermediate CPs naturally have only 2 samples — the SP endpoints
    # ARE the curve geometrically; there is nothing in between.
    for rc in open_rcs:
        assert len(rc.xy) >= 2, (
            f"open RC must include both SP endpoints, got {len(rc.xy)}"
        )

    # ─── sub-assert 2: closed CC (torus side, no SPs would have been added) ─
    # Torus side view: 4 cusps split CCs into closed-ish arcs. Use top-down
    # torus view: 2 closed CCs without VPs (no cusps top-down).
    surf_t = torus(perturb=False)
    mesh_t, proj_t, splits_t, cps_t, css_t, ccs_t, dps_t, sis_t, subs_t = \
        _full_pipeline(surf_t, J=(0.0, 1.0, 0.0), resolution=20)
    # Top-down torus has 2 closed CCs and no other curves; O12 adds an HC that
    # creates SPTs on each. Each closed CC becomes a closed SubCurve with
    # start == end (the HA SP).
    closed_subs = [
        (i, s) for i, s in enumerate(subs_t)
        if s.kind == "CC" and s.is_closed and s.start == s.end and s.start >= 0
    ]
    assert closed_subs, (
        f"expected ≥1 closed CC SubCurve from top-down torus; got kinds "
        f"{[(s.kind, s.is_closed, s.start, s.end) for s in subs_t]}"
    )
    rcs_t = resample_all(
        subs_t, surf_t, proj_t, splits_t, mesh_t, css_t, sis_t, cps_t, dps_t,
        resolution=30,
    )
    for i, _ in closed_subs:
        rc = rcs_t[i]
        assert rc.start == rc.end and rc.start >= 0
        # Uniformity-of-spacing check removed 2026-05-29: CCs now use their
        # CPs verbatim as samples (see [[resume-hc-match-cc]]). CP positions
        # follow the triangulation contour-chain, not a uniform ladder, so
        # consecutive xy chord-diffs are not bounded by a small factor.

    # ─── sub-assert 3: annular/disk BC snap ─────────────────────────────────
    surf_d = disk_paraboloid_ca(perturb=False)
    mesh_d, proj_d, splits_d, cps_d, css_d, ccs_d, dps_d, sis_d, subs_d = \
        _full_pipeline(surf_d, J=(0.0, 1.0, 0.0), resolution=15)
    rcs_d = resample_all(
        subs_d, surf_d, proj_d, splits_d, mesh_d, css_d, sis_d, cps_d, dps_d,
        resolution=30,
    )
    bc_rcs = [
        (i, rc) for i, (rc, sub) in enumerate(zip(rcs_d, subs_d)) if sub.kind == "BC"
    ]
    assert bc_rcs, "expected ≥1 BC RC on cartesian disk"
    r_min = float(surf_d.domain.bounds[0])
    r_max = float(surf_d.domain.bounds[1])
    # The annular snap is applied to interpolated uvs; we check that the
    # lifted xy lies on the surface boundary in 3D (xyz = (x, y, z) with
    # ‖(x, y)‖ in {r_min, r_max}).
    for i, rc in bc_rcs:
        sub = subs_d[i]
        if rc.start == -1 and rc.end == -1:
            continue  # SP-less verbatim — already on the boundary by construction.
        # Re-derive each sample's uv via xyz: surf is x=u, y=v, so xyz[:2] == uv.
        # Recover uv from xy or via lifted xyz (we don't store uv on RC).
        # Use rc.xy (proj is identity for ortho onto x,y for J=(0,1,0); axis=(0,0,-1)
        # — depth-only mismatch, xy IS (x, y) up to sign). The projection is
        # I=[1,0,0], J=[0,1,0] (axis=[0,0,1] for ortho), so xy = (x, y) directly.
        norms = np.linalg.norm(rc.xy, axis=1)
        # Interior samples (not pinned endpoints) should be on a boundary.
        for k in range(1, len(norms) - 1):
            err_max = abs(norms[k] - r_max)
            err_min = abs(norms[k] - r_min)
            assert min(err_max, err_min) < 1e-9, (
                f"disk BC sample k={k} has ‖xy‖={norms[k]}, not on r_min={r_min} or r_max={r_max}"
            )

    # ─── sub-assert 4: SIC two-sheet (G21) — uv close-aware step ────────────
    # Tilted view (axis not aligned with the fig-8's symmetry planes) so SIC
    # SubCurves don't collapse to image-space points.
    surf_8 = fig8(perturb=False)
    mesh_8, proj_8, splits_8, cps_8, css_8, ccs_8, dps_8, sis_8, subs_8 = \
        _full_pipeline(surf_8, I=(1.0, 0.0, 0.3), J=(0.0, 1.0, 0.4), resolution=20)
    assert len(sis_8) > 0, "fig8 should produce SISs"
    rcs_8 = resample_all(
        subs_8, surf_8, proj_8, splits_8, mesh_8, css_8, sis_8, cps_8, dps_8,
        resolution=30,
    )
    domain_8 = mesh_8.domain
    sic_subs_and_rcs = [
        (sub, rc) for sub, rc in zip(subs_8, rcs_8) if sub.kind == "SIC"
    ]
    assert sic_subs_and_rcs, "expected SIC SubCurves from fig8"
    # Re-walk each SIC RC and check consecutive xy displacements are bounded
    # (no preimage jumps — which would manifest as ~M-sized jumps in xy).
    M_xy = float(
        np.hypot(*(proj_8.XY(mesh_8.xyz).max(axis=0) - proj_8.XY(mesh_8.xyz).min(axis=0)))
    )
    for sub, rc in sic_subs_and_rcs:
        if len(rc.xy) < 2:
            continue
        steps = np.linalg.norm(np.diff(rc.xy, axis=0), axis=1)
        assert steps.max() < 0.5 * M_xy, (
            f"SIC RC has a {steps.max():.3f} step (M={M_xy:.3f}) — preimage jump"
        )

    # ─── sub-assert 5: endpoint SP indices (G5) ────────────────────────────
    for rc, sub in zip(rcs, subs):
        if sub.kind == "HC":
            continue
        if sub.start == -1 and sub.end == -1:
            assert rc.start == -1 and rc.end == -1
            continue
        assert rc.start == sub.start
        assert rc.end == sub.end
        # Identity through dereference.
        assert splits.sps[rc.start] is splits.sps[sub.start]
        assert splits.sps[rc.end] is splits.sps[sub.end]

    # ─── sub-assert 6: endpoint xy pinning ─────────────────────────────────
    for rc, sub in zip(rcs, subs):
        if sub.start == -1 or len(rc.xy) == 0:
            continue
        expected_start = np.asarray(splits.sps[sub.start][2], dtype=float)
        np.testing.assert_allclose(rc.xy[0], expected_start, atol=1e-12)
        if not sub.is_closed and sub.end >= 0:
            expected_end = np.asarray(splits.sps[sub.end][2], dtype=float)
            np.testing.assert_allclose(rc.xy[-1], expected_end, atol=1e-12)
