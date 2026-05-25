"""Tests for surface_play.visibility (Layer O, steps O15-O16)."""

import numpy as np

from surface_play.curves import ResampledCurve
from surface_play.visibility import (
    LPInfeasibleError,
    bfs_visibility,
    break_dtype,
    compute_projection_breaks,
    lp_refine_visibility,
)


def _make_rc(kind, xy, depth, dir_=None, start=0, end=0, vc_in=0, vc_out=0):
    return ResampledCurve(
        kind=kind, start=start, end=end,
        depth=np.asarray(depth, dtype=float),
        xy=np.asarray(xy, dtype=float),
        dir=np.asarray(dir_, dtype=float) if dir_ is not None else None,
        vc_in=vc_in, vc_out=vc_out,
    )


def test_compute_projection_breaks():
    # ─── sub-assert 1: non-overlapping projected curves → no BKs ───────────
    rc_a = _make_rc("BC", [[0, 0], [1, 0]], [0.0, 0.0],
                    dir_=[[0, 1], [0, 1]])
    rc_b = _make_rc("BC", [[0, 2], [1, 2]], [0.0, 0.0],
                    dir_=[[0, 1], [0, 1]])
    bks = compute_projection_breaks([rc_a, rc_b], surface=None, projection=None)
    assert len(bks) == 0, f"non-overlapping → expected 0 BKs, got {len(bks)}"

    # ─── sub-assert 2: CC in front of an inner BC → BKs with |delta_v|=2 ─
    # Convention: axis = I × J points toward viewer, so larger projection.Z
    # means closer to viewer. CC depth = +1 (front, occluder).
    # Inner BC depth = 0 (back, occluded).
    cc = _make_rc(
        "CC",
        xy=[[0.0, 0.0], [1.0, 0.0]],
        depth=[1.0, 1.0],
        dir_=[[0.0, 1.0], [0.0, 1.0]],  # surface is "above" the CC
    )
    bc_inner = _make_rc(
        "BC",
        xy=[[0.5, -0.5], [0.5, 0.5]],
        depth=[0.0, 0.0],
        dir_=[[1.0, 0.0], [1.0, 0.0]],
    )
    bks = compute_projection_breaks([cc, bc_inner], surface=None, projection=None)
    assert len(bks) == 1, f"single crossing → expected 1 BK, got {len(bks)}"
    bk = bks[0]
    assert int(bk["rc_idx"]) == 1, (
        f"BK should be on the OCCLUDED rc_idx=1 (inner BC, deeper), got {int(bk['rc_idx'])}"
    )
    assert abs(int(bk["delta_v"])) == 2, (
        f"CC-as-occluder → |delta_v|=2, got {abs(int(bk['delta_v']))}"
    )

    # ─── sub-assert 3: real fixture — fig8 ortho ─────────────────────────
    # Build a tiny realistic pipeline; just check BK count > 0 (fig8 has
    # known self-occlusion) and delta_v ∈ {-2, -1, +1, +2}.
    from surface_play.contour import (
        build_contour_curves,
        build_contour_segments,
        find_contour_points,
        find_vps,
    )
    from surface_play.curves import build_bcs, resample_all
    from surface_play.helpers import build_helper_curves
    from surface_play.intersections import (
        build_sics,
        build_sis_pairs,
        find_double_points,
        find_triple_points,
        tp_dtype,
    )
    from surface_play.mesh import build_mesh
    from surface_play.projection import Projection
    from surface_play.splitting import (
        SplitArrays, assemble_subcurves,
        split_at_cdps, split_bcs_at_bcps, split_bcs_at_bdps,
        split_bcs_at_corners, split_ccs_at_vps, split_sics_at_tps,
    )
    from surface_play.test_fixtures import fig8, perturb_axis

    surf = fig8(perturb=False)
    mesh = build_mesh(surf.domain, surf, resolution=20, jitter=True, seed=42)
    bcs = build_bcs(mesh)
    proj = Projection(surf,
                      I=perturb_axis([1.0, 0.0, 0.0], seed=1),
                      J=perturb_axis([0.0, 0.0, 1.0], seed=2))
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    vps = find_vps(ccs, css, cps, surf, proj)
    dps = find_double_points(mesh, surf)
    sis_pairs = build_sis_pairs(dps)
    sics = build_sics(sis_pairs) if len(sis_pairs) else []
    tps = (find_triple_points(sis_pairs, dps, mesh, surf)
           if len(sis_pairs) >= 3 else np.empty(0, dtype=tp_dtype))
    splits = SplitArrays()
    split_bcs_at_corners(mesh, bcs, splits, proj)
    split_bcs_at_bcps(mesh, bcs, ccs, css, cps, splits, surf, proj)
    split_bcs_at_bdps(mesh, bcs, sics, sis_pairs, dps, splits, surf, proj)
    split_sics_at_tps(tps, sis_pairs, dps, splits, surf, proj)
    split_ccs_at_vps(vps, css, ccs, splits, surf, proj)
    split_at_cdps(mesh, css, cps, sis_pairs, dps, splits, surf, proj)
    hcs = build_helper_curves(bcs, ccs, sics, css, sis_pairs, cps, dps,
                              mesh, splits, proj, surf, mesh.domain)
    subs = assemble_subcurves(bcs, ccs, sics, hcs, mesh, css, sis_pairs, splits)
    rcs = resample_all(subs, surf, proj, splits, mesh, css, sis_pairs, cps, dps,
                       resolution=30)
    bks = compute_projection_breaks(rcs, surf, proj)
    # fig8 has self-occlusion → expect non-trivial BK count.
    assert len(bks) > 0, "fig8 ortho should produce projection breaks"
    assert set(int(b["delta_v"]) for b in bks) <= {-2, -1, 1, 2}, (
        f"delta_v out of expected set: {sorted(set(int(b['delta_v']) for b in bks))}"
    )

    # ─── sub-assert 4: same-RC self-crossing IS emitted (A) ─────────────
    # A BC that crosses itself in projection: e.g. an X shape on a single RC.
    # Sample order forms a path that crosses itself.
    self_x = _make_rc(
        "BC",
        xy=[[-1.0, -1.0], [1.0, 1.0], [-1.0, 1.0], [1.0, -1.0]],
        depth=[0.0, 1.0, 0.0, 1.0],
        dir_=[[0.0, 1.0]] * 4,
    )
    bks = compute_projection_breaks([self_x], surface=None, projection=None)
    assert len(bks) >= 1, (
        f"same-RC self-crossing (A) → ≥1 BK, got {len(bks)}"
    )
    assert all(int(b["rc_idx"]) == 0 for b in bks), (
        "all BKs should be on rc_idx=0 (only RC)"
    )

    # ─── sub-assert 5: HC crossing BC → no BK (HC doesn't occlude) ───────
    # HC at depth=1 (front, would be the occluder by depth alone).
    # BC at depth=0 (back). HC.kind not in {BC, CC} → no BK emitted.
    bc_h = _make_rc(
        "BC",
        xy=[[0.0, -1.0], [0.0, 1.0]],
        depth=[0.0, 0.0],
        dir_=[[1.0, 0.0], [1.0, 0.0]],
    )
    hc = _make_rc(
        "HC",
        xy=[[-1.0, 0.0], [1.0, 0.0]],
        depth=[1.0, 1.0],
        dir_=None,
    )
    bks = compute_projection_breaks([bc_h, hc], surface=None, projection=None)
    assert len(bks) == 0, (
        f"HC-as-occluder should produce 0 BKs, got {len(bks)} "
        f"(delta_v={[int(b['delta_v']) for b in bks]})"
    )


# ── O16: bfs_visibility ─────────────────────────────────────────────────────

def test_bfs_visibility():
    # ─── sub-assert 1: paraboloid ortho → anchor vis=0, all vis ≤ 0 ──────
    from surface_play.contour import (
        build_contour_curves, build_contour_segments,
        find_contour_points, find_vps,
    )
    from surface_play.curves import build_bcs, resample_all
    from surface_play.helpers import build_helper_curves
    from surface_play.intersections import (
        build_sics, build_sis_pairs, find_double_points,
        find_triple_points, tp_dtype,
    )
    from surface_play.mesh import build_mesh
    from surface_play.projection import Projection
    from surface_play.splitting import (
        SplitArrays, assemble_subcurves,
        split_at_cdps, split_bcs_at_bcps, split_bcs_at_bdps,
        split_bcs_at_corners, split_ccs_at_vps, split_sics_at_tps,
    )
    from surface_play.test_fixtures import paraboloid, perturb_axis

    def _pipeline(surf, I=(1, 0, 0), J=(0, 1, 0), res=15):
        mesh = build_mesh(surf.domain, surf, resolution=res, jitter=True, seed=42)
        bcs = build_bcs(mesh)
        proj = Projection(surf,
                          I=perturb_axis(I, seed=1),
                          J=perturb_axis(J, seed=2))
        cps = find_contour_points(mesh, proj)
        css = build_contour_segments(cps, mesh)
        ccs = build_contour_curves(css, cps)
        vps = find_vps(ccs, css, cps, surf, proj)
        dps = find_double_points(mesh, surf)
        sis_pairs = build_sis_pairs(dps)
        sics = build_sics(sis_pairs) if len(sis_pairs) else []
        tps = (find_triple_points(sis_pairs, dps, mesh, surf)
               if len(sis_pairs) >= 3 else np.empty(0, dtype=tp_dtype))
        splits = SplitArrays()
        split_bcs_at_corners(mesh, bcs, splits, proj)
        split_bcs_at_bcps(mesh, bcs, ccs, css, cps, splits, surf, proj)
        split_bcs_at_bdps(mesh, bcs, sics, sis_pairs, dps, splits, surf, proj)
        split_sics_at_tps(tps, sis_pairs, dps, splits, surf, proj)
        split_ccs_at_vps(vps, css, ccs, splits, surf, proj)
        split_at_cdps(mesh, css, cps, sis_pairs, dps, splits, surf, proj)
        hcs = build_helper_curves(bcs, ccs, sics, css, sis_pairs, cps, dps,
                                  mesh, splits, proj, surf, mesh.domain)
        subs = assemble_subcurves(bcs, ccs, sics, hcs, mesh, css, sis_pairs, splits)
        rcs = resample_all(subs, surf, proj, splits, mesh, css, sis_pairs, cps, dps,
                           resolution=30)
        return rcs, splits, surf, proj

    rcs, splits, surf, proj = _pipeline(paraboloid(perturb=False), J=(0, 0, 1))
    breaks = compute_projection_breaks(rcs, surf, proj)
    vis = bfs_visibility(rcs, breaks, splits, anchor_mode="leftmost")
    assert len(vis) == len(rcs)
    for rc in rcs:
        v = vis[id(rc)]
        assert (v <= 0).all(), (
            f"paraboloid: vis should be ≤0 everywhere, got max={int(v.max())}"
        )

    # ─── sub-assert 2: single closed RC, no breaks → all vis=0 ──────────
    rc_loop = _make_rc(
        "BC",
        xy=[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]],
        depth=[0.0] * 5,
        dir_=[[0, 1]] * 5,
        start=-1, end=-1,  # SP-less closed
    )
    vis_loop = bfs_visibility([rc_loop], np.empty(0, dtype=break_dtype),
                              splits=None, anchor_mode="leftmost")
    np.testing.assert_array_equal(vis_loop[id(rc_loop)], np.zeros(5, dtype=np.int32))

    # ─── sub-assert 3: synthetic break → after-break samples vis=-2 ─────
    rc_one = _make_rc(
        "BC",
        xy=[[0, 0], [1, 0], [2, 0], [3, 0]],
        depth=[0.0] * 4,
        dir_=[[0, 1]] * 4,
        start=-1, end=-1,
    )
    breaks_one = np.array([(0, 1, 0.5, (1.5, 0.0), -2)], dtype=break_dtype)
    v = bfs_visibility([rc_one], breaks_one, splits=None, anchor_mode="leftmost")
    arr = v[id(rc_one)]
    # Anchor is the leftmost sample = sample 0 (x=0). Break is between samples 1 and 2.
    assert arr[0] == 0 and arr[1] == 0, f"pre-break should be 0, got {arr.tolist()}"
    assert arr[2] == -2 and arr[3] == -2, (
        f"post-break should be -2, got {arr.tolist()}"
    )

    # ─── sub-assert 4: multi-curve through SP — propagation across SP ────
    # Two RCs sharing SP idx=0 at one endpoint. RC_a vc_out=-1; RC_b vc_in=0.
    # Anchor on RC_a sample 0; expect RC_b's first sample = vis_a_end - (-1) + 0.
    rc_a = _make_rc(
        "CC",
        xy=[[0, 0], [1, 0]],
        depth=[0.0, 0.0],
        dir_=[[0, 1], [0, 1]],
        start=-1, end=0, vc_out=-1,
    )
    rc_b = _make_rc(
        "CC",
        xy=[[1, 0], [2, 0]],
        depth=[0.0, 0.0],
        dir_=[[0, 1], [0, 1]],
        start=0, end=-1, vc_in=0,
    )
    v2 = bfs_visibility([rc_a, rc_b], np.empty(0, dtype=break_dtype),
                        splits=None, anchor_mode="leftmost")
    arr_a = v2[id(rc_a)]
    arr_b = v2[id(rc_b)]
    # RC_a: anchor sample 0 vis=0, sample 1 (end) vis=0 (no breaks).
    # RC_b start: vis_a_end - vc_out_a + vc_in_b = 0 - (-1) + 0 = +1.
    assert arr_a[0] == 0 and arr_a[1] == 0, f"rc_a vis: {arr_a.tolist()}"
    assert arr_b[0] == 1 and arr_b[1] == 1, f"rc_b vis: {arr_b.tolist()}"

    # ─── sub-assert 5: determinism (G23) — repeated runs identical ──────
    v_run_1 = bfs_visibility(rcs, breaks, splits, anchor_mode="leftmost")
    v_run_2 = bfs_visibility(rcs, breaks, splits, anchor_mode="leftmost")
    assert set(v_run_1.keys()) == set(v_run_2.keys())
    for k in v_run_1:
        np.testing.assert_array_equal(v_run_1[k], v_run_2[k])


# ── O17: lp_refine_visibility ───────────────────────────────────────────────

def test_lp_refine_visibility():
    # ─── sub-assert 1: no-correction case (consistent input → LP unchanged) ─
    # Single 3-sample open RC, one break with delta_v=-1 between samples 0-1.
    # Anchor at leftmost (sample 0); BFS gives [0, -1, -1].
    rc = _make_rc(
        "BC",
        xy=[[0, 0], [1, 0], [2, 0]],
        depth=[0.0, 0.0, 0.0],
        dir_=[[0, 1]] * 3,
        start=-1, end=-1,
    )
    bks = np.array([(0, 0, 0.5, (0.5, 0.0), -1)], dtype=break_dtype)
    v_bfs = bfs_visibility([rc], bks, splits=None, anchor_mode="leftmost")
    v_lp = lp_refine_visibility([rc], bks, splits=None, vis_bfs=v_bfs,
                                anchors="leftmost")
    np.testing.assert_array_equal(v_lp[id(rc)], v_bfs[id(rc)])

    # ─── sub-assert 2: single bad break → LP corrects to all-zero ─────────
    # 3-sample open RC, break with delta_v=+3 (impossible since vis ≤ 0).
    # LP must slack vc to 0, giving s=3, vis = [0, 0, 0].
    rc2 = _make_rc(
        "BC",
        xy=[[0, 0], [1, 0], [2, 0]],
        depth=[0.0, 0.0, 0.0],
        dir_=[[0, 1]] * 3,
        start=-1, end=-1,
    )
    bks2 = np.array([(0, 0, 0.5, (0.5, 0.0), 3)], dtype=break_dtype)
    v_lp2 = lp_refine_visibility([rc2], bks2, splits=None, vis_bfs=None,
                                 anchors="leftmost")
    np.testing.assert_array_equal(v_lp2[id(rc2)], np.zeros(3, dtype=np.int32))

    # ─── sub-assert 3: anchors="extremes" → 4 anchor samples all vis=0 ──
    # Three small RCs whose leftmost/rightmost/topmost/bottommost are
    # distinct samples.
    rc_a = _make_rc("BC", xy=[[-5, 0], [-4, 0]], depth=[0, 0],
                    dir_=[[0, 1]] * 2, start=-1, end=-1)
    rc_b = _make_rc("BC", xy=[[5, 0], [4, 0]], depth=[0, 0],
                    dir_=[[0, 1]] * 2, start=-1, end=-1)
    rc_c = _make_rc("BC", xy=[[0, -5], [0, 5]], depth=[0, 0],
                    dir_=[[1, 0]] * 2, start=-1, end=-1)
    v_lp3 = lp_refine_visibility(
        [rc_a, rc_b, rc_c], np.empty(0, dtype=break_dtype),
        splits=None, vis_bfs=None, anchors="extremes",
    )
    # Leftmost = rc_a[0] (x=-5); rightmost = rc_b[0] (x=5);
    # bottommost = rc_c[0] (y=-5); topmost = rc_c[1] (y=5).
    assert int(v_lp3[id(rc_a)][0]) == 0
    assert int(v_lp3[id(rc_b)][0]) == 0
    assert int(v_lp3[id(rc_c)][0]) == 0
    assert int(v_lp3[id(rc_c)][1]) == 0

    # ─── sub-assert 4: infeasibility raises LPInfeasibleError ─────────────
    # Two RCs sharing an SP at each end → forms a cycle. Choose vc constants
    # so the cycle sum ≠ 0, making SP coupling impossible to satisfy.
    # RC_x: start=0, end=1, vc_in=0, vc_out=0
    # RC_y: start=0, end=1, vc_in=-1, vc_out=0  (different vc_in)
    # SP 0: ref RC_x (start), coupling row for RC_y start:
    #   vis_y[0] - vis_x[0] = (-1) - 0 = -1
    # SP 1: ref RC_x (end), coupling row for RC_y end:
    #   vis_y[N-1] - vis_x[N-1] = 0 - 0 = 0
    # Combined with propagation (no breaks → vis const within each RC):
    #   vis_y[0] = vis_y[N-1], vis_x[0] = vis_x[N-1]
    # So vis_y - vis_x = -1 AND vis_y - vis_x = 0 → infeasible.
    rc_x = _make_rc("BC", xy=[[0, 0], [1, 0]], depth=[0, 0],
                    dir_=[[0, 1]] * 2, start=0, end=1, vc_in=0, vc_out=0)
    rc_y = _make_rc("BC", xy=[[0, 1], [1, 1]], depth=[0, 0],
                    dir_=[[0, 1]] * 2, start=0, end=1, vc_in=-1, vc_out=0)
    try:
        lp_refine_visibility(
            [rc_x, rc_y], np.empty(0, dtype=break_dtype),
            splits=None, vis_bfs=None, anchors="leftmost",
        )
        raised = False
    except LPInfeasibleError:
        raised = True
    assert raised, "expected LPInfeasibleError for contradictory SP coupling"

    # ─── sub-assert 5: paraboloid ortho → LP matches BFS ────────────────
    # Reuse `_pipeline` from test_bfs_visibility scope is not accessible;
    # rebuild minimally inline.
    from surface_play.contour import (
        build_contour_curves, build_contour_segments,
        find_contour_points, find_vps,
    )
    from surface_play.curves import build_bcs, resample_all
    from surface_play.helpers import build_helper_curves
    from surface_play.intersections import (
        build_sics, build_sis_pairs, find_double_points,
        find_triple_points, tp_dtype,
    )
    from surface_play.mesh import build_mesh
    from surface_play.projection import Projection
    from surface_play.splitting import (
        SplitArrays, assemble_subcurves,
        split_at_cdps, split_bcs_at_bcps, split_bcs_at_bdps,
        split_bcs_at_corners, split_ccs_at_vps, split_sics_at_tps,
    )
    from surface_play.test_fixtures import paraboloid, perturb_axis

    surf = paraboloid(perturb=False)
    mesh = build_mesh(surf.domain, surf, resolution=15, jitter=True, seed=42)
    bcs = build_bcs(mesh)
    proj = Projection(surf,
                      I=perturb_axis([1.0, 0.0, 0.0], seed=1),
                      J=perturb_axis([0.0, 0.0, 1.0], seed=2))
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    vps = find_vps(ccs, css, cps, surf, proj)
    dps = find_double_points(mesh, surf)
    sis_pairs = build_sis_pairs(dps)
    sics = build_sics(sis_pairs) if len(sis_pairs) else []
    tps = (find_triple_points(sis_pairs, dps, mesh, surf)
           if len(sis_pairs) >= 3 else np.empty(0, dtype=tp_dtype))
    splits = SplitArrays()
    split_bcs_at_corners(mesh, bcs, splits, proj)
    split_bcs_at_bcps(mesh, bcs, ccs, css, cps, splits, surf, proj)
    split_bcs_at_bdps(mesh, bcs, sics, sis_pairs, dps, splits, surf, proj)
    split_sics_at_tps(tps, sis_pairs, dps, splits, surf, proj)
    split_ccs_at_vps(vps, css, ccs, splits, surf, proj)
    split_at_cdps(mesh, css, cps, sis_pairs, dps, splits, surf, proj)
    hcs = build_helper_curves(bcs, ccs, sics, css, sis_pairs, cps, dps,
                              mesh, splits, proj, surf, mesh.domain)
    subs = assemble_subcurves(bcs, ccs, sics, hcs, mesh, css, sis_pairs, splits)
    rcs5 = resample_all(subs, surf, proj, splits, mesh, css, sis_pairs, cps, dps,
                        resolution=30)
    bks5 = compute_projection_breaks(rcs5, surf, proj)
    v_bfs5 = bfs_visibility(rcs5, bks5, splits, anchor_mode="leftmost")
    v_lp5 = lp_refine_visibility(rcs5, bks5, splits, vis_bfs=v_bfs5,
                                 anchors="leftmost")
    for rc in rcs5:
        np.testing.assert_array_equal(v_lp5[id(rc)], v_bfs5[id(rc)])
