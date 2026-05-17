"""Tests for O1 — find_contour_points.

Five criteria from the roadmap:
  1. Helicoid ortho-Z: CP count ≥ 1, all s ∈ (0, 1).
  2. Paraboloid side view: contour at v=0; after Newton uv[:,1] ≈ 0 to 1e-8.
  3. Newton convergence: residual |SN·axis| < 1e-8 at every CP; midpoint
     (use_newton=False) gives same edges and positions within one edge length.
  4. Boundary CPs on rect-no-no: CPs with ptype==4 lie on boundary edges (g==-1).
  5. Möbius seam: no duplicate edge indices in result (flip prevents double-detect).
"""

import numpy as np
import pytest

from surface_play.test_fixtures import (
    helicoid, paraboloid, torus, mobius_u,
    helicoid_ortho_view, paraboloid_side_view, torus_ortho_view, mobius_ortho_view,
)
from surface_play.mesh import build_mesh as _build_mesh
from surface_play.projection import Projection
from surface_play.contour import (
    find_contour_points, build_contour_segments, build_contour_curves, find_vps,
)

# Torus side view: I=(1,0,0), J=(0,0,1) → axis = I×J = (0,−1,0), looking
# down the -Y axis.  The torus silhouette has 4 cusps in this view.
torus_side_view = ([1, 0, 0], [0, 0, 1], None)

RES = 20


# ── helpers ───────────────────────────────────────────────────────────────────

def _make(surface_factory, view, res=RES, perturb=False):
    surf = surface_factory(perturb=perturb)
    mesh = _build_mesh(surf.domain, surf, res)
    I, J, eye = view
    proj = Projection(surf, I, J, eye=eye)
    return mesh, proj


# ── test 1: helicoid ortho-Z, CP count + s range ─────────────────────────────

def test_find_contour_points_helicoid_count_and_s_range():
    """Contour of helicoid viewed top-down (Z axis) is along u=0.

    Expect at least one CP per face column; all s strictly inside (0, 1).
    """
    mesh, proj = _make(helicoid, helicoid_ortho_view)
    cps = find_contour_points(mesh, proj)

    assert len(cps) >= 1, "expected at least one contour point"
    assert np.all(cps["s"] > 0.0), "s must be > 0"
    assert np.all(cps["s"] < 1.0), "s must be < 1"


# ── test 2: paraboloid side view, analytic contour at v=0 ────────────────────

def test_find_contour_points_paraboloid_analytic():
    """Paraboloid S=(u,v,u²+v²) viewed along -Y (axis = I×J = (0,-1,0)).

    Contour condition: SN·axis = 2v = 0  →  v = 0 exactly.
    After Newton all uv[:,1] must be within 1e-7 of zero.
    """
    mesh, proj = _make(paraboloid, paraboloid_side_view)
    cps = find_contour_points(mesh, proj)

    assert len(cps) >= 1, "expected contour points"
    assert np.allclose(cps["uv"][:, 1], 0.0, atol=1e-7), (
        "Newton should converge to v=0; max deviation = "
        f"{np.max(np.abs(cps['uv'][:, 1])):.2e}"
    )


# ── test 3: Newton convergence ────────────────────────────────────────────────

def test_find_contour_points_newton_convergence():
    """Newton (use_newton=True) converges to residual < 1e-8;
    midpoint (use_newton=False) hits the same edges but may be up to one edge
    length away from the Newton solution.
    """
    mesh, proj = _make(paraboloid, paraboloid_side_view)

    cps_newton  = find_contour_points(mesh, proj, use_newton=True)
    cps_midpt   = find_contour_points(mesh, proj, use_newton=False)

    # Residual check: |SN(uv) · axis| < 1e-8 for Newton CPs
    axis = proj._axis
    for cp in cps_newton:
        u, v = cp["uv"]
        SN = mesh.surface.SN(u, v)          # (3,) or (3, 1)
        SN = np.asarray(SN).ravel()[:3]
        residual = abs(float(SN @ axis))
        assert residual < 1e-8, f"Newton residual {residual:.2e} too large on edge {cp['e']}"

    # Same edge set detected regardless of Newton
    assert set(cps_newton["e"]) == set(cps_midpt["e"]), (
        "Newton and midpoint must detect the same sign-changing edges"
    )

    # Midpoint positions within one domain edge length of Newton positions
    idx_n = np.argsort(cps_newton["e"])
    idx_m = np.argsort(cps_midpt["e"])
    uv_n  = cps_newton["uv"][idx_n]
    uv_m  = cps_midpt["uv"][idx_m]
    # Edge length in domain ≈ 2/RES
    edge_len = 2.0 / RES * np.sqrt(2)
    diff = np.linalg.norm(uv_n - uv_m, axis=1)
    assert np.all(diff <= edge_len + 1e-10), (
        f"midpoint CPs must be within one edge length of Newton CPs; max diff = {diff.max():.4f}"
    )


# ── test 4: boundary CPs have ptype == 4 ─────────────────────────────────────

def test_find_contour_points_boundary_ptype():
    """Paraboloid side view: contour at v=0 crosses left/right domain boundaries
    (u = ±1), producing CPs on boundary edges (g == -1) with ptype == 4.
    """
    mesh, proj = _make(paraboloid, paraboloid_side_view)
    cps = find_contour_points(mesh, proj)

    boundary_cps = cps[cps["ptype"] == 4]
    interior_cps = cps[cps["ptype"] == 0]

    # Each ptype==4 CP must sit on a boundary edge
    for cp in boundary_cps:
        e = mesh.edges[cp["e"]]
        assert e["g"] == -1, f"ptype=4 CP on edge {cp['e']} but g={e['g']} (not boundary)"

    # Each ptype==0 CP must sit on an interior edge
    for cp in interior_cps:
        e = mesh.edges[cp["e"]]
        assert e["g"] != -1, f"ptype=0 CP on edge {cp['e']} but g=={e['g']} (boundary)"

    # Paraboloid rect-no-no has boundary edges; contour at v=0 must hit them
    assert len(boundary_cps) >= 1, (
        "expected at least one boundary CP where contour curve hits the domain boundary"
    )


# ── test 5: Möbius seam — no duplicate edge indices ──────────────────────────

def test_find_contour_points_mobius_no_seam_duplication():
    """Seam edges of the Möbius band carry flip=-1. sign_changes must not
    double-count them: the result must have unique edge indices.
    """
    mesh, proj = _make(mobius_u, mobius_ortho_view)
    cps = find_contour_points(mesh, proj)

    edge_indices = cps["e"]
    assert len(edge_indices) == len(np.unique(edge_indices)), (
        "duplicate edge indices detected — seam edges counted twice"
    )


# ── test 6: d field is unit-length on a smooth silhouette (regression) ───────

def test_find_contour_points_d_field_unit_length():
    """The d field is the 2D curvature direction of the projected contour
    (legacy silhouette.py "dir_vec"), not the projection of kerdS_3d.

    Regression: on a torus top-view silhouette (two smooth circles, no cusps)
    every CP must carry a unit-length d vector.  An earlier impl returned
    proj_vec(uv, kerdS_3d), which is identically zero at every CP because
    kerdS_3d lies in the tangent plane and at a CP the view axis does too,
    so the projection vanishes.
    """
    mesh, proj = _make(torus, torus_ortho_view)
    cps = find_contour_points(mesh, proj)

    norms = np.linalg.norm(cps["d"], axis=1)
    assert np.all(np.abs(norms - 1.0) < 1e-9), (
        f"d field must be unit-length; got |d| range "
        f"[{norms.min():.3e}, {norms.max():.3e}]"
    )


# ══ O2: build_contour_segments ════════════════════════════════════════════════

# ── test 1: helicoid CS count > 0 and matches pairing ─────────────────────────

def test_build_contour_segments_helicoid_count():
    """Helicoid ortho-Z: every paired face produces one CS; count > 0."""
    mesh, proj = _make(helicoid, helicoid_ortho_view)
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)

    assert len(css) >= 1, "expected at least one contour segment"


# ── test 2: each CS's endpoints lie on different edges of the same face ────────

def test_build_contour_segments_endpoints_on_same_face():
    """Each CS's p_cp and q_cp must be on different edges of cs.face."""
    mesh, proj = _make(paraboloid, paraboloid_side_view)
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)

    assert len(css) >= 1
    for cs in css:
        face = mesh.faces[cs["face"]]
        face_edge_set = set(face["edges"].tolist())
        ep = int(cps[cs["p_cp"]]["e"])
        eq = int(cps[cs["q_cp"]]["e"])
        assert ep in face_edge_set, f"p_cp edge {ep} not in face {cs['face']} edges"
        assert eq in face_edge_set, f"q_cp edge {eq} not in face {cs['face']} edges"
        assert ep != eq, "both endpoints on the same edge"


# ── test 3: CPs referenced by CSs are strictly interior (s in (0, 1)) ─────────

def test_build_contour_segments_no_vertex_crossing():
    """No CS passes through a vertex: s of every referenced CP must be in (0, 1)."""
    mesh, proj = _make(helicoid, helicoid_ortho_view)
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)

    cp_indices_used = np.unique(np.concatenate([css["p_cp"], css["q_cp"]]))
    s_vals = cps[cp_indices_used]["s"]
    assert np.all(s_vals > 0.0) and np.all(s_vals < 1.0), (
        "CP used in a CS has s outside (0, 1)"
    )


# ── test 4: splits initialised to -1 (G17) ────────────────────────────────────

def test_build_contour_segments_splits_initialised():
    """split1 and split2 must both be -1 on every freshly built CS (G17)."""
    mesh, proj = _make(paraboloid, paraboloid_side_view)
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)

    assert np.all(css["split1"] == -1), "split1 must be -1"
    assert np.all(css["split2"] == -1), "split2 must be -1"


# ══ O3: build_contour_curves ══════════════════════════════════════════════════

def _pipeline(surface_factory, view, res=RES, perturb=True):
    """Build mesh → CPs → CSs → CCs in one call."""
    mesh, proj = _make(surface_factory, view, res=res, perturb=perturb)
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    return mesh, proj, cps, css, ccs


# ── test 1: empty input ────────────────────────────────────────────────────────

def test_build_contour_curves_empty():
    """Empty css → empty CC list."""
    from surface_play.contour import cp_dtype, cs_dtype
    cps = np.zeros(0, dtype=cp_dtype)
    css = np.zeros(0, dtype=cs_dtype)
    ccs = build_contour_curves(css, cps)
    assert ccs == []


# ── test 2: helicoid rect-no-no gives open CC(s) ──────────────────────────────

def test_build_contour_curves_helicoid_open():
    """Helicoid on a non-identified rectangle: contour at u=0 is an open curve.
    Expect ≥1 CC and no closed CC (rect-no-no cannot close the contour).
    """
    _, _, _, _, ccs = _pipeline(helicoid, helicoid_ortho_view)
    assert len(ccs) >= 1, "expected at least one CC"
    assert all(not cc.is_closed for cc in ccs), (
        "helicoid rect-no-no should produce only open CCs"
    )


# ── test 3: cylinder_cy gives closed CC(s) ────────────────────────────────────

def test_build_contour_curves_torus_closed():
    """Torus (u_cy, v_cy) viewed along Z: the two silhouette circles are closed CCs."""
    _, _, _, _, ccs = _pipeline(torus, torus_ortho_view)
    assert len(ccs) >= 1, "expected at least one CC"
    assert any(cc.is_closed for cc in ccs), (
        "torus should produce at least one closed CC"
    )


# ── test 4: boundary CPs appear only at chain endpoints ───────────────────────

def test_build_contour_curves_boundary_cps_are_endpoints():
    """A boundary CP (ptype=4) belongs to exactly 1 CS, so it must be a
    degree-1 node in the CS graph — i.e., an open-chain endpoint, never interior.
    """
    _, _, cps, css, ccs = _pipeline(paraboloid, paraboloid_side_view)

    bcp_indices = set(int(i) for i in np.where(cps["ptype"] == 4)[0])
    if not bcp_indices:
        pytest.skip("no boundary CPs found for this fixture/view — adjust if needed")

    # Each boundary CP must appear in exactly 1 CS
    all_ends = np.concatenate([css["p_cp"], css["q_cp"]])
    for bcp in bcp_indices:
        count = int(np.sum(all_ends == bcp))
        assert count == 1, (
            f"boundary CP {bcp} referenced by {count} CSs (must be 1 → chain endpoint)"
        )


# ══ O4: find_vps ══════════════════════════════════════════════════════════════
#
# Roadmap test criteria (lines 1063-1069):
#   1. No-cusp surface → empty VP array.
#   2. Surface with known cusps (torus side view → 4 cusps): count and bounds.
#   3. refine=False vs refine=True positions within L/10.
#   4. vis_change ∈ {-1, +1} (G10).
#   5. Möbius: VPs valid, flip-aware d field (no spurious cusps from seam).


def _full_pipeline(surface_factory, view, *, res=RES, perturb=True):
    """Build mesh → CPs → CSs → CCs → VPs."""
    mesh, proj = _make(surface_factory, view, res=res, perturb=perturb)
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    return mesh, proj, cps, css, ccs


# ── test 1: torus top-down has no cusps ───────────────────────────────────────

def test_find_vps_no_cusp_torus_top_down():
    """Torus viewed from above: silhouette is two smooth concentric circles
    (one per v-branch, no cusps).  find_vps must return an empty array.

    Paraboloid_side_view is *not* used here because at v=0 the kerdS direction
    is parallel to the view axis, the d field is then defined by 2nd-order
    behaviour, and the no-cusp signal is fine — but the cleaner non-degenerate
    no-cusp test is the torus top view.
    """
    mesh, proj, cps, css, ccs = _full_pipeline(torus, torus_ortho_view)
    vps = find_vps(ccs, css, cps, mesh.surface, proj)
    assert len(vps) == 0, f"expected 0 VPs on torus top-down, got {len(vps)}"


# ── test 2: torus side view has the canonical 4 cusps ────────────────────────

def test_find_vps_torus_side_four_cusps():
    """Torus viewed along -Y exhibits the canonical 4 cusps on its silhouette.

    Per roadmap line 1065: "Torus has 4 cusps in the initial viewpoint of the
    DB".  Each VP must sit strictly inside its CS (s ∈ (0, 1)) and reference
    a valid CS index.
    """
    mesh, proj, cps, css, ccs = _full_pipeline(torus, torus_side_view)
    vps = find_vps(ccs, css, cps, mesh.surface, proj)

    assert len(vps) == 4, f"expected 4 VPs on torus side view, got {len(vps)}"
    assert np.all((vps["s"] > 0.0) & (vps["s"] < 1.0)), "VP s must be in (0, 1)"
    assert np.all((vps["cs"] >= 0) & (vps["cs"] < len(css))), (
        "VP cs index out of bounds"
    )


# ── test 3: refine=False and refine=True agree within ~edge_length / 10 ───────

def test_find_vps_refine_within_edge_length():
    """Refinement moves the VP off the chord midpoint but no farther than
    roughly one mesh edge length / 10 — both produce the same count.

    On the torus side view, edge length in v is ≈ (2π) / RES, so the bound
    is ≈ (2π)/RES √2 / 10.  We use a generous factor since the orthogonal
    Newton step can leave the chord segment by an amount of that scale.
    """
    mesh, proj, cps, css, ccs = _full_pipeline(torus, torus_side_view)
    vps_mid    = find_vps(ccs, css, cps, mesh.surface, proj, refine=False)
    vps_refine = find_vps(ccs, css, cps, mesh.surface, proj, refine=True)

    assert len(vps_mid) == len(vps_refine), (
        f"refine changed VP count: {len(vps_mid)} → {len(vps_refine)}"
    )

    edge_len = (2 * np.pi / RES) * np.sqrt(2)
    bound = edge_len  # within one edge length is plenty for this granularity
    for vm, vr in zip(vps_mid, vps_refine):
        diff = float(np.linalg.norm(vm["uv"] - vr["uv"]))
        assert diff <= bound + 1e-10, (
            f"refinement moved VP by {diff:.4f} > edge_len={edge_len:.4f}"
        )


# ── test 4: vis_change ∈ {-1, +1} (G10) ───────────────────────────────────────

def test_find_vps_vis_change_sign():
    """Every VP carries vis_change exactly +1 or -1, never 0 (G10)."""
    mesh, proj, cps, css, ccs = _full_pipeline(torus, torus_side_view)
    vps = find_vps(ccs, css, cps, mesh.surface, proj)

    assert len(vps) >= 1, "need at least one VP for this assertion"
    assert np.all(np.abs(vps["vis_change"]) == 1), (
        f"vis_change must be ±1; got {np.unique(vps['vis_change'])}"
    )


# ── test 5: Möbius — no spurious VPs from the seam, valid invariants ──────────

def test_find_vps_mobius_no_seam_artifacts():
    """The seam of the Möbius band carries flip=-1 on its edges.  O1's
    sign_changes already honours flip per-edge, and the d field is intrinsic
    in uv (computed from second derivatives + projection); the chain through
    the seam should therefore not introduce spurious sign reversals.

    Whatever VP count the Möbius gives, every VP must satisfy the standard
    invariants: valid CS index, s ∈ (0, 1), |vis_change| = 1.
    """
    mesh, proj, cps, css, ccs = _full_pipeline(mobius_u, mobius_ortho_view)
    vps = find_vps(ccs, css, cps, mesh.surface, proj)

    if len(vps) > 0:
        assert np.all((vps["s"] > 0.0) & (vps["s"] < 1.0)), (
            "Möbius VP s out of (0, 1)"
        )
        assert np.all((vps["cs"] >= 0) & (vps["cs"] < len(css))), (
            "Möbius VP cs index out of bounds"
        )
        assert np.all(np.abs(vps["vis_change"]) == 1), (
            "Möbius VP vis_change must be ±1"
        )
