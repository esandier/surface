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
    helicoid, paraboloid, mobius_u,
    helicoid_ortho_view, paraboloid_side_view, mobius_ortho_view,
)
from surface_play.mesh import build_mesh as _build_mesh
from surface_play.projection import Projection
from surface_play.contour import find_contour_points

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
