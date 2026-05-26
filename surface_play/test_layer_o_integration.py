"""Layer O integration test.

Runs the full per-viewpoint pipeline (build_mesh → contours → splits → helpers
→ assemble_subcurves → resample_all → compute_projection_breaks →
bfs_visibility → lp_refine_visibility) on every fixture × multiple random
axes, and asserts:

  1. LP feasibility (no LPInfeasibleError on any trial).
  2. BFS values are all ≤ 0 (visibility constraint).
  3. LP values are all ≤ 0.
  4. BFS == LP on every RC sample (LP repair count is 0).

If any of (1)-(4) fail on any trial, the test fails with the offending
fixture/trial reported.
"""

import numpy as np
import pytest

from surface_play.contour import (build_contour_curves, build_contour_segments,
    find_contour_points, find_vps)
from surface_play.curves import build_bcs, resample_all
from surface_play.helpers import build_helper_curves
from surface_play.intersections import (build_sics, build_sis_pairs,
    find_double_points, find_triple_points, tp_dtype)
from surface_play.mesh import build_mesh
from surface_play.projection import Projection
from surface_play.splitting import (SplitArrays, assemble_subcurves, split_at_cdps,
    split_bcs_at_bcps, split_bcs_at_bdps, split_bcs_at_corners,
    split_ccs_at_vps, split_sics_at_tps)
from surface_play.test_fixtures import (paraboloid, torus, mobius_u, mobius_v,
    disk_paraboloid_po, disk_paraboloid_ca)
from surface_play.visibility import (bfs_visibility, compute_projection_breaks,
    lp_refine_visibility, LPInfeasibleError)


FIXTURES = [
    ("paraboloid",      paraboloid),
    ("torus",           torus),
    ("mobius_u",        mobius_u),
    ("mobius_v",        mobius_v),
    ("disk_para_po",    disk_paraboloid_po),
    ("disk_para_ca",    disk_paraboloid_ca),
]


def _random_axis(rng):
    a = rng.standard_normal(3); a /= np.linalg.norm(a)
    helper = np.array([1.0, 0.0, 0.0])
    if abs(a @ helper) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    I = helper - (helper @ a) * a; I /= np.linalg.norm(I)
    J = np.cross(a, I)
    return a, I, J


def _run_pipeline(surf, I, J, resolution):
    mesh = build_mesh(surf.domain, surf, resolution=resolution, jitter=True, seed=42)
    proj = Projection(surf, I=I.tolist(), J=J.tolist())
    bcs = build_bcs(mesh)
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
                              mesh, splits, proj, surf, surf.domain)
    subs = assemble_subcurves(bcs, ccs, sics, hcs, mesh, css, sis_pairs, splits)
    rcs = resample_all(subs, surf, proj, splits, mesh, css, sis_pairs, cps, dps)
    breaks = compute_projection_breaks(rcs, surf, proj)
    bfs = bfs_visibility(rcs, breaks, splits)
    lp = lp_refine_visibility(rcs, breaks, splits, vis_bfs=bfs)
    return rcs, bfs, lp


# Parametrize over (fixture, trial) so failures pinpoint the exact case.
# N_TRIALS=4 keeps the suite fast (~30 s); the standalone diag_sweep.py
# script handles 10+ trials when needed.
N_TRIALS = 4
SEED = 2026
RESOLUTION = 100


def _generate_cases():
    rng = np.random.default_rng(SEED)
    out = []
    for fixture_name, factory in FIXTURES:
        for t in range(N_TRIALS):
            a, I, J = _random_axis(rng)
            out.append((fixture_name, factory, t, a, I, J))
    return out


_CASES = _generate_cases()


@pytest.mark.parametrize(
    "fixture_name,factory,trial,axis,I,J",
    _CASES,
    ids=[f"{c[0]}_t{c[2]}" for c in _CASES],
)
def test_layer_o_pipeline(fixture_name, factory, trial, axis, I, J):
    """Layer O pipeline is feasible and BFS == LP on every random viewpoint.

    Per Layer O sign-off gate (roadmap line 1513-1518): the full per-viewpoint
    pipeline must complete without LPInfeasibleError, produce all-nonpositive
    visibilities (BFS and LP), and BFS must agree with LP (no cycle-deficit
    repair needed). Failure on any trial isolates the offending fixture +
    axis for diagnosis.
    """
    surf = factory(perturb=False)
    try:
        rcs, bfs, lp = _run_pipeline(surf, I, J, RESOLUTION)
    except LPInfeasibleError as e:
        pytest.fail(
            f"{fixture_name} trial {trial}: LP INFEASIBLE at axis={axis.round(4).tolist()} "
            f"({e})"
        )

    # Aggregate stats across all RCs.
    bfs_flat = np.concatenate([bfs[id(rc)] for rc in rcs]) if rcs else np.array([0])
    lp_flat = np.concatenate([lp[id(rc)] for rc in rcs]) if rcs else np.array([0])
    bfs_max = int(bfs_flat.max()); bfs_min = int(bfs_flat.min())
    lp_max = int(lp_flat.max());  lp_min = int(lp_flat.min())
    n_diff = sum(1 for rc in rcs if not np.array_equal(bfs[id(rc)], lp[id(rc)]))

    assert bfs_max <= 0, (
        f"{fixture_name} t{trial}: BFS has positive value {bfs_max} "
        f"(range=[{bfs_min},{bfs_max}]) at axis={axis.round(4).tolist()}"
    )
    assert lp_max <= 0, (
        f"{fixture_name} t{trial}: LP has positive value {lp_max} "
        f"(range=[{lp_min},{lp_max}]) at axis={axis.round(4).tolist()}"
    )
    assert n_diff == 0, (
        f"{fixture_name} t{trial}: BFS != LP on {n_diff} RC(s) — LP had to "
        f"repair a cycle deficit. axis={axis.round(4).tolist()}"
    )
