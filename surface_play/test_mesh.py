import numpy as np
import pytest

from surface_play.domain import Domain
from surface_play.mesh import _generate_rect_mesh


def test_generate_rect_mesh():
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))
    resolution = 10
    uv, tris = _generate_rect_mesh(domain, resolution)

    N = len(uv)
    M = len(tris)
    n_bnd = 4 * resolution

    # 1. N ~ resolution^2, M ~ 2N, Euler chi = 1 (disk topology)
    assert N > resolution ** 2 * 0.5, f"Too few vertices: {N}"

    edges = set()
    for tri in tris:
        for k in range(3):
            edges.add(tuple(sorted([int(tri[k]), int(tri[(k + 1) % 3])])))
    E = len(edges)
    chi = N - E + M
    assert chi == 1, f"Euler characteristic should be 1 (disk), got {chi}"

    # Euler implies exactly M = 2N - n_bnd - 2 for a disk triangulation
    assert M == 2 * N - n_bnd - 2, f"M={M} should equal 2N-n_bnd-2={2*N-n_bnd-2}"

    # 2. 4*resolution boundary vertices — distinct and exactly on boundary
    u_min, u_max, v_min, v_max = domain.bounds
    bverts = uv[:n_bnd]
    assert len(np.unique(bverts, axis=0)) == n_bnd, "Boundary vertices not all distinct"

    on_boundary = (
        np.isclose(bverts[:, 0], u_min)
        | np.isclose(bverts[:, 0], u_max)
        | np.isclose(bverts[:, 1], v_min)
        | np.isclose(bverts[:, 1], v_max)
    )
    assert on_boundary.all(), "Some boundary vertices not on the rectangle boundary"

    # 3. Triangle areas: min > 0 and sum == domain area
    v0 = uv[tris[:, 0]]
    v1 = uv[tris[:, 1]]
    v2 = uv[tris[:, 2]]
    d1 = v1 - v0
    d2 = v2 - v0
    areas = 0.5 * np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0])
    assert areas.min() > 0, "Degenerate triangle found (area == 0)"
    domain_area = (u_max - u_min) * (v_max - v_min)
    assert abs(areas.sum() - domain_area) < 1e-10, (
        f"Area sum {areas.sum()} != domain area {domain_area}"
    )

    # 4. Dtypes
    assert uv.dtype == np.float64, f"uv dtype should be float64, got {uv.dtype}"
    assert tris.dtype == np.int32, f"tris dtype should be int32, got {tris.dtype}"
