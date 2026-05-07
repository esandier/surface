import math

import numpy as np
import pytest

from surface_play.domain import Domain
from surface_play.mesh import _generate_disk_mesh, _generate_rect_mesh


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


def _euler_chi(uv, tris):
    edges = set()
    for tri in tris:
        for k in range(3):
            edges.add(tuple(sorted([int(tri[k]), int(tri[(k + 1) % 3])])))
    return len(uv) - len(edges) + len(tris)


def _tri_areas(uv, tris):
    v0, v1, v2 = uv[tris[:, 0]], uv[tris[:, 1]], uv[tris[:, 2]]
    d1, d2 = v1 - v0, v2 - v0
    return 0.5 * np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0])


def test_generate_disk_mesh():
    # 1. Unit disk: chi = 1
    domain_disk = Domain(type="disk", bounds=(0.0, 1.0, 0.0, 2 * math.pi))
    uv_d, tris_d = _generate_disk_mesh(domain_disk, resolution=20)

    assert _euler_chi(uv_d, tris_d) == 1, "Unit disk: Euler chi should be 1"

    # 3. Outer boundary vertices within 1e-12 of r_max=1
    n_outer_d = max(3, round(math.pi * 20))
    outer_d = uv_d[:n_outer_d]
    r_outer_d = np.sqrt(outer_d[:, 0] ** 2 + outer_d[:, 1] ** 2)
    assert np.all(np.abs(r_outer_d - 1.0) < 1e-12), "Outer ring not on r=1"

    # 3. Areas: min > 0, sum == polygon area (not π — boundary is a polygon)
    areas_d = _tri_areas(uv_d, tris_d)
    assert areas_d.min() > 0, "Degenerate triangle in disk mesh"
    # Regular n-gon inscribed in unit circle: A = (n/2)*sin(2π/n)
    poly_area_d = 0.5 * n_outer_d * math.sin(2 * math.pi / n_outer_d)
    assert abs(areas_d.sum() - poly_area_d) < 1e-10, (
        f"Disk area sum {areas_d.sum()} != polygon area {poly_area_d}"
    )

    # 4. Dtypes
    assert uv_d.dtype == np.float64
    assert tris_d.dtype == np.int32

    # 2. Unit annulus r_min=0.3, r_max=1: chi = 0
    domain_ann = Domain(type="annulus", bounds=(0.3, 1.0, 0.0, 2 * math.pi))
    uv_a, tris_a = _generate_disk_mesh(domain_ann, resolution=20)

    assert _euler_chi(uv_a, tris_a) == 0, "Annulus: Euler chi should be 0"

    # No triangles inside inner ring
    centroids = (uv_a[tris_a[:, 0]] + uv_a[tris_a[:, 1]] + uv_a[tris_a[:, 2]]) / 3
    r_centroids = np.sqrt(centroids[:, 0] ** 2 + centroids[:, 1] ** 2)
    assert r_centroids.min() >= 0.3 - 1e-10, "Triangles found inside inner ring"

    # 3. Boundary vertices: outer on r=1, inner on r=0.3
    n_outer_a = max(3, round(math.pi * 20))
    n_inner_a = max(3, round(math.pi * 0.3 * 20))
    outer_a = uv_a[:n_outer_a]
    inner_a = uv_a[n_outer_a : n_outer_a + n_inner_a]
    assert np.all(np.abs(np.sqrt(outer_a[:, 0] ** 2 + outer_a[:, 1] ** 2) - 1.0) < 1e-12)
    assert np.all(np.abs(np.sqrt(inner_a[:, 0] ** 2 + inner_a[:, 1] ** 2) - 0.3) < 1e-12)

    # 3. Areas: min > 0, sum == polygon annulus area
    areas_a = _tri_areas(uv_a, tris_a)
    assert areas_a.min() > 0
    poly_outer_a = 0.5 * n_outer_a * math.sin(2 * math.pi / n_outer_a)
    poly_inner_a = 0.5 * n_inner_a * 0.3 ** 2 * math.sin(2 * math.pi / n_inner_a)
    poly_area_a = poly_outer_a - poly_inner_a
    assert abs(areas_a.sum() - poly_area_a) < 1e-10, (
        f"Annulus area sum {areas_a.sum()} != polygon annulus area {poly_area_a}"
    )

    # 4. Dtypes
    assert uv_a.dtype == np.float64
    assert tris_a.dtype == np.int32
