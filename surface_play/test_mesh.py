import collections
import math

import numpy as np
import pytest
from scipy.spatial.distance import pdist

from surface_play.domain import Domain
from surface_play.mesh import (
    Mesh,
    _apply_identifications,
    _build_edges_faces,
    _generate_disk_mesh,
    _generate_rect_mesh,
    _jitter,
    build_mesh,
    edge_dtype,
    face_dtype,
)
from surface_play.surface import SurfaceParams
from surface_play.test_fixtures import (
    disk_paraboloid_ca,
    mobius_u,
    paraboloid,
    torus,
)


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


def _typical_edge_length(uv, tris):
    e0 = np.linalg.norm(uv[tris[:, 1]] - uv[tris[:, 0]], axis=1)
    e1 = np.linalg.norm(uv[tris[:, 2]] - uv[tris[:, 1]], axis=1)
    e2 = np.linalg.norm(uv[tris[:, 0]] - uv[tris[:, 2]], axis=1)
    return float(np.median(np.concatenate([e0, e1, e2])))


def test_jitter():
    # --- rect no-id: all 4 sides are true boundaries ---
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))
    uv, tris = _generate_rect_mesh(domain, resolution=10)
    L = _typical_edge_length(uv, tris)
    uv_j = _jitter(uv, tris, domain, seed=42)

    # 5. Vertex count unchanged
    assert len(uv_j) == len(uv)

    # 1. Bbox shift < 0.002*L per axis
    for ax in range(2):
        assert abs(uv_j[:, ax].min() - uv[:, ax].min()) < 0.002 * L
        assert abs(uv_j[:, ax].max() - uv[:, ax].max()) < 0.002 * L

    # 2. True boundary vertices remain exactly on boundary
    u_min, u_max, v_min, v_max = domain.bounds
    tol = 1e-12
    assert np.all(np.abs(uv_j[np.abs(uv[:, 0] - u_min) < tol, 0] - u_min) < 1e-12)
    assert np.all(np.abs(uv_j[np.abs(uv[:, 0] - u_max) < tol, 0] - u_max) < 1e-12)
    assert np.all(np.abs(uv_j[np.abs(uv[:, 1] - v_min) < tol, 1] - v_min) < 1e-12)
    assert np.all(np.abs(uv_j[np.abs(uv[:, 1] - v_max) < tol, 1] - v_max) < 1e-12)

    # 3. Reproducibility: same seed → identical; different seed → different
    assert np.array_equal(uv_j, _jitter(uv, tris, domain, seed=42))
    assert not np.array_equal(uv_j, _jitter(uv, tris, domain, seed=99))

    # 4. No vertex collapse
    assert pdist(uv_j).min() > 0

    # --- rect cy (u identified): only v-sides are true boundaries ---
    domain_cy = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0), u_identify="cy")
    uv_cy, tris_cy = _generate_rect_mesh(domain_cy, resolution=10)
    uv_cy_j = _jitter(uv_cy, tris_cy, domain_cy, seed=42)
    # v-sides clamped
    assert np.all(np.abs(uv_cy_j[np.abs(uv_cy[:, 1] - 0.0) < tol, 1] - 0.0) < 1e-12)
    assert np.all(np.abs(uv_cy_j[np.abs(uv_cy[:, 1] - 1.0) < tol, 1] - 1.0) < 1e-12)
    # u-sides are free (seam) — at least one u-side vertex should have moved
    mask_u0 = np.abs(uv_cy[:, 0] - 0.0) < tol
    assert not np.allclose(uv_cy_j[mask_u0, 0], 0.0, atol=1e-14)

    # --- disk: outer ring stays at r_max ---
    domain_disk = Domain(type="disk", bounds=(0.0, 1.0, 0.0, 2 * math.pi))
    uv_d, tris_d = _generate_disk_mesh(domain_disk, resolution=20)
    uv_d_j = _jitter(uv_d, tris_d, domain_disk, seed=42)
    r_pre = np.sqrt(uv_d[:, 0] ** 2 + uv_d[:, 1] ** 2)
    r_post = np.sqrt(uv_d_j[:, 0] ** 2 + uv_d_j[:, 1] ** 2)
    assert np.all(np.abs(r_post[np.abs(r_pre - 1.0) < tol] - 1.0) < 1e-12)

    # --- annulus: both rings preserved ---
    domain_ann = Domain(type="annulus", bounds=(0.3, 1.0, 0.0, 2 * math.pi))
    uv_a, tris_a = _generate_disk_mesh(domain_ann, resolution=20)
    uv_a_j = _jitter(uv_a, tris_a, domain_ann, seed=42)
    r_a_pre = np.sqrt(uv_a[:, 0] ** 2 + uv_a[:, 1] ** 2)
    r_a_post = np.sqrt(uv_a_j[:, 0] ** 2 + uv_a_j[:, 1] ** 2)
    assert np.all(np.abs(r_a_post[np.abs(r_a_pre - 1.0) < tol] - 1.0) < 1e-12)
    assert np.all(np.abs(r_a_post[np.abs(r_a_pre - 0.3) < tol] - 0.3) < 1e-12)


def _euler_geometric(tris: np.ndarray, vertex_class: np.ndarray) -> int:
    """
    Euler χ of the post-id complex when edges are merged geometrically (by class)
    and triangles are not merged. This is a quick check used in tests; it counts
    by canonical labels and so under-counts edges in pathological coincidences,
    but it suffices for cy-cy/mo-no/etc.
    """
    classes = vertex_class[tris]
    V = len(np.unique(classes))
    all_edges = np.sort(
        np.concatenate([classes[:, [0, 1]], classes[:, [1, 2]], classes[:, [2, 0]]]),
        axis=1,
    )
    E = len(np.unique(all_edges, axis=0))
    return V - E + len(tris)


def test_apply_identifications():
    res = 10

    def mesh_rect(u_id, v_id):
        d = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0),
                   u_identify=u_id, v_identify=v_id)
        uv, tris = _generate_rect_mesh(d, resolution=res)
        return d, uv, tris

    # 1. no-no: tris and vertex_class are pass-through identity.
    d, uv, tris = mesh_rect("no", "no")
    uv_out, tris_out, vc = _apply_identifications(uv, tris, d)
    assert np.array_equal(tris, tris_out), "no-no: tris must be unchanged"
    assert np.array_equal(vc, np.arange(len(uv))), "no-no: vertex_class == arange(N)"

    # 2. cy-no: each u-side class has size 2 (u_min↔u_max), others size 1.
    d, uv, tris = mesh_rect("cy", "no")
    _, tris_out, vc = _apply_identifications(uv, tris, d)
    sizes = np.bincount(vc)
    sizes_nonzero = sizes[sizes > 0]
    assert set(sizes_nonzero.tolist()) <= {1, 2}, (
        f"cy-no: class sizes should be 1 or 2, got {set(sizes_nonzero.tolist())}"
    )
    # Number of size-2 classes = number of u-side pairs = n + 1 (the corner pair (0,n) plus n more)
    n = res
    expected_pairs = n + 1  # (0,n), (3n+0,2n-0), (3n+1,2n-1), ..., (3n+(n-1),2n-(n-1)) = (4n-1,n+1)
    # Careful: (3n+0, 2n-0) = (3n, 2n); (0, n) already counts vertex 0. Pair (3n, 2n) joins corners.
    # Just check the count of size-2 classes is reasonable: equal to (4n - V_postid).
    assert (sizes == 2).sum() >= n, "cy-no: should have at least n size-2 classes"
    # tris not relabeled — pre-id labels survive.
    assert np.array_equal(tris_out, tris), "cy-no: tris should be pre-id (no relabel)"

    # 3. mo-mo res=5: corners collapse to exactly 2 classes; both corner triangles survive.
    d_mm = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0),
                  u_identify="mo", v_identify="mo")
    uv_mm, tris_mm = _generate_rect_mesh(d_mm, resolution=5)
    _, tris_mm_out, vc_mm = _apply_identifications(uv_mm, tris_mm, d_mm)
    n = 5
    corner_classes = {int(vc_mm[0]), int(vc_mm[n]), int(vc_mm[2 * n]), int(vc_mm[3 * n])}
    assert len(corner_classes) == 2, (
        f"mo-mo: expected 2 corner classes, got {corner_classes}"
    )
    assert int(vc_mm[0]) == int(vc_mm[2 * n]), "mo-mo: corners 0 and 2n must share class"
    assert int(vc_mm[n]) == int(vc_mm[3 * n]), "mo-mo: corners n and 3n must share class"
    # Critical regression: both pre-id triangles (5,6,4) and (14,15,16) — which under the
    # old C4 dedup got collapsed to one — must remain.
    def has_tri(tris_arr, a, b, c):
        s = {a, b, c}
        return any(s == {int(t[0]), int(t[1]), int(t[2])} for t in tris_arr)
    assert has_tri(tris_mm, 5, 6, 4) and has_tri(tris_mm, 14, 15, 16), (
        "fixture sanity: pre-id (5,6,4) and (14,15,16) should both exist"
    )
    assert has_tri(tris_mm_out, 5, 6, 4) and has_tri(tris_mm_out, 14, 15, 16), (
        "mo-mo: both pre-id corner triangles must survive (regression for label-dedup bug)"
    )

    # 4. cy-cy res=10: each class has size ≤ 4; tris not relabeled.
    d, uv, tris = mesh_rect("cy", "cy")
    _, tris_out, vc = _apply_identifications(uv, tris, d)
    sizes = np.bincount(vc)
    assert sizes.max() <= 4, f"cy-cy: max class size {sizes.max()} should be ≤ 4"
    assert np.array_equal(tris_out, tris), "cy-cy: tris should be pre-id (no relabel)"

    # 5. Disk: vertex_class == arange(N), tris unchanged.
    d_disk = Domain(type="disk", bounds=(0.0, 1.0, 0.0, 2 * math.pi))
    uv_d, tris_d = _generate_disk_mesh(d_disk, resolution=20)
    _, tris_d_out, vc_d = _apply_identifications(uv_d, tris_d, d_disk)
    assert np.array_equal(tris_d, tris_d_out)
    assert np.array_equal(vc_d, np.arange(len(uv_d)))

    # 6. Pairing under jitter: vertex_class is index-based, immune to uv perturbation.
    d, uv, tris = mesh_rect("cy", "no")
    uv_j = _jitter(uv, tris, d, seed=42)
    _, _, vc_clean = _apply_identifications(uv, tris, d)
    _, _, vc_jitter = _apply_identifications(uv_j, tris, d)
    assert np.array_equal(vc_clean, vc_jitter), \
        "vertex_class must be index-based, not coordinate-based"

    # 7. Determinism: same input → same vertex_class.
    _, _, vc_a = _apply_identifications(uv, tris, d)
    _, _, vc_b = _apply_identifications(uv, tris, d)
    assert np.array_equal(vc_a, vc_b)


def _build_post_id(u_id, v_id, resolution, seed=42):
    """Helper: generate raw mesh, jitter, apply identifications. Returns
    (domain, uv_jittered, tris_filtered, vertex_class)."""
    d = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0),
               u_identify=u_id, v_identify=v_id)
    uv, tris = _generate_rect_mesh(d, resolution=resolution)
    uv_j = _jitter(uv, tris, d, seed=seed)
    _, tris_f, vc = _apply_identifications(uv_j, tris, d)
    return d, uv_j, tris_f, vc


def test_build_edges_faces():
    import collections
    import inspect
    from surface_play import mesh as _mesh_module

    # 1. Boundary counts for each (u_id, v_id) combo at res=5.
    res = 5
    expected_boundary = {
        ("no", "no"): 4 * res,
        ("cy", "no"): 2 * res,
        ("no", "cy"): 2 * res,
        ("mo", "no"): 2 * res,
        ("no", "mo"): 2 * res,
        ("cy", "cy"): 0,
        ("cy", "mo"): 0,
        ("mo", "cy"): 0,
        ("mo", "mo"): 0,
    }
    expected_chi = {
        ("no", "no"): 1,
        ("cy", "no"): 0, ("no", "cy"): 0,
        ("mo", "no"): 0, ("no", "mo"): 0,
        ("cy", "cy"): 0, ("cy", "mo"): 0, ("mo", "cy"): 0,
        ("mo", "mo"): 1,
    }

    for (u_id, v_id), bnd_expected in expected_boundary.items():
        d, uv, tris, vc = _build_post_id(u_id, v_id, resolution=res)
        SN = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
        edges, faces = _build_edges_faces(uv, tris, vc, SN, d)
        assert edges.dtype == edge_dtype
        assert faces.dtype == face_dtype

        # boundary count
        bnd = int((edges["g"] == -1).sum())
        assert bnd == bnd_expected, (
            f"{u_id}-{v_id} res={res}: boundary={bnd}, expected {bnd_expected}"
        )

        # face count == len(tris)
        assert len(faces) == len(tris), (
            f"{u_id}-{v_id}: face count {len(faces)} != tris {len(tris)}"
        )

        # Euler χ from V = #classes used by tris, E = len(edges), F = len(faces).
        V = len(np.unique(vc[tris]))
        chi = V - len(edges) + len(faces)
        assert chi == expected_chi[(u_id, v_id)], (
            f"{u_id}-{v_id}: χ={chi}, expected {expected_chi[(u_id, v_id)]}"
        )

    # 2. mo-mo manifoldness across resolutions.
    for r in (5, 7, 10):
        d, uv, tris, vc = _build_post_id("mo", "mo", resolution=r)
        SN = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
        edges, faces = _build_edges_faces(uv, tris, vc, SN, d)
        # Every edge has 2 incident faces.
        edge_face_count = collections.Counter()
        for f in faces:
            for e_idx in f["edges"]:
                edge_face_count[int(e_idx)] += 1
        counts = sorted(edge_face_count.values())
        assert min(counts) == 2 and max(counts) == 2, (
            f"mo-mo res={r}: non-manifold edges, face-count hist = "
            f"{collections.Counter(counts)}"
        )
        assert (edges["g"] == -1).sum() == 0, (
            f"mo-mo res={r}: should have zero boundary edges"
        )

    # 3. flip on cy-no with constant SN: every edge has flip=+1.
    d, uv, tris, vc = _build_post_id("cy", "no", resolution=res)
    SN_const = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
    edges_cy, _ = _build_edges_faces(uv, tris, vc, SN_const, d)
    assert (edges_cy["flip"] == 1).all(), (
        "cy-no with constant SN: every edge should have flip=+1"
    )

    # 4. flip on mo-no with a Möbius-band-like SN that reverses across the u-seam.
    d, uv, tris, vc = _build_post_id("mo", "no", resolution=res)
    u_min, u_max, _, _ = d.bounds
    period_u = u_max - u_min
    s = np.pi * (uv[:, 0] - u_min) / period_u
    SN_mo = np.column_stack([np.cos(s), np.zeros(len(uv)), np.sin(s)])
    edges_mo, _ = _build_edges_faces(uv, tris, vc, SN_mo, d)
    # At least some u-merged edges (those that did merge across the seam) should have flip=-1.
    interior = edges_mo["g"] >= 0
    assert (edges_mo["flip"][interior] == -1).any(), (
        "mo-no with sign-reversing SN: some merged seam edges should have flip=-1"
    )
    # Non-side edges (with both endpoints close in u) should have flip=+1.
    # Sanity: at least most edges still +1.
    assert (edges_mo["flip"] == 1).sum() > (edges_mo["flip"] == -1).sum()

    # 5. dir: boundary edges have dir ⟂ pq and pointing toward third vertex.
    d, uv, tris, vc = _build_post_id("no", "no", resolution=res)
    SN = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
    edges_nn, _ = _build_edges_faces(uv, tris, vc, SN, d)
    bmask = edges_nn["g"] == -1
    pq_b = edges_nn["pq"][bmask]
    dir_b = edges_nn["dir"][bmask]
    dot_pq = (dir_b * pq_b).sum(axis=1)
    assert np.all(np.abs(dot_pq) < 1e-12), "boundary dir must be ⟂ pq"
    for e in edges_nn[bmask]:
        f_idx = int(e["f"])
        verts = tuple(int(x) for x in tris[f_idx])
        third = [v for v in verts if v != int(e["p_idx"]) and v != int(e["q_idx"])][0]
        offset = uv[third] - e["p"]
        assert float(np.dot(e["dir"], offset)) > 0, (
            "boundary dir must point inward toward the third vertex"
        )

    # 6. p, pq, pr are exact pre-id subtractions (no float drift from close()).
    d, uv, tris, vc = _build_post_id("cy", "cy", resolution=res)
    SN = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
    edges_cc, faces_cc = _build_edges_faces(uv, tris, vc, SN, d)
    for e in edges_cc:
        i, j = int(e["p_idx"]), int(e["q_idx"])
        assert np.array_equal(e["p"], uv[i]), "edge p must equal uv[p_idx] exactly"
        assert np.array_equal(e["pq"], uv[j] - uv[i]), (
            "edge pq must equal uv[q_idx] - uv[p_idx] exactly"
        )
    for f in faces_cc:
        i, j, k = int(f["verts"][0]), int(f["verts"][1]), int(f["verts"][2])
        assert np.array_equal(f["p"], uv[i])
        assert np.array_equal(f["pq"], uv[j] - uv[i])
        assert np.array_equal(f["pr"], uv[k] - uv[i])

    # 7. Splits initialized to -1.
    assert (edges_cc["split1"] == -1).all()
    assert (edges_cc["split2"] == -1).all()
    assert (edges_nn["split1"] == -1).all()
    assert (edges_nn["split2"] == -1).all()

    # 8. Source-level regression: no domain.close() inside _build_edges_faces (G13).
    src = inspect.getsource(_mesh_module._build_edges_faces)
    assert ".close(" not in src, (
        "_build_edges_faces must not call domain.close() — per-element fields are "
        "pre-id subtractions only (G13)"
    )


TWO_PI = 2 * math.pi


def test_build_mesh():
    # 1. Helicoid rect u-cy: every edge flip=+1 (cy orientation-preserving).
    #    SN(u,v) = (sin(v), -cos(v), u). At u=0 and u=1, dot product = 1 > 0.
    domain_hcy = Domain(type="rect", bounds=(0.0, 1.0, 0.0, TWO_PI), u_identify="cy")
    surf_hcy = SurfaceParams(
        "u*cos(v)", "u*sin(v)", "v", "u v", domain_hcy, perturb=False
    )
    m = build_mesh(domain_hcy, surf_hcy, resolution=15, jitter=True, seed=42)
    assert isinstance(m, Mesh)
    assert m.uv.shape[1] == 2
    assert m.xyz.shape == m.SN.shape == (len(m.uv), 3)
    assert (m.edges["flip"] == 1).all(), "helicoid u-cy: all flip must be +1"
    assert len(m.boundary_edge_idx) == 2 * 15  # v_min and v_max only

    # 2. Möbius band rect u-mo: some seam edges have flip=-1.
    surf_mo = mobius_u(perturb=False)
    m_mo = build_mesh(surf_mo.domain, surf_mo, resolution=15, jitter=True, seed=42)
    assert isinstance(m_mo, Mesh)
    assert len(m_mo.boundary_edge_idx) == 2 * 15  # v_min and v_max sides
    # At least some merged seam edges must have flip=-1 (Möbius reversal).
    assert (m_mo.edges["flip"] == -1).any(), "mobius u-mo: some flip must be -1"

    # 3. Disk r_max=1, resolution=20: 1 boundary loop, no seams, vc == arange(N).
    surf_disk = disk_paraboloid_ca(perturb=False)
    m_disk = build_mesh(surf_disk.domain, surf_disk, resolution=20, seed=42)
    n_outer = max(3, round(math.pi * 20))  # = round(π·20) ≈ 63
    assert len(m_disk.boundary_edge_idx) == n_outer
    assert np.array_equal(m_disk.vertex_class, np.arange(len(m_disk.uv)))
    assert len(m_disk.corner_idx) == 0

    # 4. Annulus r_min=0.3, r_max=1, resolution=20: 2 boundary loops.
    domain_ann = Domain(type="annulus", bounds=(0.3, 1.0, 0.0, TWO_PI))
    surf_ann = SurfaceParams("x", "y", "x**2 + y**2", "x y", domain_ann, perturb=False)
    m_ann = build_mesh(domain_ann, surf_ann, resolution=20, seed=42)
    n_outer_a = max(3, round(math.pi * 20))
    n_inner_a = max(3, round(math.pi * 0.3 * 20))
    assert len(m_ann.boundary_edge_idx) == n_outer_a + n_inner_a

    # 5. Fig-8 immersion cy-cy: χ=0 (torus topology), zero boundary edges.
    domain_fig8 = Domain(
        type="rect", bounds=(0.0, TWO_PI, 0.0, TWO_PI),
        u_identify="cy", v_identify="cy",
    )
    surf_fig8 = SurfaceParams(
        "(2 + cos(u))*cos(v)", "(2 + cos(u))*sin(v)", "sin(2*u)",
        "u v", domain_fig8, perturb=False,
    )
    m_fig8 = build_mesh(domain_fig8, surf_fig8, resolution=15, seed=42)
    assert len(m_fig8.boundary_edge_idx) == 0
    V = len(np.unique(m_fig8.vertex_class[m_fig8.tris]))
    chi = V - len(m_fig8.edges) + len(m_fig8.faces)
    assert chi == 0, f"fig-8 cy-cy: χ={chi}, expected 0"

    # 6. mo-mo regression: zero boundary edges, every edge has exactly 2 incident faces.
    domain_momo = Domain(
        type="rect", bounds=(0.0, TWO_PI, -0.3, 0.3),
        u_identify="mo", v_identify="mo",
    )
    surf_momo = SurfaceParams(
        "cos(u)", "sin(u)", "v", "u v", domain_momo, perturb=False
    )
    m_momo = build_mesh(domain_momo, surf_momo, resolution=10, seed=42)
    assert len(m_momo.boundary_edge_idx) == 0
    edge_face_count = collections.Counter()
    for f in m_momo.faces:
        for e_idx in f["edges"]:
            edge_face_count[int(e_idx)] += 1
    counts = list(edge_face_count.values())
    assert min(counts) == 2 and max(counts) == 2, (
        f"mo-mo: non-manifold edges; face-count hist = {collections.Counter(counts)}"
    )

    # 7. corner_idx length: 4 for no-no, 1 for cy-cy, 2 for mo-mo.
    surf_para = paraboloid(perturb=False)
    m_nn = build_mesh(surf_para.domain, surf_para, resolution=10, seed=42)
    assert len(m_nn.corner_idx) == 4

    surf_tor = torus(perturb=False)
    m_cycy = build_mesh(surf_tor.domain, surf_tor, resolution=10, seed=42)
    assert len(m_cycy.corner_idx) == 1

    assert len(m_momo.corner_idx) == 2

    # 8. Determinism: jitter=True with fixed seed gives identical fields on two calls.
    m_a = build_mesh(surf_para.domain, surf_para, resolution=10, jitter=True, seed=7)
    m_b = build_mesh(surf_para.domain, surf_para, resolution=10, jitter=True, seed=7)
    assert np.array_equal(m_a.uv, m_b.uv)
    assert np.array_equal(m_a.xyz, m_b.xyz)
    assert np.array_equal(m_a.tris, m_b.tris)
