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


def test_apply_identifications():
    res = 10

    def mesh_rect(u_id, v_id):
        d = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0),
                   u_identify=u_id, v_identify=v_id)
        uv, tris = _generate_rect_mesh(d, resolution=res)
        return d, uv, tris

    # 1. no-no: pass-through (no compaction, tris unchanged).
    d, uv, tris = mesh_rect("no", "no")
    uv_out, tris_out, on_u, on_v, corners = _apply_identifications(uv, tris, d)
    assert np.array_equal(uv, uv_out), "no-no: uv unchanged"
    assert np.array_equal(tris, tris_out), "no-no: tris unchanged"
    assert not on_u.any() and not on_v.any()
    assert len(corners) == 4, "no-no: 4 distinct corners"

    # 2. cy-no: u_min/u_max merge. K = N - (n + 1) since each of n interior u-pairs and
    #    the corner pair (0,n) collapses one duplicate. (Corner (2n, 3n) also collapses.)
    d, uv, tris = mesh_rect("cy", "no")
    uv_c, tris_c, on_u, on_v, corners = _apply_identifications(uv, tris, d)
    n = res
    # Compacted uv size: original N minus (n + 1) merged duplicates
    # (n interior u-side pairs + 1 corner pair collapses, plus the (2n, 3n) corner pair).
    assert len(uv_c) == len(uv) - (n + 1), (
        f"cy-no res={res}: compacted size {len(uv_c)} != {len(uv) - (n + 1)}"
    )
    # tris indices must all be < len(uv_c).
    assert tris_c.max() < len(uv_c), "cy-no: tris indices must be in compacted range"
    # u-seam membership: 2 * (n + 1) - 2 = 2n vertices on a u-side -> after merge half remain
    # but corner survives once -> total on_u = 2n (n+1 left-side classes merged with right).
    # Just check at least some are on_u and none on_v.
    assert on_u.sum() > 0 and not on_v.any()
    assert len(corners) == 2, "cy-no: top-left/bottom-left corners survive as 2 distinct"

    # 3. cy-cy res=10: corners all merge to one class.
    d, uv, tris = mesh_rect("cy", "cy")
    uv_c, tris_c, on_u, on_v, corners = _apply_identifications(uv, tris, d)
    assert len(corners) == 1, "cy-cy: all 4 corners merge to 1 class"
    assert tris_c.max() < len(uv_c)
    # All boundary classes are on at least one seam.
    # Every original boundary vertex appears in either on_u or on_v of its class.
    assert (on_u | on_v).sum() > 0

    # 4. mo-mo rejected.
    d_mm = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0),
                  u_identify="mo", v_identify="mo")
    uv_mm, tris_mm = _generate_rect_mesh(d_mm, resolution=5)
    with pytest.raises(ValueError, match="mo, mo"):
        _apply_identifications(uv_mm, tris_mm, d_mm)

    # 5. Disk: pass-through, no corners.
    d_disk = Domain(type="disk", bounds=(0.0, 1.0, 0.0, 2 * math.pi))
    uv_d, tris_d = _generate_disk_mesh(d_disk, resolution=20)
    uv_d_out, tris_d_out, on_u_d, on_v_d, corners_d = _apply_identifications(
        uv_d, tris_d, d_disk
    )
    assert np.array_equal(uv_d, uv_d_out)
    assert np.array_equal(tris_d, tris_d_out)
    assert len(corners_d) == 0

    # 6. Pairing is index-based, immune to uv jitter.
    d, uv, tris = mesh_rect("cy", "no")
    uv_j = _jitter(uv, tris, d, seed=42)
    _, tris_clean, _, _, _ = _apply_identifications(uv, tris, d)
    _, tris_jit, _, _, _ = _apply_identifications(uv_j, tris, d)
    assert np.array_equal(tris_clean, tris_jit), \
        "identification must be index-based, not coordinate-based"

    # 7. Determinism: same input → same compaction.
    _, tris_a, _, _, _ = _apply_identifications(uv, tris, d)
    _, tris_b, _, _, _ = _apply_identifications(uv, tris, d)
    assert np.array_equal(tris_a, tris_b)


def _build_post_id(u_id, v_id, resolution, seed=42):
    """Helper: generate raw mesh, jitter, apply identifications.

    Returns ``(domain, uv_compacted, tris_canonical, on_u_seam, on_v_seam,
    uv_pre, tris_pre)`` — the trailing pair are the pre-identification
    `_jitter`-ed uv and the original `_generate_rect_mesh` tris (needed by
    `_build_edges_faces` for pre-id pq/pr subtractions, per roadmap line 552).
    """
    d = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0),
               u_identify=u_id, v_identify=v_id)
    uv, tris = _generate_rect_mesh(d, resolution=resolution)
    uv_j = _jitter(uv, tris, d, seed=seed)
    uv_c, tris_c, on_u, on_v, _ = _apply_identifications(uv_j, tris, d)
    return d, uv_c, tris_c, on_u, on_v, uv_j, tris


def test_build_edges_faces():
    import collections

    # 1. Boundary counts for each supported (u_id, v_id) combo at res=5.
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
    }
    expected_chi = {
        ("no", "no"): 1,
        ("cy", "no"): 0, ("no", "cy"): 0,
        ("mo", "no"): 0, ("no", "mo"): 0,
        ("cy", "cy"): 0, ("cy", "mo"): 0, ("mo", "cy"): 0,
    }

    for (u_id, v_id), bnd_expected in expected_boundary.items():
        d, uv, tris, on_u, on_v, uv_pre, tris_pre = _build_post_id(u_id, v_id, resolution=res)
        SN = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
        edges, faces = _build_edges_faces(uv, tris, uv_pre, tris_pre, SN, on_u, on_v, d)
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

        # Euler χ = V - E + F where V is the post-compaction vertex count used by tris.
        V = len(np.unique(tris))
        chi = V - len(edges) + len(faces)
        assert chi == expected_chi[(u_id, v_id)], (
            f"{u_id}-{v_id}: χ={chi}, expected {expected_chi[(u_id, v_id)]}"
        )

    # 2. flip on cy-no with constant SN: every edge has flip=+1.
    d, uv, tris, on_u, on_v, uv_pre, tris_pre = _build_post_id("cy", "no", resolution=res)
    SN_const = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
    edges_cy, _ = _build_edges_faces(uv, tris, uv_pre, tris_pre, SN_const, on_u, on_v, d)
    assert (edges_cy["flip"] == 1).all(), (
        "cy-no with constant SN: every edge should have flip=+1"
    )

    # 3. flip on mo-no: with SN that flips sign across the mo-seam (mimicking
    # the Möbius non-orientability), edges crossing identification copies have
    # flip=-1.  Set SN at vertices with u > 0.5 to +z, u <= 0.5 to -z — so
    # canonical merged seam vertices on opposite u-copies disagree on SN sign.
    d, uv, tris, on_u, on_v, uv_pre, tris_pre = _build_post_id("mo", "no", resolution=res)
    SN = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
    SN[:, 2] = np.where(uv[:, 0] > 0.5, 1.0, -1.0)
    edges_mo, _ = _build_edges_faces(uv, tris, uv_pre, tris_pre, SN, on_u, on_v, d)
    assert (edges_mo["flip"] == -1).any(), (
        "mo-no: SN flipping across the seam should produce some flip=-1 edges"
    )
    assert (edges_mo["flip"] == 1).sum() > (edges_mo["flip"] == -1).sum()

    # 4. dir: boundary edges have dir ⟂ pq and pointing toward third vertex.
    d, uv, tris, on_u, on_v, uv_pre, tris_pre = _build_post_id("no", "no", resolution=res)
    SN = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
    edges_nn, _ = _build_edges_faces(uv, tris, uv_pre, tris_pre, SN, on_u, on_v, d)
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

    # 5. p, pq, pr are pre-identification subtractions (roadmap line 552).
    #    Geometric contract: every edge/face pq, pr is the *short* vector
    #    between its endpoints, bounded above by a few mesh edge lengths.
    #    Compacted subtraction `uv[j] - uv[i]` agrees on non-seam edges and
    #    differs by exactly ±period (modulo jitter noise) on seam-crossing
    #    edges.
    d, uv, tris, on_u, on_v, uv_pre, tris_pre = _build_post_id("cy", "cy", resolution=res)
    SN = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
    edges_cc, faces_cc = _build_edges_faces(uv, tris, uv_pre, tris_pre, SN, on_u, on_v, d)
    typical = max(d.period_u, d.period_v) / res  # ~ regular-grid edge length

    for e in edges_cc:
        i, j = int(e["p_idx"]), int(e["q_idx"])
        assert np.array_equal(e["p"], uv[i]), "edge p must equal uv[p_idx] exactly"
        # The pre-id pq must be a short vector: well under the domain extent.
        assert float(np.linalg.norm(e["pq"])) < 5.0 * typical, (
            f"edge {i}-{j}: pre-id pq must be short, got |pq|="
            f"{np.linalg.norm(e['pq']):.4f}, expected ≤ {5*typical:.4f}"
        )
    for f in faces_cc:
        i, j, k = int(f["verts"][0]), int(f["verts"][1]), int(f["verts"][2])
        assert np.array_equal(f["p"], uv[i])
        assert float(np.linalg.norm(f["pq"])) < 5.0 * typical
        assert float(np.linalg.norm(f["pr"])) < 5.0 * typical

    # Seam-crossing edges *must* exist on cy-cy and *must* differ from the
    # naive compacted subtraction (otherwise the fix is a no-op).
    seam_diff_count = 0
    for e in edges_cc:
        i, j = int(e["p_idx"]), int(e["q_idx"])
        naive = uv[j] - uv[i]
        if not np.allclose(e["pq"], naive, atol=1e-10):
            seam_diff_count += 1
    assert seam_diff_count > 0, (
        "expected at least one seam-crossing edge to use pre-id pq distinct "
        "from `uv[q_idx] - uv[p_idx]` on a cy-cy mesh"
    )

    # 6. Splits initialized to -1.
    assert (edges_cc["split1"] == -1).all()
    assert (edges_cc["split2"] == -1).all()
    assert (edges_nn["split1"] == -1).all()
    assert (edges_nn["split2"] == -1).all()

    # 7. Manifoldness for closed surfaces: every edge has 2 incident faces.
    for combo in (("cy", "cy"), ("cy", "mo"), ("mo", "cy")):
        d, uv, tris, on_u, on_v, uv_pre, tris_pre = _build_post_id(*combo, resolution=res)
        SN = np.tile(np.array([0.0, 0.0, 1.0]), (len(uv), 1))
        edges, faces = _build_edges_faces(uv, tris, uv_pre, tris_pre, SN, on_u, on_v, d)
        edge_face_count = collections.Counter()
        for f in faces:
            for e_idx in f["edges"]:
                edge_face_count[int(e_idx)] += 1
        counts = sorted(edge_face_count.values())
        assert min(counts) == 2 and max(counts) == 2, (
            f"{combo}: non-manifold edges, face-count hist = "
            f"{collections.Counter(counts)}"
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

    # 3. Disk r_max=1, resolution=20: 1 boundary loop, no seams, no corners.
    surf_disk = disk_paraboloid_ca(perturb=False)
    m_disk = build_mesh(surf_disk.domain, surf_disk, resolution=20, seed=42)
    n_outer = max(3, round(math.pi * 20))  # = round(π·20) ≈ 63
    assert len(m_disk.boundary_edge_idx) == n_outer
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
    V = len(np.unique(m_fig8.tris))
    chi = V - len(m_fig8.edges) + len(m_fig8.faces)
    assert chi == 0, f"fig-8 cy-cy: χ={chi}, expected 0"

    # 6. mo-mo rejected.
    domain_momo = Domain(
        type="rect", bounds=(0.0, TWO_PI, -0.3, 0.3),
        u_identify="mo", v_identify="mo",
    )
    surf_momo = SurfaceParams(
        "cos(u)", "sin(u)", "v", "u v", domain_momo, perturb=False
    )
    with pytest.raises(ValueError, match="mo, mo"):
        build_mesh(domain_momo, surf_momo, resolution=10, seed=42)

    # 7. corner_idx length: 4 for no-no, 1 for cy-cy.
    surf_para = paraboloid(perturb=False)
    m_nn = build_mesh(surf_para.domain, surf_para, resolution=10, seed=42)
    assert len(m_nn.corner_idx) == 4

    surf_tor = torus(perturb=False)
    m_cycy = build_mesh(surf_tor.domain, surf_tor, resolution=10, seed=42)
    assert len(m_cycy.corner_idx) == 1

    # 8. Determinism: jitter=True with fixed seed gives identical fields on two calls.
    m_a = build_mesh(surf_para.domain, surf_para, resolution=10, jitter=True, seed=7)
    m_b = build_mesh(surf_para.domain, surf_para, resolution=10, jitter=True, seed=7)
    assert np.array_equal(m_a.uv, m_b.uv)
    assert np.array_equal(m_a.xyz, m_b.xyz)
    assert np.array_equal(m_a.tris, m_b.tris)


def test_antipodal_disk():
    """Antipodal boundary gluing closes the unit disk into ℝP²."""
    TWO_PI = 2 * math.pi

    # Domain validation: antipodal is disk/annulus only — but ANY outer radius
    # (no r_max=1 requirement; the involution is σ_R(z) = -R²/z̄).
    with pytest.raises(ValueError, match="only valid on a disk"):
        Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0),
               boundary_identify="antipodal")
    # Annulus is allowed: outer glued, inner boundary kept (Möbius band, χ=0).
    dom_ann = Domain(type="annulus", bounds=(0.3, 1.0, 0.0, TWO_PI),
                     boundary_identify="antipodal")
    surf_ann = SurfaceParams("u", "v", "u*u + v*v", "u v", dom_ann, perturb=False)
    m_ann = build_mesh(dom_ann, surf_ann, resolution=16, jitter=True, seed=1)
    chi_ann = len(m_ann.uv) - len(m_ann.edges) + len(m_ann.faces)
    assert chi_ann == 0, f"antipodal annulus (Möbius band): χ={chi_ann}, expected 0"
    assert len(m_ann.boundary_edge_idx) > 0, "antipodal annulus keeps inner boundary"
    assert (m_ann.edges["flip"] == -1).any()
    # Non-unit outer radius closes to ℝP² just the same.
    dom_r2 = Domain(type="disk", bounds=(0.0, 2.0, 0.0, TWO_PI),
                    boundary_identify="antipodal")
    m_r2 = build_mesh(dom_r2, SurfaceParams("u", "v", "u+v", "u v", dom_r2,
                                            perturb=False), resolution=16, seed=1)
    assert len(m_r2.uv) - len(m_r2.edges) + len(m_r2.faces) == 1
    assert len(m_r2.boundary_edge_idx) == 0

    dom = Domain(type="disk", bounds=(0.0, 1.0, 0.0, TWO_PI),
                 boundary_identify="antipodal")
    assert dom.is_antipodal

    # σ_R involution (R=1): σ(z) = -z/|z|²; localize picks the closer of {q, σ(q)}.
    assert np.allclose(dom._sigma(np.array([[0.5, 0.0]])), [[-2.0, 0.0]])  # interior→exterior
    # A near-boundary point: localize keeps the SAME-side rep (q here), not σ(q).
    p_i = np.array([0.1, 0.2]); q_i = np.array([0.15, 0.18])
    assert np.allclose(dom.localize(p_i, q_i), [0.15, 0.18])
    # interpolate maps an outside result back inside the disk.
    res = dom.interpolate(np.array([0.9, 0.0]), np.array([-0.9, 0.0]), 1.0)
    assert np.hypot(*res) <= 1.0 + 1e-9

    surf = SurfaceParams("u", "v", "u*u + v*v", "u v", dom, perturb=False)
    m = build_mesh(dom, surf, resolution=16, jitter=True, seed=1)

    # ℝP²: χ = V − E + F = 1; no boundary; some orientation-reversing seam edges.
    V = len(m.uv)
    chi = V - len(m.edges) + len(m.faces)
    assert chi == 1, f"antipodal disk: χ={chi}, expected 1 (ℝP²)"
    assert len(m.boundary_edge_idx) == 0, "antipodal disk: boundary must be glued"
    assert (m.edges["flip"] == -1).any(), "antipodal seam edges must carry flip=-1"
    # Every edge is interior (shared by 2 faces) on a closed surface.
    assert (m.edges["g"] >= 0).all()
    assert np.isfinite(m.xyz).all()

    # Determinism.
    m2 = build_mesh(dom, surf, resolution=16, jitter=True, seed=1)
    assert np.array_equal(m.uv, m2.uv) and np.array_equal(m.tris, m2.tris)
