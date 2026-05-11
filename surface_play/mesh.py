import math

import numpy as np
import triangle

from surface_play.domain import Domain


edge_dtype = np.dtype([
    ("p_idx", "i4"), ("q_idx", "i4"),
    ("p", "2f8"),
    ("pq", "2f8"),
    ("f", "i4"), ("g", "i4"),
    ("dir", "2f8"),
    ("flip", "i1"),
    ("split1", "i4"), ("split2", "i4"),
])

face_dtype = np.dtype([
    ("verts", "3i4"),
    ("edges", "3i4"),
    ("p", "2f8"),
    ("pq", "2f8"),
    ("pr", "2f8"),
])


def _boundary_edge_count(tris: np.ndarray) -> int:
    """Count edges that appear in exactly one triangle (boundary edges)."""
    all_edges = np.sort(
        np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1
    )
    _, counts = np.unique(all_edges, axis=0, return_counts=True)
    return int((counts == 1).sum())


def _generate_rect_mesh(domain: Domain, resolution: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns:
      uv:   (N, 2) float64 — vertex (u, v) coordinates.
      tris: (M, 3) int32   — triangle vertex indices.
    No identifications applied yet (paired vertices have distinct indices). No jitter.
    """
    u_min, u_max, v_min, v_max = domain.bounds
    du = (u_max - u_min) / resolution
    dv = (v_max - v_min) / resolution
    n_bv = 4 * resolution

    bverts = np.empty((n_bv, 2), dtype=np.float64)
    # Bottom: u from u_min..u_max-du, v=v_min
    for i in range(resolution):
        bverts[i] = [u_min + i * du, v_min]
    # Right: u=u_max, v from v_min..v_max-dv
    for i in range(resolution):
        bverts[resolution + i] = [u_max, v_min + i * dv]
    # Top: u from u_max..u_min+du, v=v_max
    for i in range(resolution):
        bverts[2 * resolution + i] = [u_max - i * du, v_max]
    # Left: u=u_min, v from v_max..v_min+dv
    for i in range(resolution):
        bverts[3 * resolution + i] = [u_min, v_max - i * dv]

    idx = np.arange(n_bv, dtype=np.int32)
    bsegs = np.column_stack([idx, np.roll(idx, -1)])

    rect_area = (u_max - u_min) * (v_max - v_min)
    area = rect_area / resolution ** 2

    pslg = {"vertices": bverts, "segments": bsegs}
    result = triangle.triangulate(pslg, opts=f"pYq30a{area:.17g}")

    uv = np.asarray(result["vertices"], dtype=np.float64)
    tris = np.asarray(result["triangles"], dtype=np.int32)
    return uv, tris


SIDE_NONE = 0
SIDE_U_MIN = 1
SIDE_U_MAX = 2
SIDE_V_MIN = 3
SIDE_V_MAX = 4


def _build_edges_faces(
    uv: np.ndarray,
    tris: np.ndarray,
    vertex_class: np.ndarray,
    SN_per_vertex: np.ndarray,
    domain: Domain,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build edges and faces structured arrays. `tris` carries PRE-id vertex labels
    (from C4). Edges merge across identified sides by the geometric rule of G13 —
    two pre-id edges merge iff both endpoints lie on a single identified side and
    their endpoints correspond under that side's pairing. Faces are never merged.
    Per-element geometric fields are direct pre-id subtractions (no seam lifting).
    """
    N = len(uv)
    M = len(tris)

    is_rect = domain.type == "rect"
    u_identified = is_rect and domain.u_identify in ("cy", "mo")
    v_identified = is_rect and domain.v_identify in ("cy", "mo")

    on_u_min = np.zeros(N, dtype=bool)
    on_u_max = np.zeros(N, dtype=bool)
    on_v_min = np.zeros(N, dtype=bool)
    on_v_max = np.zeros(N, dtype=bool)
    if is_rect:
        n = _boundary_edge_count(tris) // 4
        if n > 0:
            on_v_min[0:n] = True
            on_u_max[n:2 * n] = True
            on_v_max[2 * n:3 * n] = True
            on_u_min[3 * n:4 * n] = True
            on_u_min[0] = True
            on_v_min[n] = True
            on_u_max[2 * n] = True
            on_v_max[3 * n] = True

    def edge_side(a: int, b: int) -> int:
        if u_identified:
            if on_u_min[a] and on_u_min[b]:
                return SIDE_U_MIN
            if on_u_max[a] and on_u_max[b]:
                return SIDE_U_MAX
        if v_identified:
            if on_v_min[a] and on_v_min[b]:
                return SIDE_V_MIN
            if on_v_max[a] and on_v_max[b]:
                return SIDE_V_MAX
        return SIDE_NONE

    records: dict[tuple[int, int], list[tuple[int, int, int, int, int]]] = {}
    face_edge_idx = np.empty((M, 3), dtype=np.int32)

    def key_of(a: int, b: int, side: int) -> tuple[int, int]:
        if side == SIDE_NONE:
            return (a, b) if a < b else (b, a)
        ca = int(vertex_class[a])
        cb = int(vertex_class[b])
        return (ca, cb) if ca < cb else (cb, ca)

    for f_idx in range(M):
        i = int(tris[f_idx, 0])
        j = int(tris[f_idx, 1])
        k = int(tris[f_idx, 2])
        local_edges = ((i, j, k), (j, k, i), (k, i, j))
        for local_pos, (a, b, third) in enumerate(local_edges):
            side = edge_side(a, b)
            key = key_of(a, b, side)
            rec_list = records.get(key)
            if rec_list is None:
                rec_list = []
                records[key] = rec_list
            rec_list.append((f_idx, third, a, b, side))

    key_to_idx = {key: idx for idx, key in enumerate(records.keys())}

    for f_idx in range(M):
        i = int(tris[f_idx, 0])
        j = int(tris[f_idx, 1])
        k = int(tris[f_idx, 2])
        local_edges = ((i, j, k), (j, k, i), (k, i, j))
        for local_pos, (a, b, _third) in enumerate(local_edges):
            side = edge_side(a, b)
            key = key_of(a, b, side)
            face_edge_idx[f_idx, local_pos] = key_to_idx[key]

    n_edges = len(records)
    edges = np.zeros(n_edges, dtype=edge_dtype)

    canonical_sides = (SIDE_U_MIN, SIDE_V_MIN)

    for key, recs in records.items():
        e_idx = key_to_idx[key]
        if len(recs) > 2:
            raise ValueError(
                f"Edge key {key} has {len(recs)} face references (expected 1 or 2). "
                f"Records: {recs}"
            )

        if len(recs) == 2:
            r0, r1 = recs
            if r0[4] in canonical_sides:
                canonical_rec, partner_rec = r0, r1
            elif r1[4] in canonical_sides:
                canonical_rec, partner_rec = r1, r0
            else:
                canonical_rec, partner_rec = r0, r1
        else:
            canonical_rec, partner_rec = recs[0], None

        f_idx, third, a, b, side = canonical_rec
        if a < b:
            p_idx, q_idx = a, b
        else:
            p_idx, q_idx = b, a

        p = uv[p_idx]
        pq = uv[q_idx] - p
        edges[e_idx]["p_idx"] = p_idx
        edges[e_idx]["q_idx"] = q_idx
        edges[e_idx]["p"] = p
        edges[e_idx]["pq"] = pq

        if partner_rec is None:
            edges[e_idx]["f"] = f_idx
            edges[e_idx]["g"] = -1
            d = np.array([-pq[1], pq[0]], dtype=np.float64)
            offset = uv[third] - p
            if float(np.dot(d, offset)) < 0.0:
                d = -d
            nrm = float(np.linalg.norm(d))
            if nrm > 0.0:
                d = d / nrm
            edges[e_idx]["dir"] = d
        else:
            edges[e_idx]["f"] = f_idx
            edges[e_idx]["g"] = partner_rec[0]
            edges[e_idx]["dir"] = (0.0, 0.0)

        if partner_rec is not None and side != SIDE_NONE:
            pa, pb = partner_rec[2], partner_rec[3]
            if int(vertex_class[pa]) == int(vertex_class[p_idx]):
                p_partner = pa
            else:
                p_partner = pb
            edges[e_idx]["flip"] = (
                1 if float(np.dot(SN_per_vertex[p_idx], SN_per_vertex[p_partner])) > 0.0
                else -1
            )
        else:
            edges[e_idx]["flip"] = (
                1 if float(np.dot(SN_per_vertex[p_idx], SN_per_vertex[q_idx])) > 0.0
                else -1
            )

        edges[e_idx]["split1"] = -1
        edges[e_idx]["split2"] = -1

    faces = np.zeros(M, dtype=face_dtype)
    for f_idx in range(M):
        i = int(tris[f_idx, 0])
        j = int(tris[f_idx, 1])
        k = int(tris[f_idx, 2])
        p = uv[i]
        faces[f_idx]["verts"] = (i, j, k)
        faces[f_idx]["edges"] = face_edge_idx[f_idx]
        faces[f_idx]["p"] = p
        faces[f_idx]["pq"] = uv[j] - p
        faces[f_idx]["pr"] = uv[k] - p

    return edges, faces


def _jitter(
    uv: np.ndarray,
    tris: np.ndarray,
    domain: Domain,
    *,
    seed: int | None = None,
) -> np.ndarray:
    """
    Returns a jittered copy of uv. Independent jitter per vertex by ±0.1% of typical
    edge length. Reprojects true boundary vertices onto the domain boundary. Run BEFORE
    `_apply_identifications` (C4): paired vertices receive independent jitter — the
    discarded one is dropped at identification, so coherence is unnecessary.
    RNG seeded by `seed` (default = id(uv) for repro).
    """
    if seed is None:
        seed = id(uv)

    e0 = np.linalg.norm(uv[tris[:, 1]] - uv[tris[:, 0]], axis=1)
    e1 = np.linalg.norm(uv[tris[:, 2]] - uv[tris[:, 1]], axis=1)
    e2 = np.linalg.norm(uv[tris[:, 0]] - uv[tris[:, 2]], axis=1)
    L = float(np.median(np.concatenate([e0, e1, e2])))

    rng = np.random.default_rng(seed)
    delta = rng.uniform(-1.0, 1.0, uv.shape) * 0.001 * L
    uv_new = uv + delta

    tol = 1e-12
    if domain.type == "rect":
        u_min, u_max, v_min, v_max = domain.bounds
        if domain.u_identify == "no":
            uv_new[np.abs(uv[:, 0] - u_min) < tol, 0] = u_min
            uv_new[np.abs(uv[:, 0] - u_max) < tol, 0] = u_max
        if domain.v_identify == "no":
            uv_new[np.abs(uv[:, 1] - v_min) < tol, 1] = v_min
            uv_new[np.abs(uv[:, 1] - v_max) < tol, 1] = v_max
    elif domain.type in ("disk", "annulus"):
        r_min, r_max = domain.bounds[0], domain.bounds[1]
        r = np.sqrt(uv[:, 0] ** 2 + uv[:, 1] ** 2)
        ring_radii = [r_max] + ([r_min] if r_min > 0 else [])
        for target_r in ring_radii:
            mask = np.abs(r - target_r) < tol
            r_new = np.sqrt(uv_new[mask, 0] ** 2 + uv_new[mask, 1] ** 2)
            uv_new[mask] = uv_new[mask] / r_new[:, None] * target_r

    return uv_new


def _apply_identifications(
    uv_jittered: np.ndarray,
    tris: np.ndarray,
    domain: Domain,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute vertex equivalence classes induced by the side-identification rules.
    `tris` is NOT relabeled — pre-id indices are preserved (G13). Identification of
    edges and faces is a geometric matter handled in C5, not a label rewrite here.

    Returns (uv_unchanged, tris_filtered, vertex_class):
      - tris_filtered drops any triangle whose three vertices map to fewer than three
        distinct equivalence classes (geometrically degenerate post-id). Pre-id labels
        are preserved on surviving rows.
      - vertex_class[i] is the canonical (smallest pre-id index) representative of i's
        class. Vertices not on identified sides keep vertex_class[i] == i.
    """
    N = len(uv_jittered)

    if domain.type != "rect" or (
        domain.u_identify == "no" and domain.v_identify == "no"
    ):
        return uv_jittered, tris.copy(), np.arange(N, dtype=np.int32)

    n = _boundary_edge_count(tris) // 4

    pairs: list[tuple[int, int]] = []

    if domain.u_identify == "cy":
        pairs.append((0, n))
        for i in range(n):
            pairs.append((3 * n + i, 2 * n - i))
    elif domain.u_identify == "mo":
        pairs.append((0, 2 * n))
        for i in range(n):
            pairs.append((3 * n + i, n + i))

    if domain.v_identify == "cy":
        pairs.append((n, 2 * n))
        for i in range(n):
            pairs.append((i, 3 * n - i))
    elif domain.v_identify == "mo":
        pairs.append((n, 3 * n))
        for i in range(n):
            pairs.append((i, 2 * n + i))

    parent = np.arange(N, dtype=np.int32)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            if ra > rb:
                ra, rb = rb, ra
            parent[rb] = ra

    vertex_class = np.array([find(i) for i in range(N)], dtype=np.int32)

    classes_per_tri = vertex_class[tris]
    distinct = (
        (classes_per_tri[:, 0] != classes_per_tri[:, 1])
        & (classes_per_tri[:, 1] != classes_per_tri[:, 2])
        & (classes_per_tri[:, 0] != classes_per_tri[:, 2])
    )
    tris_filtered = tris[distinct].astype(np.int32)

    return uv_jittered, tris_filtered, vertex_class


def _generate_disk_mesh(domain: Domain, resolution: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Disk (r_min == 0) or annulus (r_min > 0). Returns (uv, tris) in cartesian (u, v).
    """
    r_min, r_max = domain.bounds[0], domain.bounds[1]
    bbox_diag = domain.bbox_diag()  # = 2 * r_max

    n_outer = max(3, round(2 * math.pi * r_max * resolution / bbox_diag))
    theta_outer = np.linspace(0, 2 * math.pi, n_outer, endpoint=False)
    outer_verts = np.column_stack(
        [r_max * np.cos(theta_outer), r_max * np.sin(theta_outer)]
    )
    idx_o = np.arange(n_outer, dtype=np.int32)
    outer_segs = np.column_stack([idx_o, np.roll(idx_o, -1)])

    if r_min > 0:
        n_inner = max(3, round(2 * math.pi * r_min * resolution / bbox_diag))
        theta_inner = np.linspace(0, 2 * math.pi, n_inner, endpoint=False)
        inner_verts = np.column_stack(
            [r_min * np.cos(theta_inner), r_min * np.sin(theta_inner)]
        )
        idx_i = n_outer + np.arange(n_inner, dtype=np.int32)
        inner_segs = np.column_stack([idx_i, np.roll(idx_i, -1)])

        all_verts = np.vstack([outer_verts, inner_verts])
        all_segs = np.vstack([outer_segs, inner_segs])
        domain_area = math.pi * (r_max ** 2 - r_min ** 2)
        pslg = {
            "vertices": all_verts,
            "segments": all_segs,
            "holes": np.array([[0.0, 0.0]]),
        }
    else:
        all_verts = outer_verts
        all_segs = outer_segs
        domain_area = math.pi * r_max ** 2
        pslg = {"vertices": all_verts, "segments": all_segs}

    area = domain_area / resolution ** 2
    result = triangle.triangulate(pslg, opts=f"pYq30a{area:.17g}")

    uv = np.asarray(result["vertices"], dtype=np.float64)
    tris = np.asarray(result["triangles"], dtype=np.int32)
    return uv, tris
