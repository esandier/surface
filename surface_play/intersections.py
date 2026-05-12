"""intersections.py — P5: sweep_segments; C8: candidate_pairs (BVH broad phase);
C9: find_double_points; C10: build_sis_pairs.

Shared 2D segment-segment intersection primitive used by the domain sweep
(self-intersection preimages, CS×SIS, etc.) and the view-plane sweep.
See Modular_rewrite_roadmap.md §P5; G3 (close-aware) is enforced when
`domain` is a rectangular domain with identification.

C8: `candidate_pairs(e_bbox, f_bbox)` returns overlapping AABB pairs via a
stack-based BVH with axis-cycling partitioning, Numba-JIT compiled. See G14.

C9: `find_double_points(mesh, surface)` returns DP records (edge-vs-face hits
in 3D), with sibling-pair consumption (G2) and seam-aware vertex-class skip
(G13). Uses C5 per-element (p, pq, pr) for uv recovery (G15).

C10: `build_sis_pairs(dps)` returns SIS records (DP-DP edges of the
self-intersection graph) with vectorized A1/A2 face-sharing test and flip
detection. See G17 for split sentinel.
"""

import warnings
from typing import Optional

import numpy as np
from numba import njit


intersect_dtype = np.dtype([
    ("uv",  "f8", 2),
    ("a",   "i4"),
    ("b",   "i4"),
    ("t_a", "f8"),
    ("t_b", "f8"),
])


dp_dtype = np.dtype([
    ("xyz",         "f8", 3),
    ("uv1",         "f8", 2),
    ("uv2",         "f8", 2),
    ("E1",          "i4"),
    ("E2",          "i4"),
    ("F2",          "i4"),
    ("A1",          "i4", 2),
    ("A2",          "i4", 2),
    ("type",        "U2"),
    ("ptype",       "i4"),
    ("on_boundary", "?"),
])


sis_dtype = np.dtype([
    ("p_dp",   "i4"),
    ("q_dp",   "i4"),
    ("flip",   "i1"),
    ("split1", "i4"),
    ("split2", "i4"),
])


def _half_period_adjust(diff: np.ndarray, period: float) -> np.ndarray:
    return (diff + 0.5 * period) % period - 0.5 * period


def _close_pair_array(p0: np.ndarray, p1: np.ndarray, domain) -> np.ndarray:
    """Vectorized `domain.close(p0, p1)` for arrays of shape (..., 2)."""
    p1c = np.array(p1, dtype=float, copy=True)
    if domain is None or getattr(domain, "type", None) != "rect":
        return p1c
    d = p1c - p0
    if domain.u_identify in ("cy", "mo"):
        d[..., 0] = _half_period_adjust(d[..., 0], domain.period_u)
    if domain.v_identify in ("cy", "mo"):
        d[..., 1] = _half_period_adjust(d[..., 1], domain.period_v)
    return p0 + d


def sweep_segments(
    seg_a_uv0: np.ndarray,
    seg_a_uv1: np.ndarray,
    seg_b_uv0: Optional[np.ndarray],
    seg_b_uv1: Optional[np.ndarray],
    domain,                                # Optional[Domain]
    *,
    self_sweep: bool = False,
    tol: float = 1e-9,
) -> np.ndarray:
    """Find 2D segment-segment intersections; structured array of `intersect_dtype`.

    If `domain` is a rectangular Domain with identification, segment endpoints
    are anchored close-aware (G3): each segment's q is shifted toward its own p,
    and at pairing time b's anchor is shifted into a's frame via half-period
    adjustment. If `domain` is None or non-identified, the sweep is flat 2D.

    self_sweep=True requires seg_b_*=None; A is swept against itself with the
    `i < j` filter. Acceptance is strict: t_a, t_b ∈ (tol, 1-tol).
    """
    seg_a_uv0 = np.asarray(seg_a_uv0, dtype=float).reshape(-1, 2)
    seg_a_uv1 = np.asarray(seg_a_uv1, dtype=float).reshape(-1, 2)

    if self_sweep:
        if seg_b_uv0 is not None or seg_b_uv1 is not None:
            raise ValueError("self_sweep=True requires seg_b_uv0 and seg_b_uv1 to be None")
        seg_b_uv0 = seg_a_uv0
        seg_b_uv1 = seg_a_uv1
    else:
        if seg_b_uv0 is None or seg_b_uv1 is None:
            raise ValueError("cross-sweep requires seg_b_uv0 and seg_b_uv1")
        seg_b_uv0 = np.asarray(seg_b_uv0, dtype=float).reshape(-1, 2)
        seg_b_uv1 = np.asarray(seg_b_uv1, dtype=float).reshape(-1, 2)

    n = seg_a_uv0.shape[0]
    m = seg_b_uv0.shape[0]
    if n == 0 or m == 0:
        return np.empty(0, dtype=intersect_dtype)

    seg_a_uv1_loc = _close_pair_array(seg_a_uv0, seg_a_uv1, domain)
    seg_b_uv1_loc = _close_pair_array(seg_b_uv0, seg_b_uv1, domain)
    a_dir = seg_a_uv1_loc - seg_a_uv0          # (n, 2)
    b_dir = seg_b_uv1_loc - seg_b_uv0          # (m, 2)

    use_close = (domain is not None and getattr(domain, "type", None) == "rect")
    u_id = getattr(domain, "u_identify", "no") if domain is not None else "no"
    v_id = getattr(domain, "v_identify", "no") if domain is not None else "no"

    out_uv: list[np.ndarray] = []
    out_a: list[np.ndarray] = []
    out_b: list[np.ndarray] = []
    out_ta: list[np.ndarray] = []
    out_tb: list[np.ndarray] = []
    batch = max(1, 16384 // max(m, 1))

    for i0 in range(0, n, batch):
        i1 = min(i0 + batch, n)
        a_anchor = seg_a_uv0[i0:i1, None, :]            # (B, 1, 2)
        a_d = a_dir[i0:i1, None, :]                     # (B, 1, 2)
        b_anchor = seg_b_uv0[None, :, :]                # (1, m, 2)
        b_d = b_dir[None, :, :]                         # (1, m, 2)

        diff = b_anchor - a_anchor                      # (B, m, 2), broadcast copy
        diff = np.broadcast_to(diff, (i1 - i0, m, 2)).copy()
        if use_close:
            if u_id in ("cy", "mo"):
                diff[..., 0] = _half_period_adjust(diff[..., 0], domain.period_u)
            if v_id in ("cy", "mo"):
                diff[..., 1] = _half_period_adjust(diff[..., 1], domain.period_v)
        b_anchor_local = a_anchor + diff                # (B, m, 2)
        b_end_local = b_anchor_local + b_d              # (B, m, 2)

        a_end = a_anchor + a_d
        a_lo = np.minimum(a_anchor, a_end)
        a_hi = np.maximum(a_anchor, a_end)
        b_lo = np.minimum(b_anchor_local, b_end_local)
        b_hi = np.maximum(b_anchor_local, b_end_local)
        overlap = (
            (a_lo[..., 0] <= b_hi[..., 0] + tol)
            & (b_lo[..., 0] <= a_hi[..., 0] + tol)
            & (a_lo[..., 1] <= b_hi[..., 1] + tol)
            & (b_lo[..., 1] <= a_hi[..., 1] + tol)
        )

        if self_sweep:
            i_idx = np.arange(i0, i1)[:, None]
            j_idx = np.arange(m)[None, :]
            overlap = overlap & (i_idx < j_idx)

        if not overlap.any():
            continue

        ii_w, jj_w = np.where(overlap)
        p0 = a_anchor[ii_w, 0, :]                        # (K, 2)
        u_vec = a_d[ii_w, 0, :]                          # (K, 2)
        q0 = b_anchor_local[ii_w, jj_w, :]               # (K, 2)
        v_vec = b_d[0, jj_w, :]                          # (K, 2)

        denom = u_vec[:, 0] * v_vec[:, 1] - u_vec[:, 1] * v_vec[:, 0]
        ok = np.abs(denom) > 1e-14
        diff0 = q0 - p0
        denom_safe = np.where(ok, denom, 1.0)
        t_a = (diff0[:, 0] * v_vec[:, 1] - diff0[:, 1] * v_vec[:, 0]) / denom_safe
        t_b = (diff0[:, 0] * u_vec[:, 1] - diff0[:, 1] * u_vec[:, 0]) / denom_safe

        hit = ok & (t_a > tol) & (t_a < 1.0 - tol) & (t_b > tol) & (t_b < 1.0 - tol)
        if not hit.any():
            continue

        sel = np.where(hit)[0]
        out_uv.append(p0[sel] + t_a[sel, None] * u_vec[sel])
        out_a.append((ii_w[sel] + i0).astype(np.int32))
        out_b.append(jj_w[sel].astype(np.int32))
        out_ta.append(t_a[sel])
        out_tb.append(t_b[sel])

    if not out_uv:
        return np.empty(0, dtype=intersect_dtype)

    uv_all = np.concatenate(out_uv)
    a_all = np.concatenate(out_a)
    b_all = np.concatenate(out_b)
    ta_all = np.concatenate(out_ta)
    tb_all = np.concatenate(out_tb)

    res = np.empty(uv_all.shape[0], dtype=intersect_dtype)
    res["uv"] = uv_all
    res["a"] = a_all
    res["b"] = b_all
    res["t_a"] = ta_all
    res["t_b"] = tb_all
    return res


@njit(cache=True)
def _bvh_kernel(e_bbox, f_bbox):
    n_e = e_bbox.shape[0]
    n_f = f_bbox.shape[0]

    max_res = n_e * 15
    out_e = np.empty(max_res, dtype=np.int32)
    out_f = np.empty(max_res, dtype=np.int32)
    count = 0

    e_indices = np.arange(n_e, dtype=np.int32)
    f_indices = np.arange(n_f, dtype=np.int32)

    stack = np.empty((64, 5), dtype=np.int32)
    stack[0] = [0, n_e, 0, n_f, 0]
    stack_ptr = 1

    while stack_ptr > 0:
        stack_ptr -= 1
        e_start, e_end, f_start, f_end, depth = stack[stack_ptr]

        n_curr_e = e_end - e_start
        n_curr_f = f_end - f_start

        if n_curr_e * n_curr_f < 5000 or depth > 15:
            for i in range(e_start, e_end):
                ie = e_indices[i]
                eb = e_bbox[ie]
                for j in range(f_start, f_end):
                    iface = f_indices[j]
                    fb = f_bbox[iface]
                    if (eb[0] <= fb[3] and eb[3] >= fb[0] and
                        eb[1] <= fb[4] and eb[4] >= fb[1] and
                        eb[2] <= fb[5] and eb[5] >= fb[2]):
                        if count < max_res:
                            out_e[count], out_f[count] = ie, iface
                            count += 1
            continue

        axis = depth % 3

        f_coords = f_bbox[f_indices[f_start:f_end], axis]
        mid = (np.min(f_coords) + np.max(f_coords)) * 0.5

        f_split = f_start
        for i in range(f_start, f_end):
            if f_bbox[f_indices[i], axis] <= mid:
                tmp = f_indices[f_split]
                f_indices[f_split] = f_indices[i]
                f_indices[i] = tmp
                f_split += 1

        e_split = e_start
        for i in range(e_start, e_end):
            if e_bbox[e_indices[i], axis] <= mid:
                tmp = e_indices[e_split]
                e_indices[e_split] = e_indices[i]
                e_indices[i] = tmp
                e_split += 1

        if e_split > e_start and f_split > f_start:
            stack[stack_ptr] = [e_start, e_split, f_start, f_split, depth + 1]
            stack_ptr += 1
        if e_split < e_end and f_split < f_end:
            stack[stack_ptr] = [e_split, e_end, f_split, f_end, depth + 1]
            stack_ptr += 1

    return out_e[:count], out_f[:count]


def candidate_pairs(e_bbox: np.ndarray, f_bbox: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Stack-based BVH AABB-overlap broad phase (C8 / G14).

    Parameters
    ----------
    e_bbox, f_bbox : np.ndarray, shape (N, 6)
        Per-item AABBs as [xmin, ymin, zmin, xmax, ymax, zmax].

    Returns
    -------
    (e_idx, f_idx) : tuple of np.ndarray (int64)
        Deduplicated overlapping AABB pairs.
    """
    e_bbox = np.ascontiguousarray(e_bbox, dtype=np.float64)
    f_bbox = np.ascontiguousarray(f_bbox, dtype=np.float64)
    if e_bbox.shape[0] == 0 or f_bbox.shape[0] == 0:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64))

    res_e, res_f = _bvh_kernel(e_bbox, f_bbox)

    if len(res_e) == 0:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64))

    combined = res_e.astype(np.int64) << 32 | res_f.astype(np.int64)
    sort_idx = np.argsort(combined)
    combined = combined[sort_idx]
    mask = np.ones(len(combined), dtype=np.bool_)
    mask[1:] = combined[1:] != combined[:-1]

    return (res_e[sort_idx][mask].astype(np.int64),
            res_f[sort_idx][mask].astype(np.int64))


def moller_trumbore_batch(
    origins: np.ndarray,
    dirs: np.ndarray,
    v0: np.ndarray,
    v1: np.ndarray,
    v2: np.ndarray,
    *,
    eps: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized Möller-Trumbore ray/triangle intersection.

    Each row is one (ray, triangle) test. Ray = origins[i] + t * dirs[i].
    Returns (t, u, v, hit) where hit is the boolean mask of valid intersections
    with t ∈ (0, 1), u, v ≥ 0, u + v ≤ 1.
    """
    edge1 = v1 - v0
    edge2 = v2 - v0
    h = np.cross(dirs, edge2)
    a = np.einsum("ij,ij->i", edge1, h)
    det_ok = np.abs(a) > eps
    a_safe = np.where(det_ok, a, 1.0)
    f = 1.0 / a_safe
    s = origins - v0
    u = f * np.einsum("ij,ij->i", s, h)
    q = np.cross(s, edge1)
    v = f * np.einsum("ij,ij->i", dirs, q)
    t = f * np.einsum("ij,ij->i", edge2, q)
    hit = (
        det_ok
        & (u >= 0.0)
        & (v >= 0.0)
        & (u + v <= 1.0)
        & (t > 0.0)
        & (t < 1.0)
    )
    return t, u, v, hit


def find_double_points(mesh, surface, *, delta: float = 1e-6) -> np.ndarray:
    """Detect 3D self-intersections of `mesh` as edge-vs-face hits (C9).

    Pipeline: 3D AABBs → C8 BVH broad phase → seam-aware vertex-class skip
    (G13) → batched Möller-Trumbore → per-hit Appendix classification with
    sibling-pair consumption (G2). uv1/uv2 are recovered from C5 per-element
    (p, pq, pr) — already in fundamental domain (G15).

    `surface` is accepted for API parity with downstream callers but unused
    in C9 (xyz already populated on `mesh`).
    """
    del surface  # unused; mesh.xyz is precomputed in C6.

    edges = mesh.edges
    faces = mesh.faces
    tris = mesh.tris
    xyz = mesh.xyz
    vertex_class = mesh.vertex_class

    n_e = len(edges)
    n_f = len(faces)
    if n_e == 0 or n_f == 0:
        return np.empty(0, dtype=dp_dtype)

    ep = edges["p_idx"]
    eq = edges["q_idx"]
    ef = edges["f"]
    eg = edges["g"]

    # 1. Per-edge / per-face 3D AABBs.
    e_pts_p = xyz[ep]
    e_pts_q = xyz[eq]
    e_bbox = np.hstack([np.minimum(e_pts_p, e_pts_q), np.maximum(e_pts_p, e_pts_q)])
    f_pts = xyz[tris]                                  # (n_f, 3, 3)
    f_bbox = np.hstack([f_pts.min(axis=1), f_pts.max(axis=1)])

    # 2. C8 BVH broad phase.
    ei, fi = candidate_pairs(e_bbox, f_bbox)
    if len(ei) == 0:
        return np.empty(0, dtype=dp_dtype)

    # 3. Vertex-class skip (G13).
    cls_ep = vertex_class[ep[ei]]                      # (K,)
    cls_eq = vertex_class[eq[ei]]
    cls_fv = vertex_class[tris[fi]]                    # (K, 3)
    shared = (
        (cls_ep[:, None] == cls_fv).any(axis=1)
        | (cls_eq[:, None] == cls_fv).any(axis=1)
    )
    keep = ~shared
    ei = ei[keep]
    fi = fi[keep]
    if len(ei) == 0:
        return np.empty(0, dtype=dp_dtype)

    # 4. Batched Möller-Trumbore.
    origins = xyz[ep[ei]]
    dirs = xyz[eq[ei]] - origins
    v0 = xyz[tris[fi, 0]]
    v1 = xyz[tris[fi, 1]]
    v2 = xyz[tris[fi, 2]]
    t, u, v, hit = moller_trumbore_batch(origins, dirs, v0, v1, v2)
    if not np.any(hit):
        return np.empty(0, dtype=dp_dtype)

    ei = ei[hit]; fi = fi[hit]
    t = t[hit]; u = u[hit]; v = v[hit]
    v0 = v0[hit]; v1 = v1[hit]; v2 = v2[hit]
    bary = np.column_stack([1.0 - u - v, u, v])
    min_bary = bary.min(axis=1)
    xyz_hit = v0 + u[:, None] * (v1 - v0) + v[:, None] * (v2 - v0)

    N = len(ei)
    hit_index = {(int(ei[k]), int(fi[k])): k for k in range(N)}

    def _other_face_of_edge(e_idx: int, f_not: int) -> int:
        f0, f1 = int(ef[e_idx]), int(eg[e_idx])
        if f0 == f_not:
            return f1
        if f1 == f_not:
            return f0
        return -1

    def _vertex_opposite_to_edge_in_face(f_idx: int, e_idx: int) -> int:
        tri = tris[f_idx]
        pe = int(vertex_class[int(ep[e_idx])])
        qe = int(vertex_class[int(eq[e_idx])])
        for j in range(3):
            cls_v = int(vertex_class[int(tri[j])])
            if cls_v != pe and cls_v != qe:
                return j
        return -1

    def _ef_record(e_idx: int, f_idx: int, t_k: float, u_k: float, v_k: float,
                   xyz_k: np.ndarray) -> tuple:
        e_p = edges[e_idx]["p"]
        e_pq = edges[e_idx]["pq"]
        f_p = faces[f_idx]["p"]
        f_pq = faces[f_idx]["pq"]
        f_pr = faces[f_idx]["pr"]
        uv1 = e_p + t_k * e_pq
        uv2 = f_p + u_k * f_pq + v_k * f_pr
        on_bnd = (int(eg[e_idx]) == -1)
        A1 = (int(ef[e_idx]), int(eg[e_idx]))
        A2 = (int(f_idx), -1)
        return (xyz_k, uv1, uv2, int(e_idx), -1, int(f_idx),
                A1, A2, "EF", 2 if on_bnd else 0, on_bnd)

    def _ee_record(e1: int, e2: int, f2: int, t_k: float, u_k: float, v_k: float,
                   xyz_k: np.ndarray) -> tuple:
        e1_p = edges[e1]["p"]
        e1_pq = edges[e1]["pq"]
        f2_p = faces[f2]["p"]
        f2_pq = faces[f2]["pq"]
        f2_pr = faces[f2]["pr"]
        uv1 = e1_p + t_k * e1_pq
        uv2 = f2_p + u_k * f2_pq + v_k * f2_pr
        on_bnd = (int(eg[e1]) == -1) or (int(eg[e2]) == -1)
        A1 = (int(ef[e1]), int(eg[e1]))
        A2 = (int(ef[e2]), int(eg[e2]))
        return (xyz_k, uv1, uv2, int(e1), int(e2), -1,
                A1, A2, "EE", 2 if on_bnd else 0, on_bnd)

    consumed = np.zeros(N, dtype=bool)
    records: list[tuple] = []
    ee_keys: set[tuple[int, int]] = set()

    # 5. Per-hit Appendix classification (G2): sort by descending min_bary,
    # process serially, mark sibling pairs consumed.
    order = np.argsort(-min_bary)
    for k in order:
        if consumed[k]:
            continue
        consumed[k] = True
        e1, f2 = int(ei[k]), int(fi[k])
        if min_bary[k] >= delta:
            records.append(_ef_record(e1, f2, t[k], u[k], v[k], xyz_hit[k]))
            continue

        # Hit lies on edge of F2 opposite the smallest-bary vertex.
        opp = int(np.argmin(bary[k]))
        e2 = int(faces[f2]["edges"][(opp + 1) % 3])
        if e2 == e1:
            continue

        f1, f1p = int(ef[e1]), int(eg[e1])
        f2p = _other_face_of_edge(e2, f2)

        sib_e1_f2p = hit_index.get((e1, f2p)) if f2p >= 0 else None
        sib_e2_f1 = hit_index.get((e2, f1)) if f1 >= 0 else None
        sib_e2_f1p = hit_index.get((e2, f1p)) if f1p >= 0 else None

        if sib_e2_f1 is not None:
            q_idx, q_f = sib_e2_f1, f1
        elif sib_e2_f1p is not None:
            q_idx, q_f = sib_e2_f1p, f1p
        else:
            q_idx, q_f = None, -1

        if q_idx is None:
            key = (min(e1, e2), max(e1, e2))
            if key not in ee_keys:
                ee_keys.add(key)
                records.append(_ee_record(e1, e2, f2, t[k], u[k], v[k], xyz_hit[k]))
        else:
            j_opp_e1 = _vertex_opposite_to_edge_in_face(q_f, e1)
            if j_opp_e1 >= 0 and bary[q_idx][j_opp_e1] >= delta:
                records.append(_ef_record(
                    e2, q_f, t[q_idx], u[q_idx], v[q_idx], xyz_hit[q_idx]
                ))
            else:
                key = (min(e1, e2), max(e1, e2))
                if key not in ee_keys:
                    ee_keys.add(key)
                    records.append(_ee_record(e1, e2, f2, t[k], u[k], v[k], xyz_hit[k]))

        for s in (sib_e1_f2p, sib_e2_f1, sib_e2_f1p):
            if s is not None:
                consumed[s] = True

    if not records:
        return np.empty(0, dtype=dp_dtype)

    out = np.empty(len(records), dtype=dp_dtype)
    for i, rec in enumerate(records):
        (out["xyz"][i], out["uv1"][i], out["uv2"][i], out["E1"][i],
         out["E2"][i], out["F2"][i], out["A1"][i], out["A2"][i],
         out["type"][i], out["ptype"][i], out["on_boundary"][i]) = rec
    return out


def build_sis_pairs(dps: np.ndarray) -> np.ndarray:
    """Pair DPs into SIS records by A1/A2 face-sharing (C10).

    Two DPs `i < j` are an SIS edge iff their `A1`/`A2` face-pair fields share
    faces in one of two orientations:
      - unflipped: `A1[i] ∩ A1[j] ≠ ∅` and `A2[i] ∩ A2[j] ≠ ∅` → flip=+1
      - flipped:   `A1[i] ∩ A2[j] ≠ ∅` and `A2[i] ∩ A1[j] ≠ ∅` → flip=-1
    If both match (pathological), warn and default to flip=+1.
    `split1`, `split2` are initialized to -1 (G17 — populated in Layer O).
    """
    n = len(dps)
    if n < 2:
        return np.empty(0, dtype=sis_dtype)

    i_arr, j_arr = np.triu_indices(n, k=1)
    A1i = dps["A1"][i_arr]
    A2i = dps["A2"][i_arr]
    A1j = dps["A1"][j_arr]
    A2j = dps["A2"][j_arr]

    def _shares(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return ((a[:, :, None] == b[:, None, :]) & (a[:, :, None] >= 0)).any(axis=(1, 2))

    unflipped = _shares(A1i, A1j) & _shares(A2i, A2j)
    flipped = _shares(A1i, A2j) & _shares(A2i, A1j)

    ambiguous = unflipped & flipped
    if ambiguous.any():
        amb_idx = np.flatnonzero(ambiguous)
        for k in amb_idx:
            warnings.warn(
                f"ambiguous A1/A2 match for SIS pair "
                f"({int(i_arr[k])}, {int(j_arr[k])}); defaulting to flip=+1",
                stacklevel=2,
            )

    keep = unflipped | flipped
    n_out = int(keep.sum())
    out = np.empty(n_out, dtype=sis_dtype)
    if n_out == 0:
        return out
    out["p_dp"] = i_arr[keep].astype(np.int32)
    out["q_dp"] = j_arr[keep].astype(np.int32)
    out["flip"] = np.where(unflipped[keep], 1, -1).astype(np.int8)
    out["split1"] = -1
    out["split2"] = -1
    return out
