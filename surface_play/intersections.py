"""intersections.py — P5: sweep_segments; C8: candidate_pairs (BVH broad phase);
C9: find_double_points; C10: build_sis_pairs; C11: build_sics; C12:
find_triple_points.

Shared 2D segment-segment intersection primitive used by the domain sweep
(self-intersection preimages, CS×SIS, etc.) and the view-plane sweep.
See Modular_rewrite_roadmap.md §P5; G3 (close-aware) is enforced when
`domain` is a rectangular domain with identification.

C8: `candidate_pairs(e_bbox, f_bbox)` returns overlapping AABB pairs via a
stack-based spatial-split BVH with axis-cycling partitioning, Numba-JIT
compiled. Straddling AABBs recurse into both children (unlike object-split).
See G14.

C9: `find_double_points(mesh, surface)` returns DP records (edge-vs-face hits
in 3D), with sibling-pair consumption (G2) and seam-aware vertex-class skip
(G13). Uses C5 per-element (p, pq, pr) for uv recovery (G15).

C10: `build_sis_pairs(dps)` returns SIS records (DP-DP edges of the
self-intersection graph) with vectorized A1/A2 face-sharing test and flip
detection. See G17 for split sentinel.

C11: `build_sics(sis_pairs)` chains SIS records into SelfIntersectingCurve
objects by delegating to P4 `make_lines` on the DP-index graph.

C12: `find_triple_points(sis_pairs, dps, mesh, surface)` reconstructs the
two preimage segments of every SIS and sweeps them in self-mode (G3
close-aware). Three preimages meeting at one 3D point produce three domain
hits whose `f_other` values pair up as (F2,F3), (F1,F3), (F1,F2); 3D
verification (`‖S(P_i)-S(P_j)‖ < xyz_tol`) discards spurious interlocks.
"""

import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np
from numba import njit

from surface_play.curves import make_lines


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


tp_dtype = np.dtype([
    ("xyz",         "f8", 3),
    ("sis_indices", "i4", 3),
    ("uv",          "f8", (3, 2)),
    ("faces",       "i4", 3),
])


_tp_seg_dtype = np.dtype([
    ("uv0",     "f8", 2),
    ("uv1",     "f8", 2),
    ("f_here",  "i4"),
    ("f_other", "i4"),
    ("sis_idx", "i4"),
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
    """Spatial-split BVH: AABBs straddling a partition appear in both children.

    Unlike object-split, an item with `min <= mid` AND `max > mid` is recursed
    into both halves. This avoids the legacy bug where wide AABBs were placed
    in one child only and missed overlap pairs in the other (see
    `bvh_object_split_limitation.md`).

    Guards: bail to brute-force in-node if (1) the split would inflate total
    work (no-progress: parent items mostly straddle), or (2) the workspace /
    output buffer runs out. Output is raw pairs (with duplicates from
    spatial overlap of leaves); `candidate_pairs` dedupes downstream.
    """
    n_e = e_bbox.shape[0]
    n_f = f_bbox.shape[0]

    max_depth = 15
    leaf_threshold = 5000

    # Workspace: indices into original arrays, growing append-only as items
    # are duplicated into children. 32× the input is generous for typical
    # geometries; overflow falls back to brute force in-node (still correct).
    e_capacity = n_e * 32 + 64
    f_capacity = n_f * 32 + 64
    e_workspace = np.empty(e_capacity, dtype=np.int32)
    f_workspace = np.empty(f_capacity, dtype=np.int32)

    for i in range(n_e):
        e_workspace[i] = i
    for i in range(n_f):
        f_workspace[i] = i
    e_next = n_e
    f_next = n_f

    stack_capacity = 4096
    stack = np.empty((stack_capacity, 5), dtype=np.int32)
    stack[0, 0] = 0
    stack[0, 1] = n_e
    stack[0, 2] = 0
    stack[0, 3] = n_f
    stack[0, 4] = 0
    stack_ptr = 1

    # Generous output buffer: spatial split can repeat true pairs across
    # straddling leaves; dedup happens in `candidate_pairs`.
    max_res = (n_e + n_f) * 50 + 4096
    out_e = np.empty(max_res, dtype=np.int32)
    out_f = np.empty(max_res, dtype=np.int32)
    count = 0

    while stack_ptr > 0:
        stack_ptr -= 1
        e_start = stack[stack_ptr, 0]
        e_end = stack[stack_ptr, 1]
        f_start = stack[stack_ptr, 2]
        f_end = stack[stack_ptr, 3]
        depth = stack[stack_ptr, 4]
        n_curr_e = e_end - e_start
        n_curr_f = f_end - f_start

        # Leaf or fallback: brute-force this node.
        if n_curr_e * n_curr_f < leaf_threshold or depth >= max_depth:
            for i in range(e_start, e_end):
                ie = e_workspace[i]
                eb0 = e_bbox[ie, 0]; eb1 = e_bbox[ie, 1]; eb2 = e_bbox[ie, 2]
                eb3 = e_bbox[ie, 3]; eb4 = e_bbox[ie, 4]; eb5 = e_bbox[ie, 5]
                for j in range(f_start, f_end):
                    iface = f_workspace[j]
                    if (eb0 <= f_bbox[iface, 3] and eb3 >= f_bbox[iface, 0] and
                        eb1 <= f_bbox[iface, 4] and eb4 >= f_bbox[iface, 1] and
                        eb2 <= f_bbox[iface, 5] and eb5 >= f_bbox[iface, 2]):
                        if count < max_res:
                            out_e[count] = ie
                            out_f[count] = iface
                            count += 1
            continue

        axis = depth % 3
        axis_max = axis + 3

        # Midpoint over f's full extent on this axis.
        f_min_a = f_bbox[f_workspace[f_start], axis]
        f_max_a = f_bbox[f_workspace[f_start], axis_max]
        for j in range(f_start + 1, f_end):
            idx = f_workspace[j]
            v_min = f_bbox[idx, axis]
            v_max = f_bbox[idx, axis_max]
            if v_min < f_min_a: f_min_a = v_min
            if v_max > f_max_a: f_max_a = v_max
        mid = (f_min_a + f_max_a) * 0.5

        # Need worst-case 2× current size on each side (every item straddles).
        # If workspace can't accommodate, fall back to brute-force this node.
        if (e_next + 2 * n_curr_e > e_capacity or
            f_next + 2 * n_curr_f > f_capacity or
            stack_ptr + 2 > stack_capacity):
            for i in range(e_start, e_end):
                ie = e_workspace[i]
                eb0 = e_bbox[ie, 0]; eb1 = e_bbox[ie, 1]; eb2 = e_bbox[ie, 2]
                eb3 = e_bbox[ie, 3]; eb4 = e_bbox[ie, 4]; eb5 = e_bbox[ie, 5]
                for j in range(f_start, f_end):
                    iface = f_workspace[j]
                    if (eb0 <= f_bbox[iface, 3] and eb3 >= f_bbox[iface, 0] and
                        eb1 <= f_bbox[iface, 4] and eb4 >= f_bbox[iface, 1] and
                        eb2 <= f_bbox[iface, 5] and eb5 >= f_bbox[iface, 2]):
                        if count < max_res:
                            out_e[count] = ie
                            out_f[count] = iface
                            count += 1
            continue

        e_left_start = e_next
        for i in range(e_start, e_end):
            idx = e_workspace[i]
            if e_bbox[idx, axis] <= mid:
                e_workspace[e_next] = idx
                e_next += 1
        e_left_end = e_next
        e_right_start = e_next
        for i in range(e_start, e_end):
            idx = e_workspace[i]
            if e_bbox[idx, axis_max] > mid:
                e_workspace[e_next] = idx
                e_next += 1
        e_right_end = e_next

        f_left_start = f_next
        for i in range(f_start, f_end):
            idx = f_workspace[i]
            if f_bbox[idx, axis] <= mid:
                f_workspace[f_next] = idx
                f_next += 1
        f_left_end = f_next
        f_right_start = f_next
        for i in range(f_start, f_end):
            idx = f_workspace[i]
            if f_bbox[idx, axis_max] > mid:
                f_workspace[f_next] = idx
                f_next += 1
        f_right_end = f_next

        e_left_n = e_left_end - e_left_start
        e_right_n = e_right_end - e_right_start
        f_left_n = f_left_end - f_left_start
        f_right_n = f_right_end - f_right_start

        # No-progress guard: if total child work doesn't reduce parent work
        # by ≥25%, bail. Prevents pathological recursion when items straddle.
        parent_work = n_curr_e * n_curr_f
        child_work = e_left_n * f_left_n + e_right_n * f_right_n
        if child_work * 4 >= parent_work * 3:
            # Rewind appended workspace; brute-force this node.
            e_next = e_left_start
            f_next = f_left_start
            for i in range(e_start, e_end):
                ie = e_workspace[i]
                eb0 = e_bbox[ie, 0]; eb1 = e_bbox[ie, 1]; eb2 = e_bbox[ie, 2]
                eb3 = e_bbox[ie, 3]; eb4 = e_bbox[ie, 4]; eb5 = e_bbox[ie, 5]
                for j in range(f_start, f_end):
                    iface = f_workspace[j]
                    if (eb0 <= f_bbox[iface, 3] and eb3 >= f_bbox[iface, 0] and
                        eb1 <= f_bbox[iface, 4] and eb4 >= f_bbox[iface, 1] and
                        eb2 <= f_bbox[iface, 5] and eb5 >= f_bbox[iface, 2]):
                        if count < max_res:
                            out_e[count] = ie
                            out_f[count] = iface
                            count += 1
            continue

        if e_left_n > 0 and f_left_n > 0:
            stack[stack_ptr, 0] = e_left_start
            stack[stack_ptr, 1] = e_left_end
            stack[stack_ptr, 2] = f_left_start
            stack[stack_ptr, 3] = f_left_end
            stack[stack_ptr, 4] = depth + 1
            stack_ptr += 1
        if e_right_n > 0 and f_right_n > 0:
            stack[stack_ptr, 0] = e_right_start
            stack[stack_ptr, 1] = e_right_end
            stack[stack_ptr, 2] = f_right_start
            stack[stack_ptr, 3] = f_right_end
            stack[stack_ptr, 4] = depth + 1
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

    # 3. Skip pairs sharing a vertex. Vertex labels are already canonical (compacted
    # mesh), so plain label equality is sufficient — no vertex_class lookup needed.
    cls_ep = ep[ei]                                    # (K,)
    cls_eq = eq[ei]
    cls_fv = tris[fi]                                  # (K, 3)
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
        pe = int(ep[e_idx])
        qe = int(eq[e_idx])
        for j in range(3):
            v = int(tri[j])
            if v != pe and v != qe:
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


@dataclass
class SelfIntersectingCurve:
    sis_indices: np.ndarray   # 1D int array; signs encode reversal as in make_lines output
    is_closed: bool


def _first_shared_face(a: np.ndarray, b: np.ndarray) -> int:
    """Lowest-index face shared by length-2 A-sets `a` and `b` (-1 if none).

    `-1` entries are sentinel padding (see C9 dp_dtype) and ignored.
    """
    a_set = {int(x) for x in a if int(x) >= 0}
    b_set = {int(y) for y in b if int(y) >= 0}
    common = a_set & b_set
    return min(common) if common else -1


def find_triple_points(
    sis_pairs: np.ndarray,
    dps: np.ndarray,
    mesh,
    surface,
    *,
    xyz_tol: float = 1e-3,
) -> np.ndarray:
    """Detect triple points: three SIS preimages meeting at one 3D point (C12).

    Pipeline: build 2 preimage segments per SIS (using `flip`, `uv1/uv2`,
    `A1/A2`) → G3 close-aware self-sweep → keep hits where both preimages
    share a face (`f_here[a] == f_here[b]`) → group by `(f_other, f_other)`
    pair-key, look for triples (F1, F2, F3) where all three sorted pair-keys
    (F2,F3), (F1,F3), (F1,F2) are present → 3D verify via `surface.S`.

    Returns a structured array of `tp_dtype` (empty if no TPs).
    """
    K = len(sis_pairs)
    if K < 3:
        return np.empty(0, dtype=tp_dtype)

    # 1. Build per-SIS preimage segments (two per SIS).
    segs = np.empty(2 * K, dtype=_tp_seg_dtype)
    n_seg = 0
    for k in range(K):
        p = int(sis_pairs["p_dp"][k])
        q = int(sis_pairs["q_dp"][k])
        f = int(sis_pairs["flip"][k])

        A_p_uv = dps["uv1"][p]
        A_q_uv = dps["uv1"][q] if f == 1 else dps["uv2"][q]
        B_p_uv = dps["uv2"][p]
        B_q_uv = dps["uv2"][q] if f == 1 else dps["uv1"][q]

        A_p_set = dps["A1"][p]
        A_q_set = dps["A1"][q] if f == 1 else dps["A2"][q]
        B_p_set = dps["A2"][p]
        B_q_set = dps["A2"][q] if f == 1 else dps["A1"][q]

        f_here_A = _first_shared_face(A_p_set, A_q_set)
        f_here_B = _first_shared_face(B_p_set, B_q_set)
        if f_here_A < 0 or f_here_B < 0:
            continue

        segs[n_seg] = (A_p_uv, A_q_uv, f_here_A, f_here_B, k)
        segs[n_seg + 1] = (B_p_uv, B_q_uv, f_here_B, f_here_A, k)
        n_seg += 2

    segs = segs[:n_seg]
    if n_seg < 3:
        return np.empty(0, dtype=tp_dtype)

    # 2. Close-aware self-sweep on all preimage segments.
    domain = getattr(mesh, "domain", None)
    hits = sweep_segments(
        segs["uv0"], segs["uv1"], None, None, domain, self_sweep=True
    )
    if len(hits) == 0:
        return np.empty(0, dtype=tp_dtype)

    # 3. Keep hits where both preimages share a face.
    same_face = segs["f_here"][hits["a"]] == segs["f_here"][hits["b"]]
    hits = hits[same_face]
    if len(hits) == 0:
        return np.empty(0, dtype=tp_dtype)

    fo_a = segs["f_other"][hits["a"]]
    fo_b = segs["f_other"][hits["b"]]
    pair_lo = np.minimum(fo_a, fo_b).astype(np.int64)
    pair_hi = np.maximum(fo_a, fo_b).astype(np.int64)

    pair_map: dict[tuple[int, int], list[int]] = {}
    face_partners: dict[int, set[int]] = {}
    for i in range(len(hits)):
        key = (int(pair_lo[i]), int(pair_hi[i]))
        pair_map.setdefault(key, []).append(i)
        face_partners.setdefault(key[0], set()).add(key[1])
        face_partners.setdefault(key[1], set()).add(key[0])

    # 4. Iterate ordered triples (F1, F2, F3) and verify in 3D.
    out: list[np.ndarray] = []
    seen: set[tuple[int, int, int]] = set()
    for F1 in sorted(face_partners.keys()):
        plist = sorted(face_partners[F1])
        for ia in range(len(plist)):
            F2 = plist[ia]
            if F2 <= F1:
                continue
            for ib in range(ia + 1, len(plist)):
                F3 = plist[ib]
                if F3 <= F2:
                    continue
                if (F2, F3) not in pair_map:
                    continue
                triple = (F1, F2, F3)
                if triple in seen:
                    continue
                seen.add(triple)

                cands_F1 = pair_map[(F2, F3)]   # hit on F1
                cands_F2 = pair_map[(F1, F3)]   # hit on F2
                cands_F3 = pair_map[(F1, F2)]   # hit on F3

                rec = _try_emit_tp(
                    cands_F1, cands_F2, cands_F3,
                    F1, F2, F3, hits, segs, surface, xyz_tol,
                )
                if rec is not None:
                    out.append(rec)

    if not out:
        return np.empty(0, dtype=tp_dtype)

    result = np.empty(len(out), dtype=tp_dtype)
    for i, r in enumerate(out):
        result[i] = r
    # Stable order: by sorted sis_indices tuple, then by faces.
    order = np.lexsort((
        result["faces"][:, 2], result["faces"][:, 1], result["faces"][:, 0],
        result["sis_indices"][:, 2],
        result["sis_indices"][:, 1],
        result["sis_indices"][:, 0],
    ))
    return result[order]


def _try_emit_tp(
    cands_F1, cands_F2, cands_F3,
    F1, F2, F3, hits, segs, surface, xyz_tol,
):
    """Iterate hit combinations for face triple (F1<F2<F3); emit first that
    passes 3D verification. Returns a 0-d tp_dtype record or None.
    """
    for h1 in cands_F1:
        for h2 in cands_F2:
            for h3 in cands_F3:
                P1 = hits["uv"][h1]
                P2 = hits["uv"][h2]
                P3 = hits["uv"][h3]
                S1 = np.asarray(surface.S(float(P1[0]), float(P1[1])),
                                dtype=float).reshape(3)
                S2 = np.asarray(surface.S(float(P2[0]), float(P2[1])),
                                dtype=float).reshape(3)
                S3 = np.asarray(surface.S(float(P3[0]), float(P3[1])),
                                dtype=float).reshape(3)
                d12 = float(np.linalg.norm(S1 - S2))
                d13 = float(np.linalg.norm(S1 - S3))
                d23 = float(np.linalg.norm(S2 - S3))
                if max(d12, d13, d23) >= xyz_tol:
                    continue

                sis_set = {
                    int(segs["sis_idx"][hits["a"][h1]]),
                    int(segs["sis_idx"][hits["b"][h1]]),
                    int(segs["sis_idx"][hits["a"][h2]]),
                    int(segs["sis_idx"][hits["b"][h2]]),
                    int(segs["sis_idx"][hits["a"][h3]]),
                    int(segs["sis_idx"][hits["b"][h3]]),
                }
                if len(sis_set) != 3:
                    continue

                rec = np.zeros((), dtype=tp_dtype)
                rec["xyz"] = (S1 + S2 + S3) / 3.0
                rec["sis_indices"] = sorted(sis_set)
                rec["uv"][0] = P1
                rec["uv"][1] = P2
                rec["uv"][2] = P3
                rec["faces"] = (F1, F2, F3)
                return rec
    return None


def build_sics(sis_pairs: np.ndarray) -> list[SelfIntersectingCurve]:
    """Chain SIS records into SICs via P4 `make_lines` on the DP-index graph (C11).

    Nodes = DP indices, edges = SIS records. Each chain returned by `make_lines`
    is a 1D array of signed 1-indexed SIS ids: `+k` means SIS k-1 traversed
    `p_dp → q_dp`, `-k` means the reverse. Closed chains duplicate the first
    token at the end (`chain[0] == chain[-1]`).
    """
    if len(sis_pairs) == 0:
        return []

    segments = np.column_stack([sis_pairs["p_dp"], sis_pairs["q_dp"]]).astype(np.intp)
    chains = make_lines(segments)

    result = []
    for chain in chains:
        is_closed = bool(chain[0] == chain[-1])
        result.append(SelfIntersectingCurve(
            sis_indices=np.asarray(chain, dtype=np.intp),
            is_closed=is_closed,
        ))
    return result
