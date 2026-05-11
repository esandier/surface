"""intersections.py — P5: sweep_segments; C8: candidate_pairs (BVH broad phase).

Shared 2D segment-segment intersection primitive used by the domain sweep
(self-intersection preimages, CS×SIS, etc.) and the view-plane sweep.
See Modular_rewrite_roadmap.md §P5; G3 (close-aware) is enforced when
`domain` is a rectangular domain with identification.

C8: `candidate_pairs(e_bbox, f_bbox)` returns overlapping AABB pairs via a
stack-based BVH with axis-cycling partitioning, Numba-JIT compiled. See G14.
"""

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
