"""intersections.py — P5: sweep_segments

Shared 2D segment-segment intersection primitive used by the domain sweep
(self-intersection preimages, CS×SIS, etc.) and the view-plane sweep.
See Modular_rewrite_roadmap.md §P5; G3 (close-aware) is enforced when
`domain` is a rectangular domain with identification.
"""

from typing import Optional

import numpy as np


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
