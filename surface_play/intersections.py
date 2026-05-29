"""intersections.py — P5: sweep_segments; C8/C9: edge-vs-face broad phase +
find_double_points; C10: build_sis_pairs; C11: build_sics; C12:
find_triple_points.

Shared 2D segment-segment intersection primitive used by the domain sweep
(self-intersection preimages, CS×SIS, etc.) and the view-plane sweep.
See Modular_rewrite_roadmap.md §P5; G3 (close-aware) is enforced when
`domain` is a rectangular domain with identification.

C8 (broad phase): `edge_face_candidates(e_bbox, f_bbox, ep, eq, tris)` returns
non-adjacent overlapping edge/face AABB pairs via a uniform spatial grid
(cell side = max element extent ⇒ bounded insertion), with the shared-vertex
skip (G13) fused into the Numba kernel so the ~1M 1-ring adjacency overlaps of
a connected mesh never materialise. (Superseded the spatial-split BVH broad
phase — see [[construction-perf-2026-05-29]].)

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


@njit(cache=True)
def _sweep_segments_kernel(
    a_uv0,           # (n, 2)
    a_d,             # (n, 2) — pre-computed close-aware direction
    b_uv0,           # (m, 2)
    b_d,             # (m, 2)
    use_close,       # bool
    u_period_mode,   # int: 0 = no wrap, 1 = wrap with period_u
    v_period_mode,   # int: 0 = no wrap, 1 = wrap with period_v
    period_u,        # float (ignored if u_period_mode == 0)
    period_v,        # float
    self_sweep,      # bool
    tol,             # float
):
    """Numba kernel: O(n*m) brute-force pairwise segment intersection.

    Mirrors `sweep_segments` semantics exactly:
      - AABB overlap filter (close-aware diff if `use_close`),
      - 2x2 linear solve for (t_a, t_b), strict interior acceptance,
      - returns arrays of intersections in unspecified order.
    """
    n = a_uv0.shape[0]
    m = b_uv0.shape[0]

    cap = 64
    out_uv = np.empty((cap, 2), dtype=np.float64)
    out_a = np.empty(cap, dtype=np.int32)
    out_b = np.empty(cap, dtype=np.int32)
    out_ta = np.empty(cap, dtype=np.float64)
    out_tb = np.empty(cap, dtype=np.float64)
    k = 0

    for i in range(n):
        a0u = a_uv0[i, 0]
        a0v = a_uv0[i, 1]
        adu = a_d[i, 0]
        adv = a_d[i, 1]
        a1u = a0u + adu
        a1v = a0v + adv
        if adu >= 0.0:
            a_lo_u = a0u; a_hi_u = a1u
        else:
            a_lo_u = a1u; a_hi_u = a0u
        if adv >= 0.0:
            a_lo_v = a0v; a_hi_v = a1v
        else:
            a_lo_v = a1v; a_hi_v = a0v

        j_start = i + 1 if self_sweep else 0
        for j in range(j_start, m):
            b0u = b_uv0[j, 0]
            b0v = b_uv0[j, 1]
            bdu = b_d[j, 0]
            bdv = b_d[j, 1]

            diff_u = b0u - a0u
            diff_v = b0v - a0v
            if use_close:
                if u_period_mode == 1:
                    diff_u = (diff_u + 0.5 * period_u) % period_u - 0.5 * period_u
                if v_period_mode == 1:
                    diff_v = (diff_v + 0.5 * period_v) % period_v - 0.5 * period_v
            b_anc_u = a0u + diff_u
            b_anc_v = a0v + diff_v
            b_end_u = b_anc_u + bdu
            b_end_v = b_anc_v + bdv
            if bdu >= 0.0:
                b_lo_u = b_anc_u; b_hi_u = b_end_u
            else:
                b_lo_u = b_end_u; b_hi_u = b_anc_u
            if bdv >= 0.0:
                b_lo_v = b_anc_v; b_hi_v = b_end_v
            else:
                b_lo_v = b_end_v; b_hi_v = b_anc_v

            # AABB overlap (close-aware in B's local frame).
            if a_lo_u > b_hi_u + tol:
                continue
            if b_lo_u > a_hi_u + tol:
                continue
            if a_lo_v > b_hi_v + tol:
                continue
            if b_lo_v > a_hi_v + tol:
                continue

            denom = adu * bdv - adv * bdu
            if -1e-14 < denom < 1e-14:
                continue

            d0_u = b_anc_u - a0u
            d0_v = b_anc_v - a0v
            t_a = (d0_u * bdv - d0_v * bdu) / denom
            t_b = (d0_u * adv - d0_v * adu) / denom

            if t_a <= tol or t_a >= 1.0 - tol:
                continue
            if t_b <= tol or t_b >= 1.0 - tol:
                continue

            if k >= cap:
                new_cap = cap * 2
                new_uv = np.empty((new_cap, 2), dtype=np.float64)
                new_uv[:cap] = out_uv
                out_uv = new_uv
                new_a = np.empty(new_cap, dtype=np.int32)
                new_a[:cap] = out_a
                out_a = new_a
                new_b = np.empty(new_cap, dtype=np.int32)
                new_b[:cap] = out_b
                out_b = new_b
                new_ta = np.empty(new_cap, dtype=np.float64)
                new_ta[:cap] = out_ta
                out_ta = new_ta
                new_tb = np.empty(new_cap, dtype=np.float64)
                new_tb[:cap] = out_tb
                out_tb = new_tb
                cap = new_cap

            out_uv[k, 0] = a0u + t_a * adu
            out_uv[k, 1] = a0v + t_a * adv
            out_a[k] = i
            out_b[k] = j
            out_ta[k] = t_a
            out_tb[k] = t_b
            k += 1

    return out_uv[:k], out_a[:k], out_b[:k], out_ta[:k], out_tb[:k]


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
    a_d = np.ascontiguousarray(seg_a_uv1_loc - seg_a_uv0)
    b_d = np.ascontiguousarray(seg_b_uv1_loc - seg_b_uv0)
    a_uv0_c = np.ascontiguousarray(seg_a_uv0, dtype=np.float64)
    b_uv0_c = np.ascontiguousarray(seg_b_uv0, dtype=np.float64)

    use_close = bool(domain is not None and getattr(domain, "type", None) == "rect")
    u_id = getattr(domain, "u_identify", "no") if domain is not None else "no"
    v_id = getattr(domain, "v_identify", "no") if domain is not None else "no"
    u_mode = 1 if (use_close and u_id in ("cy", "mo")) else 0
    v_mode = 1 if (use_close and v_id in ("cy", "mo")) else 0
    p_u = float(getattr(domain, "period_u", 0.0)) if u_mode == 1 else 0.0
    p_v = float(getattr(domain, "period_v", 0.0)) if v_mode == 1 else 0.0

    uv_all, a_all, b_all, ta_all, tb_all = _sweep_segments_kernel(
        a_uv0_c, a_d, b_uv0_c, b_d,
        use_close, u_mode, v_mode, p_u, p_v,
        bool(self_sweep), float(tol),
    )

    if uv_all.shape[0] == 0:
        return np.empty(0, dtype=intersect_dtype)

    res = np.empty(uv_all.shape[0], dtype=intersect_dtype)
    res["uv"] = uv_all
    res["a"] = a_all
    res["b"] = b_all
    res["t_a"] = ta_all
    res["t_b"] = tb_all
    return res


@njit(cache=True)
def _grid_ef_kernel(e_bbox, f_bbox, ep, eq, tris,
                    origin, inv_h, nx, ny, nz, cap):
    """Uniform-grid edge-vs-face AABB broad phase with FUSED shared-vertex skip.

    Bins faces into a cell grid (cell side = max element AABB extent, so every
    element spans ≤2 cells per axis → bounded insertion, no missed overlaps),
    then for each edge tests only the faces in the cells its AABB touches.

    A pair sharing a mesh vertex can never be a transverse self-intersection
    (it's incident geometry), so it is skipped inline — never emitted. On a
    non-self-intersecting surface this drops the emitted-pair count from ~1M
    (all 1-ring adjacency) to ~0. Equivalent to the post-broad-phase
    shared-vertex filter in `find_double_points`, but fused so the adjacency
    pairs cost nothing downstream.

    Returns `(out_e, out_f, total)`. Arrays are filled up to `min(total, cap)`;
    if `total > cap` the caller retries with a larger buffer (correctness).
    """
    n_e = e_bbox.shape[0]
    n_f = f_bbox.shape[0]
    ncells = nx * ny * nz

    # CSR bin of faces: count per cell, prefix-sum, then fill.
    f_count = np.zeros(ncells + 1, dtype=np.int64)
    for jf in range(n_f):
        cx0 = int((f_bbox[jf, 0] - origin[0]) * inv_h[0])
        cy0 = int((f_bbox[jf, 1] - origin[1]) * inv_h[1])
        cz0 = int((f_bbox[jf, 2] - origin[2]) * inv_h[2])
        cx1 = int((f_bbox[jf, 3] - origin[0]) * inv_h[0])
        cy1 = int((f_bbox[jf, 4] - origin[1]) * inv_h[1])
        cz1 = int((f_bbox[jf, 5] - origin[2]) * inv_h[2])
        for cx in range(cx0, cx1 + 1):
            for cy in range(cy0, cy1 + 1):
                for cz in range(cz0, cz1 + 1):
                    f_count[(cx * ny + cy) * nz + cz + 1] += 1
    for c in range(ncells):
        f_count[c + 1] += f_count[c]
    f_items = np.empty(f_count[ncells], dtype=np.int32)
    f_fill = f_count[:ncells].copy()
    for jf in range(n_f):
        cx0 = int((f_bbox[jf, 0] - origin[0]) * inv_h[0])
        cy0 = int((f_bbox[jf, 1] - origin[1]) * inv_h[1])
        cz0 = int((f_bbox[jf, 2] - origin[2]) * inv_h[2])
        cx1 = int((f_bbox[jf, 3] - origin[0]) * inv_h[0])
        cy1 = int((f_bbox[jf, 4] - origin[1]) * inv_h[1])
        cz1 = int((f_bbox[jf, 5] - origin[2]) * inv_h[2])
        for cx in range(cx0, cx1 + 1):
            for cy in range(cy0, cy1 + 1):
                for cz in range(cz0, cz1 + 1):
                    c = (cx * ny + cy) * nz + cz
                    f_items[f_fill[c]] = jf
                    f_fill[c] += 1

    out_e = np.empty(cap, dtype=np.int32)
    out_f = np.empty(cap, dtype=np.int32)
    total = 0
    for ie in range(n_e):
        eb0 = e_bbox[ie, 0]; eb1 = e_bbox[ie, 1]; eb2 = e_bbox[ie, 2]
        eb3 = e_bbox[ie, 3]; eb4 = e_bbox[ie, 4]; eb5 = e_bbox[ie, 5]
        pa = ep[ie]; qa = eq[ie]
        cx0 = int((eb0 - origin[0]) * inv_h[0])
        cy0 = int((eb1 - origin[1]) * inv_h[1])
        cz0 = int((eb2 - origin[2]) * inv_h[2])
        cx1 = int((eb3 - origin[0]) * inv_h[0])
        cy1 = int((eb4 - origin[1]) * inv_h[1])
        cz1 = int((eb5 - origin[2]) * inv_h[2])
        for cx in range(cx0, cx1 + 1):
            for cy in range(cy0, cy1 + 1):
                for cz in range(cz0, cz1 + 1):
                    c = (cx * ny + cy) * nz + cz
                    for p in range(f_count[c], f_count[c + 1]):
                        jf = f_items[p]
                        v0 = tris[jf, 0]; v1 = tris[jf, 1]; v2 = tris[jf, 2]
                        if (pa == v0 or pa == v1 or pa == v2 or
                                qa == v0 or qa == v1 or qa == v2):
                            continue  # shared vertex → incident, not a DP
                        if (eb0 <= f_bbox[jf, 3] and eb3 >= f_bbox[jf, 0] and
                                eb1 <= f_bbox[jf, 4] and eb4 >= f_bbox[jf, 1] and
                                eb2 <= f_bbox[jf, 5] and eb5 >= f_bbox[jf, 2]):
                            if total < cap:
                                out_e[total] = ie
                                out_f[total] = jf
                            total += 1
    return out_e, out_f, total


def edge_face_candidates(
    e_bbox: np.ndarray, f_bbox: np.ndarray,
    ep: np.ndarray, eq: np.ndarray, tris: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Edge-vs-face broad phase for self-intersection (C9 steps 2+3 fused).

    Uniform-grid AABB overlap with shared-vertex pairs skipped inline. Returns
    the deduplicated, combined-key-sorted `(e_idx, f_idx)` of non-adjacent
    overlapping pairs — identical to a plain AABB-overlap broad phase followed
    by dropping pairs that share a mesh vertex, but far cheaper because the ~1M
    adjacency overlaps of a connected mesh are never materialised.
    """
    e_bbox = np.ascontiguousarray(e_bbox, dtype=np.float64)
    f_bbox = np.ascontiguousarray(f_bbox, dtype=np.float64)
    n_e, n_f = e_bbox.shape[0], f_bbox.shape[0]
    if n_e == 0 or n_f == 0:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64))

    ep = np.ascontiguousarray(ep, dtype=np.int32)
    eq = np.ascontiguousarray(eq, dtype=np.int32)
    tris = np.ascontiguousarray(tris, dtype=np.int32)

    allb = np.vstack([e_bbox, f_bbox])
    lo = allb[:, :3].min(axis=0)
    hi = allb[:, 3:].max(axis=0)
    span = hi - lo
    # Cell side = max per-axis element extent ⇒ each AABB spans ≤2 cells/axis.
    h = np.maximum(allb[:, 3:] - allb[:, :3], 0.0).max(axis=0)
    h = np.where(h <= 0.0, np.where(span > 0.0, span, 1.0), h)
    inv_h = np.ascontiguousarray(1.0 / h)
    n = np.maximum(np.floor(span * inv_h).astype(np.int64) + 1, 1)
    nx, ny, nz = int(n[0]), int(n[1]), int(n[2])
    origin = np.ascontiguousarray(lo)

    # Generous initial buffer; fused emission ≈ #true near-intersections, so
    # overflow (→ retry) effectively never fires on real geometry.
    cap = max((n_e + n_f) * 4, 1 << 16)
    res_e, res_f, total = _grid_ef_kernel(
        e_bbox, f_bbox, ep, eq, tris, origin, inv_h, nx, ny, nz, cap
    )
    if total > cap:  # rare: re-run with exact capacity
        cap = total
        res_e, res_f, total = _grid_ef_kernel(
            e_bbox, f_bbox, ep, eq, tris, origin, inv_h, nx, ny, nz, cap
        )
    if total == 0:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64))
    res_e = res_e[:total]
    res_f = res_f[:total]

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

    # 2+3. Broad phase + shared-vertex skip, fused. Uniform-grid AABB overlap
    # that never emits a pair sharing a mesh vertex (canonical labels, so plain
    # equality suffices — no vertex_class lookup). Equivalent to a plain
    # AABB-overlap broad phase followed by the post-filter, but avoids
    # materialising the ~1M 1-ring adjacency overlaps of a connected mesh
    # (G13/G14).
    ei, fi = edge_face_candidates(e_bbox, f_bbox, ep, eq, tris)
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
