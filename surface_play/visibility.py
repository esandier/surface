"""visibility.py — Layer O, steps O15-O17.

O15 — `compute_projection_breaks`: detect xy-plane crossings between
ResampledCurve segments; emit one break per crossing where the occluder is
a BC or CC.

O16 — `bfs_visibility`: propagate visibility from an anchor (vis=0) across
all RCs, through within-RC breaks and across SP-shared endpoints.
"""

from __future__ import annotations

from collections import deque
from typing import TYPE_CHECKING, Literal

import numpy as np

from surface_play.intersections import sweep_segments

if TYPE_CHECKING:
    from surface_play.curves import ResampledCurve
    from surface_play.projection import Projection
    from surface_play.splitting import SplitArrays
    from surface_play.surface import SurfaceParams


break_dtype = np.dtype([
    ("rc_idx",     "i4"),
    ("sample_idx", "i4"),
    ("t",          "f8"),
    ("xy",         "2f8"),
    ("delta_v",    "i1"),
])


def compute_projection_breaks(
    rcs: "list[ResampledCurve]",
    surface: "SurfaceParams",
    projection: "Projection",
) -> np.ndarray:
    """Detect xy-plane crossings between ResampledCurve segments → BKs.

    Spec §"Projection intersections" (lines 470-485). Per-hit:
    - depth on each side via linear interp of `rc.depth`; greater depth = OCCLUDED.
    - Occluder kind determines |delta_v|: 2 for CC, 1 for BC, else skip.
    - Sign: `-` if `T · dir > 0` where T = occluded tangent, dir = occluder's
      view-plane normal toward the surface (from `rc.dir`, populated by O14).
    - HC / SIC as occluder: never emit (they don't occlude).
    - Same-RC self-crossings included (an RC can occlude itself when its
      projection loops).
    """
    del surface  # not needed; depth already on rcs

    # Flatten segments: one row per polyline edge.
    p0s: list[np.ndarray] = []
    p1s: list[np.ndarray] = []
    seg_rc_idx: list[int] = []
    seg_sample_idx: list[int] = []
    for ri, rc in enumerate(rcs):
        n = len(rc.xy)
        for si in range(n - 1):
            p0s.append(rc.xy[si])
            p1s.append(rc.xy[si + 1])
            seg_rc_idx.append(ri)
            seg_sample_idx.append(si)

    if not p0s:
        return np.empty(0, dtype=break_dtype)

    p0 = np.asarray(p0s, dtype=float)
    p1 = np.asarray(p1s, dtype=float)
    seg_rc_idx_arr = np.asarray(seg_rc_idx, dtype=np.int32)
    seg_sample_idx_arr = np.asarray(seg_sample_idx, dtype=np.int32)

    hits = sweep_segments(p0, p1, None, None, None, self_sweep=True)
    if len(hits) == 0:
        return np.empty(0, dtype=break_dtype)

    out: list[tuple] = []
    for h in hits:
        a = int(h["a"]); b = int(h["b"])
        ta = float(h["t_a"]); tb = float(h["t_b"])
        xy_hit = np.asarray(h["uv"], dtype=float).copy()

        ri_a = int(seg_rc_idx_arr[a]); sa = int(seg_sample_idx_arr[a])
        ri_b = int(seg_rc_idx_arr[b]); sb = int(seg_sample_idx_arr[b])
        rc_a = rcs[ri_a]; rc_b = rcs[ri_b]

        depth_a = (1.0 - ta) * float(rc_a.depth[sa]) + ta * float(rc_a.depth[sa + 1])
        depth_b = (1.0 - tb) * float(rc_b.depth[sb]) + tb * float(rc_b.depth[sb + 1])

        # axis = I × J points TOWARD viewer, so projection.Z returns larger
        # values for points CLOSER to the viewer. The OCCLUDER is the one
        # with greater depth (closer); the OCCLUDED is the smaller-depth side.
        if depth_a <= depth_b:
            occl_rc, occl_ri, occl_si, t_occ = rc_a, ri_a, sa, ta
            ocer_rc, ocer_si, t_ocer = rc_b, sb, tb
        else:
            occl_rc, occl_ri, occl_si, t_occ = rc_b, ri_b, sb, tb
            ocer_rc, ocer_si, t_ocer = rc_a, sa, ta

        # Occluder kind → magnitude.
        if ocer_rc.kind == "CC":
            mag = 2
        elif ocer_rc.kind == "BC":
            mag = 1
        else:
            continue  # HC / SIC don't occlude

        if ocer_rc.dir is None:
            continue  # defensive: BC/CC should have dir from O14

        # T = analytic image-space tangent of the occluded RC at the crossing.
        if occl_rc.tan is not None:
            T = (
                (1.0 - t_occ) * np.asarray(occl_rc.tan[occl_si], dtype=float)
                + t_occ * np.asarray(occl_rc.tan[occl_si + 1], dtype=float)
            )
        else:
            T = occl_rc.xy[occl_si + 1] - occl_rc.xy[occl_si]

        # dir_ocer: in-image normal to the occluder's image tangent.
        # Sign chosen so that dir_ocer · old_dir > 0, where old_dir is the
        # lifted-uv-inward direction projected to image (rc.dir from O14).
        # The old_dir has the right 3D sign (uses surface orientation) but
        # is not perpendicular to the image tangent when the surface tilts
        # in the depth direction at the boundary. Decomposing old_dir into
        # along-tangent + perpendicular components and keeping just the
        # perpendicular gives an in-image-perpendicular direction with the
        # correct sign convention.
        old_dir = (
            (1.0 - t_ocer) * np.asarray(ocer_rc.dir[ocer_si], dtype=float)
            + t_ocer * np.asarray(ocer_rc.dir[ocer_si + 1], dtype=float)
        )
        if ocer_rc.tan is not None:
            ocer_tan = (
                (1.0 - t_ocer) * np.asarray(ocer_rc.tan[ocer_si], dtype=float)
                + t_ocer * np.asarray(ocer_rc.tan[ocer_si + 1], dtype=float)
            )
            t_norm2 = float(ocer_tan @ ocer_tan)
            if t_norm2 > 0.0:
                # Project old_dir onto plane perpendicular to ocer_tan.
                proj_coeff = float(old_dir @ ocer_tan) / t_norm2
                dir_ocer = old_dir - proj_coeff * ocer_tan
            else:
                dir_ocer = old_dir
        else:
            dir_ocer = old_dir

        sign = -1 if float(T @ dir_ocer) > 0.0 else 1
        delta_v = int(sign * mag)

        out.append((
            occl_ri, occl_si, t_occ,
            (float(xy_hit[0]), float(xy_hit[1])),
            delta_v,
        ))

    if not out:
        return np.empty(0, dtype=break_dtype)
    return np.array(out, dtype=break_dtype)


# ── O16 ───────────────────────────────────────────────────────────────────────

def _pick_anchors(rcs: "list[ResampledCurve]",
                  mode: Literal["leftmost", "extremes"]
                  ) -> list[tuple[int, int]]:
    """Return list of (rc_idx, sample_idx). G23 deterministic tie-break by
    (rc_idx, sample_idx).
    """
    if not rcs:
        return []
    if mode == "leftmost":
        best = None  # (x, rc_idx, sample_idx)
        for ri, rc in enumerate(rcs):
            for si in range(len(rc.xy)):
                key = (float(rc.xy[si, 0]), ri, si)
                if best is None or key < best:
                    best = key
        return [(best[1], best[2])]
    if mode == "extremes":
        cands = []
        for ri, rc in enumerate(rcs):
            for si in range(len(rc.xy)):
                cands.append((float(rc.xy[si, 0]), float(rc.xy[si, 1]), ri, si))
        leftmost = min(cands, key=lambda c: (c[0], c[2], c[3]))
        rightmost = min(cands, key=lambda c: (-c[0], c[2], c[3]))
        bottom = min(cands, key=lambda c: (c[1], c[2], c[3]))
        top = min(cands, key=lambda c: (-c[1], c[2], c[3]))
        seen: list[tuple[int, int]] = []
        for c in (leftmost, rightmost, bottom, top):
            key = (c[2], c[3])
            if key not in seen:
                seen.append(key)
        return seen
    raise ValueError(f"unknown anchor_mode {mode!r}")


def bfs_visibility(
    rcs: "list[ResampledCurve]",
    breaks: np.ndarray,
    splits: "SplitArrays",
    *,
    anchor_mode: Literal["leftmost", "extremes"] = "leftmost",
) -> dict[int, np.ndarray]:
    """Propagate visibility from anchor(s) (vis=0) through breaks and SP joins.

    Spec §"Anchor and BFS propagation" (lines 487-497). Returns `{id(rc):
    vis_per_sample}`. Unreachable RCs (islands) come back as all-zeros; O17
    LP will rescue them.
    """
    del splits  # not needed: SP-sharing is keyed by integer SP index on rc

    n = len(rcs)
    if n == 0:
        return {}

    # SP index → list of (rc_idx, "start"|"end").
    sp_to_eps: dict[int, list[tuple[int, str]]] = {}
    for ri, rc in enumerate(rcs):
        if rc.start >= 0:
            sp_to_eps.setdefault(int(rc.start), []).append((ri, "start"))
        # Closed RC (start == end): registering both would double-count.
        if rc.end >= 0 and rc.end != rc.start:
            sp_to_eps.setdefault(int(rc.end), []).append((ri, "end"))

    # Per-RC cumulative segment delta from breaks.
    seg_delta = [np.zeros(max(len(rc.xy) - 1, 0), dtype=np.int32) for rc in rcs]
    for bk in breaks:
        ri = int(bk["rc_idx"])
        si = int(bk["sample_idx"])
        if 0 <= ri < n and 0 <= si < len(seg_delta[ri]):
            seg_delta[ri][si] += int(bk["delta_v"])

    rc_vis: list[np.ndarray | None] = [None] * n

    def _fill(ri: int, si: int, v0: int) -> None:
        rc = rcs[ri]
        N = len(rc.xy)
        arr = np.zeros(N, dtype=np.int32)
        arr[si] = v0
        for k in range(si, N - 1):
            arr[k + 1] = arr[k] + int(seg_delta[ri][k])
        for k in range(si, 0, -1):
            arr[k - 1] = arr[k] - int(seg_delta[ri][k - 1])
        rc_vis[ri] = arr

    anchors = _pick_anchors(rcs, anchor_mode)

    queue: deque[int] = deque()
    for ri, si in anchors:
        if rc_vis[ri] is None:
            # If the anchor sample lies AT an SP (sample 0 or last), pin the
            # SP's shared visibility to 0 — not the sample value — because an
            # SP IS the silhouette boundary, so its shared vis is 0 by
            # definition. The rc's near-SP sample value then equals
            # `0 + delta_self` = vc_in (or vc_out), placing this rc on its
            # correct side of the silhouette. Interior anchors keep vis=0 at
            # the sample (no SP convention applies).
            rc = rcs[ri]
            N = len(rc.xy)
            if si == 0 and rc.start >= 0:
                v0 = int(rc.vc_in)
            elif si == N - 1 and rc.end >= 0:
                v0 = int(rc.vc_out)
            else:
                v0 = 0
            _fill(ri, si, v0)
            queue.append(ri)

    while queue:
        ri = queue.popleft()
        rc = rcs[ri]
        for which, sp_idx, sample_idx, delta_self in (
            ("start", rc.start, 0, int(rc.vc_in)),
            ("end", rc.end, len(rc.xy) - 1, int(rc.vc_out)),
        ):
            if sp_idx < 0:
                continue
            v_self = int(rc_vis[ri][sample_idx])
            for (nri, nwhich) in sp_to_eps.get(int(sp_idx), []):
                if nri == ri and nwhich == which:
                    continue
                if rc_vis[nri] is not None:
                    continue
                neighbor = rcs[nri]
                if nwhich == "start":
                    delta_n = int(neighbor.vc_in)
                    n_si = 0
                else:
                    delta_n = int(neighbor.vc_out)
                    n_si = len(neighbor.xy) - 1
                v_n = v_self - delta_self + delta_n
                _fill(nri, n_si, v_n)
                queue.append(nri)

    out: dict[int, np.ndarray] = {}
    for ri in range(n):
        if rc_vis[ri] is None:
            rc_vis[ri] = np.zeros(len(rcs[ri].xy), dtype=np.int32)
        out[id(rcs[ri])] = rc_vis[ri]
    return out


# ── O17 ───────────────────────────────────────────────────────────────────────

class LPInfeasibleError(RuntimeError):
    """G7: the LP visibility system has no feasible solution."""


def lp_refine_visibility(
    rcs: "list[ResampledCurve]",
    breaks: np.ndarray,
    splits: "SplitArrays",
    vis_bfs: dict[int, np.ndarray] | None = None,
    *,
    anchors: Literal["leftmost", "extremes"] = "leftmost",
) -> dict[int, np.ndarray]:
    """Refine visibility via a sparse L1-slack LP (spec §"LP pass", G24).

    Variables:
      - vis[k] per sample (upper-bounded at 0, free below).
      - vc_i per break (free).
      - s_i ≥ 0 (L1 slack for |vc_i − vc_i⁰|).

    Equality rows:
      - per-segment propagation: vis[k+1] − vis[k] − Σvc_i = 0
      - SP coupling: vis_o − vis_ref = δ_o − δ_ref (δ = vc_in/vc_out of RC)
      - anchor: vis = 0

    Inequalities: vc_i − s_i ≤ vc_i⁰ and −vc_i − s_i ≤ −vc_i⁰.
    Objective: Σ s_i. Solver: scipy.optimize.linprog(method="highs").
    Result: vis rounded to int. `vis_bfs` accepted for API symmetry; unused.
    On status=2: raise `LPInfeasibleError`.
    """
    from scipy.optimize import linprog
    from scipy.sparse import csr_matrix

    del splits, vis_bfs  # not needed in LP construction

    n = len(rcs)
    if n == 0:
        return {}

    # Variable layout.
    rc_offsets: list[int] = []
    offset = 0
    for rc in rcs:
        rc_offsets.append(offset)
        offset += len(rc.xy)
    N_vis = offset
    N_breaks = len(breaks)
    N_vars = N_vis + 2 * N_breaks
    # col(vis[ri][k]) = rc_offsets[ri] + k
    # col(vc_i)       = N_vis + i
    # col(s_i)        = N_vis + N_breaks + i

    # Breaks indexed by (rc, segment).
    bks_in_seg: dict[tuple[int, int], list[int]] = {}
    for bi, bk in enumerate(breaks):
        bks_in_seg.setdefault((int(bk["rc_idx"]), int(bk["sample_idx"])), []).append(bi)

    eq_rows: list[int] = []
    eq_cols: list[int] = []
    eq_data: list[float] = []
    eq_b: list[float] = []
    row = 0

    # Per-segment propagation.
    for ri, rc in enumerate(rcs):
        N = len(rc.xy)
        for k in range(N - 1):
            eq_rows.append(row); eq_cols.append(rc_offsets[ri] + k + 1); eq_data.append(1.0)
            eq_rows.append(row); eq_cols.append(rc_offsets[ri] + k);     eq_data.append(-1.0)
            for bi in bks_in_seg.get((ri, k), []):
                eq_rows.append(row); eq_cols.append(N_vis + bi);          eq_data.append(-1.0)
            eq_b.append(0.0)
            row += 1

    # SP coupling.
    sp_to_eps: dict[int, list[tuple[int, str]]] = {}
    for ri, rc in enumerate(rcs):
        if rc.start >= 0:
            sp_to_eps.setdefault(int(rc.start), []).append((ri, "start"))
        if rc.end >= 0 and rc.end != rc.start:
            sp_to_eps.setdefault(int(rc.end), []).append((ri, "end"))

    for sp_idx, eps in sp_to_eps.items():
        if len(eps) < 2:
            continue
        ref_ri, ref_which = eps[0]
        ref_si = 0 if ref_which == "start" else len(rcs[ref_ri].xy) - 1
        ref_delta = (int(rcs[ref_ri].vc_in) if ref_which == "start"
                     else int(rcs[ref_ri].vc_out))
        for (ri_o, which_o) in eps[1:]:
            si_o = 0 if which_o == "start" else len(rcs[ri_o].xy) - 1
            delta_o = (int(rcs[ri_o].vc_in) if which_o == "start"
                       else int(rcs[ri_o].vc_out))
            eq_rows.append(row); eq_cols.append(rc_offsets[ri_o] + si_o);  eq_data.append(1.0)
            eq_rows.append(row); eq_cols.append(rc_offsets[ref_ri] + ref_si); eq_data.append(-1.0)
            eq_b.append(float(delta_o - ref_delta))
            row += 1

    # Anchors. When the anchor sample sits at an SP endpoint of its RC, pin to
    # the SP's "shared" value (vc_in / vc_out), not 0 — matches bfs_visibility.
    anchor_pairs = _pick_anchors(rcs, anchors)
    for (ri, si) in anchor_pairs:
        rc = rcs[ri]
        N = len(rc.xy)
        if si == 0 and rc.start >= 0:
            v0 = float(rc.vc_in)
        elif si == N - 1 and rc.end >= 0:
            v0 = float(rc.vc_out)
        else:
            v0 = 0.0
        eq_rows.append(row); eq_cols.append(rc_offsets[ri] + si); eq_data.append(1.0)
        eq_b.append(v0)
        row += 1

    A_eq = csr_matrix(
        (np.asarray(eq_data, dtype=float),
         (np.asarray(eq_rows, dtype=int), np.asarray(eq_cols, dtype=int))),
        shape=(row, N_vars),
    ) if row > 0 else None
    b_eq = np.asarray(eq_b, dtype=float) if row > 0 else None

    # Inequalities.
    ub_rows: list[int] = []
    ub_cols: list[int] = []
    ub_data: list[float] = []
    ub_b: list[float] = []
    for bi in range(N_breaks):
        vc0 = float(breaks[bi]["delta_v"])
        # vc_i - s_i <= vc0
        r = 2 * bi
        ub_rows.append(r); ub_cols.append(N_vis + bi);              ub_data.append(1.0)
        ub_rows.append(r); ub_cols.append(N_vis + N_breaks + bi);   ub_data.append(-1.0)
        ub_b.append(vc0)
        # -vc_i - s_i <= -vc0
        r = 2 * bi + 1
        ub_rows.append(r); ub_cols.append(N_vis + bi);              ub_data.append(-1.0)
        ub_rows.append(r); ub_cols.append(N_vis + N_breaks + bi);   ub_data.append(-1.0)
        ub_b.append(-vc0)
    A_ub = csr_matrix(
        (np.asarray(ub_data, dtype=float),
         (np.asarray(ub_rows, dtype=int), np.asarray(ub_cols, dtype=int))),
        shape=(2 * N_breaks, N_vars),
    ) if N_breaks > 0 else None
    b_ub = np.asarray(ub_b, dtype=float) if N_breaks > 0 else None

    # Bounds.
    bounds: list[tuple[float | None, float | None]] = []
    bounds.extend([(None, 0.0)] * N_vis)
    bounds.extend([(None, None)] * N_breaks)
    bounds.extend([(0.0, None)] * N_breaks)

    # Objective: minimise Σ s_i.
    c = np.zeros(N_vars, dtype=float)
    c[N_vis + N_breaks : N_vis + 2 * N_breaks] = 1.0

    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=bounds, method="highs")

    if res.status == 2:
        raise LPInfeasibleError(f"LP infeasible: {res.message}")
    if not res.success:
        raise RuntimeError(f"LP failed: status={res.status} message={res.message}")

    vis_flat = np.round(res.x[:N_vis]).astype(np.int32)
    out: dict[int, np.ndarray] = {}
    for ri in range(n):
        N = len(rcs[ri].xy)
        out[id(rcs[ri])] = vis_flat[rc_offsets[ri] : rc_offsets[ri] + N]
    return out
