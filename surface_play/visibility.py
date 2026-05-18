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

        # Greater depth = farther from viewer = OCCLUDED.
        if depth_a >= depth_b:
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

        T = occl_rc.xy[occl_si + 1] - occl_rc.xy[occl_si]
        dir_ocer = (
            (1.0 - t_ocer) * np.asarray(ocer_rc.dir[ocer_si], dtype=float)
            + t_ocer * np.asarray(ocer_rc.dir[ocer_si + 1], dtype=float)
        )

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
            _fill(ri, si, 0)
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
