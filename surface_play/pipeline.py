"""pipeline.py — W1: Construction/Outline orchestrator.

Pure assembly layer between Django views and the Layer C/O algorithms.
Wires `build_construction` (Layer C) and the full O1-O17 chain (Layer O)
into the two artifacts the frontend consumes:

  * `SurfaceInit`     — projection-independent (mesh + threejs payload),
                        cached per (record_pk, structural-hash).
  * `OutlineResult`   — projection-dependent visibility-split polylines.

No new algorithms; no silhouette imports (W1 sign-off gate, test 7).

Spec: `Modular_rewrite_roadmap.md` lines 1548-1610.
"""

from __future__ import annotations

import functools
import hashlib
from dataclasses import dataclass
from typing import Literal

import numpy as np

from surface_play import settings as _settings
from surface_play.construction import ConstructionResult, build_construction
from surface_play.contour import (
    build_contour_curves,
    build_contour_segments,
    find_contour_points,
    find_vps,
)
from surface_play.curves import resample_all
from surface_play.domain import Domain
from surface_play.helpers import build_helper_curves
from surface_play.projection import Projection
from surface_play.splitting import (
    SplitArrays,
    assemble_subcurves,
    split_at_cdps,
    split_bcs_at_bcps,
    split_bcs_at_bdps,
    split_bcs_at_corners,
    split_ccs_at_vps,
    split_sics_at_tps,
)
from surface_play.surface import SurfaceParams
from surface_play.visibility import (
    LPInfeasibleError,
    bfs_visibility,
    compute_projection_breaks,
    lp_refine_visibility,
)


# ── Dataclasses ──────────────────────────────────────────────────────────────


@dataclass
class SurfaceInit:
    record_pk: int
    cache_key: str
    domain: Domain
    surface: SurfaceParams
    construction: ConstructionResult
    threejs: dict


@dataclass
class OutlineResult:
    lines_by_visibility: dict[int, list[list[tuple[float, float]]]]
    si_lines_by_visibility: dict[int, list[list[tuple[float, float]]]]
    origin: tuple[float, float]


# ── mesh_to_threejs ──────────────────────────────────────────────────────────


def mesh_to_threejs(construction: ConstructionResult) -> dict:
    """Projection-independent three.js mesh payload.

    Returns one row per compacted vertex (post-C4 identification) — the
    `Mesh.uv/.xyz/.SN` arrays are already deduplicated, so `tris` indexes
    them directly with no further remap.
    """
    mesh = construction.mesh
    return {
        "vertices": mesh.xyz.tolist(),
        "faces": mesh.tris.tolist(),
        "normals": mesh.SN.tolist(),
        "uvs": mesh.uv.tolist(),
    }


# ── build_surface_init ───────────────────────────────────────────────────────


def _record_signature(record) -> tuple:
    """Hashable structural signature of a SurfaceRecord-like object.

    Duck-typed: any object exposing the fields below works (lets tests
    pass a SimpleNamespace without standing up the Django ORM).
    """
    return (
        str(record.X),
        str(record.Y),
        str(record.Z),
        str(record.parameter_names),
        float(record.u_min), float(record.u_max),
        float(record.v_min), float(record.v_max),
        str(record.u_identify), str(record.v_identify),
        str(record.domain_type), str(record.coord_type),
        float(record.r_min), float(record.r_max),
        str(record.output_type),
    )


def _build_domain_and_surface(sig: tuple) -> tuple[Domain, SurfaceParams]:
    (X, Y, Z, param_names,
     u_min, u_max, v_min, v_max,
     u_identify, v_identify,
     domain_type, coord_type,
     r_min, r_max,
     output_type) = sig

    if domain_type == "rect":
        bounds = (u_min, u_max, v_min, v_max)
        domain = Domain(
            type="rect", bounds=bounds,
            u_identify=u_identify, v_identify=v_identify,
            coord_type=coord_type,
        )
    else:
        # disk: bounds are (r_min, r_max, 0, 2π) per Domain contract.
        bounds = (r_min, r_max, 0.0, 2.0 * np.pi)
        domain = Domain(
            type="disk", bounds=bounds, coord_type=coord_type,
        )

    surface = SurfaceParams(
        X=X, Y=Y, Z=Z,
        parameter_names=param_names,
        domain=domain,
        output_type=output_type,
    )
    return domain, surface


@functools.lru_cache(maxsize=_settings.SURFACE_CACHE_SIZE)
def _cached_surface_init(
    record_pk: int,
    sig: tuple,
    resolution: int,
    jitter: bool,
    seed: int | None,
) -> SurfaceInit:
    domain, surface = _build_domain_and_surface(sig)
    construction = build_construction(
        domain, surface, resolution, jitter=jitter, seed=seed,
    )
    threejs = mesh_to_threejs(construction)
    cache_key = hashlib.sha256(
        repr(sig + (resolution,)).encode("utf-8")
    ).hexdigest()
    return SurfaceInit(
        record_pk=int(record_pk),
        cache_key=cache_key,
        domain=domain,
        surface=surface,
        construction=construction,
        threejs=threejs,
    )


def build_surface_init(
    record,
    *,
    resolution: int | None = None,
    jitter: bool = True,
    seed: int | None = None,
) -> SurfaceInit:
    """Build (and cache) the projection-independent half of the pipeline.

    Cache key: `(record.pk, structural-record-sig, resolution, jitter, seed)`.
    Repeated calls with identical args return the same `SurfaceInit`
    object (is-equality), enabling cheap cache hits across POST requests.
    """
    if resolution is None:
        resolution = _settings.GRID_RESOLUTION
    sig = _record_signature(record)
    return _cached_surface_init(int(record.pk), sig, int(resolution), bool(jitter), seed)


# ── build_outline ────────────────────────────────────────────────────────────


def _run_outline_pipeline(
    construction: ConstructionResult,
    projection: Projection,
    surface: SurfaceParams,
    *,
    newton_cusp: bool,
    canvas_resolution: int | None,
    project_resampled: bool,
    propagation: str,
):
    """Run O1→O17 and return (rcs, breaks, splits, vis_by_id).

    Structurally identical to `probe_visibility._build_outline` minus the
    diagnostic bundle. `vis_by_id` is the chosen visibility map keyed by
    `id(rc)` (BFS or LP per `propagation`).
    """
    mesh = construction.mesh
    bcs = construction.bcs
    dps = construction.dps
    sis_pairs = construction.sis_pairs
    sics = construction.sics
    tps = construction.tps

    cps = find_contour_points(mesh, projection, use_newton=newton_cusp)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    vps = find_vps(ccs, css, cps, surface, projection)

    splits = SplitArrays()
    split_bcs_at_corners(mesh, bcs, splits, projection)
    split_bcs_at_bcps(mesh, bcs, ccs, css, cps, splits, surface, projection)
    split_bcs_at_bdps(mesh, bcs, sics, sis_pairs, dps, splits, surface, projection)
    split_sics_at_tps(tps, sis_pairs, dps, splits, surface, projection)
    split_ccs_at_vps(vps, css, ccs, splits, surface, projection)
    split_at_cdps(mesh, css, cps, sis_pairs, dps, splits, surface, projection)

    hcs = build_helper_curves(
        bcs, ccs, sics, css, sis_pairs, cps, dps,
        mesh, splits, projection, surface, surface.domain,
    )
    subs = assemble_subcurves(bcs, ccs, sics, hcs, mesh, css, sis_pairs, splits)

    rcs = resample_all(
        subs, surface, projection, splits, mesh, css, sis_pairs, cps, dps,
        resolution=canvas_resolution,
        project_resampled=project_resampled,
    )

    breaks = compute_projection_breaks(rcs, surface, projection)
    bfs = bfs_visibility(rcs, breaks, splits)

    if propagation == "BFS":
        vis = bfs
    elif propagation in ("LP1", "LP4"):
        try:
            lp = lp_refine_visibility(rcs, breaks, splits, vis_bfs=bfs)
            vis = lp
        except LPInfeasibleError:
            vis = bfs
    else:
        raise ValueError(f"unknown propagation mode {propagation!r}")

    return rcs, breaks, vis


def _split_rc_by_visibility(
    rc, vis: np.ndarray, breaks_by_seg: dict[tuple[int, int], list],
    rc_idx: int,
) -> list[tuple[int, list[tuple[float, float]]]]:
    """Walk an RC's projected polyline, split at within-segment breaks.

    Returns list of `(vis_level, [(x, y), ...])` polyline chunks. Each
    polyline has a single, constant visibility integer. SIC RCs are
    returned as a single polyline (visibility = current `vis[0]`); the
    caller decides how to store them.

    `breaks_by_seg[(rc_idx, k)]` is a sorted list of `(t, delta_v)` for
    the segment `k → k+1`.
    """
    sx = rc.xy[:, 0]
    sy = rc.xy[:, 1]
    chunks: list[tuple[int, list[tuple[float, float]]]] = []

    if len(rc.xy) < 2:
        return chunks

    curr_pts: list[tuple[float, float]] = [(float(sx[0]), float(sy[0]))]
    curr_vis = int(vis[0])

    def _flush(next_vis: int, anchor: tuple[float, float] | None = None):
        nonlocal curr_pts, curr_vis
        if len(curr_pts) >= 2:
            chunks.append((curr_vis, curr_pts))
        curr_pts = [anchor] if anchor is not None else []
        curr_vis = next_vis

    for k in range(len(rc.xy) - 1):
        segment_breaks = breaks_by_seg.get((rc_idx, k), [])
        prev_t = 0.0
        running_vis = curr_vis
        for t, dv in segment_breaks:
            if t <= prev_t:
                running_vis += int(dv)
                continue
            xb = float(sx[k] + t * (sx[k + 1] - sx[k]))
            yb = float(sy[k] + t * (sy[k + 1] - sy[k]))
            curr_pts.append((xb, yb))
            _flush(running_vis + int(dv), anchor=(xb, yb))
            running_vis += int(dv)
            prev_t = t
        # No more breaks on this segment — append the segment endpoint.
        curr_pts.append((float(sx[k + 1]), float(sy[k + 1])))
        # curr_vis should equal running_vis at this point (segment-end
        # value); reconcile in case of float-equal multi-break above.
        curr_vis = running_vis

    if len(curr_pts) >= 2:
        chunks.append((curr_vis, curr_pts))
    return chunks


def build_outline(
    init: SurfaceInit,
    I,
    J,
    O,
    eye,
    *,
    newton_cusp: bool = True,
    canvas_resolution: int | None = None,
    project_resampled: bool = False,
    propagation: Literal["BFS", "LP1", "LP4"] = "BFS",
) -> OutlineResult:
    """Project the cached construction under `(I, J, O, eye)` and assemble
    per-visibility polylines."""
    projection = Projection(init.surface, I=I, J=J, O=O, eye=eye)

    rcs, breaks, vis = _run_outline_pipeline(
        init.construction, projection, init.surface,
        newton_cusp=newton_cusp,
        canvas_resolution=canvas_resolution,
        project_resampled=project_resampled,
        propagation=propagation,
    )

    breaks_by_seg: dict[tuple[int, int], list[tuple[float, int]]] = {}
    has_t_field = (
        breaks.dtype.names is not None and "t" in breaks.dtype.names
        if len(breaks) else False
    )
    for b in breaks:
        ri = int(b["rc_idx"])
        si = int(b["sample_idx"])
        t = float(b["t"]) if has_t_field else 0.5
        dv = int(b["delta_v"])
        breaks_by_seg.setdefault((ri, si), []).append((t, dv))
    for key in breaks_by_seg:
        breaks_by_seg[key].sort(key=lambda td: td[0])

    lines_by_visibility: dict[int, list[list[tuple[float, float]]]] = {}
    si_lines_by_visibility: dict[int, list[list[tuple[float, float]]]] = {}

    for ri, rc in enumerate(rcs):
        if rc.kind == "HC":
            continue  # G22 — invisible bridges.
        v = vis[id(rc)]
        bucket = (
            si_lines_by_visibility if rc.kind == "SIC"
            else lines_by_visibility
        )
        if rc.kind == "SIC":
            # SICs are drawn solid; do not split at projection breaks.
            xy = [(float(p[0]), float(p[1])) for p in rc.xy]
            if len(xy) >= 2:
                key = int(v[0])
                bucket.setdefault(key, []).append(xy)
            continue
        for vis_level, pts in _split_rc_by_visibility(rc, v, breaks_by_seg, ri):
            bucket.setdefault(int(vis_level), []).append(pts)

    origin_xy = projection.XY(np.zeros(3))
    origin = (float(origin_xy[0]), float(origin_xy[1]))

    return OutlineResult(
        lines_by_visibility=lines_by_visibility,
        si_lines_by_visibility=si_lines_by_visibility,
        origin=origin,
    )
