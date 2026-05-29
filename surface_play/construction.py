"""construction.py — C13: `build_construction` orchestrator.

Sequential composition of the Layer C steps into a single
`ConstructionResult` (the Construction phase output, cached upstream per
spec §"Two computation phases"):

    build_mesh → build_bcs → find_double_points → build_sis_pairs
        → build_sics → find_triple_points

Pure data-flow; no side effects. All gotchas G1-G16 are handled by the
underlying functions; this layer adds none.
"""

from dataclasses import dataclass

import numpy as np

from surface_play.curves import BoundaryCurve, build_bcs
from surface_play.domain import Domain
from surface_play.intersections import (
    SelfIntersectingCurve,
    build_sics,
    build_sis_pairs,
    find_double_points,
    find_triple_points,
)
from surface_play.mesh import Mesh, build_mesh
from surface_play.surface import SurfaceParams


@dataclass
class ConstructionResult:
    mesh: Mesh
    bcs: list[BoundaryCurve]
    dps: np.ndarray
    sis_pairs: np.ndarray
    sics: list[SelfIntersectingCurve]
    tps: np.ndarray


def build_construction(
    domain: Domain,
    surface: SurfaceParams,
    resolution: int,
    *,
    jitter: bool = True,
    seed: int | None = None,
    mesh: Mesh | None = None,
) -> ConstructionResult:
    """Run the full Construction phase and return all artifacts.

    `mesh`: reuse an already-built mesh (e.g. the one cached for the 3D canvas
    display) instead of regenerating it. When ``None`` the mesh is built here.
    """
    if mesh is None:
        mesh = build_mesh(domain, surface, resolution, jitter=jitter, seed=seed)
    bcs = build_bcs(mesh)
    dps = find_double_points(mesh, surface)
    sis_pairs = build_sis_pairs(dps)
    sics = build_sics(sis_pairs)
    tps = find_triple_points(sis_pairs, dps, mesh, surface)
    return ConstructionResult(
        mesh=mesh,
        bcs=bcs,
        dps=dps,
        sis_pairs=sis_pairs,
        sics=sics,
        tps=tps,
    )
