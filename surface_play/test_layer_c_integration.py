"""Layer C end-to-end integration tests.

Runs `build_construction` on every fixture and asserts the
`ConstructionResult` shape (BC / DP / SIC / TP counts) matches expectations.
Also checks determinism (same seed → identical results).
"""

from dataclasses import dataclass

import numpy as np
import pytest

from surface_play.construction import build_construction
from surface_play.test_fixtures import (
    disk_paraboloid_po, fig8, helicoid, mobius_u, paraboloid, torus,
)


@dataclass(frozen=True)
class Expected:
    name: str
    bc_count: int
    bc_all_closed: bool
    dp_min: int          # inclusive lower bound on DP count
    dp_max: int          # inclusive upper bound
    sic_count: int
    tp_count: int


CASES = [
    # Embedded rectangles: 1 closed BC (G16), no DPs/SICs/TPs.
    (helicoid,   Expected("helicoid",   1, True, 0, 0, 0, 0)),
    (paraboloid, Expected("paraboloid", 1, True, 0, 0, 0, 0)),
    # Closed surface: no boundary.
    (torus,      Expected("torus",      0, True, 0, 0, 0, 0)),
    # Möbius mo-no: one closed boundary (wraps twice under the twist).
    (mobius_u,   Expected("mobius_u",   1, True, 0, 0, 0, 0)),
    # Disk: one closed outer ring.
    (disk_paraboloid_po,
                 Expected("disk_po",    1, True, 0, 0, 0, 0)),
    # Fig-8 immersion: closed surface, 1 SIC, no TPs.
    (fig8,       Expected("fig8",       0, True, 50, 100, 1, 0)),
]


@pytest.mark.parametrize("factory,expected", CASES, ids=[e.name for _, e in CASES])
def test_layer_c_integration(factory, expected):
    surface = factory()
    resolution = 20 if expected.name == "fig8" else 15
    result = build_construction(surface.domain, surface, resolution=resolution,
                                seed=0)

    assert len(result.bcs) == expected.bc_count, (
        f"{expected.name}: BC count {len(result.bcs)} != {expected.bc_count}"
    )
    if expected.bc_all_closed:
        for i, bc in enumerate(result.bcs):
            assert bc.is_closed, f"{expected.name}: BC{i} not closed"

    n_dp = len(result.dps)
    assert expected.dp_min <= n_dp <= expected.dp_max, (
        f"{expected.name}: DP count {n_dp} outside "
        f"[{expected.dp_min}, {expected.dp_max}]"
    )

    assert len(result.sics) == expected.sic_count, (
        f"{expected.name}: SIC count {len(result.sics)} != {expected.sic_count}"
    )
    for i, sic in enumerate(result.sics):
        assert sic.is_closed, f"{expected.name}: SIC{i} not closed"

    assert len(result.tps) == expected.tp_count, (
        f"{expected.name}: TP count {len(result.tps)} != {expected.tp_count}"
    )


@pytest.mark.parametrize("factory", [fig8, helicoid, mobius_u])
def test_layer_c_determinism(factory):
    """Same seed → identical numeric arrays; equal-length object lists."""
    surface = factory()
    res = 20 if factory is fig8 else 15
    r1 = build_construction(surface.domain, surface, resolution=res, seed=7)
    r2 = build_construction(surface.domain, surface, resolution=res, seed=7)

    np.testing.assert_array_equal(r1.mesh.uv, r2.mesh.uv)
    np.testing.assert_array_equal(r1.mesh.tris, r2.mesh.tris)
    np.testing.assert_array_equal(r1.dps, r2.dps)
    np.testing.assert_array_equal(r1.sis_pairs, r2.sis_pairs)
    np.testing.assert_array_equal(r1.tps, r2.tps)

    assert len(r1.bcs) == len(r2.bcs)
    assert len(r1.sics) == len(r2.sics)
    for s1, s2 in zip(r1.sics, r2.sics):
        np.testing.assert_array_equal(s1.sis_indices, s2.sis_indices)
        assert s1.is_closed == s2.is_closed


def test_layer_c_different_seeds_perturb_mesh():
    """Different seeds produce different jittered meshes but same topology
    (BC / SIC / TP counts unchanged for fig-8)."""
    surface = fig8()
    r1 = build_construction(surface.domain, surface, resolution=20, seed=1)
    r2 = build_construction(surface.domain, surface, resolution=20, seed=2)

    assert not np.array_equal(r1.mesh.uv, r2.mesh.uv)
    assert len(r1.bcs) == len(r2.bcs)
    assert len(r1.sics) == len(r2.sics)
    assert len(r1.tps) == len(r2.tps)
