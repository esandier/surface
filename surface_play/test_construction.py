"""Tests for C13: `build_construction` orchestrator."""

import numpy as np

from surface_play.construction import ConstructionResult, build_construction
from surface_play.test_fixtures import (
    disk_paraboloid_po, fig8, helicoid, mobius_u,
)


def test_build_construction_helicoid():
    """Helicoid rect-no-no res=15: non-empty BCs (one outer loop), no DPs/SICs/TPs."""
    surface = helicoid()
    result = build_construction(surface.domain, surface, resolution=15)

    assert isinstance(result, ConstructionResult)
    assert len(result.bcs) >= 1
    assert len(result.dps) == 0
    assert len(result.sis_pairs) == 0
    assert result.sics == []
    assert len(result.tps) == 0


def test_build_construction_mobius():
    """Möbius mo-no res=15: exactly 1 closed BC, no DPs/SICs/TPs (embedded)."""
    surface = mobius_u()
    result = build_construction(surface.domain, surface, resolution=15)

    assert len(result.bcs) == 1
    assert result.bcs[0].is_closed is True
    assert len(result.dps) == 0
    assert len(result.sis_pairs) == 0
    assert result.sics == []
    assert len(result.tps) == 0


def test_build_construction_disk():
    """Disk paraboloid res=20: exactly 1 closed BC, no DPs/SICs/TPs."""
    surface = disk_paraboloid_po()
    result = build_construction(surface.domain, surface, resolution=20)

    assert len(result.bcs) == 1
    assert result.bcs[0].is_closed is True
    assert len(result.dps) == 0
    assert len(result.sis_pairs) == 0
    assert result.sics == []
    assert len(result.tps) == 0


def test_build_construction_fig8():
    """Fig-8 cy-cy res=20: 0 BCs (no boundary), 1 SIC with DPs and SISs, 0 TPs."""
    surface = fig8()
    result = build_construction(surface.domain, surface, resolution=20)

    assert len(result.bcs) == 0
    assert len(result.dps) > 0
    assert len(result.sis_pairs) > 0
    assert len(result.sics) == 1
    assert result.sics[0].is_closed is True
    assert len(result.tps) == 0


def test_build_construction_determinism():
    """Same seed → identical numeric arrays and equal-length lists."""
    surface = fig8()
    r1 = build_construction(surface.domain, surface, resolution=20, seed=123)
    r2 = build_construction(surface.domain, surface, resolution=20, seed=123)

    np.testing.assert_array_equal(r1.mesh.uv, r2.mesh.uv)
    np.testing.assert_array_equal(r1.mesh.tris, r2.mesh.tris)
    np.testing.assert_array_equal(r1.dps, r2.dps)
    np.testing.assert_array_equal(r1.sis_pairs, r2.sis_pairs)
    np.testing.assert_array_equal(r1.tps, r2.tps)

    assert len(r1.bcs) == len(r2.bcs)
    assert len(r1.sics) == len(r2.sics)
    for c1, c2 in zip(r1.sics, r2.sics):
        np.testing.assert_array_equal(c1.sis_indices, c2.sis_indices)
        assert c1.is_closed == c2.is_closed
