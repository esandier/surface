"""Tests for P4: make_lines, P6: sign_changes, C7: build_bcs."""

import math

import numpy as np
import pytest
from surface_play.curves import BoundaryCurve, build_bcs, make_lines, sign_changes
from surface_play.domain import Domain
from surface_play.mesh import build_mesh
from surface_play.surface import SurfaceParams
from surface_play.test_fixtures import (
    disk_paraboloid_ca,
    mobius_u,
    paraboloid,
    torus,
)

TWO_PI = 2 * math.pi


def _abs_chain(chain: np.ndarray) -> list[int]:
    """Return segment indices (0-based, unsigned) in traversal order."""
    return [abs(int(x)) - 1 for x in chain]


def _is_closed(chain: np.ndarray) -> bool:
    """A chain is closed if it ends with the same segment as it began (loop token)."""
    # Closed: last entry's abs == first entry's abs (the loop closes)
    # Actually by our convention a closed loop appends the closing edge again.
    # Simpler check: build the vertex sequence and see if start == end.
    return False  # overridden below per test


# ---------------------------------------------------------------------------
# Test 1: 4 segments forming an open path  0-1-2-3-4  (not closed)
# ---------------------------------------------------------------------------
class TestOpenPath:
    # Vertices 0,1,2,3,4 connected: (0,1),(1,2),(2,3),(3,4)
    segs = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])

    def test_one_chain(self):
        chains = make_lines(self.segs)
        assert len(chains) == 1

    def test_length_4(self):
        chains = make_lines(self.segs)
        chain = chains[0]
        # 4 edges, all visited
        assert len(_abs_chain(chain)) == 4

    def test_all_segments_covered(self):
        chains = make_lines(self.segs)
        covered = set(_abs_chain(chains[0]))
        assert covered == {0, 1, 2, 3}


# ---------------------------------------------------------------------------
# Test 2: 4 segments forming a closed quad  0-1-2-3-0
# ---------------------------------------------------------------------------
class TestClosedQuad:
    segs = np.array([[0, 1], [1, 2], [2, 3], [3, 0]])

    def test_one_chain(self):
        chains = make_lines(self.segs)
        assert len(chains) == 1

    def test_all_4_segments(self):
        chains = make_lines(self.segs)
        abs_ids = _abs_chain(chains[0])
        # closed chain has 5 entries: 4 edges + repeat of closing edge
        assert set(abs_ids) == {0, 1, 2, 3}

    def test_is_closed(self):
        chains = make_lines(self.segs)
        chain = chains[0]
        abs_ids = _abs_chain(chain)
        # first and last entry refer to the same segment (loop closure)
        assert abs_ids[0] == abs_ids[-1]


# ---------------------------------------------------------------------------
# Test 3: Two disconnected loops
# ---------------------------------------------------------------------------
class TestTwoLoops:
    # Loop A: 0-1-2-0
    # Loop B: 3-4-5-3
    segs = np.array([
        [0, 1], [1, 2], [2, 0],   # loop A
        [3, 4], [4, 5], [5, 3],   # loop B
    ])

    def test_two_chains(self):
        chains = make_lines(self.segs)
        assert len(chains) == 2

    def test_both_closed(self):
        chains = make_lines(self.segs)
        for chain in chains:
            abs_ids = _abs_chain(chain)
            assert abs_ids[0] == abs_ids[-1], "each loop should close on itself"

    def test_all_6_segments_covered(self):
        chains = make_lines(self.segs)
        covered = set()
        for c in chains:
            covered |= set(_abs_chain(c))
        assert covered == {0, 1, 2, 3, 4, 5}


# ---------------------------------------------------------------------------
# Test 4: empty input
# ---------------------------------------------------------------------------
def test_empty_input():
    assert make_lines(np.empty((0, 2), dtype=int)) == []


# ---------------------------------------------------------------------------
# P6: sign_changes
# ---------------------------------------------------------------------------
class TestSignChanges:
    def test_default_flip(self):
        vals_p = np.array([1, -1, 2, -3])
        vals_q = np.array([-1, 1, 3, -2])
        mask = sign_changes(vals_p, vals_q)
        assert mask.tolist() == [True, True, False, False]

    def test_with_flip(self):
        vals_p = np.array([1, -1, 2, -3])
        vals_q = np.array([-1, 1, 3, -2])
        flip = np.array([-1, 1, 1, 1])
        mask = sign_changes(vals_p, vals_q, flip)
        assert mask.tolist() == [False, True, False, False]

    def test_empty(self):
        mask = sign_changes(np.array([]), np.array([]))
        assert mask.shape == (0,)
        assert mask.dtype == bool


# ---------------------------------------------------------------------------
# C7: build_bcs
# ---------------------------------------------------------------------------

def _bcs(surf, resolution=8):
    """Build a mesh and return its BoundaryCurves."""
    return build_bcs(build_mesh(surf.domain, surf, resolution=resolution, seed=42))


def test_build_bcs():
    # 1. Rect no-id: 1 BC, closed (G16 — four sides form one loop).
    bcs = _bcs(paraboloid(perturb=False))
    assert len(bcs) == 1
    assert bcs[0].is_closed

    # 2. Rect cy-no: 2 BCs, both closed (top and bottom rims).
    domain_cy = Domain(type="rect", bounds=(0.0, 1.0, 0.0, TWO_PI), u_identify="cy")
    surf_cy = SurfaceParams("u*cos(v)", "u*sin(v)", "v", "u v", domain_cy, perturb=False)
    bcs = _bcs(surf_cy)
    assert len(bcs) == 2
    assert all(bc.is_closed for bc in bcs)

    # 3. Rect cy-cy: 0 BCs.
    bcs = _bcs(torus(perturb=False))
    assert len(bcs) == 0

    # 4. Rect mo-no: 1 BC, closed (mo merges the two v-side open arcs into one loop).
    bcs = _bcs(mobius_u(perturb=False))
    assert len(bcs) == 1
    assert bcs[0].is_closed

    # 5. Rect mo-mo: 0 BCs.
    domain_momo = Domain(
        type="rect", bounds=(0.0, TWO_PI, -0.3, 0.3),
        u_identify="mo", v_identify="mo",
    )
    surf_momo = SurfaceParams("cos(u)", "sin(u)", "v", "u v", domain_momo, perturb=False)
    bcs = build_bcs(build_mesh(domain_momo, surf_momo, resolution=8, seed=42))
    assert len(bcs) == 0

    # 6. Disk: 1 closed BC.
    bcs = _bcs(disk_paraboloid_ca(perturb=False), resolution=20)
    assert len(bcs) == 1
    assert bcs[0].is_closed

    # 7. Annulus: 2 closed BCs.
    domain_ann = Domain(type="annulus", bounds=(0.3, 1.0, 0.0, TWO_PI))
    surf_ann = SurfaceParams("x", "y", "x**2 + y**2", "x y", domain_ann, perturb=False)
    bcs = build_bcs(build_mesh(domain_ann, surf_ann, resolution=20, seed=42))
    assert len(bcs) == 2
    assert all(bc.is_closed for bc in bcs)
