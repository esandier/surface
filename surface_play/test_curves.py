"""Tests for P4: make_lines."""

import numpy as np
import pytest
from surface_play.curves import make_lines


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
