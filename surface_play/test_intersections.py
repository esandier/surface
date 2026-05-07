"""Tests for P5: sweep_segments."""

import math
import time

import numpy as np
import pytest

from surface_play.domain import Domain
from surface_play.intersections import intersect_dtype, sweep_segments


TWO_PI = 2 * math.pi


def test_cross_simple():
    """Two crossing diagonals on a unit square produce one hit at (0.5, 0.5)."""
    a0 = np.array([[0.0, 0.0]])
    a1 = np.array([[1.0, 1.0]])
    b0 = np.array([[0.0, 1.0]])
    b1 = np.array([[1.0, 0.0]])
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))

    hits = sweep_segments(a0, a1, b0, b1, domain)

    assert hits.dtype == intersect_dtype
    assert len(hits) == 1
    np.testing.assert_allclose(hits["uv"][0], (0.5, 0.5), atol=1e-12)
    assert int(hits["a"][0]) == 0
    assert int(hits["b"][0]) == 0
    assert abs(hits["t_a"][0] - 0.5) < 1e-12
    assert abs(hits["t_b"][0] - 0.5) < 1e-12


def test_self_same_segment_no_hit():
    """Self-mode, same segment included twice → no hit (i<j filter; collinear → singular)."""
    a0 = np.array([[0.0, 0.0], [0.0, 0.0]])
    a1 = np.array([[1.0, 1.0], [1.0, 1.0]])
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))

    hits = sweep_segments(a0, a1, None, None, domain, self_sweep=True)

    assert len(hits) == 0


def test_close_aware_seam_crossing():
    """Cyl-u rect (0, 2π): a wrapping horizontal segment vs a vertical near the seam.

    Without `domain`, the wrap is not understood: the parametric intersection
    of the literal lines lies at t_a ≈ 1.008 (just outside (0, 1)) and the
    bbox-x ranges don't even overlap → no hit. With `domain`, close-aware
    anchoring of seg_a's q to (0.10 + 2π, 0.5) places the cross at t_a ≈ 0.75,
    t_b = 0.5 → exactly one hit.
    """
    domain_cy = Domain(type="rect", bounds=(0.0, TWO_PI, 0.0, 1.0), u_identify="cy")
    a0 = np.array([[6.18, 0.5]])
    a1 = np.array([[0.10, 0.5]])
    b0 = np.array([[0.05, 0.4]])
    b1 = np.array([[0.05, 0.6]])

    hits_with = sweep_segments(a0, a1, b0, b1, domain_cy)
    assert len(hits_with) == 1
    assert int(hits_with["a"][0]) == 0
    assert int(hits_with["b"][0]) == 0
    np.testing.assert_allclose(hits_with["uv"][0], (0.05 + TWO_PI, 0.5), atol=1e-9)
    assert 0.0 < hits_with["t_a"][0] < 1.0
    assert 0.0 < hits_with["t_b"][0] < 1.0
    assert abs(hits_with["t_b"][0] - 0.5) < 1e-9

    hits_without = sweep_segments(a0, a1, b0, b1, None)
    assert len(hits_without) == 0


def test_parallel_non_overlapping_no_hit():
    """Two horizontal parallel segments at distinct y → no hit (denom = 0)."""
    a0 = np.array([[0.0, 0.0]])
    a1 = np.array([[1.0, 0.0]])
    b0 = np.array([[0.0, 0.5]])
    b1 = np.array([[1.0, 0.5]])
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))

    hits = sweep_segments(a0, a1, b0, b1, domain)

    assert len(hits) == 0


def test_empty_inputs():
    """Empty input arrays → empty result, both cross and self modes."""
    empty = np.empty((0, 2))
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))

    out_cross = sweep_segments(empty, empty, empty, empty, domain)
    assert len(out_cross) == 0
    assert out_cross.dtype == intersect_dtype

    out_self = sweep_segments(empty, empty, None, None, domain, self_sweep=True)
    assert len(out_self) == 0
    assert out_self.dtype == intersect_dtype


def test_performance_smoke():
    """1000 random short segments, self-sweep, completes in <1s with sane hit count."""
    rng = np.random.default_rng(42)
    n = 1000
    p0 = rng.uniform(0.0, 1.0, (n, 2))
    p1 = p0 + rng.uniform(-0.1, 0.1, (n, 2))
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))

    t0 = time.perf_counter()
    hits = sweep_segments(p0, p1, None, None, domain, self_sweep=True)
    elapsed = time.perf_counter() - t0

    assert elapsed < 1.0, f"sweep took {elapsed:.3f}s, expected < 1s"
    assert hits.dtype == intersect_dtype
    # Sanity: short random segments in a unit square → some hits but not pathological.
    assert 0 <= len(hits) < n * n // 4
