import math
import numpy as np
import pytest

from surface_play.domain import Domain

TWO_PI = 2 * math.pi


def cyl_rect():
    """Cylinder-identified rectangle: u ∈ [0, 2π], v ∈ [0, 1]."""
    return Domain(type="rect", bounds=(0.0, TWO_PI, 0.0, 1.0), u_identify="cy")


def plain_rect():
    return Domain(type="rect", bounds=(0.0, TWO_PI, 0.0, 1.0))


def unit_square():
    return Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))


def unit_disk():
    return Domain(type="disk", bounds=(0.0, 1.0, 0.0, TWO_PI), coord_type="po")


# ── Test 1: close wraps u across identification ──────────────────────────────

def test_close_cyl_wraps_u():
    d = cyl_rect()
    p = np.array([0.1, 0.5])
    q = np.array([6.0, 0.5])
    result = d.close(p, q)
    expected = np.array([6.0 - TWO_PI, 0.5])
    np.testing.assert_allclose(result, expected, atol=1e-12)


# ── Test 2: close leaves nearby point unchanged ───────────────────────────────

def test_close_no_wrap_needed():
    d = cyl_rect()
    p = np.array([0.1, 0.5])
    q = np.array([0.2, 0.5])
    result = d.close(p, q)
    np.testing.assert_allclose(result, np.array([0.2, 0.5]), atol=1e-12)


def test_close_unidentified_axis_unchanged():
    d = cyl_rect()  # v_identify="no"
    p = np.array([0.1, 0.1])
    q = np.array([0.2, 0.9])  # large v shift — must not be wrapped
    result = d.close(p, q)
    assert result[1] == pytest.approx(0.9)


# ── Test 3: vectorized close matches element-wise ────────────────────────────

def test_close_vectorized():
    d = cyl_rect()
    ps = np.array([[0.1, 0.5], [0.1, 0.5], [3.0, 0.0]])
    qs = np.array([[6.0, 0.5], [0.2, 0.5], [6.0, 0.7]])
    result = d.close(ps, qs)
    expected = np.stack([
        d.close(ps[i], qs[i]) for i in range(len(ps))
    ])
    np.testing.assert_allclose(result, expected, atol=1e-12)


def test_close_disk_unchanged():
    d = unit_disk()
    p = np.array([0.5, 1.0])
    q = np.array([0.8, 5.0])
    result = d.close(p, q)
    np.testing.assert_allclose(result, q, atol=1e-15)


# ── Test 4: bary_edge and bary_face interpolation ────────────────────────────

def test_bary_edge():
    d = plain_rect()
    p = np.array([1.0, 2.0])
    pq = np.array([4.0, 0.0])
    np.testing.assert_allclose(d.bary_edge(p, pq, 0.0), [1.0, 2.0])
    np.testing.assert_allclose(d.bary_edge(p, pq, 0.5), [3.0, 2.0])
    np.testing.assert_allclose(d.bary_edge(p, pq, 1.0), [5.0, 2.0])


def test_bary_face():
    d = plain_rect()
    p = np.array([0.0, 0.0])
    pq = np.array([1.0, 0.0])
    pr = np.array([0.0, 1.0])
    np.testing.assert_allclose(d.bary_face(p, pq, pr, 0.0, 0.0), [0.0, 0.0])
    np.testing.assert_allclose(d.bary_face(p, pq, pr, 1.0, 0.0), [1.0, 0.0])
    np.testing.assert_allclose(d.bary_face(p, pq, pr, 0.0, 1.0), [0.0, 1.0])
    np.testing.assert_allclose(d.bary_face(p, pq, pr, 0.5, 0.25), [0.5, 0.25])


# ── Test 5: bbox_diag and width sanity values ─────────────────────────────────

def test_bbox_unit_square():
    d = unit_square()
    assert d.bbox_diag() == pytest.approx(math.sqrt(2))
    assert d.width() == pytest.approx(math.sqrt(2))


def test_bbox_unit_disk():
    d = unit_disk()
    assert d.width() == pytest.approx(2.0)   # 2·r_max = 2·1
    assert d.bbox_diag() == pytest.approx(2.0)
    assert d.bbox_diag() > 0.0


def test_period_properties():
    d = cyl_rect()
    assert d.period_u == pytest.approx(TWO_PI)
    assert math.isnan(d.period_v)

    d2 = Domain(type="rect", bounds=(0.0, 1.0, 0.0, TWO_PI), v_identify="mo")
    assert math.isnan(d2.period_u)
    assert d2.period_v == pytest.approx(TWO_PI)
