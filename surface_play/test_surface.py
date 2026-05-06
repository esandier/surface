"""P2 tests: SurfaceParams parametrization, perturbation, lambdified callables.

Criteria from Modular_rewrite_roadmap.md §P2:
1. Helicoid sanity: S(0,0) ≈ 0; SN(0,0) parallel to (0, 1, 0).
2. Perturbation magnitude ≤ 0.001·bbox_diag over a 100-point grid.
3. Identification preservation under perturbation (cy on u/v, mo on u/v).
4. Vectorized: S(u_array, v_array) shape (3, N) and matches scalar calls.
5. cse speedup ≥ 3× on a 1000-point grid.
"""

import math
import time

import numpy as np
import pytest

from surface_play.domain import Domain
from surface_play.surface import SurfaceParams
from surface_play.test_fixtures import (
    cylinder_cy, disk_paraboloid_ca, disk_paraboloid_po,
    helicoid, mobius_u, mobius_v, paraboloid, torus, TWO_PI,
)


# ── Test 1: Helicoid sanity ──────────────────────────────────────────────────


def test_helicoid_origin_unperturbed():
    s = SurfaceParams(
        "u*cos(v)", "u*sin(v)", "v", "u v",
        Domain(type="rect", bounds=(-1.0, 1.0, -math.pi, math.pi)),
        perturb=False,
    )
    val = s.S(0.0, 0.0)
    assert val.shape == (3,)
    np.testing.assert_allclose(val, [0.0, 0.0, 0.0], atol=1e-12)

    sn = s.SN(0.0, 0.0)
    # Su(0,0) = (1, 0, 0); Sv(0,0) = (0, 0, 1); SN = (0, -1, 0) — parallel to (0, 1, 0).
    sn_norm = sn / np.linalg.norm(sn)
    target = np.array([0.0, 1.0, 0.0])
    cos_angle = abs(float(np.dot(sn_norm, target)))
    assert cos_angle == pytest.approx(1.0, abs=1e-12)


# ── Test 2: Perturbation magnitude ───────────────────────────────────────────


def test_perturbation_magnitude_below_threshold():
    s_perturbed = helicoid(perturb=True)
    s_pure = helicoid(perturb=False)

    u_min, u_max, v_min, v_max = s_perturbed.domain.bounds
    grid = 10  # 10 × 10 = 100 points
    us = np.linspace(u_min, u_max, grid)
    vs = np.linspace(v_min, v_max, grid)
    UU, VV = np.meshgrid(us, vs, indexing="ij")

    s_p = s_perturbed.S(UU.ravel(), VV.ravel())  # (3, 100)
    s_o = s_pure.S(UU.ravel(), VV.ravel())

    diff = s_p - s_o
    norms = np.linalg.norm(diff, axis=0)
    max_diff = float(norms.max())
    assert max_diff < 0.001 * s_pure.bbox_diag, (
        f"max ‖δS‖ = {max_diff:.3e} >= {0.001 * s_pure.bbox_diag:.3e}"
    )


# ── Test 3: Identification preservation under perturbation ───────────────────


def test_identification_cy_u_torus():
    s = torus(perturb=True)
    u_min, u_max, _, _ = s.domain.bounds
    vs = np.linspace(0.1, TWO_PI - 0.1, 7)
    for v in vs:
        np.testing.assert_allclose(
            s.S(u_min, float(v)), s.S(u_max, float(v)), atol=1e-12,
        )


def test_identification_cy_v_torus():
    s = torus(perturb=True)
    _, _, v_min, v_max = s.domain.bounds
    us = np.linspace(0.1, TWO_PI - 0.1, 7)
    for u in us:
        np.testing.assert_allclose(
            s.S(float(u), v_min), s.S(float(u), v_max), atol=1e-12,
        )


def test_identification_mo_u_mobius():
    s = mobius_u(perturb=True)
    u_min, u_max, v_min, v_max = s.domain.bounds
    vs = np.linspace(v_min + 0.05, v_max - 0.05, 7)
    for v in vs:
        v_mirror = v_max + v_min - float(v)
        np.testing.assert_allclose(
            s.S(u_min, float(v)), s.S(u_max, v_mirror), atol=1e-12,
        )


def test_identification_mo_v_mobius():
    s = mobius_v(perturb=True)
    u_min, u_max, v_min, v_max = s.domain.bounds
    us = np.linspace(u_min + 0.05, u_max - 0.05, 7)
    for u in us:
        u_mirror = u_max + u_min - float(u)
        np.testing.assert_allclose(
            s.S(float(u), v_min), s.S(u_mirror, v_max), atol=1e-12,
        )


# ── Test 4: Vectorized eval ──────────────────────────────────────────────────


def test_vectorized_S_shape_matches_scalar():
    s = paraboloid(perturb=False)
    N = 7
    us = np.linspace(-1.0, 1.0, N)
    vs = np.linspace(-1.0, 1.0, N)

    out = s.S(us, vs)
    assert out.shape == (3, N)

    for i in range(N):
        scalar = s.S(float(us[i]), float(vs[i]))
        assert scalar.shape == (3,)
        np.testing.assert_allclose(out[:, i], scalar, atol=1e-12)


def test_vectorized_all_seven_outputs():
    s = paraboloid(perturb=True)
    N = 5
    us = np.linspace(-0.5, 0.5, N)
    vs = np.linspace(-0.5, 0.5, N)
    for fn_name in ("S", "Su", "Sv", "Suu", "Suv", "Svv", "SN"):
        fn = getattr(s, fn_name)
        out = fn(us, vs)
        assert out.shape == (3, N), f"{fn_name} returned shape {out.shape}"


# ── Test 5: cse speedup ──────────────────────────────────────────────────────


def test_cse_speedup_at_least_3x():
    s = torus(perturb=True)
    u_min, u_max, v_min, v_max = s.domain.bounds
    us = np.linspace(u_min, u_max, 1000)
    vs = np.linspace(v_min, v_max, 1000)

    # Warm up (JIT/cache)
    for _ in range(2):
        s._eval_all(us, vs)
        s._naive_S(us, vs); s._naive_Su(us, vs); s._naive_Sv(us, vs)
        s._naive_Suu(us, vs); s._naive_Suv(us, vs); s._naive_Svv(us, vs)
        s._naive_SN(us, vs)

    REPEATS = 10

    # Combined cse'd evaluator (single call to _eval_fn inside _eval_all).
    t0 = time.perf_counter()
    for _ in range(REPEATS):
        s._eval_fn(us, vs)
    t_cse = time.perf_counter() - t0

    # Naive: 7 separate lambdified calls, no cse, no shared subexpressions.
    t0 = time.perf_counter()
    for _ in range(REPEATS):
        s._naive_S(us, vs)
        s._naive_Su(us, vs)
        s._naive_Sv(us, vs)
        s._naive_Suu(us, vs)
        s._naive_Suv(us, vs)
        s._naive_Svv(us, vs)
        s._naive_SN(us, vs)
    t_naive = time.perf_counter() - t0

    speedup = t_naive / t_cse
    print(
        f"\n[cse speedup] naive={t_naive*1000:.1f}ms  cse={t_cse*1000:.1f}ms"
        f"  -> {speedup:.2f}x"
    )
    assert speedup >= 3.0, f"cse speedup {speedup:.2f}x < 3x"


# ── coord_type='po': polar→cartesian substitution ───────────────────────────


def test_po_paraboloid_evaluates_in_cartesian():
    """After polar→cartesian substitution the eval interface takes (x, y)."""
    s = disk_paraboloid_po(perturb=False)
    # Internal _u, _v should be cartesian symbols
    assert str(s._u) == "_x"
    assert str(s._v) == "_y"

    # Z = x² + y² after substitution + simplify
    val = s.S(0.3, 0.4)
    np.testing.assert_allclose(val, [0.3, 0.4, 0.3**2 + 0.4**2], atol=1e-12)

    # Vectorized call
    xs = np.array([0.0, 0.3, -0.5, 0.5])
    ys = np.array([0.0, 0.4, 0.0, -0.5])
    out = s.S(xs, ys)
    expected = np.stack([xs, ys, xs**2 + ys**2])
    np.testing.assert_allclose(out, expected, atol=1e-12)


def test_po_and_ca_paraboloid_agree():
    """Same paraboloid written in polar vs cartesian gives identical S, Su, Sv, SN."""
    s_po = disk_paraboloid_po(perturb=False)
    s_ca = disk_paraboloid_ca(perturb=False)

    rng = np.random.default_rng(0)
    rs = rng.uniform(0.1, 0.95, 25)
    thetas = rng.uniform(0.0, TWO_PI, 25)
    xs = rs * np.cos(thetas)
    ys = rs * np.sin(thetas)

    for fn_name in ("S", "Su", "Sv", "Suu", "Suv", "Svv", "SN"):
        po_val = getattr(s_po, fn_name)(xs, ys)
        ca_val = getattr(s_ca, fn_name)(xs, ys)
        np.testing.assert_allclose(
            po_val, ca_val, atol=1e-10,
            err_msg=f"{fn_name} differs between po and ca formulations",
        )


def test_po_paraboloid_normal_orientation():
    """SN of z=x²+y² at (x, y) is parallel to (-2x, -2y, 1)."""
    s = disk_paraboloid_po(perturb=False)
    sn = s.SN(0.3, 0.4)
    target = np.array([-2 * 0.3, -2 * 0.4, 1.0])
    sn_norm = sn / np.linalg.norm(sn)
    target_norm = target / np.linalg.norm(target)
    cos_angle = abs(float(np.dot(sn_norm, target_norm)))
    assert cos_angle == pytest.approx(1.0, abs=1e-12)


def test_po_rejected_on_rect_domain():
    """coord_type='po' is only valid for disk/annulus."""
    with pytest.raises(ValueError, match="coord_type='po'"):
        SurfaceParams(
            "u*cos(v)", "u*sin(v)", "u**2", "u v",
            Domain(type="rect", bounds=(0, 1, 0, TWO_PI), coord_type="po"),
            perturb=False,
        )


# ── Bonus: output equivalence between cse and naive paths ────────────────────


def test_cse_results_match_naive():
    s = torus(perturb=True)
    us = np.linspace(0.1, TWO_PI - 0.1, 50)
    vs = np.linspace(0.1, TWO_PI - 0.1, 50)

    cse_S, cse_Su, cse_Sv, cse_Suu, cse_Suv, cse_Svv, cse_SN = s._eval_all(us, vs)

    nv_S = np.array(s._naive_S(us, vs))
    nv_Su = np.array(s._naive_Su(us, vs))
    nv_Sv = np.array(s._naive_Sv(us, vs))
    nv_Suu = np.array(s._naive_Suu(us, vs))
    nv_Suv = np.array(s._naive_Suv(us, vs))
    nv_Svv = np.array(s._naive_Svv(us, vs))
    nv_SN = np.array(s._naive_SN(us, vs))

    np.testing.assert_allclose(cse_S, nv_S, atol=1e-12)
    np.testing.assert_allclose(cse_Su, nv_Su, atol=1e-12)
    np.testing.assert_allclose(cse_Sv, nv_Sv, atol=1e-12)
    np.testing.assert_allclose(cse_Suu, nv_Suu, atol=1e-12)
    np.testing.assert_allclose(cse_Suv, nv_Suv, atol=1e-12)
    np.testing.assert_allclose(cse_Svv, nv_Svv, atol=1e-12)
    np.testing.assert_allclose(cse_SN, nv_SN, atol=1e-12)
