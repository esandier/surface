"""P3 tests: Projection (ortho/perspective).

Criteria from Modular_rewrite_roadmap.md §P3:
1. Ortho I=[1,0,0], J=[0,1,0], O=None: XY([1,2,3]) == (1, 2); Z([1,2,3]) == 3.
2. Ortho with non-default O=[5,0,0]: XY([6,2,3]) == (1, 2).
3. Persp eye=[0,0,5], O=[0,0,5]: XY([0,0,5]) == (0, 0).
4. Persp O != eye → ValueError.
5. per_vertex_viewer_dot shape (N,), matches per-vertex loop on a 50-vertex mesh.

Persp ker_param / proj_vec follow-up sub-step tests:
6. Persp ker_param converges to ortho ker_param as eye recedes along axis.
7. Persp ker_param at a known persp contour point: kerdS parallel to (S - eye).
8. Persp proj_vec matches finite-difference of XY along dir3d.
"""

import math

import numpy as np
import pytest

from surface_play.domain import Domain
from surface_play.projection import Projection
from surface_play.surface import SurfaceParams
from surface_play.test_fixtures import helicoid, paraboloid


# ── Test 1: ortho default O ──────────────────────────────────────────────────


def test_ortho_default_O():
    surf = paraboloid(perturb=False)
    proj = Projection(surf, I=[1, 0, 0], J=[0, 1, 0])
    assert proj.mode == "ortho"
    np.testing.assert_allclose(proj.XY(np.array([1.0, 2.0, 3.0])), [1.0, 2.0])
    assert proj.Z(np.array([1.0, 2.0, 3.0])) == pytest.approx(3.0)


# ── Test 2: ortho with non-default O ─────────────────────────────────────────


def test_ortho_non_default_O():
    surf = paraboloid(perturb=False)
    proj = Projection(surf, I=[1, 0, 0], J=[0, 1, 0], O=[5, 0, 0])
    np.testing.assert_allclose(proj.XY(np.array([6.0, 2.0, 3.0])), [1.0, 2.0])


# ── Test 3: persp XY(eye) ────────────────────────────────────────────────────


def test_persp_eye_to_origin():
    surf = paraboloid(perturb=False)
    proj = Projection(
        surf, I=[1, 0, 0], J=[0, 1, 0], O=[0, 0, 5], eye=[0, 0, 5],
    )
    assert proj.mode == "persp"
    np.testing.assert_allclose(proj.XY(np.array([0.0, 0.0, 5.0])), [0.0, 0.0])


# ── Test 4: persp O ≠ eye rejected ───────────────────────────────────────────


def test_persp_O_must_equal_eye():
    surf = paraboloid(perturb=False)
    with pytest.raises(ValueError):
        Projection(
            surf, I=[1, 0, 0], J=[0, 1, 0], O=[0, 0, 0], eye=[0, 0, 5],
        )


# ── Test 5: per_vertex_viewer_dot vectorized vs loop ─────────────────────────


class _FakeMesh:
    def __init__(self, SN, xyz=None):
        self.SN = SN
        self.xyz = xyz


def test_per_vertex_viewer_dot_matches_loop():
    surf = paraboloid(perturb=False)

    rng = np.random.default_rng(0)
    N = 50
    SN = rng.standard_normal((N, 3))
    xyz = rng.standard_normal((N, 3))

    # Ortho: dot with the constant (no-arg) viewer direction.
    proj = Projection(surf, I=[1, 0, 0], J=[0, 1, 0])
    mesh = _FakeMesh(SN=SN, xyz=xyz)
    out = proj.per_vertex_viewer_dot(mesh)
    assert out.shape == (N,)
    vd_const = proj.viewer_direction()
    expected = np.array([float(np.dot(SN[i], vd_const)) for i in range(N)])
    np.testing.assert_allclose(out, expected, atol=1e-12)

    # Persp: dot with per-vertex (xyz - eye), unnormalized.
    eye = np.array([0.0, 0.0, 5.0])
    proj_p = Projection(
        surf, I=[1, 0, 0], J=[0, 1, 0], O=eye.tolist(), eye=eye.tolist(),
    )
    out_p = proj_p.per_vertex_viewer_dot(mesh)
    assert out_p.shape == (N,)
    expected_p = np.array([
        float(np.dot(SN[i], proj_p.viewer_direction(xyz[i]))) for i in range(N)
    ])
    np.testing.assert_allclose(out_p, expected_p, atol=1e-12)

    # Persp without xyz raises (the design point of dropping the public `axis`).
    with pytest.raises(ValueError):
        proj_p.viewer_direction()


# ── Persp ker_param / proj_vec (P3 follow-up sub-step) ───────────────────────


def _paraboloid_no_perturb():
    """Z = u² + v² over [-1, 1]² without perturbation, for clean analytical checks."""
    return SurfaceParams(
        "u", "v", "u**2 + v**2", "u v",
        Domain(type="rect", bounds=(-1.0, 1.0, -1.0, 1.0)),
        perturb=False,
    )


def test_persp_ker_param_converges_to_ortho():
    """As eye recedes along axis, the persp Jacobian → ortho Jacobian (up to 1/z),
    so right-singular vectors agree (up to sign).

    Use helicoid: paraboloid would give a degenerate SVD (J = identity,
    σ_1 = σ_2, kernel direction arbitrary).
    """
    surf = helicoid(perturb=False)
    p = np.array([0.5, math.pi / 4])

    proj_o = Projection(surf, I=[1, 0, 0], J=[0, 1, 0])
    kp_o = proj_o.ker_param(p)

    last = None
    for h in (50.0, 500.0, 5000.0):
        eye = [0.0, 0.0, h]
        proj_p = Projection(surf, I=[1, 0, 0], J=[0, 1, 0], O=eye, eye=eye)
        kp_p = proj_p.ker_param(p)
        # SVD sign-ambiguous → compare |cos angle| → 1 monotonically with h.
        cos_abs = abs(float(np.dot(kp_o, kp_p)))
        assert cos_abs > 0.99, f"h={h}: |cos|={cos_abs}"
        if last is not None:
            assert (1.0 - cos_abs) < (1.0 - last), \
                f"convergence not monotone: h={h} |cos|={cos_abs} prev={last}"
        last = cos_abs


def test_persp_ker_param_at_paraboloid_contour():
    """Paraboloid Z=u²+v² with eye=(0,0,-h): the persp contour is the circle
    u²+v² = h. At any contour point p, kerdS should be parallel to S(p) - eye
    (the viewing ray direction)."""
    surf = _paraboloid_no_perturb()
    h = 0.25
    eye = [0.0, 0.0, -h]
    proj = Projection(surf, I=[1, 0, 0], J=[0, 1, 0], O=eye, eye=eye)

    # On the contour circle r = sqrt(h) = 0.5; pick p = (0.5, 0).
    p = np.array([0.5, 0.0])

    # Sanity: this point really is on the persp contour (SN ⊥ (S - eye)).
    S_p = surf.S(p[0], p[1])
    SN_p = surf.SN(p[0], p[1])
    d = S_p - np.array(eye)
    assert abs(float(SN_p @ d)) < 1e-12

    # kerdS should be parallel to d.
    w = proj.kerdS(p)
    cos_abs = abs(float(np.dot(w, d))) / (np.linalg.norm(w) * np.linalg.norm(d))
    assert cos_abs == pytest.approx(1.0, abs=1e-9)


def test_persp_proj_vec_finite_difference():
    """proj_vec is the linearization of XY at S(uv) applied to dir3d, so
    (XY(S(uv) + ε·dir3d) − XY(S(uv))) / ε → proj_vec(uv, dir3d) as ε → 0."""
    surf = _paraboloid_no_perturb()
    eye = [0.0, 0.0, 5.0]
    proj = Projection(surf, I=[1, 0, 0], J=[0, 1, 0], O=eye, eye=eye)

    uv = np.array([0.3, -0.2])
    rng = np.random.default_rng(1)
    for _ in range(4):
        dir3d = rng.standard_normal(3)
        analytical = proj.proj_vec(uv, dir3d)

        S_p = np.asarray(surf.S(uv[0], uv[1]), dtype=float).reshape(3)
        eps = 1e-6
        numerical = (proj.XY(S_p + eps * dir3d) - proj.XY(S_p)) / eps
        np.testing.assert_allclose(analytical, numerical, atol=1e-5, rtol=1e-5)
