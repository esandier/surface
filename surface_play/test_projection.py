"""P3 tests: Projection (ortho/perspective).

Criteria from Modular_rewrite_roadmap.md §P3:
1. Ortho I=[1,0,0], J=[0,1,0], O=None: XY([1,2,3]) == (1, 2); Z([1,2,3]) == 3.
2. Ortho with non-default O=[5,0,0]: XY([6,2,3]) == (1, 2).
3. Persp eye=[0,0,5], O=[0,0,5]: XY([0,0,5]) == (0, 0).
4. Persp O != eye → ValueError.
5. per_vertex_viewer_dot shape (N,), matches per-vertex loop on a 50-vertex mesh.
"""

import numpy as np
import pytest

from surface_play.projection import Projection
from surface_play.test_fixtures import paraboloid


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
