"""Library of small sympy surfaces used as fixtures by Layer P/C/O tests.

Per roadmap §3 "Test strategy": paraboloid, torus-cyl, Möbius band, plus
helpers (helicoid, cylinder_cy) used by P2 tests.
"""

import math

from surface_play.domain import Domain
from surface_play.surface import SurfaceParams

TWO_PI = 2 * math.pi


def helicoid(perturb: bool = True) -> SurfaceParams:
    """Helicoid on a non-identified rectangle."""
    domain = Domain(type="rect", bounds=(-1.0, 1.0, -math.pi, math.pi))
    return SurfaceParams("u*cos(v)", "u*sin(v)", "v", "u v", domain, perturb=perturb)


def paraboloid(perturb: bool = True) -> SurfaceParams:
    """Paraboloid Z = u² + v² over the unit square."""
    domain = Domain(type="rect", bounds=(-1.0, 1.0, -1.0, 1.0))
    return SurfaceParams("u", "v", "u**2 + v**2", "u v", domain, perturb=perturb)


def cylinder_cy(perturb: bool = True) -> SurfaceParams:
    """Vertical cylinder of unit radius, height 1; v_identify='cy'."""
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, TWO_PI), v_identify="cy")
    return SurfaceParams("cos(v)", "sin(v)", "u", "u v", domain, perturb=perturb)


def torus(perturb: bool = True) -> SurfaceParams:
    """Standard torus (R=2, r=1); u_identify='cy', v_identify='cy'."""
    domain = Domain(
        type="rect", bounds=(0.0, TWO_PI, 0.0, TWO_PI),
        u_identify="cy", v_identify="cy",
    )
    return SurfaceParams(
        "(2 + cos(v)) * cos(u)",
        "(2 + cos(v)) * sin(u)",
        "sin(v)",
        "u v", domain, perturb=perturb,
    )


def mobius_u(perturb: bool = True) -> SurfaceParams:
    """Möbius band with u_identify='mo' on bounds (0, 2π).

    Identification (0, v) ~ (2π, -v) holds because v_max + v_min = 0.
    """
    domain = Domain(
        type="rect", bounds=(0.0, TWO_PI, -0.3, 0.3), u_identify="mo",
    )
    return SurfaceParams(
        "(2 + v*cos(u/2)) * cos(u)",
        "(2 + v*cos(u/2)) * sin(u)",
        "v * sin(u/2)",
        "u v", domain, perturb=perturb,
    )


def mobius_v(perturb: bool = True) -> SurfaceParams:
    """Möbius band with v_identify='mo' (axes swapped vs mobius_u)."""
    domain = Domain(
        type="rect", bounds=(-0.3, 0.3, 0.0, TWO_PI), v_identify="mo",
    )
    return SurfaceParams(
        "(2 + u*cos(v/2)) * cos(v)",
        "(2 + u*cos(v/2)) * sin(v)",
        "u * sin(v/2)",
        "u v", domain, perturb=perturb,
    )


def fig8(perturb: bool = True) -> SurfaceParams:
    """Fig-8 torus immersion (roadmap line 653); cy-cy on (0, 2π)².

    Self-intersects along a closed SIC. Used as the canonical SIC fixture
    for C8-C12 tests.
    """
    domain = Domain(
        type="rect", bounds=(0.0, TWO_PI, 0.0, TWO_PI),
        u_identify="cy", v_identify="cy",
    )
    return SurfaceParams(
        "(2 + cos(u))*cos(v)",
        "(2 + cos(u))*sin(v)",
        "sin(2*u)",
        "u v", domain, perturb=perturb,
    )


def disk_paraboloid_po(perturb: bool = True) -> SurfaceParams:
    """Paraboloid Z = r² over the unit disk, polar parametrization."""
    domain = Domain(type="disk", bounds=(0.0, 1.0, 0.0, TWO_PI), coord_type="po")
    return SurfaceParams(
        "r*cos(theta)", "r*sin(theta)", "r**2",
        "r theta", domain, perturb=perturb,
    )


def disk_paraboloid_ca(perturb: bool = True) -> SurfaceParams:
    """Paraboloid Z = x²+y² over the unit disk, cartesian parametrization."""
    domain = Domain(type="disk", bounds=(0.0, 1.0, 0.0, TWO_PI), coord_type="ca")
    return SurfaceParams(
        "x", "y", "x**2 + y**2",
        "x y", domain, perturb=perturb,
    )


# ── Canonical viewpoints for Layer O tests ────────────────────────────────────
# Each entry: (I, J, eye).  eye=None → orthographic.

# Helicoid: ortho, axis = I×J = (0,0,1) — top-down Z view; contour at u=0.
helicoid_ortho_view = ([1, 0, 0], [0, 1, 0], None)

# Paraboloid: ortho, axis = I×J = (0,-1,0) — side Y view; contour at v=0.
paraboloid_side_view = ([1, 0, 0], [0, 0, 1], None)

# Torus: ortho, axis = (0,0,1) — Z view; gives a rich silhouette with cusps.
torus_ortho_view = ([1, 0, 0], [0, 1, 0], None)

# Möbius: ortho, axis = (0,0,1); used for seam/flip regression tests.
mobius_ortho_view = ([1, 0, 0], [0, 1, 0], None)

# Fig-8: perspective from Z=5; used for SIC visibility integration tests.
fig8_persp_view = ([1, 0, 0], [0, 1, 0], [0, 0, 5])

# ── Smoke: each fixture builds without error ─────────────────────────────────

import pytest


@pytest.mark.parametrize(
    "factory",
    [helicoid, paraboloid, cylinder_cy, torus, mobius_u, mobius_v, fig8,
     disk_paraboloid_po, disk_paraboloid_ca],
)
def test_fixture_builds(factory):
    s = factory(perturb=True)
    assert s.bbox_diag > 0.0
    assert s.S(0.5, 0.5).shape == (3,)
