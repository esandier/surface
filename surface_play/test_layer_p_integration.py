"""Layer P integration smoke test.

End-to-end check that the Layer P primitives compose. On a helicoid fixture,
build Domain/SurfaceParams/Projection, generate a 2-triangle manual mesh,
compute SN at vertices, then exercise sign_changes / make_lines /
sweep_segments. Asserts shapes, dtypes, and basic well-formedness — not
algorithmic correctness (covered by per-step tests).

Spec: Modular_rewrite_roadmap.md §"Layer P integration".
"""

import numpy as np
import pytest

from surface_play.curves import make_lines, sign_changes
from surface_play.domain import Domain
from surface_play.intersections import intersect_dtype, sweep_segments
from surface_play.projection import Projection
from surface_play.surface import SurfaceParams


@pytest.fixture
def helicoid_setup():
    """Helicoid S(u,v) = (v·cos u, v·sin u, u) on a rect domain, ortho projection."""
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.5, 1.5))
    surface = SurfaceParams(
        X="v*cos(u)", Y="v*sin(u)", Z="u",
        parameter_names="u v",
        domain=domain,
        perturb=False,
    )
    projection = Projection(surface, I=[1.0, 0.0, 0.0], J=[0.0, 1.0, 0.0])
    return domain, surface, projection


def test_layer_p_compose(helicoid_setup):
    domain, surface, projection = helicoid_setup

    # 2-triangle manual mesh: a square in (u, v) split along the 0-2 diagonal.
    uv = np.array([
        [0.1, 0.6],
        [0.9, 0.6],
        [0.9, 1.4],
        [0.1, 1.4],
    ], dtype=float)
    tris = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.intp)
    assert uv.shape == (4, 2)
    assert tris.shape == (2, 3)

    # SN at the four vertices.
    SN = np.stack([
        np.asarray(surface.SN(float(u), float(v)), dtype=float).reshape(3)
        for u, v in uv
    ])
    assert SN.shape == (4, 3)
    assert np.all(np.isfinite(SN))

    # Projection sanity: viewer_direction in ortho is the (3,) axis.
    view = projection.viewer_direction()
    assert view.shape == (3,)
    assert np.isclose(np.linalg.norm(view), 1.0)

    # Edges of the 2 triangles (unique, undirected).
    edges = np.array([
        [0, 1],
        [1, 2],
        [0, 2],   # shared diagonal
        [2, 3],
        [0, 3],
    ], dtype=np.intp)

    # sign_changes on a synthetic per-edge value array. Use the first two edges
    # as positive→negative crossings, the rest as same-sign (no crossing).
    vals_p = np.array([+1.0, +2.0, +1.0, +1.0, +1.0])
    vals_q = np.array([-1.0, -3.0, +1.0, +1.0, +1.0])
    mask = sign_changes(vals_p, vals_q)
    assert mask.dtype == bool
    assert mask.shape == (5,)
    assert mask.tolist() == [True, True, False, False, False]

    # With per-edge flip ∈ {-1,+1} (Möbius-style): flipping the first sign
    # cancels its crossing, leaving only edge 1.
    flip = np.array([-1, +1, +1, +1, +1])
    mask_f = sign_changes(vals_p, vals_q, flip)
    assert mask_f.tolist() == [False, True, False, False, False]

    # make_lines on the manually-constructed edge set. Topology: vertex 0 has
    # degree 3 (edges 0, 2, 4), so chains may stop at the branch point — we
    # only check that the call returns a list of int arrays without error.
    chains = make_lines(edges)
    assert isinstance(chains, list)
    for c in chains:
        assert isinstance(c, np.ndarray)
        assert c.dtype == np.intp
        assert c.ndim == 1
        for s in c:
            assert 0 < abs(int(s)) <= len(edges)

    # sweep_segments on two crossing edge sets: a horizontal pair vs. a vertical
    # pair, each pair containing one segment that crosses one in the other set.
    seg_a0 = np.array([[0.0, 0.5], [0.0, 0.0]])     # one crosses, one does not
    seg_a1 = np.array([[1.0, 0.5], [1.0, 0.0]])
    seg_b0 = np.array([[0.5, 0.0], [2.0, 0.0]])
    seg_b1 = np.array([[0.5, 1.0], [2.0, 1.0]])
    hits = sweep_segments(seg_a0, seg_a1, seg_b0, seg_b1, domain=None)
    assert hits.dtype == intersect_dtype
    assert hits.ndim == 1
    # Exactly one true crossing: a[0] = (0,0.5)-(1,0.5) crosses b[0] = (0.5,0)-(0.5,1).
    assert len(hits) == 1
    assert int(hits["a"][0]) == 0
    assert int(hits["b"][0]) == 0
    assert np.allclose(hits["uv"][0], [0.5, 0.5])

    # Empty self-sweep returns the empty structured array.
    empty = sweep_segments(
        np.empty((0, 2)), np.empty((0, 2)), None, None,
        domain=None, self_sweep=True,
    )
    assert empty.dtype == intersect_dtype
    assert len(empty) == 0
