"""W1 — pipeline.py orchestrator tests.

7 assertions per roadmap lines 1601-1608:
  1. build_surface_init returns a SurfaceInit whose pieces match standalone calls.
  2. Cache hit: same args ⇒ same Python object (is-equality).
  3. Cache miss on resolution: different resolution ⇒ different object.
  4. mesh_to_threejs payload is JSON-serializable; vertex/face counts match Mesh.
  5. build_outline on the Layer-O-verified fixture set × random axes —
     non-empty `lines_by_visibility`, vis ≤ 0 on every key, origin = (0, 0)
     when O = world origin.
  6. Debug-knob plumbing: propagation="LP1" hits lp_refine_visibility;
     project_resampled=True reaches O14 (resample_all).
  7. No silhouette import: subprocess check that the legacy module is not
     pulled in transitively.
"""

import json
import math
import subprocess
import sys
import textwrap
from types import SimpleNamespace

import numpy as np
import pytest

from surface_play import pipeline
from surface_play.construction import build_construction
from surface_play.domain import Domain
from surface_play.surface import SurfaceParams


TWO_PI = 2.0 * math.pi


# ── Fixture records — same surfaces as test_layer_o_integration.py ───────────
#
# Each entry is the kwargs of a SimpleNamespace that quacks like a
# SurfaceRecord. They mirror the Domain + SurfaceParams arguments of the
# matching test_fixtures.py factories.
FIXTURE_RECORDS = {
    "paraboloid": dict(
        X="u", Y="v", Z="u**2 + v**2",
        parameter_names="u v",
        u_min=-1.0, u_max=1.0, v_min=-1.0, v_max=1.0,
        u_identify="no", v_identify="no",
        domain_type="rect", coord_type="ca",
        r_min=0.0, r_max=1.0, output_type="ca",
    ),
    "torus": dict(
        X="(2 + cos(v)) * cos(u)",
        Y="(2 + cos(v)) * sin(u)",
        Z="sin(v)",
        parameter_names="u v",
        u_min=0.0, u_max=TWO_PI, v_min=0.0, v_max=TWO_PI,
        u_identify="cy", v_identify="cy",
        domain_type="rect", coord_type="ca",
        r_min=0.0, r_max=1.0, output_type="ca",
    ),
    "mobius_u": dict(
        X="(2 + v*cos(u/2)) * cos(u)",
        Y="(2 + v*cos(u/2)) * sin(u)",
        Z="v * sin(u/2)",
        parameter_names="u v",
        u_min=0.0, u_max=TWO_PI, v_min=-0.3, v_max=0.3,
        u_identify="mo", v_identify="no",
        domain_type="rect", coord_type="ca",
        r_min=0.0, r_max=1.0, output_type="ca",
    ),
    "mobius_v": dict(
        X="(2 + u*cos(v/2)) * cos(v)",
        Y="(2 + u*cos(v/2)) * sin(v)",
        Z="u * sin(v/2)",
        parameter_names="u v",
        u_min=-0.3, u_max=0.3, v_min=0.0, v_max=TWO_PI,
        u_identify="no", v_identify="mo",
        domain_type="rect", coord_type="ca",
        r_min=0.0, r_max=1.0, output_type="ca",
    ),
    "disk_para_po": dict(
        X="r*cos(theta)", Y="r*sin(theta)", Z="r**2",
        parameter_names="r theta",
        u_min=0.0, u_max=1.0, v_min=0.0, v_max=TWO_PI,
        u_identify="no", v_identify="no",
        domain_type="disk", coord_type="po",
        r_min=0.0, r_max=1.0, output_type="ca",
    ),
    "disk_para_ca": dict(
        X="x", Y="y", Z="x**2 + y**2",
        parameter_names="x y",
        u_min=0.0, u_max=1.0, v_min=0.0, v_max=TWO_PI,
        u_identify="no", v_identify="no",
        domain_type="disk", coord_type="ca",
        r_min=0.0, r_max=1.0, output_type="ca",
    ),
}


def _make_record(fixture: str = "paraboloid", *, pk: int = 4242,
                 **overrides) -> SimpleNamespace:
    fields = dict(FIXTURE_RECORDS[fixture])
    fields.update(overrides)
    fields["pk"] = pk
    fields["name"] = f"{fixture}_pipeline_test"
    return SimpleNamespace(**fields)


def _random_axis(rng: np.random.Generator):
    """Pick a random viewing axis + an orthonormal in-plane (I, J) frame.

    Identical convention to test_layer_o_integration._random_axis.
    """
    a = rng.standard_normal(3)
    a /= np.linalg.norm(a)
    helper = np.array([1.0, 0.0, 0.0])
    if abs(a @ helper) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    I = helper - (helper @ a) * a
    I /= np.linalg.norm(I)
    J = np.cross(a, I)
    return a, I, J


@pytest.fixture(autouse=True)
def _clear_cache():
    """Each test gets a fresh LRU cache so is-equality and miss tests don't
    spuriously hit leftovers from earlier tests in the module."""
    pipeline._cached_surface_init.cache_clear()
    pipeline._cached_mesh.cache_clear()
    yield
    pipeline._cached_surface_init.cache_clear()
    pipeline._cached_mesh.cache_clear()


# ── 1. build_surface_init parity with standalone calls ───────────────────────


def test_build_surface_init_matches_standalone():
    rec = _make_record("paraboloid")
    init = pipeline.build_surface_init(rec, resolution=18, jitter=False, seed=7)

    domain = Domain(
        type="rect",
        bounds=(rec.u_min, rec.u_max, rec.v_min, rec.v_max),
        u_identify=rec.u_identify, v_identify=rec.v_identify,
        coord_type=rec.coord_type,
    )
    surf = SurfaceParams(
        X=rec.X, Y=rec.Y, Z=rec.Z,
        parameter_names=rec.parameter_names,
        domain=domain, output_type=rec.output_type,
    )
    standalone = build_construction(domain, surf, 18, jitter=False, seed=7)

    assert init.domain.type == "rect"
    assert init.surface.bbox_diag == pytest.approx(surf.bbox_diag, rel=1e-10)
    assert init.construction.mesh.uv.shape == standalone.mesh.uv.shape
    assert init.construction.mesh.tris.shape == standalone.mesh.tris.shape
    np.testing.assert_allclose(
        init.construction.mesh.xyz, standalone.mesh.xyz, atol=1e-12,
    )


# ── 1b. build_mesh_init: mesh-only canvas path ───────────────────────────────


def test_build_mesh_init_no_construction_and_shares_mesh():
    """build_mesh_init returns mesh + threejs without running construction, and
    the mesh object is reused (is-identity) by the later full build, so the
    self-intersection-only cost is all that's paid on the first outline."""
    rec = _make_record("paraboloid")
    mi = pipeline.build_mesh_init(rec, resolution=16, seed=5)

    # No self-intersection artifacts on the mesh-only init.
    assert isinstance(mi, pipeline.MeshInit)
    assert not hasattr(mi, "construction")
    # threejs payload shape parity with the mesh.
    assert len(mi.threejs["vertices"]) == len(mi.mesh.xyz)
    assert len(mi.threejs["faces"]) == len(mi.mesh.tris)

    # The full init reuses the exact same Mesh object (shared _cached_mesh).
    si = pipeline.build_surface_init(rec, resolution=16, seed=5)
    assert mi.mesh is si.construction.mesh
    assert mi.cache_key == si.cache_key
    assert mi.threejs == si.threejs


def test_build_mesh_init_cache_hit_same_object():
    rec = _make_record("paraboloid")
    a = pipeline.build_mesh_init(rec, resolution=14, seed=1)
    b = pipeline.build_mesh_init(rec, resolution=14, seed=1)
    assert a.mesh is b.mesh


# ── 2. Cache hit returns the same Python object ──────────────────────────────


def test_cache_hit_same_object():
    rec = _make_record("paraboloid")
    a = pipeline.build_surface_init(rec, resolution=14, seed=1)
    b = pipeline.build_surface_init(rec, resolution=14, seed=1)
    assert a is b


# ── 3. Cache miss on resolution change ───────────────────────────────────────


def test_cache_miss_on_resolution():
    rec = _make_record("paraboloid")
    a = pipeline.build_surface_init(rec, resolution=14, seed=1)
    c = pipeline.build_surface_init(rec, resolution=20, seed=1)
    assert a is not c
    assert a.cache_key != c.cache_key


# ── 4. mesh_to_threejs JSON round-trip + count parity ────────────────────────


def test_mesh_to_threejs_json_roundtrip():
    rec = _make_record("paraboloid")
    init = pipeline.build_surface_init(rec, resolution=12, seed=2)
    payload = init.threejs

    encoded = json.dumps(payload)
    decoded = json.loads(encoded)

    assert set(decoded.keys()) == {"vertices", "faces", "normals", "uvs"}
    assert len(decoded["vertices"]) == init.construction.mesh.uv.shape[0]
    assert len(decoded["uvs"]) == init.construction.mesh.uv.shape[0]
    assert len(decoded["normals"]) == init.construction.mesh.uv.shape[0]
    assert len(decoded["faces"]) == init.construction.mesh.tris.shape[0]


# ── 5. build_outline parametrized over (fixture × random axis) ───────────────
#
# Modeled on test_layer_o_integration.py: same fixture set, same
# random-axis seeding, same `vis ≤ 0` invariant. Pipeline test runs at
# resolution=60 with 2 trials per fixture to keep the suite responsive
# (integration runs at res=100 × 4 trials).

_OUTLINE_RESOLUTION = 60
_OUTLINE_N_TRIALS = 2
_OUTLINE_SEED = 2026


def _outline_cases():
    rng = np.random.default_rng(_OUTLINE_SEED)
    out = []
    for name in FIXTURE_RECORDS:
        for t in range(_OUTLINE_N_TRIALS):
            _, I, J = _random_axis(rng)
            out.append((name, t, I.tolist(), J.tolist()))
    return out


@pytest.mark.parametrize(
    "fixture,trial,I,J",
    _outline_cases(),
    ids=[f"{c[0]}_t{c[1]}" for c in _outline_cases()],
)
def test_build_outline_random_axes(fixture, trial, I, J):
    """build_outline on every (Layer-O-verified fixture × random axis) trial.

    Asserts: non-empty lines, all vis keys ≤ 0, origin = projection.XY(world
    origin) with O = (0,0,0) (so origin should be (0, 0) in ortho mode).
    """
    rec = _make_record(fixture)
    init = pipeline.build_surface_init(
        rec, resolution=_OUTLINE_RESOLUTION, seed=42,
    )
    result = pipeline.build_outline(
        init, I=I, J=J, O=[0.0, 0.0, 0.0], eye=None,
    )

    assert isinstance(result, pipeline.OutlineResult)
    total = sum(len(v) for v in result.lines_by_visibility.values()) + \
            sum(len(v) for v in result.si_lines_by_visibility.values())
    assert total > 0, f"{fixture} t{trial}: no outline polylines produced"

    for bucket_name, bucket in (
        ("lines", result.lines_by_visibility),
        ("si_lines", result.si_lines_by_visibility),
    ):
        for vis_key, polys in bucket.items():
            assert isinstance(vis_key, int)
            assert vis_key <= 0, (
                f"{fixture} t{trial}: {bucket_name} has vis={vis_key} > 0 "
                f"(I={I}, J={J})"
            )
            for poly in polys:
                assert len(poly) >= 2

    # Ortho with O = world origin → projection.XY(0) = (0, 0).
    assert result.origin == (0.0, 0.0)


# ── 6. Debug-knob plumbing — LP1 + project_resampled ─────────────────────────


def test_debug_knob_plumbing_lp1_and_project_resampled(monkeypatch):
    rec = _make_record("paraboloid")
    init = pipeline.build_surface_init(rec, resolution=14, seed=3)

    import surface_play.pipeline as pipe_mod

    lp_calls = []
    original_lp = pipe_mod.lp_refine_visibility

    def lp_spy(rcs, breaks, splits, vis_bfs=None, **kw):
        lp_calls.append(True)
        return original_lp(rcs, breaks, splits, vis_bfs=vis_bfs, **kw)

    monkeypatch.setattr(pipe_mod, "lp_refine_visibility", lp_spy)

    seen_project_resampled = []
    original_resample = pipe_mod.resample_all

    def resample_spy(*args, project_resampled=None, **kw):
        seen_project_resampled.append(project_resampled)
        return original_resample(*args, project_resampled=project_resampled, **kw)

    monkeypatch.setattr(pipe_mod, "resample_all", resample_spy)

    pipeline.build_outline(
        init, I=[1.0, 0.0, 0.0], J=[0.0, 1.0, 0.0],
        O=[0.0, 0.0, 0.0], eye=None,
        propagation="LP1",
        project_resampled=True,
    )

    assert len(lp_calls) == 1, "LP refinement was not invoked under propagation='LP1'"
    assert seen_project_resampled == [True], (
        f"project_resampled did not propagate to resample_all: {seen_project_resampled}"
    )


# ── 7. No silhouette import in pipeline ──────────────────────────────────────


def test_no_silhouette_import():
    code = textwrap.dedent("""
        import sys
        import surface_play.pipeline  # noqa: F401
        leaked = [m for m in sys.modules
                  if m.endswith('surface_play.silhouette')
                  or m.endswith('surface_play.silhouette_cleaned_by_gemini')]
        assert not leaked, leaked
    """)
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, (
        f"silhouette leaked into sys.modules after importing pipeline\n"
        f"stdout: {result.stdout}\nstderr: {result.stderr}"
    )


# ── 8. Repeat build_outline on same init must not accumulate split slots ─────


def test_repeat_build_outline_on_same_init():
    """`mesh.edges` and `sis_pairs` live in the LRU-cached ConstructionResult
    and their split1/split2 fields are mutated by `_run_outline_pipeline`.
    Without a per-run reset, a second call sees stale slot indices and
    eventually overflows G17 (`SplitSlotOverflowError`). Regression for the
    Sinus pk=13 / side-view 500 hit during W3 smoke-testing.
    """
    rec = _make_record("paraboloid")
    init = pipeline.build_surface_init(rec, resolution=20, seed=11)
    I, J, O = [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]
    first = pipeline.build_outline(init, I=I, J=J, O=O, eye=None)
    second = pipeline.build_outline(init, I=I, J=J, O=O, eye=None)
    # Same input → same vis topology; chiefly we're asserting the second
    # call returned at all (would have raised SplitSlotOverflowError pre-fix).
    assert sorted(first.lines_by_visibility.keys()) == \
           sorted(second.lines_by_visibility.keys())
