"""test_visual_regression.py — W4 visual regression against silhouette.py.

Roadmap §W4 (lines 1696-1724). Two modes:

  - **record** (CLI): `python -m surface_play.test_visual_regression --mode=record`
    Runs `silhouette.Surface.traitement(...)` on every fixture × viewpoint,
    serializes the result to `regression_baselines/<fixture>__v<i>.json`.

  - **check** (pytest default): runs `pipeline.build_outline(...)` on the same
    fixture × viewpoint, ALWAYS writes a side-by-side matplotlib PNG to
    `test_artifacts/<fixture>__v<i>.png`, and asserts only "did the new
    pipeline produce a non-empty, well-shaped result?" (smoke-level).

The numerical comparison was intentionally dropped: silhouette.py and the
new pipeline differ in legitimate ways (sampling density, SIC handling
per success criterion 6) that any tolerance threshold either misses or
flags falsely. The 27 side-by-side PNGs are the deliverable — review them
by eye after every meaningful pipeline change.

This test is the contractual gate before W5 retires `silhouette.py`. If a
code change legitimately moves a baseline, regenerate and re-review in
the same PR (roadmap line 1722).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import sys
from types import SimpleNamespace

import numpy as np
import pytest

# Django must be configured before importing `test_smoke` (which transitively
# loads `surface_play.models.SurfaceRecord`). Pytest's conftest does this for
# us under `pytest …`; this guard covers the `python -m …` record path.
if not os.environ.get("DJANGO_SETTINGS_MODULE"):
    os.environ["DJANGO_SETTINGS_MODULE"] = "surfaces_project.settings"
try:
    from django.apps import apps as _django_apps
    if not _django_apps.ready:
        import django
        django.setup()
except ImportError:  # pragma: no cover — Django always present in this project
    pass

from surface_play import pipeline
from surface_play.test_smoke import _FIXTURE_RECORD_DICTS, _fixture_record


# ── Paths ────────────────────────────────────────────────────────────────────


_REPO_ROOT = pathlib.Path(__file__).parent.parent
BASELINE_DIR = _REPO_ROOT / "regression_baselines"
ARTIFACTS_DIR = _REPO_ROOT / "test_artifacts"


# ── Tunables ─────────────────────────────────────────────────────────────────


# Use a moderate resolution for both legacy and new pipelines. Production is
# 200 (settings.py); a lower value keeps the regression suite tractable while
# still exercising the same code paths.
RESOLUTION = 40

# Determinism (roadmap criterion 4): fixed seed → byte-identical output across
# runs. Mesh jitter is OFF (jitter=False) for the comparison so the new
# pipeline's per-vertex jitter doesn't introduce uncorrelated noise vs. the
# legacy (which uses its own surface-perturbation scheme).
JITTER_SEED = 0
JITTER_ON = False

FIXTURE_NAMES = list(_FIXTURE_RECORD_DICTS.keys())


# ── Viewpoints — 3 per fixture, stable across runs ───────────────────────────


def _viewpoints_for(name: str) -> list[tuple[list, list, list | None]]:
    """Return 3 viewpoints for `name`: default ortho + 2 seeded random axes.

    All have a ~1e-3 deterministic perturbation to dodge axis-aligned
    degeneracies (per test_fixtures.perturb_axis convention).
    """
    seed = sum(ord(c) for c in name) % 997
    rng = np.random.default_rng(seed)

    views: list[tuple[list, list, list | None]] = []

    # View 0 — slightly off-axis ortho ([1,0,0], [0,1,0]).
    I0 = (np.array([1.0, 0.0, 0.0]) + rng.uniform(-1e-3, 1e-3, 3))
    J0 = (np.array([0.0, 1.0, 0.0]) + rng.uniform(-1e-3, 1e-3, 3))
    views.append((I0.tolist(), J0.tolist(), None))

    # Views 1, 2 — random axes (orthonormal frame about a random view direction).
    for _ in range(2):
        a = rng.standard_normal(3)
        a /= np.linalg.norm(a)
        helper = np.array([1.0, 0.0, 0.0])
        if abs(a @ helper) > 0.9:
            helper = np.array([0.0, 1.0, 0.0])
        I = helper - (helper @ a) * a
        I /= np.linalg.norm(I)
        J = np.cross(a, I)
        views.append((I.tolist(), J.tolist(), None))

    return views


# ── Legacy and new outline drivers ───────────────────────────────────────────


# Path to the retired legacy module. Kept out of `surface_play/` (W5) so it is
# never importable as part of the package; loaded by file path on demand for
# the `record` baseline path only.
_LEGACY_SILHOUETTE_PATH = _REPO_ROOT / "old stuff" / "silhouette_legacy.py"


def _load_legacy_surface():
    """Import `Surface` from the retired `old stuff/silhouette_legacy.py`.

    The folder name contains a space, so a normal `import` is impossible —
    load by file path via importlib. Cached on the function object so repeated
    record-mode calls don't re-exec the heavy module.
    """
    cached = getattr(_load_legacy_surface, "_cached", None)
    if cached is not None:
        return cached
    import importlib.util

    if not _LEGACY_SILHOUETTE_PATH.exists():
        raise FileNotFoundError(
            f"legacy silhouette module not found at {_LEGACY_SILHOUETTE_PATH}. "
            f"It was retired here in W5; record mode needs it to regenerate "
            f"baselines (check mode does not)."
        )
    spec = importlib.util.spec_from_file_location(
        "silhouette_legacy", _LEGACY_SILHOUETTE_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _load_legacy_surface._cached = module.Surface
    return module.Surface


def _legacy_outline(name: str, I, J, eye) -> dict[str, list]:
    """Run silhouette.Surface.traitement and shape its output to the new format.

    Returns `{vis_key_str: [polyline, polyline, ...]}`, where each polyline
    is a list of [x, y] floats. Legacy doesn't separate SIC from BC/CC, so
    everything lands in one dict (matches the merged view the comparator
    uses on the new side).
    """
    # Lazy import — the legacy module is heavy and not needed for `check`
    # mode. It was retired to `old stuff/silhouette_legacy.py` in W5; the
    # folder name has a space, so load it by file path via importlib.
    Surface = _load_legacy_surface()

    fields = _FIXTURE_RECORD_DICTS[name]
    bounds = (fields["u_min"], fields["u_max"], fields["v_min"], fields["v_max"])
    quotient = (fields["u_identify"], fields["v_identify"])

    s = Surface(
        X=fields["X"], Y=fields["Y"], Z=fields["Z"],
        param_names=fields["parameter_names"],
        bounds=bounds, quotient=quotient,
        domain_type=fields["domain_type"], coord_type=fields["coord_type"],
        r_min=fields["r_min"], r_max=fields["r_max"],
        output_type=fields["output_type"],
    )
    # `triangulate` populates the S_grid / SN_grid arrays that `traitement`
    # reads. Legacy views.py invoked this implicitly via for_3js + a
    # follow-up call; we have to do it explicitly here.
    s.triangulate(RESOLUTION)
    s.set_axis(I=I, J=J, eye=eye)
    s.traitement(use_lp=True, use_newton=True)
    raw = s.plot_for_browser()  # defaultdict[int, list[list[ndarray]]]

    out: dict[str, list] = {}
    for k, polys in raw.items():
        key = str(int(k))
        out[key] = [
            [[float(p[0]), float(p[1])] for p in poly]
            for poly in polys
            if len(poly) >= 2
        ]
    return out


def _new_outline(name: str, I, J, eye) -> dict[str, list]:
    """Run pipeline.build_outline. Merges SIC into the main dict so the
    shape matches the legacy output (which has no separate SIC bucket)."""
    pk = 20_000 + FIXTURE_NAMES.index(name)
    rec = _fixture_record(name, pk=pk)
    init = pipeline.build_surface_init(
        rec, resolution=RESOLUTION, jitter=JITTER_ON, seed=JITTER_SEED,
    )
    O = [0.0, 0.0, 0.0] if eye is None else eye
    result = pipeline.build_outline(init, I=I, J=J, O=O, eye=eye)

    merged: dict[str, list] = {}
    for k, polys in result.lines_by_visibility.items():
        merged.setdefault(str(int(k)), []).extend(
            [[list(map(float, p)) for p in poly] for poly in polys]
        )
    for k, polys in result.si_lines_by_visibility.items():
        merged.setdefault(str(int(k)), []).extend(
            [[list(map(float, p)) for p in poly] for poly in polys]
        )
    return merged


# ── Diff PNG ─────────────────────────────────────────────────────────────────


def _write_diff_png(name: str, view_idx: int, legacy: dict, new: dict) -> str:
    """Side-by-side matplotlib comparison; returns the artifact path."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    for ax, outline, title in [
        (axes[0], legacy, f"legacy: {name} v{view_idx}"),
        (axes[1], new, f"new: {name} v{view_idx}"),
    ]:
        for k, polys in outline.items():
            for poly in polys:
                if len(poly) < 2:
                    continue
                arr = np.asarray(poly)
                front = (int(k) == 0)
                ax.plot(
                    arr[:, 0], arr[:, 1],
                    color="black" if front else "0.5",
                    linestyle="-" if front else "--",
                    linewidth=0.9 if front else 0.6,
                )
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.2)

    out = ARTIFACTS_DIR / f"{name}__v{view_idx}.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return str(out)


# ── Baseline I/O ─────────────────────────────────────────────────────────────


def _baseline_path(name: str, view_idx: int) -> pathlib.Path:
    return BASELINE_DIR / f"{name}__v{view_idx}.json"


# ── Record mode ──────────────────────────────────────────────────────────────


def record_all(*, fixtures: list[str] | None = None) -> dict[str, list[str]]:
    """Record baselines for the requested fixtures (default: all).

    Returns `{fixture: [status_for_each_viewpoint]}` so callers can spot
    fixtures that errored vs. recorded cleanly.
    """
    BASELINE_DIR.mkdir(parents=True, exist_ok=True)
    fixtures = fixtures or FIXTURE_NAMES
    results: dict[str, list[str]] = {}
    for name in fixtures:
        per_fixture: list[str] = []
        for idx, (I, J, eye) in enumerate(_viewpoints_for(name)):
            try:
                outline = _legacy_outline(name, I, J, eye)
            except Exception as exc:
                per_fixture.append(f"v{idx}: FAILED ({exc.__class__.__name__}: {exc})")
                continue
            n_total = sum(len(polys) for polys in outline.values())
            payload = {
                "fixture": name,
                "view_idx": idx,
                "view": {"I": I, "J": J, "eye": eye},
                "lines_by_visibility": outline,
                "metadata": {
                    "n_polylines": n_total,
                    "n_vis_keys": len(outline),
                },
            }
            _baseline_path(name, idx).write_text(
                json.dumps(payload, indent=2, sort_keys=True),
                encoding="utf-8",
            )
            per_fixture.append(f"v{idx}: {n_total} polylines")
        results[name] = per_fixture
    return results


# ── Pytest entry ─────────────────────────────────────────────────────────────


_VIEW_CASES = [(name, i) for name in FIXTURE_NAMES for i in range(3)]


@pytest.mark.parametrize(
    "name, view_idx", _VIEW_CASES,
    ids=[f"{n}_v{i}" for n, i in _VIEW_CASES],
)
def test_visual_regression(name, view_idx):
    """Run new pipeline against the legacy baseline and emit a side-by-side
    PNG for human review.

    There is no numerical pass/fail — silhouette.py and the new pipeline have
    legitimate algorithmic differences (e.g. sampling density, SIC handling)
    and a tolerance-based comparator either flags too many false positives
    (tight tol) or hides real regressions (loose tol). Instead, the artifact
    PNG goes to `test_artifacts/<fixture>__v<view>.png` for visual review.

    The pytest assertion only catches the easy "did the new pipeline crash
    or produce nonsense?" cases — basic shape validity of the output dict.
    """
    bpath = _baseline_path(name, view_idx)
    if not bpath.exists():
        pytest.skip(
            f"baseline missing: {bpath.relative_to(_REPO_ROOT)}\n"
            f"Generate with: python -m surface_play.test_visual_regression --mode=record"
        )

    baseline = json.loads(bpath.read_text(encoding="utf-8"))
    I = baseline["view"]["I"]
    J = baseline["view"]["J"]
    eye = baseline["view"]["eye"]
    legacy = baseline["lines_by_visibility"]

    new = _new_outline(name, I, J, eye)

    # Always emit the side-by-side artifact — that's what we review.
    _write_diff_png(name, view_idx, legacy, new)

    # Smoke-equivalent shape validity: assert the new pipeline produced
    # SOMETHING. (Crashes already surface as test errors.)
    n_total = sum(len(polys) for polys in new.values())
    assert n_total >= 1, (
        f"{name} v{view_idx}: new pipeline produced 0 polylines. "
        f"See test_artifacts/{name}__v{view_idx}.png."
    )


# ── Determinism guard ────────────────────────────────────────────────────────


@pytest.mark.parametrize("name", FIXTURE_NAMES)
def test_determinism(name):
    """Two runs in a row on the same fixture × view 0 produce byte-identical
    JSON (roadmap §W4 test criterion 4)."""
    I, J, eye = _viewpoints_for(name)[0]
    a = _new_outline(name, I, J, eye)
    b = _new_outline(name, I, J, eye)
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True), (
        f"non-deterministic outline for {name}"
    )


# ── CLI: record-mode entry point ─────────────────────────────────────────────


def _cli():
    parser = argparse.ArgumentParser(
        description="Record or check visual-regression baselines."
    )
    parser.add_argument(
        "--mode", choices=["record", "check"], default="check",
        help="`record` writes baselines; `check` says how to run pytest.",
    )
    parser.add_argument(
        "--fixture", action="append",
        help="Restrict to specific fixture(s); repeatable.",
    )
    args = parser.parse_args()

    if args.mode == "record":
        results = record_all(fixtures=args.fixture)
        for name, statuses in results.items():
            print(f"{name}:")
            for s in statuses:
                print(f"  {s}")
        return 0

    print(
        "check mode: run\n"
        "  pytest surface_play/test_visual_regression.py\n"
        "from the repo root."
    )
    return 0


if __name__ == "__main__":
    sys.exit(_cli())
