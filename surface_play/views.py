"""views.py — W2: Django GET / POST endpoints backed by pipeline.

GET /surfaceplay/<pk>/  — renders templates/play.html with three.js mesh data.
POST /surfaceplay/<pk>/ — returns
    {lines_by_visibility, si_lines_by_visibility, origin}
for the supplied (I, J, O, eye) and optional debug knobs.

No silhouette imports (W2 sign-off gate, test 10).
"""

from __future__ import annotations

import json
import logging
from typing import Any

import numpy as np

from django.http import (
    HttpResponse,
    HttpResponseBadRequest,
    HttpResponseNotAllowed,
    JsonResponse,
)
from django.shortcuts import get_object_or_404, render
from django.urls import reverse_lazy
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import ListView
from django.views.generic.edit import CreateView, DeleteView, UpdateView

from surface_play import pipeline
from surface_play.forms import SurfaceRecordForm
from surface_play.models import SurfaceRecord


logger = logging.getLogger(__name__)


# ── Debug-dict → build_outline kwargs ────────────────────────────────────────


_DEBUG_KEY_MAP = {
    "PROPAGATION": "propagation",
    "RESOLUTION": "canvas_resolution",
    "NEWTON_CUSP": "newton_cusp",
    "PROJECT_RESAMPLED": "project_resampled",
}


def _debug_kwargs(debug: dict[str, Any] | None) -> dict[str, Any]:
    """Validate `debug` dict and map known keys to build_outline kwargs.

    Unknown keys are logged at WARNING and ignored (preserve client compat).
    """
    if not debug:
        return {}
    out: dict[str, Any] = {}
    for k, v in debug.items():
        target = _DEBUG_KEY_MAP.get(k)
        if target is None:
            logger.warning("play_outline: ignoring unknown debug key %r", k)
            continue
        out[target] = v
    return out


# ── GET / POST endpoints ─────────────────────────────────────────────────────


def _legacy_threejs_context(threejs: dict, mesh_xyz: np.ndarray) -> dict:
    """Convert structured threejs payload back to the flat shape the current
    templates/play.html expects (`positions`, `normals`, `faces`, `center`,
    `radius`). Future template work can read `threejs` directly instead.
    """
    vertices = np.asarray(threejs["vertices"], dtype=float)
    normals = np.asarray(threejs["normals"], dtype=float)
    faces = np.asarray(threejs["faces"], dtype=np.int64)

    center = mesh_xyz.mean(axis=0)
    radius = float(np.linalg.norm(mesh_xyz - center, axis=1).max())

    return {
        "positions": vertices.ravel().tolist(),
        "normals": normals.ravel().tolist(),
        "faces": faces.ravel().tolist(),
        "center": center.tolist(),
        "radius": radius,
    }


def _play_get(request, record: SurfaceRecord) -> HttpResponse:
    init = pipeline.build_surface_init(record)
    ctx = _legacy_threejs_context(init.threejs, init.construction.mesh.xyz)
    ctx.update({
        "threejs": init.threejs,
        "surface": record,
        "pk": record.pk,
        "initial_elev": record.initial_elev,
        "initial_azim": record.initial_azim,
        "initial_inplane": record.initial_inplane,
        "initial_zoom": record.initial_zoom,
        "initial_perspective": record.initial_perspective,
        "debug_ui": True,
    })
    return render(request, "play.html", ctx)


def _vis_dict_to_jsonable(d: dict[int, list]) -> dict[str, list]:
    """JsonResponse can't serialize int keys via the default encoder — convert
    to strings while keeping order by sorted-vis."""
    return {str(int(k)): d[k] for k in sorted(d.keys())}


def _play_post(request, record: SurfaceRecord) -> HttpResponse:
    try:
        data = json.loads(request.body.decode("utf-8") or "{}")
    except json.JSONDecodeError as exc:
        return HttpResponseBadRequest(f"invalid JSON body: {exc}")

    for field in ("I", "J", "O"):
        if field not in data:
            return HttpResponseBadRequest(f"missing required field: {field}")
    I = data["I"]
    J = data["J"]
    O = data["O"]
    eye = data.get("eye")  # None → ortho
    debug = data.get("debug", {}) or {}

    init = pipeline.build_surface_init(record)
    try:
        result = pipeline.build_outline(
            init, I=I, J=J, O=O, eye=eye, **_debug_kwargs(debug),
        )
    except ValueError as exc:
        # P3 raises on persp O != eye; surface as 400.
        return HttpResponseBadRequest(str(exc))

    return JsonResponse({
        "lines_by_visibility": _vis_dict_to_jsonable(result.lines_by_visibility),
        "si_lines_by_visibility": _vis_dict_to_jsonable(result.si_lines_by_visibility),
        "origin": list(result.origin),
    })


@csrf_exempt
def play(request, pk: int) -> HttpResponse:
    """GET → render the play page; POST → return outline JSON."""
    record = get_object_or_404(SurfaceRecord, pk=pk)
    if request.method == "GET":
        return _play_get(request, record)
    if request.method == "POST":
        return _play_post(request, record)
    return HttpResponseNotAllowed(["GET", "POST"])


# ── Save endpoints ───────────────────────────────────────────────────────────


@csrf_exempt
def save_view(request, pk):
    if request.method != "POST":
        return HttpResponse(status=405)
    rec = get_object_or_404(SurfaceRecord, pk=pk)
    data = json.loads(request.body.decode())
    rec.initial_elev = float(data["elev"])
    rec.initial_azim = float(data["azim"])
    rec.initial_inplane = float(data["inplane"])
    rec.initial_zoom = float(data["zoom"])
    rec.initial_perspective = bool(data.get("perspective", False))
    rec.save(update_fields=[
        "initial_elev", "initial_azim", "initial_inplane",
        "initial_zoom", "initial_perspective",
    ])
    return HttpResponse("{}", content_type="application/json")


@csrf_exempt
def save_thumbnail(request, pk):
    if request.method != "POST":
        return HttpResponse(status=405)
    rec = get_object_or_404(SurfaceRecord, pk=pk)
    data = json.loads(request.body.decode())
    I = data.get("I")
    J = data.get("J")
    eye = data.get("eye")
    try:
        from surface_play.thumbnail import render_thumbnail
        rec.thumbnail = render_thumbnail(rec, I=I, J=J, eye=eye)
        rec.save(update_fields=["thumbnail"])
    except Exception as e:
        return HttpResponse(
            json.dumps({"error": str(e)}),
            content_type="application/json",
            status=500,
        )
    return HttpResponse("{}", content_type="application/json")


# ── CRUD views (unchanged) ───────────────────────────────────────────────────


class SurfaceRecordListView(ListView):
    model = SurfaceRecord
    context_object_name = "surfaces"


class SurfaceRecordCreateView(CreateView):
    model = SurfaceRecord
    form_class = SurfaceRecordForm
    context_object_name = "surface"


class SurfaceRecordUpdateView(UpdateView):
    model = SurfaceRecord
    form_class = SurfaceRecordForm
    context_object_name = "surface"


class SurfaceRecordDeleteView(DeleteView):
    model = SurfaceRecord
    success_url = reverse_lazy("surfaces")
    context_object_name = "surface"


# ── Back-compat alias for urls.py (legacy class-based dispatch) ──────────────


class SurfacePlayView:
    """Compatibility shim — the urls.py used to wire `SurfacePlayView.as_view()`.
    The class-based form is gone; expose `as_view()` that returns the function
    view so the existing URL conf keeps working without edits."""

    @staticmethod
    def as_view():
        return play
