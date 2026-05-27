"""W3 — `render_thumbnail` tests (roadmap lines 1683-1690).

Seven assertions cover: valid SVG output, expected viewBox, non-trivial
content, LRU sharing with `pipeline._cached_surface_init`, TextField
round-trip, ortho/persp `O` derivation, and the silhouette-import gate.
"""

from __future__ import annotations

import subprocess
import sys
from unittest.mock import patch

from django.test import TestCase

from surface_play import pipeline
from surface_play.models import SurfaceRecord
from surface_play.thumbnail import render_thumbnail


def _make_paraboloid_record() -> SurfaceRecord:
    return SurfaceRecord.objects.create(
        name="paraboloid_thumbnail_test",
        X="u", Y="v", Z="u**2 + v**2",
        parameter_names="u v",
        u_min=-1.0, u_max=1.0, v_min=-1.0, v_max=1.0,
        u_identify="no", v_identify="no",
        domain_type="rect", coord_type="ca",
        r_min=0.0, r_max=1.0,
        output_type="ca",
    )


class RenderThumbnailTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.record = _make_paraboloid_record()

    def setUp(self):
        pipeline._cached_surface_init.cache_clear()
        # Generic tilted ortho view — paraboloid + axis-aligned hits
        # boundary-edge degeneracies that aren't this test's concern.
        self.I = [1.0, 0.0, 0.3]
        self.J = [0.0, 1.0, 0.4]

    # ── 1. Valid SVG ─────────────────────────────────────────────────────────

    def test_returns_valid_svg_string(self):
        svg = render_thumbnail(self.record, I=self.I, J=self.J)
        self.assertIsInstance(svg, str)
        self.assertTrue(svg.lstrip().startswith("<svg"))
        self.assertTrue(svg.rstrip().endswith("</svg>"))

    # ── 2. viewBox matches our fixed convention ──────────────────────────────

    def test_viewbox_matches_convention(self):
        svg = render_thumbnail(self.record, I=self.I, J=self.J)
        self.assertIn('viewBox="-1.1 -1.1 2.2 2.2"', svg)

    # ── 3. Non-trivial content ───────────────────────────────────────────────

    def test_non_trivial_content(self):
        svg = render_thumbnail(self.record, I=self.I, J=self.J)
        # A real paraboloid view emits well over a dozen polylines —
        # the threshold below catches "empty <svg/>" without being flaky.
        self.assertGreater(len(svg), 400)
        self.assertIn("<polyline", svg)

    # ── 4. LRU shared with views ─────────────────────────────────────────────

    def test_cache_shared_with_build_surface_init(self):
        info_before = pipeline._cached_surface_init.cache_info()
        pipeline.build_surface_init(self.record)
        info_after_warmup = pipeline._cached_surface_init.cache_info()
        self.assertEqual(
            info_after_warmup.misses, info_before.misses + 1,
            "first build_surface_init call should miss",
        )
        render_thumbnail(self.record, I=self.I, J=self.J)
        info_after_thumb = pipeline._cached_surface_init.cache_info()
        self.assertGreaterEqual(
            info_after_thumb.hits, info_after_warmup.hits + 1,
            "render_thumbnail should reuse the cached SurfaceInit",
        )
        self.assertEqual(
            info_after_thumb.misses, info_after_warmup.misses,
            "render_thumbnail must not recompute construction",
        )

    # ── 5. TextField round-trip ──────────────────────────────────────────────

    def test_storage_round_trip(self):
        svg = render_thumbnail(self.record, I=self.I, J=self.J)
        self.record.thumbnail = svg
        self.record.save(update_fields=["thumbnail"])
        reloaded = SurfaceRecord.objects.get(pk=self.record.pk)
        self.assertEqual(reloaded.thumbnail, svg)

    # ── 6. O derivation: ortho → (0,0,0); persp → eye ────────────────────────

    def test_O_derivation_ortho_vs_persp(self):
        # Spy on build_outline and short-circuit it: this test is about the O
        # `render_thumbnail` derives, not about full-pipeline correctness on
        # the chosen view (which has its own integration coverage).
        captured = []
        stub = pipeline.OutlineResult(
            lines_by_visibility={},
            si_lines_by_visibility={},
            origin=(0.0, 0.0),
        )

        def _spy(init, I, J, O, eye, **kw):
            captured.append({"O": list(O), "eye": None if eye is None else list(eye)})
            return stub

        with patch.object(pipeline, "build_outline", side_effect=_spy):
            # Ortho: eye omitted → O should be [0, 0, 0].
            render_thumbnail(self.record, I=self.I, J=self.J)
            self.assertEqual(captured[-1]["O"], [0.0, 0.0, 0.0])
            self.assertIsNone(captured[-1]["eye"])

            # Persp: eye supplied → O must equal eye (P3 invariant).
            eye = [0.5, -0.3, 3.0]
            render_thumbnail(self.record, I=self.I, J=self.J, eye=eye)
            self.assertEqual(captured[-1]["O"], eye)
            self.assertEqual(captured[-1]["eye"], eye)

    # ── 7. Silhouette-import gate ────────────────────────────────────────────

    def test_thumbnail_module_does_not_pull_in_silhouette(self):
        script = (
            "import surface_play.thumbnail, sys; "
            "leaked = [m for m in sys.modules if 'silhouette' in m]; "
            "assert not leaked, leaked"
        )
        result = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True, text=True,
        )
        self.assertEqual(
            result.returncode, 0,
            f"silhouette leaked into sys.modules: "
            f"stdout={result.stdout!r} stderr={result.stderr!r}",
        )
