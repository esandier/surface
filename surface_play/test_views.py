"""W2 — Django GET/POST endpoint tests via the test client.

10 assertions per roadmap lines 1646-1655.
"""

import json
import math
import re
from pathlib import Path

import numpy as np
from django.test import Client, TestCase
from django.urls import reverse

from surface_play import pipeline
from surface_play.models import SurfaceRecord


TWO_PI = 2.0 * math.pi


def _make_paraboloid_record() -> SurfaceRecord:
    return SurfaceRecord.objects.create(
        name="paraboloid_views_test",
        X="u", Y="v", Z="u**2 + v**2",
        parameter_names="u v",
        u_min=-1.0, u_max=1.0, v_min=-1.0, v_max=1.0,
        u_identify="no", v_identify="no",
        domain_type="rect", coord_type="ca",
        r_min=0.0, r_max=1.0,
        output_type="ca",
    )


class PlayEndpointTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.record = _make_paraboloid_record()

    def setUp(self):
        # Fresh LRU on every test so cache-observable tests aren't poisoned.
        pipeline._cached_surface_init.cache_clear()
        self.client = Client()
        self.url = reverse("surface-play", kwargs={"pk": self.record.pk})

    # ── 1. GET ───────────────────────────────────────────────────────────────

    def test_get_renders_three_js_payload(self):
        response = self.client.get(self.url)
        self.assertEqual(response.status_code, 200)
        body = response.content.decode("utf-8")
        # The legacy template consumes positions/normals/faces/center/radius.
        for key in ("positions", "normals", "faces", "center", "radius"):
            self.assertIn(key, body)

    # ── 2. POST ortho minimal ────────────────────────────────────────────────

    def test_post_ortho_returns_outline(self):
        body = {
            "I": [1.0, 0.0, 0.0],
            "J": [0.0, 0.0, 1.0],  # paraboloid_side_view
            "O": [0.0, 0.0, 0.0],
            "eye": None,
            "debug": {},
        }
        response = self.client.post(
            self.url, data=json.dumps(body),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        self.assertEqual(
            set(data.keys()),
            {"lines_by_visibility", "si_lines_by_visibility", "origin"},
        )
        self.assertTrue(len(data["lines_by_visibility"]) > 0)
        # origin = projection.XY(0) with O=(0,0,0) ⇒ (0, 0).
        self.assertEqual(data["origin"], [0.0, 0.0])

    # ── 3. POST persp ────────────────────────────────────────────────────────

    def test_post_persp_returns_outline(self):
        eye = [0.0, 0.0, 5.0]
        body = {
            "I": [1.0, 0.0, 0.0],
            "J": [0.0, 1.0, 0.0],
            "O": eye,  # required: O == eye in persp
            "eye": eye,
            "debug": {},
        }
        response = self.client.post(
            self.url, data=json.dumps(body),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        self.assertIn("origin", data)
        # Persp XY(eye) = (0, 0).
        self.assertEqual(data["origin"], [0.0, 0.0])

    # ── 4. Persp O != eye rejected ───────────────────────────────────────────

    def test_post_persp_rejects_O_not_eye(self):
        body = {
            "I": [1.0, 0.0, 0.0],
            "J": [0.0, 1.0, 0.0],
            "O": [0.0, 0.0, 0.0],
            "eye": [0.0, 0.0, 5.0],
            "debug": {},
        }
        response = self.client.post(
            self.url, data=json.dumps(body),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 400)

    # ── 5. Ortho O translates output by a constant 2D shift ──────────────────

    def test_ortho_O_translates_polylines_by_constant_shift(self):
        # Generic random axis (not aligned with the paraboloid's symmetry plane).
        # Canonical / axis-aligned views hit boundary-edge X-coordinate
        # degeneracies in the leftmost-anchor pick and produce different
        # absolute vis labels — geometry still translates, but labels and
        # polyline counts diverge. Random generic axes are reliable.
        rng = np.random.default_rng(7)
        a = rng.standard_normal(3); a /= np.linalg.norm(a)
        helper = np.array([1.0, 0.0, 0.0])
        if abs(a @ helper) > 0.9:
            helper = np.array([0.0, 1.0, 0.0])
        I_arr = helper - (helper @ a) * a; I_arr /= np.linalg.norm(I_arr)
        J_arr = np.cross(a, I_arr)
        I = I_arr.tolist(); J = J_arr.tolist()
        O1 = [0.0, 0.0, 0.0]
        O2 = [0.3, 0.7, -0.2]

        def _post(O):
            r = self.client.post(
                self.url,
                data=json.dumps({"I": I, "J": J, "O": O, "eye": None, "debug": {}}),
                content_type="application/json",
            )
            self.assertEqual(r.status_code, 200, r.content)
            return json.loads(r.content)

        a = _post(O1)
        b = _post(O2)

        # Expected shift: XY(xyz; O2) − XY(xyz; O1) = −P·(O2 − O1),
        # where P = (I·, J·). Same for the world origin.
        I_arr = np.array(I); J_arr = np.array(J)
        dO = np.array(O2) - np.array(O1)
        expected_shift = np.array([-I_arr @ dO, -J_arr @ dO])

        np.testing.assert_allclose(
            np.array(b["origin"]) - np.array(a["origin"]),
            expected_shift, atol=1e-9,
        )

        # Same polyline counts per visibility key.
        self.assertEqual(
            sorted(a["lines_by_visibility"].keys()),
            sorted(b["lines_by_visibility"].keys()),
        )
        for vis_key, polys_a in a["lines_by_visibility"].items():
            polys_b = b["lines_by_visibility"][vis_key]
            self.assertEqual(len(polys_a), len(polys_b))
            for pa, pb in zip(polys_a, polys_b):
                arr_a = np.asarray(pa)
                arr_b = np.asarray(pb)
                self.assertEqual(arr_a.shape, arr_b.shape)
                np.testing.assert_allclose(
                    arr_b - arr_a,
                    np.broadcast_to(expected_shift, arr_a.shape),
                    atol=1e-9,
                )

    # ── 6. PROPAGATION=LP4 debug knob ────────────────────────────────────────

    def test_post_debug_propagation_lp4(self):
        # Tilted generic view — paraboloid + axis-aligned views can hit
        # BFS/LP edge cases that are not the assertion's concern.
        body = {
            "I": [1.0, 0.0, 0.3], "J": [0.0, 1.0, 0.4],
            "O": [0.0, 0.0, 0.0], "eye": None,
            "debug": {"PROPAGATION": "LP4"},
        }
        response = self.client.post(
            self.url, data=json.dumps(body),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        # All vis ≤ 0 under LP refinement (paraboloid is in the verified set).
        for k in data["lines_by_visibility"].keys():
            self.assertLessEqual(int(k), 0)

    # ── 7. Unknown debug key — 200 + warning logged ──────────────────────────

    def test_post_unknown_debug_key_warns_but_succeeds(self):
        body = {
            "I": [1.0, 0.0, 0.0], "J": [0.0, 0.0, 1.0],
            "O": [0.0, 0.0, 0.0], "eye": None,
            "debug": {"BOGUS": 42},
        }
        with self.assertLogs("surface_play.views", level="WARNING") as cm:
            response = self.client.post(
                self.url, data=json.dumps(body),
                content_type="application/json",
            )
        self.assertEqual(response.status_code, 200)
        self.assertTrue(
            any("BOGUS" in msg for msg in cm.output),
            f"expected warning mentioning BOGUS, got {cm.output}",
        )

    # ── 8. Caching observable across requests ────────────────────────────────

    def test_cache_observable_across_requests(self):
        before_misses = pipeline._cached_surface_init.cache_info().misses
        # First POST: cold cache.
        self.client.post(
            self.url,
            data=json.dumps({
                "I": [1, 0, 0], "J": [0, 0, 1], "O": [0, 0, 0],
                "eye": None, "debug": {},
            }),
            content_type="application/json",
        )
        info_after_first = pipeline._cached_surface_init.cache_info()
        self.assertEqual(info_after_first.misses, before_misses + 1)

        # Second POST: should hit the LRU.
        self.client.post(
            self.url,
            data=json.dumps({
                "I": [1, 0, 0], "J": [0, 0, 1], "O": [0, 0, 0],
                "eye": None, "debug": {},
            }),
            content_type="application/json",
        )
        info_after_second = pipeline._cached_surface_init.cache_info()
        self.assertGreaterEqual(info_after_second.hits, 1)
        self.assertEqual(info_after_second.misses, info_after_first.misses)

    # ── 9. 404 on missing pk ─────────────────────────────────────────────────

    def test_404_on_missing_pk(self):
        bad_url = reverse("surface-play", kwargs={"pk": 999_999_999})
        r_get = self.client.get(bad_url)
        self.assertEqual(r_get.status_code, 404)
        r_post = self.client.post(
            bad_url,
            data=json.dumps({"I": [1, 0, 0], "J": [0, 1, 0], "O": [0, 0, 0]}),
            content_type="application/json",
        )
        self.assertEqual(r_post.status_code, 404)

    # ── 10. Static check — no silhouette imports in views.py ─────────────────

    def test_views_module_does_not_import_silhouette(self):
        views_path = Path(__file__).parent / "views.py"
        src = views_path.read_text(encoding="utf-8")
        # Match `import surface_play.silhouette[...]` or
        # `from surface_play.silhouette[...] import`. Bare mentions in
        # comments / strings are fine.
        pattern = re.compile(
            r"^\s*(from\s+(?:\.|surface_play\.)?silhouette[\w.]*\s+import"
            r"|import\s+surface_play\.silhouette[\w.]*)",
            re.MULTILINE,
        )
        self.assertIsNone(
            pattern.search(src),
            "views.py still imports silhouette — W2 sign-off gate failed",
        )
