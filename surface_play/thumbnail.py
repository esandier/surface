"""thumbnail.py — W3: SVG thumbnail rendering on `pipeline.build_outline`.

Pure migration of the legacy `silhouette.Surface`-based renderer onto the
modular pipeline. Shares the `_cached_surface_init` LRU with the Django
views, so a save-thumbnail POST that follows a play POST for the same
record reuses the cached `SurfaceInit`.

Spec: Modular_rewrite_roadmap.md lines 1661-1692.
"""

from __future__ import annotations

from surface_play import pipeline


def render_thumbnail(record, I, J, eye=None) -> str:
    """Render the supplied viewpoint to an SVG string.

    Output goes into ``SurfaceRecord.thumbnail`` (a TextField); the surface
    list template embeds it inline as HTML.
    """
    init = pipeline.build_surface_init(record)
    O = list(eye) if eye is not None else [0.0, 0.0, 0.0]
    outline = pipeline.build_outline(init, I=I, J=J, O=O, eye=eye)

    blank = (
        '<svg viewBox="-1.1 -1.1 2.2 2.2" xmlns="http://www.w3.org/2000/svg"'
        ' style="width:100%;height:100%;"></svg>'
    )

    all_pts: list[tuple[float, float]] = []
    for bucket in (outline.lines_by_visibility, outline.si_lines_by_visibility):
        for polylines in bucket.values():
            for poly in polylines:
                all_pts.extend(poly)
    if not all_pts:
        return blank

    xs = [p[0] for p in all_pts]
    ys = [p[1] for p in all_pts]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    cx = 0.5 * (x_min + x_max)
    cy = 0.5 * (y_min + y_max)
    r = 0.5 * max(x_max - x_min, y_max - y_min)
    if r == 0.0:
        return blank

    def _emit(bucket, stroke_width: str) -> list[str]:
        out: list[str] = []
        # vis == 0 ⇒ front-sheet visible; vis < 0 ⇒ hidden (dashed). Test
        # against absolute 0, not the bucket max — a SIC bucket may contain
        # only negative keys and must still draw dashed.
        for v in sorted(bucket.keys()):
            visible = v == 0
            opacity = "1" if visible else "0.35"
            dash = ' stroke-dasharray="0.05 0.05"' if not visible else ""
            for poly in bucket[v]:
                pts = " ".join(
                    f"{(p[0] - cx) / r:.4f},{-(p[1] - cy) / r:.4f}"
                    for p in poly
                )
                if pts:
                    out.append(
                        f'<polyline points="{pts}" fill="none" stroke="steelblue"'
                        f' stroke-width="{stroke_width}" stroke-opacity="{opacity}"{dash}/>'
                    )
        return out

    parts = _emit(outline.lines_by_visibility, "0.04")
    parts.extend(_emit(outline.si_lines_by_visibility, "0.06"))

    return (
        '<svg viewBox="-1.1 -1.1 2.2 2.2" xmlns="http://www.w3.org/2000/svg"'
        ' style="width:100%;height:100%;">'
        + "".join(parts)
        + "</svg>"
    )
