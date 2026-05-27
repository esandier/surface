"""Module-level constants for surface_play (independent of Django settings).

Consumer functions accept matching kwargs that default to these globals,
so callers can override per-call without mutating module state.
"""

PROJECT_RESAMPLED: bool = False
RESOLUTION: int = 30

# W1 — pipeline settings.
# `GRID_RESOLUTION` is the construction-phase mesh resolution
# (build_construction). `CANVAS_RESOLUTION` is the outline-phase resampling
# resolution (resample_all); kept as an alias of `RESOLUTION` for backward
# compatibility with O14 callers that still default to `RESOLUTION`.
GRID_RESOLUTION: int = 30
CANVAS_RESOLUTION: int = RESOLUTION
SURFACE_CACHE_SIZE: int = 16

# Chain-step buffer used in helpers.py `_cc_samples` / `_sic_samples` to
# trim CP / DP samples within this many chain steps of any cusp / VP /
# triple-point endpoint of a SubCurve. Mirrors the legacy `trim` heuristic
# (silhouette.py ~line 957). Avoids placing HAs on top of an existing
# cusp where the front-sheet normal evaluation is degenerate.
HA_CUSP_TRIM: int = 5
