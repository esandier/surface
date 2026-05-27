"""Module-level constants for surface_play (independent of Django settings).

Consumer functions accept matching kwargs that default to these globals,
so callers can override per-call without mutating module state.
"""

PROJECT_RESAMPLED: bool = False
# Single resolution knob driving both mesh density (build_construction) and
# outline sampling density (resample_all). Per-call overrides via the
# `resolution=` kwarg on build_surface_init / resample_all.
RESOLUTION: int = 200
SURFACE_CACHE_SIZE: int = 16

# Chain-step buffer used in helpers.py `_cc_samples` / `_sic_samples` to
# trim CP / DP samples within this many chain steps of any cusp / VP /
# triple-point endpoint of a SubCurve. Mirrors the legacy `trim` heuristic
# (silhouette.py ~line 957). Avoids placing HAs on top of an existing
# cusp where the front-sheet normal evaluation is degenerate.
HA_CUSP_TRIM: int = 5
