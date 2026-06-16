"""Module-level constants for surface_play (independent of Django settings).

Consumer functions accept matching kwargs that default to these globals,
so callers can override per-call without mutating module state.
"""

PROJECT_RESAMPLED: bool = False
# Resolution knob driving the construction mesh density (build_construction)
# and outline sampling density (resample_all). Per-call overrides via the
# `resolution=` kwarg on build_surface_init / resample_all.
RESOLUTION: int = 200
# Separate, coarser resolution for the 3D *display* mesh served to the canvas
# (build_mesh_init). The colored surface the user rotates does not need the
# full construction density — a coarser mesh ships a much smaller payload so
# the canvas appears sooner. The outline pipeline keeps RESOLUTION for full
# precision; canvas and outline are independent representations of S(u, v).
# NB: distinct from the roadmap's CANVAS_RESOLUTION (which there means the
# outline *sampling* density); this is the display-mesh density only.
DISPLAY_RESOLUTION: int = 80
SURFACE_CACHE_SIZE: int = 16

# Chain-step buffer used in helpers.py `_cc_samples` / `_sic_samples` to
# trim CP / DP samples within this many chain steps of any cusp / VP /
# triple-point endpoint of a SubCurve. Mirrors the legacy `trim` heuristic
# (silhouette.py ~line 957). Avoids placing HAs on top of an existing
# cusp where the front-sheet normal evaluation is degenerate.
HA_CUSP_TRIM: int = 5

# ── Resampling densification / Newton refinement ─────────────────────────────
# Centralized here as the single source of truth (e.g. for the debug panel);
# consumed in curves.resample_all.

# Densification oversampling factor, shared by the BC build-polyline and the HC
# arclength table. Both build a dense polyline at spacing `ell / DENSIFY_SUBDIV`
# (ell = M / resolution is the global coarse sample spacing, M = mesh xy-bbox
# diagonal) so cumulative image-arclength stays accurate through projection
# folds. Equivalently the dense point count is `resolution · DENSIFY_SUBDIV ·
# L / M` for a curve of image-length L. Higher = finer arclength tables (slower).
DENSIFY_SUBDIV: int = 10
# (CC Newton refinement runs to convergence — see curves._newton_cc_refine —
# so there is no iteration-count setting.)
