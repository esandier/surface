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
# all are consumed in curves.resample_all.

# BC build-polyline: subdivisions per boundary mesh-edge segment, so cumulative
# xy-arclength tracks the true projected curve through a projection fold
# (curves._densify_bc_polyline). One batched surface eval over the dense uv.
BC_DENSIFY_NSUB: int = 24
# HC (helper-curve) straight uv-line: number of dense pre-samples used to build
# an accurate arclength table before picking sample positions (resample_all HC
# branch).
HC_DENSIFY_N: int = 200
# (CC Newton refinement runs to convergence — see curves._newton_cc_refine —
# so there is no iteration-count setting.)
