# Surface Outline App — Modular Rewrite Roadmap

**Status:** Draft, Layer P only (2026-05-04). Subsequent layers added after each sign-off.

**Spec (the bible):** [Reorganization_of_the_surface_app.md](Reorganization_of_the_surface_app.md). All step "Spec" references point to sections of that document.
**Reference only (do not import; do not read without authorization):** [surface_play/silhouette.py](surface_play/silhouette.py).
**Frontend:** [templates/play.html](templates/play.html) — preserve UX; minor additions allowed for debug views.
**Historical:** [ROADMAP.md](ROADMAP.md) (previous attempt — modularize-first, deferred SICs) and [old stuff/roadmap_that_failed.md](old%20stuff/roadmap_that_failed.md) (earlier attempt — bolted SIC onto silhouette.py).

---

## 0. Goals and constraints

### Goal
Replace the monolithic `silhouette.py` with a modular package under [surface_play/](surface_play/) that implements the spec end-to-end, including self-intersections.

### Non-goals
- No frontend redesign. Debug-only additions to `play.html` are allowed and clearly gated.
- No new model fields on `SurfaceRecord` unless a step explicitly justifies one.
- No swap of computational backends (sympy → JAX, triangle → meshzoo, etc.) inside this roadmap. Optimization later.
- No retirement of `silhouette.py` until the full pipeline reaches feature parity (Layer W).
- No taking into account the file `The four open bugs from `Problems of the surface app as of april.txt`, which is stale.


### Success criteria (final, end of Layer W)
1. [surface_play/views.py](surface_play/views.py) imports nothing from `silhouette.py`.
2. Vectorized parametrization eval (`_eval_vec` shape) live.
3. All Layer P/C/O unit tests pass; smoke tests on Helicoid, torus, paraboloid, Möbius band.
4. No frontend regression on the existing surfaces in the database.
5. More consistently good visual output from complex surfaces as `Onde radiale` or `vagues`.
6. Proper visual output for fig-8 immersion which is not properly dealt with by the current `silhouette.py`.

---

## 1. Architecture

### Module layout (flat, in [surface_play/](surface_play/))

```
surface_play/
  domain.py          P1   Domain (close, bary, area, bbox, generators)
  surface.py         P2   SurfaceParams (S, SN, derivatives, perturbation, lambdify+cse, _eval_vec)
  projection.py      P3   Projection (ortho/persp, XY, Z, ker_param, kerdS, per_vertex_viewer_dot)
  curves.py          P4,P6,O14  make_lines, sign_changes, ResampledCurve, resample_all
  intersections.py   P5,C8-C12  sweep_segments, build_bvh+query_bvh, find_double_points, build_sis_pairs, build_sics, find_triple_points
  mesh.py            C1-C7   _generate_rect_mesh, _generate_disk_mesh, _apply_identifications, _jitter, _build_edges_faces, build_mesh, build_bcs
  contour.py         O1-O4   find_contour_points, build_contour_segments, build_contour_curves, find_vps
  splitting.py       O5-O11,O13   SP/SPT primitives, corner/BCP/BDP/TP/VP/CDP splits, SubCurve assembly
  helpers.py         O12   build_helper_curves (with HA splits)
  visibility.py      O15-O17   compute_projection_breaks, bfs_visibility, lp_refine_visibility
  pipeline.py        W1   build_surface_init, build_outline, mesh_to_threejs
  views.py           W2   rewrite (uses pipeline, no silhouette import)
  thumbnail.py       W3   migrate to pipeline
  test_*.py          per-step unit tests
  silhouette.py      kept on disk, not imported
```

### Two computation phases

| Phase | Triggered by | Computes | Cached? |
|---|---|---|---|
| **Construction** (Layer C) | Surface/domain params or grid resolution change | Mesh, BCs, DPs, SISs, SICs, TPs | yes (per pk + cache_key) |
| **Outline** (Layer O) | Viewpoint or projection mode change | CPs, CCs, VPs, splits, resampling, helpers, projection breaks, visibility | no |

### views.py contract (preserved)
- `GET /play/<pk>/` returns the page with three.js mesh data (replaces `Surface.for_3js()`).
- `POST /play/<pk>/` returns `{lines_by_visibility, si_lines_by_visibility, origin}` for the current axis/eye, with optional `debug` knobs.
- The frontend sends `I`, `J`, `O`, `eye`, and `debug` exactly as today.

### Settings knobs (from spec §"Debug panel")
`GRID_RESOLUTION`, `CANVAS_RESOLUTION`, `NEWTON_CUSP`, `PROJECT_RESAMPLED`, `PROPAGATION ∈ {BFS, LP1, LP4}`. All pass through the `debug` POST dict to `build_outline(...)` kwargs. No module reads `settings.py` directly except the Django app config.

---

## 2. Catalog of known gotchas (referenced by ID throughout)

These come from prior-attempt notes. Each is referenced by its ID in the steps that need to handle it.

| ID | Gotcha | Detail |
|---|---|---|
| **G1** | Perturbation phase | Spec §"Perturbation" formula multiplies by `sin((u−u_min)·freq_u·π/d_u)·sin(...)` which vanishes at `u∈{u_min, u_max}` for `freq_u=2`. On cyl identification this puts seam vertices on the SIS axis. **Fix:** add fixed phases inside the sin: `sin((u−u_min)·freq_u·π/d_u + φ_u)`. `φ_u = 0` if `u_identify=='mo'` (mo seams need pert=0 for orientation reversal); else `φ_u = 0.543367` (arbitrary, breaks accidental zero coincidences). Same for v. |
| **G2** | Sibling-pair consumption | Spec §Appendix is mandatory. After Möller-Trumbore classifies an `(E1, F2)` hit as either a face or near-edge intersection, the sibling pairs `(E1, F2')`, `(E2, F1)`, `(E2, F1')` must be marked consumed so they're not re-emitted as duplicate DPs. Process highest-min_bary first. Without this, SICs fragment into many tiny chains. |
| **G3** | Close-aware sweeps | Every 2D-domain sweep (intersections.py:sweep_segments and all callers) must use `domain.close(p, q)` to lift the second endpoint into the same fundamental domain as the first. Otherwise identified-rectangle seam segments produce false intersections (large `q-p` vectors point across the domain). Spec §"2D line sweep" is explicit. |
| **G4** | Mesh generation and perturbation | (a) `triangle` is invoked with `-Y` (no boundary refinement), so enough initial boundary points must be given to it before mesh generation. (b) **Order:** raw mesh → **jitter** (independent per vertex; paired vertices need not be coherent — each retains its own jittered uv) → **compute vertex equivalence classes** (C4: union-find on side-identification pairs; `tris` is **not** rewritten — pre-id labels are preserved) → **build edges/faces** (C5: per-face `p, pq, pr` are direct subtractions of pre-id uvs — no `close()` needed since a triangle never spans the seam; edges merge across identified sides by the geometric rule of G13). Jitter `±0.1%` of typical edge length applied to all vertices. **True boundary vertices** (rect outer sides in the no-identification case; disk `r=r_min`/`r=r_max` ring) reprojected onto the boundary after jitter (clamp coord for rect; renormalize radius for disk). |
| **G5** | SplitPoint identity | A SP shared between two curves must be **the same Python object** in both, so visibility propagation can traverse via `id(sp)` lookup. Don't `dataclass(frozen=True, eq=True)` and rely on equality — store `dict[int(id(sp)), ...]`. |
| **G6** | HC dedup at same target sample | Two helper curves can argmin to the same sample of the same target curve, although no example comes to mind. The split must be applied once and the SP shared by both HCs. |
| **G7** | LP infeasibility = bug or numerical problem in computing visibility changes | `lp_refine_visibility` may report infeasible. **Do not relax the LP** to recover. Instead, surface the failure: it indicates δv values from `split_all` violate the global constraint, which is a bug or numerical problem in computing visibility changes. Raise `LPInfeasibleError`. |
| **G8** | ker_param via SVD | Implement `ker_param(p)` as the right-singular vector (smallest singular value) of the 2×2 Jacobian of the projected surface parametrization `(u, v) ↦ XY(S(u, v))`. Ortho: `J = [[I·Su, I·Sv], [J·Su, J·Sv]]`. Persp: chain rule through the perspective divide gives `J[r, c] = basis_r·Sc − (basis_r·d / z)·(n·Sc)` with `d = S(p) − eye`, `z = n·(S(p) − eye)`, `n = (I × J)/‖·‖`. Drop the uniform `1/z` for the SVD (right-singular vectors are scale-invariant); keep it in `proj_vec`. |
| **G9** | Newton fallback | When Newton on `s ∈ (0, 1)` for `inner(SN, axis) = 0` escapes the interval, fall back to `s = 0.5`. Don't loop, don't extrapolate. |
| **G10** | δv ∈ {−1, 0, +1} | All split visibility-change formulas yield ±1 or 0 per side. Magnitudes ±2 only appear at projection breaks (CC occluder). Anything else means a sign convention got dropped. |
| **G11** | cse() before lambdify | Apply `sympy.cse()` to the joint expressions for `[S, Su, Sv, Suu, Suv, Svv, SN]` so they share subexpressions in a single lambdified callable. ~10× speedup on the construction phase. Spec §"Implementation note" line 147. |
| **G12** | Pure vectorized eval | Lambdified callables are wrapped as `def _eval_vec(fn, u, v): return np.array(fn(u, v))` so a single call evaluates over arrays. No Python `for` over points. Spec lines 147-150. |
| **G13** | Edge/face identity is geometric, not by vertex labels | Two pre-id edges merge into one post-id edge **iff** both endpoints lie on a single identified side of the rectangle (e.g., both on u-min) **and** their endpoints pair under that side's C1 identification rule (so their partners lie together on the opposite identified side). All other edges keep their own identity. **Triangles never merge** — each pre-id triangle is its own 2-cell. The label-based rule "two cells with the same canonical vertex set are the same cell" is **wrong**: under mo-mo at corners, two geometrically distinct pre-id triangles in opposite corners can have all three vertices pairwise identified, producing identical canonical triples for different cells. Treating them as duplicates removes a valid 2-cell and tears the manifold (verified: rect mo-mo at res=5 with the old dedup produced 1 edge with 3 faces and 2 edges with 1 face). Because each post-id edge and each face stays inside a single fundamental-domain copy, per-element geometric fields (`p, pq, pr, dir`) are direct subtractions of the pre-id uvs — **no `close()` is needed for them**. `close()` is still required for cross-element navigation (curve traversal across the seam, DP uv interpolation per G15). The spec §"Data stored" "computed before identification" prescription is now satisfied **literally**, by keeping pre-id vertex indices in `tris` and reading `uv` straight from those. |
| **G14** | BVH broad phase for DP search | Möller-Trumbore over all (edge, face) pairs is `O(E·F)` — intractable at `res ≥ 100`. **AABB-broadcast pre-filter is also rejected** (it's `O(E·F)` memory). Required: a stack-based BVH over edge AABBs vs face AABBs, partitioning on cycling axes (`depth % 3`), brute-force at leaves. Numba-JIT compiled (cf. legacy [old stuff/intersections_prev.py](old%20stuff/intersections_prev.py) `_numba_bvh_final_boss`). Output deduped via combined-key sort. Spec line 207. |
| **G15** | DP uv interpolation uses `close()` | When recovering domain coordinates of a DP from edge/face barycentrics, the implied vertex uv coords must be lifted by `domain.close(p, q)` before interpolation, otherwise identification seams produce wildly wrong DP uvs. Spec line 209. |
| **G16** | BCs from no-id rect are ONE closed loop | `make_lines(boundary_edges)` on an unidentified rectangle yields a single closed BC traversing all four sides — corners do **not** split it. The four open arcs from corner-to-corner emerge later in Layer O via SPTs of type `cn`. Don't try to "fix" it inside `build_bcs`. Spec lines 47-48 + our P4 decision. |
| **G17** | `split1, split2` slot count | Per spec line 245, every BE, CS, SIS carries two split slots `split1, split2` (each either -1 or an index into the global SPTs array). Layer C populates segment dtypes with `split1 = split2 = -1`; Layer O fills them. **Risk:** in pathological surfaces a segment could attract >2 splits (e.g., a SIS with a BDP at one end, a TP in the middle, and a CDP near the other end). The new code must `assert split2 == -1` before assigning a third — surface as a clear failure rather than silently dropping the third split. If it ever fires in practice we'll widen the dtype. |
| **G18** | BCP `bvis_chge` formula | Visibility change of a BC at a contour-on-boundary point (BCP): copy verbatim the `bvis_chge` formula from spec lines 258-280. It uses `kerdS` (oriented toward upper fold), `Tb` (boundary tangent in 3D), `Np` (projected fold normal). The orientation logic of `Np` and the sign cases are easy to get wrong — write the formula once in `splitting.py` and unit-test against synthetic configurations before using it in real surfaces. |
| **G19** | HC distance in domain, original coords | Helper-curve component distances are computed in the **domain** (not view-plane), and on identified rectangles **the original (pre-merge) domain rectangle coordinates** are used (spec line 403). For two seam-side candidate endpoints, this means avoiding the period-jump that `close()` would introduce — we want `d(qi, qj)` as raw rectangle-coord distance, not minimum-period distance. |
| **G20** | Resampling spacing per spec | Spec lines 79-80 + 442: `d = min(L/10, M / GRID_RESOLUTION)` where `L` is the smallest **view-plane** length of curves originating/ending at the SP, and `M` is an approximate diameter of the surface in projected space. Sample at arclengths `0, d, 2d, ...` until `L/2 − d` on each half-curve. This guarantees ≥4 points per half-curve. **Don't** measure curve length in domain space; the projection-space metric is what BK detection cares about. |
| **G21** | SIC two-sheet uv tracking on resample | A SIC has one chain (in xyz) but two domain preimages (uv1 and uv2 per DP). When resampling a SIC SubCurve we need both preimages tracked. Walk consecutive DPs using each SIS's `flip` to decide which preimage continues into which side; between consecutive samples on each preimage, use `domain.close()` for the uv linear interpolation (verified on fig-8 cyl seam in prior attempt). |
| **G22** | BKs are sidecar, not splits | Spec line 467: visibility changes from projection-space intersections (BKs) are stored on curves but the curves are **not** split at BKs. BKs are sample-anchored sidecar data: `{sample_idx, t, delta_v, xy}`. Don't accidentally insert BKs as new SPTs in `split1, split2` — they belong to a separate `breaks` array. |
| **G23** | BFS anchor determinism | The leftmost-x sample (anchor for BFS) may tie. On ties, pick the lexicographically smallest `(rc_id, sample_idx)`. Otherwise different runs of the same input produce different visibility distributions. |
| **G24** | LP variables and slack | LP1/LP4 (spec §15.2): variables = per-sample `vis[k]` (≤ 0) + per-break `vc_i` (free) + per-break slack `s_i ≥ 0`. Equality rows = segment propagation `vis[k+1] − vis[k] − Σ vc_i = 0` (with wrap for closed) + SP coupling (visibilities at shared SPs match modulo δv) + anchor `vis = 0`. Inequality rows = `vc ∓ s ≤ ±vc⁰` for L1 slack on objective `Σ s_i`. Use scipy `linprog(method="highs")` on a sparse CSR matrix. Round `vis` to int at exit; if status=2 (infeasible), raise `LPInfeasibleError` (G7). |

Layer O is the only consumer of G18-G24.

---

## 3. Test strategy

### Per-step
- Each step ships its own `test_<module>.py` covering inputs, outputs, edge cases, and the gotcha(s) it touches. Run as `pytest surface_play/test_<module>.py`.
- Test fixtures: a small library of sympy surfaces (paraboloid, torus-cyl, fig-8 immersion-cyl, Möbius band, disk-of-functions). Defined in `test_fixtures.py` (created in P2). Avoid sphere which cannot be parametrized regularly in its entirety.

### Per-layer
- A `test_<layer>_integration.py` file that runs the full layer end-to-end on 2-3 fixtures.
- Layer P integration test: build a Möbius band mesh's worth of primitive calls (close, bary, lambdify, sweep) and assert no exceptions + sane shapes.

### Intermediate views
- For Layer P: matplotlib **probe scripts** (similar to existing `diag_*.py` at repo root) — one per primitive, runnable directly. Not part of pytest.
- For Layer C: a Django dev-only debug view at `/debug/<pk>/construction/` that renders mesh + DPs + SICs as overlays on the existing 3D canvas. Added in step C-final.
- For Layer O: extension of the same debug view to overlay CCs, VPs, splits, helpers, breaks. Added incrementally per O-step.

### Visual regression (final)
- Deferred to Layer W. Compare `pipeline.build_outline(...)` JSON to `silhouette.Surface.traitement(...)` JSON on each fixture surface. Tolerances per channel.

---

## 4. Sequencing

- **Layer P (primitives, 6 steps).** No surface logic. Each primitive is reused by C and O.
- **Layer C (construction, 14 steps).** Mesh and view-independent topology.
- **Layer O (outline, 17 steps).** Per-viewpoint geometry, splits, helpers, resampling, visibility.
- **Layer W (wiring, ~5 steps).** pipeline.py, views.py rewrite, thumbnail migration, smoke + visual regression.

Each layer ends with an integration test; **work on the next layer does not start until the current layer's integration test passes and the user has signed off**.

---

## Layer P — Primitives (6 steps)

Every step has the same shape: **Spec → Files → API → Algorithm → Test criterion → Gotchas → Stop & verify**.

---

### P1 — `domain.py`: Domain class

**Spec:** §"Domain of the parametrization" (lines 90-116), §"Domain API" (line 113).
**Files:** create [surface_play/domain.py](surface_play/domain.py), create [surface_play/test_domain.py](surface_play/test_domain.py).
**API:**
```python
class Domain:
    type: Literal["rect", "disk", "annulus"]
    bounds: tuple[float, float, float, float]  # (u_min, u_max, v_min, v_max) for rect; (r_min, r_max, 0, 2π) for disk/annulus
    u_identify: Literal["no", "cy", "mo"]
    v_identify: Literal["no", "cy", "mo"]
    coord_type: Literal["ca", "po"]          # cartesian or polar (disk only)

    @property
    def period_u(self) -> float: ...         # u_max − u_min (for rect with cy/mo on u-axis); else NaN
    @property
    def period_v(self) -> float: ...         # v_max − v_min (for rect with cy/mo on v-axis); else NaN

    def close(self, p: np.ndarray, q: np.ndarray) -> np.ndarray: ...   # supports (2,) or (N, 2) inputs
    def bary_edge(self, p: np.ndarray, pq: np.ndarray, s: float) -> np.ndarray: ...     # p + s·pq
    def bary_face(self, p: np.ndarray, pq: np.ndarray, pr: np.ndarray, s: float, t: float) -> np.ndarray: ...
    def bbox_diag(self) -> float: ...
    def width(self) -> float: ...           # rect diagonal or 2·r_max
```
**Algorithm:** `close(p, q)`: for `rect` with cyl/mo on a given axis, shift `q[axis]` by an integer multiple of the period to minimize `|q[axis] − p[axis]|`. Return shifted `q`. For `disk`/`annulus`, return `q` unchanged. Vectorized: accept `p, q` of shape `(2,)` or `(N, 2)`.
**Test criterion** (`test_domain.py`):
1. `close((0.1, 0.5), (6.0, 0.5))` on `u_identify="cy"`, `bounds=(0, 2π, 0, 1)` → `(6.0 − 2π, 0.5) ≈ (-0.283, 0.5)`.
2. `close((0.1, 0.5), (0.2, 0.5))` returns `(0.2, 0.5)` unchanged.
3. Vectorized `close` on `(N, 2)` array gives same result as element-wise.
4. `bary_edge` and `bary_face` give expected interpolations.
5. `bbox_diag`, `width` sanity values for a unit square and a unit disk.
**Gotchas:** G3.
**Stop & verify:** `pytest surface_play/test_domain.py` ⇒ 5+ green. No external imports beyond numpy.

---

### P2 — `surface.py`: parametrization, perturbation, lambdified callables

**Spec:** §"Parametrization of the surface" (lines 118-152), perturbation formula lines 124-139.
**Files:** create [surface_play/surface.py](surface_play/surface.py), [surface_play/test_surface.py](surface_play/test_surface.py), [surface_play/test_fixtures.py](surface_play/test_fixtures.py).
**API:**
```python
class SurfaceParams:
    def __init__(self, X: str, Y: str, Z: str, parameter_names: str, domain: Domain,
                 perturb: bool = True, output_type: Literal["ca", "cy"] = "ca"): ...

    # all of these accept (u, v) scalars OR (N,) arrays and return (3,) or (3, N)
    S: Callable
    Su: Callable; Sv: Callable
    Suu: Callable; Suv: Callable; Svv: Callable
    SN: Callable                # cross product (Su × Sv)

    bbox_diag: float            # diameter of the surface (for perturbation magnitude)
```
**Algorithm:**
1. Parse the sympy expressions. If `coord_type == "po"`, substitute polar→cartesian in (u, v). If `output_type == "cy"`, wrap the result in cylindrical→cartesian.
2. Compute symbolic Su, Sv, Suu, Suv, Svv, SN.
3. **Perturbation (G1):** if `perturb`, add to `S`:
   ```
   freq_u = 2 if u_identify != "mo" else 1
   freq_v = 2 if v_identify != "mo" else 1
   φ_u = 0 if u_identify == "mo" else 0.543367
   φ_v = 0 if v_identify == "mo" else 0.172935
   pert = sin((u - u_min) · freq_u · π / d_u + φ_u) · sin((v - v_min) · freq_v · π / d_v + φ_v) · 0.0005 · (dX·0.346, dY·0.632, dZ·0.693)
   ```
4. **cse + single lambdify (G11):** apply `sympy.cse([S, Su, Sv, Suu, Suv, Svv, SN])` then lambdify as one callable returning a tuple. Wrap as `_eval_vec` (G12) for vectorized evaluation.
**Test criterion:**
1. Helicoid `(u cos(v), u sin(v), v)`: `S(0, 0)` ≈ `(0, 0, 0)`; `SN(0, 0)` parallel to `(0, 1, 0)`.
2. Perturbation magnitude: `‖S_perturbed − S_unperturbed‖ < 0.001 · bbox_diag` over a 100-point grid.
3. **Identification preservation under perturbation** (per G1):
   - `u_identify="cy"`: `S(u_min, v) ≈ S(u_max, v)` for sample v's, to 1e-12.
   - `u_identify="mo"`: `S(u_min, v) ≈ S(u_max, v_max+v_min−v)` for sample v's, to 1e-12.
   - same for v-axis.
4. Vectorized: `S(u_array, v_array)` returns shape `(3, N)` and matches scalar calls.
5. `cse` actually applied: timing single eval of all 7 functions vs separate evals shows ≥3× speedup on a 1000-point grid.
**Gotchas:** G1, G11, G12.
**Stop & verify:** `pytest surface_play/test_surface.py` ⇒ green; speedup assertion in test 5 reported.

---

### P3 — `projection.py`: Projection class

**Spec:** §"Projection of the surface" (lines 155-159), `kerdS` formula in spec line 273.
**Files:** create [surface_play/projection.py](surface_play/projection.py), [surface_play/test_projection.py](surface_play/test_projection.py).
**API:**
```python
class Projection:
    def __init__(self, surface: SurfaceParams, I: list[float], J: list[float],
                 O: list[float] | None = None,
                 eye: list[float] | None = None): ...
    # mode: "ortho" if eye is None else "persp"
    # O: world-space point that maps to view-plane (0, 0). Defaults to [0,0,0] in ortho;
    #    must equal eye in persp (frontend always sends O = eye in persp mode).

    def XY(self, xyz: np.ndarray) -> np.ndarray: ...  # 3D → 2D view plane (anchored at O)
    def Z(self, xyz: np.ndarray) -> float: ...        # depth scalar along image-plane normal
    def viewer_direction(self, xyz: np.ndarray | None = None) -> np.ndarray: ...
        # Ortho: xyz optional; no-arg returns the constant (3,) unit axis (= I × J normalized).
        # Persp: xyz required (raises ValueError otherwise); returns `xyz - eye` (unnormalized
        # — every consumer is sign-only on dot products with SN, which is itself unnormalized).
        # No public `axis` attribute: forces every persp call site to supply xyz, eliminating
        # the silent-bug class where ortho code reaches for `proj.axis` and breaks under persp.
    def per_vertex_viewer_dot(self, mesh) -> np.ndarray: ...         # SN · viewer_direction at every vertex
    def ker_param(self, p: np.ndarray) -> np.ndarray: ...            # 2D kernel direction in (u, v) — see G8
    def kerdS(self, uv: np.ndarray) -> np.ndarray: ...               # 3D image of ker_param via dS
    def proj_vec(self, uv: np.ndarray, dir3d: np.ndarray) -> np.ndarray: ...  # 3D vector → 2D view plane
```
**Algorithm:**
- `O` defaults: `O = [0, 0, 0]` in ortho if not supplied; assert `O == eye` in persp (frontend contract from [templates/play.html:330](templates/play.html#L330)).
- Ortho mode: image-plane normal `n = I × J` (normalized, internal `_axis`); `XY(xyz) = (I·(xyz − O), J·(xyz − O))`; `Z(xyz) = n·(xyz − O)`. `viewer_direction()` returns `n`.
- Persp mode: `viewer_direction(xyz) = xyz − eye` (unnormalized — sign-only consumers, matching SN). XY uses standard perspective divide along `n`: `XY(xyz) = (I·d, J·d)/(n·d)` with `d = xyz − eye`; eye projects to `(0, 0)` by definition (zero-vector special case). `Z(xyz) = n·(xyz − eye)` is signed depth along the image-plane normal.
- `ker_param(p)` (G8): build the 2×2 Jacobian of `(XY ∘ S)` at `p`. SVD; right-singular vector for the smallest singular value is `ker_param`. Ortho: `J = [[I·Su, I·Sv], [J·Su, J·Sv]]`; `O` drops out (constant translation). Persp: chain-rule through the perspective divide gives `J[r, c] = basis_r·Sc − (basis_r·d / z)·(n·Sc)` with `d = S(p) − eye`, `z = n·d`; the uniform `1/z` prefactor is dropped for SVD (irrelevant to right-singular vectors). At a persp contour point this is rank-deficient and the kernel direction's image under dS is parallel to `d` (along the viewing ray). `proj_vec` in persp uses the same Jacobian applied to `dir3d`, but keeps the `1/z` (linear-map scale carries information).
- `per_vertex_viewer_dot` is vectorized (`np.einsum`).

**Test criterion:**
1. Ortho with `I=[1,0,0], J=[0,1,0], O=None`: `XY([1, 2, 3]) == (1, 2)`; `Z([1, 2, 3]) == 3`.
2. **Ortho with non-default `O`:** `I=[1,0,0], J=[0,1,0], O=[5, 0, 0]`: `XY([6, 2, 3]) == (1, 2)` (the `−O` shift is applied).
3. Persp with `eye=[0,0,5]`, `O=[0,0,5]`, `I=[1,0,0], J=[0,1,0]`: `XY([0, 0, 5]) == (0, 0)` (eye projects to view-plane origin).
4. **Persp `O ≠ eye` rejected:** constructor raises `ValueError` if persp mode is given `O != eye`.
5. `per_vertex_viewer_dot` shape `(N,)`, matches per-vertex loop on a 50-vertex mesh.
6. Persp `ker_param` converges to ortho `ker_param` as `eye` recedes along axis (regression on the chain-rule formula).
7. Persp `ker_param` at a known persp contour point (paraboloid `Z=u²+v²`, `eye=(0,0,−h)`, contour circle `u²+v²=h`): `kerdS` is parallel to `S(p) − eye`.
8. Persp `proj_vec` matches the finite-difference linearization of `XY` along `dir3d` (regression on the `1/z` scale).
**Gotchas:** G8.
**Stop & verify:** `pytest surface_play/test_projection.py` ⇒ 8 green.

---

### P4 — `curves.py`: `make_lines`

**Spec:** §"Curves" (lines 174-181), §"make_lines" (line 8).
**Files:** create [surface_play/curves.py](surface_play/curves.py), [surface_play/test_curves.py](surface_play/test_curves.py).
**API:**
```python
def make_lines(segments: np.ndarray) -> list[np.ndarray]:
    """
    segments: (N, 2) int array — each row is a pair of endpoint indices.
    Returns: list of 1D int arrays. Each array is a chain of segment indices, with negative
             values denoting reversed traversal. Closed chains end where they began.
    """
```
**Algorithm:** Vectorized argsort-based adjacency. Build half-edges `[(p_i, +i), (q_i, -(i+1))]`; sort by endpoint index; consecutive pairs of equal endpoints are adjacencies. Branch points (degree ≥ 3) block chaining. Then sequential traversal from a degree-1  vertex, or from any segment for a closed loop.

**Test criterion:**
1. 4 segments forming a square (open path, no closure): one chain, length 4.
2. 4 segments forming a closed quad: one closed chain.
3. Two disconnected loops: 2 closed chains.
4. Empty input: empty list.
**Gotchas:** The routine will be used for segments data where the vertices have degree at most 2. Then the resulting curve are either closed, or open with endpoints having degree 1.
**Stop & verify:** `pytest surface_play/test_curves.py::test_make_lines` ⇒ 4 green.

---

### P5 — `intersections.py`: `sweep_segments`

**Spec:** §"2D line sweep for intersections" (line 10).
**Files:** create [surface_play/intersections.py](surface_play/intersections.py), [surface_play/test_intersections.py](surface_play/test_intersections.py).
**API (matches legacy [old stuff/intersections_prev.py](old%20stuff/intersections_prev.py)):**
```python
intersect_dtype = np.dtype([
    ("uv",  "f8", 2),    # 2D intersection point in seg_a's anchor frame
    ("a",   "i4"),       # segment index in seg_a_*
    ("b",   "i4"),       # segment index in seg_b_*
    ("t_a", "f8"),       # parametric coord on a, in (0, 1)
    ("t_b", "f8"),       # parametric coord on b, in (0, 1)
])

def sweep_segments(seg_a_uv0: np.ndarray, seg_a_uv1: np.ndarray,
                   seg_b_uv0: np.ndarray | None, seg_b_uv1: np.ndarray | None,
                   domain: Domain, *,
                   self_sweep: bool = False, tol: float = 1e-9) -> np.ndarray:
    """
    Returns a structured array of dtype intersect_dtype.
    If self_sweep: seg_b_* must be None; sweep A against itself with i < j filter.
    domain.close() / Domain.period_{u,v} are used to anchor each b-segment in a's frame.
    """
```
**Algorithm:** Bounding-box sweep on x-axis (sort segments by min-x; scan a moving window of overlapping segments by x-interval). For each candidate pair, **first apply close-aware anchor shift** (G3): if domain is given, replace `q_a` with `domain.close(p_a, q_a)`, and at hit time compare `p_b` and `q_b` after shifting them into `p_a`'s frame. Then run a standard parametric segment-intersection test (solve the 2×2 linear system); accept iff both t's lie in `(0, 1)` strictly (skip endpoints — a separate pass handles those).
**Test criterion:**
1. Cross mode: `[(0,0)→(1,1)]` vs `[(0,1)→(1,0)]` → one hit at `(0.5, 0.5)`.
2. Self mode `self_sweep=True`: input includes the same segment twice → no hit (i<j filter).
3. Close-aware: rect `(0, 2π)` cyl-u; segs `[(0.1, 0.5)→(0.2, 0.5)]` vs `[(6.0, 0.5)→(6.18, 0.5)]` → hit detected (after close shifts the second to `(6.0−2π, 0.5)→(6.18−2π, 0.5)`). Without `domain`, no hit.
4. Parallel non-overlapping segments → no hit.
5. Empty inputs → `[]`.
6. Performance smoke: 1000 random segments, self-sweep, completes in <1s and reports a sane hit count.
**Gotchas:** G3.
**Stop & verify:** `pytest surface_play/test_intersections.py::test_sweep_segments` ⇒ 6 green.

---

### P6 — `curves.py`: `sign_changes`

**Spec:** §"Detection of sign changing on a collection of segments" (line 11).
**Files:** extend [surface_play/curves.py](surface_play/curves.py); add tests to [surface_play/test_curves.py](surface_play/test_curves.py).
**API:**
```python
def sign_changes(vals_p: np.ndarray, vals_q: np.ndarray,
                 flip: np.ndarray | None = None) -> np.ndarray:
    """
    Returns boolean mask shape (N,): True where vals_p[i] * vals_q[i] * flip[i] < 0.
    flip defaults to all +1. flip ∈ {-1, +1} per segment (Möbius mesh edges have -1).
    """
```
**Algorithm:** one-liner: `(vals_p * vals_q * (flip if flip is not None else 1)) < 0`. Centralized so CP detection (§Contour points), VP detection (§Cusp points), and DP edge candidates all share one definition.
**Test criterion:**
1. `vals_p=[1,-1,2,-3], vals_q=[-1,1,3,-2]` → `[True, True, False, False]`.
2. With `flip=[-1,1,1,1]` on the same input → `[False, True, False, False]` (the first sign is flipped).
3. Empty arrays → empty mask.
**Gotchas:** none new.
**Stop & verify:** `pytest surface_play/test_curves.py::test_sign_changes` ⇒ 3 green.

---

### Layer P integration

**Files:** create [surface_play/test_layer_p_integration.py](surface_play/test_layer_p_integration.py).
**Goal:** end-to-end smoke that the primitives compose. On a helicoid  fixture:
1. Build `Domain` (rect, no identification), `SurfaceParams`, `Projection`.
2. Generate a tiny manual triangulation (e.g., 2 triangles forming a square in (u, v)).
3. Compute `SN` at the four vertices.
4. Run `sign_changes` on a synthetic per-edge value array.
5. Run `make_lines` on those edges (manually constructed).
6. Run `sweep_segments` on two crossing edge sets.
7. Assert all functions return without error and produce shapes/types per their contract.

**Layer P probe scripts (matplotlib, optional but recommended):**
- `probe_p1_close.py` — visualize `close(p, q)` on a rect-cyl domain: scatter random `q`'s, draw arrows to `close(p, q)`.
- `probe_p2_perturbation.py` — plot perturbation magnitude over (u, v) for cyl, mo, no — verify identification preserved.
- `probe_p4_make_lines.py` — render a tangled segment set and the chains `make_lines` produces.
- `probe_p5_sweep.py` — render two random segment families and the hits found.

These live at repo root alongside existing `diag_*.py`. Not part of pytest; runnable by hand for sanity-checking.

---

## Layer P sign-off gate

Before Layer C is drafted:
1. All 6 step tests + integration test green.
2. User has visually inspected at least probe_p1, probe_p2, probe_p4, probe_p5.
3. User confirms that the `Domain`, `SurfaceParams`, `Projection`, `make_lines`, `sweep_segments`, `sign_changes` APIs match what they want (signatures may change before Layer C locks them in).

---

## Layer C — Construction (14 steps)

View-independent topology: mesh assembly (rect + disk + identifications + jitter + edge/face arrays), boundary curves, BVH primitive, double points, self-intersecting segments, self-intersecting curves, triple points, plus an orchestrator and matplotlib probe scripts. Browser-side debug overlay deferred to Layer W (or earlier if the user wants).

The legacy file [old stuff/intersections_prev.py](old%20stuff/intersections_prev.py) is the algorithmic reference for steps C8-C12. Use it like the spec treats `silhouette.py` — read-only reference, copy logic where correct, improve where the new architecture warrants.

All construction outputs are cacheable (keyed on surface + grid resolution), so they re-run only when those change.

Each step has the same shape: **Spec → Files → API → Algorithm → Test criterion → Gotchas → Stop & verify.**

---

### C1 — `mesh.py`: `_generate_rect_mesh`

**Spec:** §"Domain of the parametrization" + §"Mesh" (lines 162-171).
**Files:** create [surface_play/mesh.py](surface_play/mesh.py), create [surface_play/test_mesh.py](surface_play/test_mesh.py).
**API:**
```python
def _generate_rect_mesh(domain: Domain, resolution: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns:
      uv:   (N, 2) float64 — vertex (u, v) coordinates.
      tris: (M, 3) int32   — triangle vertex indices.
    No identifications applied yet (paired vertices have distinct indices). No jitter.
    """
```
**Algorithm:** Place `resolution + 1` boundary vertices uniformly along each side of the rectangle (corners shared, total `4·resolution` distinct boundary vertices). Build PSLG with `4·resolution` constrained boundary segments. Call `triangle.triangulate(pslg, opts="pYq30a{area}")` where `area ≈ rect_area / resolution²`. Return raw vertex array (with original boundary coordinates) and triangle index array.
**Test criterion:**
1. Unit square at `resolution=10`: `N` grows as ~resolution²; `M` ≈ 2·N; Euler χ = 1 (disk topology).
2. All `4·resolution` boundary vertices distinct, lying exactly on the boundary.
3. No interior triangles missing (assert min triangle area > 0; sum of areas ≈ domain area).
4. `tris` is `int32`, `uv` is `float64`.
**Gotchas:** none new (G4(a) is satisfied: pre-placed boundary points feed `-Y`).
**Stop & verify:** `pytest surface_play/test_mesh.py::test_generate_rect_mesh` ⇒ green.

---

### C2 — `mesh.py`: `_generate_disk_mesh`

**Spec:** §"User specification" (line 99), §"Mesh" (line 162).
**Files:** extend [surface_play/mesh.py](surface_play/mesh.py); add tests to [surface_play/test_mesh.py](surface_play/test_mesh.py).
**API:**
```python
def _generate_disk_mesh(domain: Domain, resolution: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Disk (r_min == 0) or annulus (r_min > 0). Returns (uv, tris) in cartesian (u, v).
    """
```
**Algorithm:** Place `2π·r_max·resolution / (rect_diag of bbox)` vertices uniformly on the outer ring. For annulus, also place vertices on the inner ring with proportional density. PSLG: outer ring as one constrained loop, inner ring (if any) as a second constrained loop with `holes` parameter pointing into the inner disk. Call `triangle.triangulate(opts="pYq30a{area}")`.
**Test criterion:**
1. Unit disk `r_max=1`, `resolution=20`: Euler χ = 1.
2. Unit annulus `r_min=0.3, r_max=1`: Euler χ = 0; no triangles inside the inner ring (test by area).
3. Boundary vertices lie within `1e-12` of `r=r_max` (and `r=r_min` for annulus).
4. `tris` is `int32`, `uv` is `float64`.
**Gotchas:** G4(a).
**Stop & verify:** `pytest surface_play/test_mesh.py::test_generate_disk_mesh` ⇒ green.

---

### C3 — `mesh.py`: `_jitter` (pre-identification)

**Spec:** spec §"Mesh", G4(b).
**Files:** extend [surface_play/mesh.py](surface_play/mesh.py); add tests.
**API:**
```python
def _jitter(uv: np.ndarray, tris: np.ndarray, domain: Domain, *,
            seed: int | None = None) -> np.ndarray:
    """
    Returns a jittered copy of uv. Independent jitter per vertex by ±0.1% of typical
    edge length. Reprojects true boundary vertices onto the domain boundary. Run BEFORE
    `_apply_identifications` (C4): paired vertices receive independent jitter — the
    discarded one is dropped at identification, so coherence is unnecessary.
    RNG seeded by `seed` (default = id(uv) for repro).
    """
```
**Algorithm:**
1. Compute typical edge length `L` from the mesh: median of `‖uv[tris[:,1]] − uv[tris[:,0]]‖` etc. across all edges.
2. Identify true boundary vertices: for rect no-id, all 4 sides; for rect with cy on u, only top/bottom v-sides remain true boundary; etc. For disk, the `r=r_max` (and `r_min` for annulus) ring. Vertices on identified-only sides (will become seam after C4) are NOT true boundary; they get free jitter.
3. Sample `delta = uniform(-1, 1, (N, 2)) · 0.001 · L`.
4. Apply: `uv_new = uv + delta`.
5. Reproject true boundary vertices: clamp to nearest side for rect; renormalize `‖(u, v)‖` to `r_min` or `r_max` for disk/annulus.
**Test criterion:**
1. Bbox of post-jitter `uv` differs from pre-jitter by `< 0.002 · L` per axis.
2. True boundary vertices remain exactly on the boundary post-jitter (`< 1e-12` deviation).
3. Reproducibility: same seed → identical `uv_new`.
4. No two vertices collapse onto each other (min distance > 0).
5. Vertex count unchanged (jitter doesn't add or remove).
**Gotchas:** G4(b).
**Stop & verify:** `pytest surface_play/test_mesh.py::test_jitter` ⇒ green.

---

### C4 — `mesh.py`: `_apply_identifications` (vertex equivalence only)

**Spec:** §"Mesh" lines 162-176, §"Data stored — Identification of edges and faces" lines 111-119, §"User specification" lines 95-99.
**Files:** extend [surface_play/mesh.py](surface_play/mesh.py); add tests.
**API:**
```python
def _apply_identifications(uv_jittered: np.ndarray, tris: np.ndarray, domain: Domain
                           ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    For rect with cy/mo on either axis: compute the vertex equivalence classes induced
    by the side-identification rules. Does NOT rewrite indices in `tris` — pre-id labels
    are preserved (G13). Edge/face identification is the responsibility of C5 and is
    determined geometrically, not by canonical labels.

    Returns (uv_unchanged, tris_filtered, vertex_class), where:
      - uv_unchanged: input uv (pass-through).
      - tris_filtered: input tris with any triangle removed whose three vertices map to
        fewer than three distinct equivalence classes (such a triangle is geometrically
        degenerate post-id — two corners coincide). For non-mo identifications this is
        almost always a no-op; included for safety on pathological corner configurations.
        Pre-id labels are preserved on the surviving rows.
      - vertex_class: (N,) int32. `vertex_class[i]` is the canonical (smallest pre-id
        index) representative of `i`'s equivalence class. For vertices not on any
        identified side, `vertex_class[i] == i`.

    Vertex pairing follows C1's deterministic boundary indexing (NOT current jittered uvs):
    cy uses same coordinate on the partner side; mo uses mirrored.
    """
```
**Algorithm:**
1. If not a rect or no identifications: return `(uv, tris, np.arange(N))`.
2. Build pair list from C1 boundary indexing (same logic as before):
   - u-cy: `(0, n)` and `(3n+i, 2n−i)` for i in 0..n−1
   - u-mo: `(0, 2n)` and `(3n+i, n+i)` for i in 0..n−1
   - v-cy: `(n, 2n)` and `(i, 3n−i)` for i in 0..n−1
   - v-mo: `(n, 3n)` and `(i, 2n+i)` for i in 0..n−1
3. Union-find with path halving; smaller index wins as canonical. Compute `vertex_class[i] = find(i)`.
4. Filter degenerate triangles: keep row `t` iff `len(set(vertex_class[t])) == 3`.
5. **Do NOT relabel `tris`**, do NOT dedup by sorted-tuple. Edge/face identity is geometric (G13).

**Why no dedup and no relabel (regression of the old C4 mo-mo bug):** Under mo-mo at corners, opposite-corner pre-id triangles can have all three vertices pairwise identified, yielding identical canonical triples for distinct cells. The old `np.unique(sorted(new_tris))` step silently deleted one such triangle, leaving a non-manifold complex (rect mo-mo res=5: 2 boundary edges, 1 edge with 3 faces). The fix is to abandon label-based identity for edges and faces entirely; C5 builds them by the geometric rule.

**Test criterion:**
1. Rect no-no: `tris_filtered == tris`; `vertex_class == arange(N)`.
2. Rect cy-no `resolution=10`: pairs (0, n), (3n+i, 2n−i) merge into the same class; all other vertices stay singleton. Each class has size 1 (interior) or 2 (a u-side pair).
3. Rect mo-mo `resolution=5`: corner classes are exactly `{0, 2n} = {0, 10}` and `{n, 3n} = {5, 15}`; the four corners collapse to 2 classes (not 1 — verifying mo-mo is RP² with 2 corner classes). No degenerate triangles removed at this resolution (regression for the bug: pre-id triangle (5, 6, 4) maps to canonical classes {5, 6, 4} which are 3 distinct classes; pre-id triangle (14, 15, 16) also maps to {4, 5, 6} — same canonical labels, but **both must remain** in `tris_filtered`).
4. Rect cy-cy `resolution=10`: each class has size at most 4 (interior corners merge 4 boundary vertices); `tris_filtered == tris` (no degeneracies).
5. Disk: `vertex_class == arange(N)`, `tris_filtered == tris`.
6. Pairing under jitter: jitter perturbs uvs but does NOT change `vertex_class` — pairing is by C1's deterministic indexing.
7. Determinism: `vertex_class` is reproducible across runs for the same `(domain, resolution)`.

**Gotchas:** G4 order (jitter before C4); G13 (vertex equivalence only — do not touch `tris` indices).
**Stop & verify:** `pytest surface_play/test_mesh.py::test_apply_identifications` ⇒ green.

---

### C5 — `mesh.py`: `_build_edges_faces`

**Spec:** §"Data stored" (lines 103-119), §"Mesh" lines 162-176.
**Files:** extend [surface_play/mesh.py](surface_play/mesh.py); add tests.
**API:**
```python
edge_dtype = np.dtype([
    ("p_idx", "i4"), ("q_idx", "i4"),  # pre-id vertex indices of the canonical edge (G13)
    ("p",  "2f8"),                     # uv[p_idx] — pre-id, direct, no close()
    ("pq", "2f8"),                     # uv[q_idx] − uv[p_idx] — pre-id subtraction, no close()
    ("f", "i4"), ("g", "i4"),          # adjacent face indices into the faces array; g = -1 if boundary
    ("dir", "2f8"),                    # inward normal in uv (boundary edges only); zeros otherwise
    ("flip", "i1"),                    # +1 or -1: orientation across the edge (mo-aware)
    ("split1", "i4"), ("split2", "i4"), # SPT indices (G17); -1 = no split. Populated in Layer O.
])

face_dtype = np.dtype([
    ("verts", "3i4"),                  # pre-id vertex indices of this triangle (in tris winding order)
    ("edges", "3i4"),                  # edge indices into the edges array
    ("p",  "2f8"),                     # uv[verts[0]] — pre-id
    ("pq", "2f8"),                     # uv[verts[1]] − uv[verts[0]] — pre-id, no close()
    ("pr", "2f8"),                     # uv[verts[2]] − uv[verts[0]] — pre-id, no close()
])

def _build_edges_faces(uv: np.ndarray, tris: np.ndarray,
                       vertex_class: np.ndarray,
                       SN_per_vertex: np.ndarray, domain: Domain
                       ) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns (edges, faces). `tris` carries PRE-id vertex labels (from C4). Edges merge
    across identified sides by the geometric rule of G13. Faces are never merged.
    `vertex_class` is consulted to identify which pre-id edges are paired (both endpoints
    on identified sides, with C1 pairings matching). `SN_per_vertex[i]` is the surface
    normal at pre-id vertex `i`.
    """
```
**Algorithm:**

Notation: a pre-id boundary vertex's *side* is one of `{u_min, u_max, v_min, v_max, interior}`, determined from its C1 index range (using `_boundary_edge_count(tris)` to recover `n`). A pre-id edge's *side* is `S` if **both** endpoints have side `S`; otherwise the edge is `interior` (this includes edges with one endpoint on the boundary and one in the interior).

1. **Classify each pre-id vertex** by side from its C1 boundary index. Interior vertices (index ≥ `4n`) get side `interior`.
2. **Iterate over every (face, local-edge) pair**, producing a record `(a, b, third, face_idx)` where `(a, b)` is the unordered pre-id edge and `third` is the third vertex of the face. Determine the edge's side from step 1.
3. **Compute each edge's canonical key:**
   - If the edge's side `S` is `interior`, OR `S` is on an unidentified axis (e.g., u-min/u-max with `u_identify == "no"`): key = `tuple(sorted([a, b]))` — no merging.
   - Else (`S` is on an identified axis): map both `a` and `b` to their canonical (`vertex_class[a]`, `vertex_class[b]`), then key = `tuple(sorted([vertex_class[a], vertex_class[b]]))`. This collapses the pre-id edge on side `S` with its partner on the paired side. **Sanity check:** assert each endpoint actually pairs to a vertex on the partner side; if a pre-id "side edge" has an endpoint whose canonical lives elsewhere, the pairing is malformed.
4. **For each canonical key**, collect the list of (face_idx, third, oriented endpoints). Length 1 ⇒ boundary; length 2 ⇒ interior. Length ≥ 3 ⇒ raise — indicates a malformed pairing.
5. **Pick a canonical pre-id edge per key:** when two pre-id edges merge, the canonical is the one on the lower-index identified side (u-min over u-max; v-min over v-max). The face containing this canonical edge becomes `f`; the other becomes `g`. `(p_idx, q_idx)` are the canonical pre-id endpoints in the canonical face's winding order around the edge.
6. **Edge geometric fields:** `p = uv[p_idx]`, `pq = uv[q_idx] − uv[p_idx]`. Pre-id subtraction; no `close()`.
7. **Boundary `dir`:** for an edge with `g = −1`, let `r = uv[third]` from the single adjacent face `f`. Inward normal in uv: orthogonal to `pq`, oriented so `dot(dir, r − p) > 0`. (`r` and `p` are both in `f`'s fundamental-domain copy, so no `close()` needed.) Interior edges: `dir = (0, 0)`.
8. **`flip`:**
   - For unmerged edges (interior or unidentified-side): `flip = +1 if dot(SN[p_idx], SN[q_idx]) > 0 else −1`. Will be `+1` everywhere unless the surface SN is itself near-degenerate.
   - For merged side edges: let `(a, b)` be the canonical pre-id endpoints (on side `S`) and `(a', b')` their partners (on the paired side, with `vertex_class[a'] = vertex_class[a]`, etc.). `flip = +1 if dot(SN[a], SN[a']) > 0 else −1`. Captures mo orientation reversal at the seam (cy: same orientation ⇒ `+1`; mo: reversed ⇒ `−1` on a Möbius-band SN field).
9. **Faces:** one face per row of `tris`. `verts = tris[f_idx]` (pre-id labels). `p = uv[verts[0]]`, `pq = uv[verts[1]] − uv[verts[0]]`, `pr = uv[verts[2]] − uv[verts[0]]`. `edges` = the 3 canonical edge indices via the dict from step 4.
10. **Split slots (G17):** `split1 = split2 = -1` for every edge.

**Test criterion:**
1. **Edge counts:** boundary count (`g == −1`) at resolution `R` matches the topological expectation per identification:
   - no-no: `4R`     (cyl: `2R`, torus: 0)
   - cy-no / mo-no: `2R`
   - cy-cy / cy-mo / mo-cy / mo-mo: `0`  ← critical regression for the old mo-mo bug
2. **Face count == len(tris):** every pre-id triangle survives as its own post-id face.
3. **Euler χ matches topology** for every (u_id, v_id) combination at `resolution=5` and `resolution=10`: no-no → 1; cy-no/mo-no → 0 (cylinder/Möbius band); cy-cy → 0 (torus); cy-mo/mo-cy → 0 (Klein bottle); **mo-mo → 1 (RP²)**.
4. **Manifoldness regression (mo-mo):** at `resolution ∈ {5, 7, 10}`, every edge in the post-id mesh has exactly 2 incident faces (`f, g ≥ 0`) and zero boundary edges. Histogram of face-counts-per-edge is `{2: E}` only.
5. **`flip`:** on rect cy-no with a smooth SN (e.g., `SN = (0,0,1)` constant): every edge has `flip = +1`. On rect mo-no with a Möbius-band-like SN that reverses sign across the u-seam (e.g., `SN(u) = (cos(π·(u−u_min)/period_u), 0, sin(...))`): every merged u-side edge has `flip = −1`, every non-side edge has `flip = +1`.
6. **`dir`:** for every boundary edge, `dot(dir, pq) ≈ 0` (machine precision) and `dot(dir, uv[third] − uv[p_idx]) > 0`.
7. **`p` and `pq` are pre-id subtractions:** for every edge, `e.p + e.pq == uv[e.q_idx]` exactly (no float drift from `close()`). Same for face `pq`, `pr`.
8. **Splits initialized:** all `edges["split1"] == edges["split2"] == -1`.
9. **No `close()` calls in the body of `_build_edges_faces`:** static-source assertion (`"close(" not in inspect.getsource(_build_edges_faces)`) — keeps G13 honest.

**Gotchas:** G13, G17.
**Stop & verify:** `pytest surface_play/test_mesh.py::test_build_edges_faces` ⇒ green.

---

### C6 — `mesh.py`: `Mesh` dataclass + `build_mesh` orchestrator

**Spec:** consolidates spec §"Mesh" (lines 162-176).
**Files:** extend [surface_play/mesh.py](surface_play/mesh.py); add integration test.
**API:**
```python
@dataclass
class Mesh:
    domain: Domain
    surface: SurfaceParams
    uv: np.ndarray                 # (N, 2) — pre-id, jittered
    tris: np.ndarray               # (M, 3) — pre-id vertex labels (G13)
    vertex_class: np.ndarray       # (N,)  — canonical class for each pre-id vertex
    edges: np.ndarray              # structured array, edge_dtype
    faces: np.ndarray              # structured array, face_dtype
    SN: np.ndarray                 # (N, 3) — per-vertex surface normal (at pre-id uv)
    xyz: np.ndarray                # (N, 3) — per-vertex 3D position (at pre-id uv)
    boundary_edge_idx: np.ndarray  # indices into edges where g == -1
    corner_idx: np.ndarray         # indices of corner vertex classes (rect only)

def build_mesh(domain: Domain, surface: SurfaceParams, resolution: int,
               *, jitter: bool = True, seed: int | None = None) -> Mesh:
    ...
```
**Algorithm:**
1. Generate raw mesh via C1 or C2 dispatch on `domain.type` → `(uv_raw, tris_raw)`.
2. **Jitter (C3)** → `uv` (jittered, pre-id). Skip if `jitter=False`.
3. **Apply identifications (C4)** → `(uv_unchanged, tris_filtered, vertex_class)`. `tris` retains pre-id labels (G13).
4. Evaluate per-vertex `xyz = surface.S(uv[:,0], uv[:,1]).T` and `SN` at the pre-id uvs via the cse'd lambdified callable (G11, G12). Each pre-id vertex (including both members of an identified pair) gets its own evaluation; mo orientation reversal between paired vertices is captured by `edge.flip` in step 5.
5. **Build edges/faces (C5)** with `(uv, tris, vertex_class, SN, domain)` — geometric merge rule (G13). No `close()` in per-element fields.
6. Identify corner vertex classes:
   - rect no-id: 4 distinct corner vertices (classes are singletons).
   - rect cy-no / mo-no etc.: corners merge into fewer classes per the C4 vertex pairing.
   - rect mo-mo: exactly 2 corner classes (regression — verifies the topology is RP² with 2 distinct corner points).
7. Wrap in `Mesh`.
**Test criterion (integration):**
1. Helicoid rect-cy on `(0, 1) × (0, 2π)`, `resolution=15`: builds cleanly, edges["flip"] == +1 everywhere on u-seam (cy = orientation preserving).
2. Möbius band rect-mo-no, `resolution=15`: edges["flip"] == −1 on the merged u-seam edges; `boundary_edge_idx` covers exactly the v-min and v-max sides (count `= 2·15`).
3. Disk `r_max=1`, `resolution=20`: 1 boundary loop, no seams.
4. Annulus `r_min=0.3, r_max=1`, `resolution=20`: 2 boundary loops (inner + outer).
5. Fig-8 immersion `(2 + cos(u))cos(v), (2 + cos(u))sin(v), sin(u·2)` on `(0, 2π) × (0, 2π)` cy-cy: mesh builds cleanly with χ = 0 (torus topology); to be reused as the SIC fixture in C8-C10.
6. **mo-mo regression:** rect mo-mo `resolution=10` with any smooth surface: `len(edges)` agrees with χ = 1 (RP²) given `V = N − |paired_classes|/2`-style counting; zero boundary edges; every edge has 2 adjacent faces (manifold).
**Gotchas:** G4, G13.
**Stop & verify:** `pytest surface_play/test_mesh.py::test_build_mesh` ⇒ green.

---

### C7 — `curves.py`: `build_bcs`

**Spec:** §"Boundary curves (BC)" (lines 183-189).
**Files:** extend [surface_play/curves.py](surface_play/curves.py); extend [surface_play/test_curves.py](surface_play/test_curves.py).
**API:**
```python
@dataclass
class BoundaryCurve:
    edge_indices: np.ndarray   # 1D int array; signs encode reversal as in make_lines output
    is_closed: bool

def build_bcs(mesh: Mesh) -> list[BoundaryCurve]:
    ...
```
**Algorithm:** Take `mesh.edges[mesh.boundary_edge_idx]`; assemble a `(K, 2)` array of vertex indices, **mapped through `mesh.vertex_class`** so that mo-identified corners are recognized as the same chain node (e.g., the v-min and v-max boundary segments of a rect-mo-no must connect at the shared corner class). Pass to `make_lines`. Wrap each chain as `BoundaryCurve` with `is_closed = (chain[0] == chain[-1])` (in segment-graph terms: the chain returns to its starting vertex class).
**Test criterion:**
1. Rect no-id: 1 BC, **closed** (per G16; the four-arcs-from-corners come later).
2. Rect cy-no: 2 BCs, both closed (top and bottom rims).
3. Rect cy-cy: 0 BCs.
4. Rect mo-no: 1 BC, closed (the two v-edges merge into one loop because mo links them).
5. Rect mo-mo: 0 BCs.
6. Disk: 1 closed BC.
7. Annulus: 2 closed BCs.
**Gotchas:** G16.
**Stop & verify:** `pytest surface_play/test_curves.py::test_build_bcs` ⇒ 7 green.

---

### C8 — `intersections.py`: BVH broad phase

**Spec:** G14 (this roadmap), legacy `_numba_bvh_final_boss` in [old stuff/intersections_prev.py](old%20stuff/intersections_prev.py).
**Files:** extend [surface_play/intersections.py](surface_play/intersections.py); extend [surface_play/test_intersections.py](surface_play/test_intersections.py).
**API:**
```python
def candidate_pairs(e_bbox: np.ndarray, f_bbox: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    e_bbox, f_bbox: shape (N, 6) — [xmin, ymin, zmin, xmax, ymax, zmax] per item.
    Returns (e_idx, f_idx) int64 arrays of overlapping AABB pairs. Deduped.
    Internal: stack-based BVH with axis-cycling (depth % 3) partition; brute-force at
    leaves where n_e * n_f < 5000 or depth > 15. Numba-JIT compiled.
    """
```
**Algorithm:** Port from legacy `_numba_bvh_final_boss` + `_candidate_pairs` verbatim, with naming aligned to the new module conventions. Pre-allocate `out_e`, `out_f` of size `15 · n_e`; manual stack of depth 64. Dedup output via combined-key sort + uniqueness mask (legacy line 224-231).
**Test criterion:**
1. Empty inputs → empty outputs.
2. Two non-overlapping AABBs → no pair returned.
3. Two overlapping AABBs → one pair `(0, 0)`.
4. **Brute-force equivalence:** on `(N=200, M=300)` random AABBs with overlap density ~5%, BVH output equals brute-force output (set equality of pairs).
5. **Performance:** `(N=10_000, M=10_000)` random AABBs runs in `< 1s` after JIT warmup; brute-force baseline takes `> 30s`. Assert ratio.
6. Determinism: same input → same output (no RNG inside BVH).
**Gotchas:** G14.
**Stop & verify:** `pytest surface_play/test_intersections.py::test_candidate_pairs` ⇒ 6 green.

---

### C9 — `intersections.py`: `find_double_points`

**Spec:** §"Double points" (lines 201-209), §Appendix (lines 491-503), G2 + G14 + G15. Algorithmic reference: legacy `find_double_points` in [old stuff/intersections_prev.py](old%20stuff/intersections_prev.py).
**Files:** extend [surface_play/intersections.py](surface_play/intersections.py); extend [surface_play/test_intersections.py](surface_play/test_intersections.py).
**API (matches legacy):**
```python
dp_dtype = np.dtype([
    ("xyz",         "f8", 3),
    ("uv1",         "f8", 2),    # domain coords on sheet 1
    ("uv2",         "f8", 2),    # domain coords on sheet 2
    ("E1",          "i4"),       # edge involved on sheet 1
    ("E2",          "i4"),       # edge on sheet 2 (-1 if "EF" type)
    ("F2",          "i4"),       # face on sheet 2 (-1 if "EE" type)
    ("A1",          "i4", 2),    # 1 or 2 face indices adjacent to E1; -1 padded
    ("A2",          "i4", 2),    # 1 or 2 face indices on sheet 2; -1 padded
    ("type",        "U2"),       # "EF" or "EE"
    ("ptype",       "i4"),       # 0 (interior) or 2 (boundary)
    ("on_boundary", "?"),        # bool
])

def find_double_points(mesh: Mesh, surface: SurfaceParams,
                       *, delta: float = 1e-6) -> np.ndarray:
    ...
```
**Algorithm:**
1. Build per-edge and per-face 3D AABBs from `mesh.xyz`.
2. **C8 BVH** generates candidate `(edge_idx, face_idx)` pairs.
3. Skip pairs sharing a vertex **class** (vectorized: any of `mesh.vertex_class[ep] == mesh.vertex_class[fv[0/1/2]]` or `mesh.vertex_class[eq] == mesh.vertex_class[fv[0/1/2]]`). Under the new C4 contract, pre-id labels are preserved in `tris`, so an edge on the u-min seam and a face on the u-max seam will not share a raw pre-id index — they share a `vertex_class`. Comparing by class avoids the spurious-DP-at-seam false positives.
4. **Batched Möller-Trumbore** over surviving candidates (legacy `moller_trumbore_batch`): vectorized solve for `(s, u, v)`. Filter by `s ∈ (0, 1)`, `u, v ≥ 0`, `u + v ≤ 1`.
5. **Per-hit Appendix classification (G2):** sort hits by descending `min_bary`. Process serially with a `consumed` mask and an `ee_keys` set. Three outcomes per legacy:
   - `min_bary ≥ delta` → emit `"EF"` record.
   - `min_bary < delta` and sibling `(E2, F1)` or `(E2, F1')` exists with `min_bary ≥ delta` at the opposite-vertex slot → emit second `"EF"` for that sibling.
   - Otherwise → emit `"EE"` record (deduped by `(min(E1, E2), max(E1, E2))`).
   - Mark `(E1, F2)`, `(E1, F2')`, `(E2, F1)`, `(E2, F1')` siblings consumed.
6. For each emitted DP, recover `uv1, uv2` (G15):
   - On the edge side: `uv1 = e.p + s · e.pq` where `e` is the edge and `s` is the MT edge parameter. Because C5 builds `(p, pq)` from pre-id labels of one face's fundamental-domain copy of the edge (G13), the result already lies in the canonical face's fundamental-domain — no `close()` call needed.
   - On the face side: `uv2 = f.p + u · f.pq + v · f.pr` (or analogous on the EE-paired edge). Same reason.
7. Populate `A1, A2` from `mesh.edges[E1].f, .g` (with `-1` padding) and analogously for E2/F2.
**Test criterion:**
1. Helicoid rect-no-no: empty DP array.
2. Möbius band: empty DP array.
3. **Fig-8 immersion** rect-cy-cy res=20: DP count consistent across runs (anchor a value once observed and lock as regression).
4. **Sibling-pair consumption** (G2 regression): no two DPs share `(E1, F2)` or `(min(E1,E2), max(E1,E2))`.
5. **No fragmentation** at fig-8 cyl res ∈ {20, 30, 40}: DP count grows smoothly with resolution (no factor-of-2 jumps).
6. `uv1, uv2` lie within domain bounds (with tolerance for jitter).
7. **Performance:** fig-8 cyl res=30 completes in `< 5s` after JIT warmup.
**Gotchas:** G2, G3, G14, G15.
**Stop & verify:** `pytest surface_play/test_intersections.py::test_find_double_points` ⇒ 7 green.

---

### C10 — `intersections.py`: `build_sis_pairs`

**Spec:** §"Self-intersecting segments" (lines 211-213).
**Files:** extend [surface_play/intersections.py](surface_play/intersections.py).
**API (vectorized, structured-array form per spec):**
```python
sis_dtype = np.dtype([
    ("p_dp", "i4"), ("q_dp", "i4"),  # DP indices (endpoints)
    ("flip", "i1"),                  # +1 if A1∩A1' and A2∩A2' match, -1 if crossed
    ("split1", "i4"), ("split2", "i4"), # SPT indices — spec line 245; -1 = no split. Populated in Layer O.
])

def build_sis_pairs(dps: np.ndarray) -> np.ndarray:
    ...
```
**Algorithm:** Vectorized over all `(i, j)` with `i < j`. For each pair, compute:
- `unflipped = (A1[i] ∩ A1[j] ≠ ∅) and (A2[i] ∩ A2[j] ≠ ∅)`
- `flipped   = (A1[i] ∩ A2[j] ≠ ∅) and (A2[i] ∩ A1[j] ≠ ∅)`
Set membership uses the legacy trick (legacy line 347-348): `(a[:, :, None] == b[:, None, :]) & (a >= 0)).any(axis=(1, 2))` over the fixed-size `(2,)` A1/A2 arrays with -1 padding.
- If exactly `unflipped` → emit `flip=+1`.
- If exactly `flipped` → emit `flip=-1`.
- If both true → log warning, emit `flip=+1` (defensive default).
- If neither → no pair.
- Initialize `split1 = split2 = -1` for every emitted SIS (G17). Populated in Layer O.
**Test criterion:**
1. Empty `dps` → empty array.
2. Two DPs with `A1[0] = A1[1]` and `A2[0] = A2[1]` → one pair, `flip=+1`.
3. Two DPs with `A1[0] = A2[1]` and `A2[0] = A1[1]` (orientation-reversing): emit `flip=-1`.
4. Fig-8 immersion cyl-cyl res=20: SIS count consistent with DP count (= `dp_count − 1` for a single closed SIC).
5. No SIS connects DPs that share no faces.
6. **Splits initialized:** all `sis["split1"] == sis["split2"] == -1` (G17 — Layer O populates).
**Gotchas:** G17.
**Stop & verify:** `pytest surface_play/test_intersections.py::test_build_sis_pairs` ⇒ 6 green.

---

### C11 — `intersections.py`: `build_sics`

**Spec:** §"Self-intersecting curves" (lines 191-200).
**Files:** extend [surface_play/intersections.py](surface_play/intersections.py).
**API:**
```python
@dataclass
class SelfIntersectingCurve:
    sis_indices: np.ndarray   # 1D int array, sign-encoded as in make_lines output
    is_closed: bool

def build_sics(sis_pairs: np.ndarray) -> list[SelfIntersectingCurve]:
    ...
```
**Algorithm:** Build `(K, 2)` array of DP-index pairs from `sis_pairs[:, ('p_dp', 'q_dp')]`; pass to `make_lines`. Wrap each chain.
**Test criterion:**
1. Empty input → empty list.
2. Fig-8 immersion cyl-cyl at res=20: exactly 1 SIC, closed.
3. Fig-8 immersion cyl-cyl at res=30 and res=40: still exactly 1 SIC each (regression for the prior fragmentation bug).
4. SIS count in the chain matches the SIS array length (no leftovers in single-SIC case).
**Gotchas:** G16-style — the chain reflects the SIS topology, not any post-processing splits.
**Stop & verify:** `pytest surface_play/test_intersections.py::test_build_sics` ⇒ 4 green.

---

### C12 — `intersections.py`: `find_triple_points`

**Spec:** §"Triple points" (lines 299-336).
**Files:** extend [surface_play/intersections.py](surface_play/intersections.py).
**API:**
```python
tp_dtype = np.dtype([
    ("xyz", "3f8"),        # 3D triple-point position
    ("sis_indices", "3i4"), # the three SISs that meet at this TP (indices into sis_pairs)
    ("uv", "(3,2)f8"),     # domain preimage points P1 ∈ F1, P2 ∈ F2, P3 ∈ F3
    ("faces", "3i4"),      # F1, F2, F3
])

def find_triple_points(sis_pairs: np.ndarray, dps: np.ndarray,
                       mesh: Mesh, surface: SurfaceParams,
                       *, xyz_tol: float = 1e-3) -> np.ndarray:
    ...
```
Note: SICs are NOT a parameter. SIS-level data (`sis_pairs` with its `flip` flag, plus DP `uv1, uv2` and `A1, A2`) is enough to reconstruct each SIS's two preimage segments. Per spec line 322, splits are recorded per-SIS (not per-SIC), so the TP dtype only needs `sis_indices`. SIC membership of each SIS is derivable downstream if needed (`sics[k].sis_indices` lookup).

**Algorithm:**
1. For each SIS in `sis_pairs`, build its two preimage segments using `dps[p_dp].uv{1,2}`, `dps[q_dp].uv{1,2}`, and the SIS's `flip`:
   - `flip == +1`: preimage A = `(uv1[p_dp], uv1[q_dp])`, preimage B = `(uv2[p_dp], uv2[q_dp])`.
   - `flip == -1`: preimage A = `(uv1[p_dp], uv2[q_dp])`, preimage B = `(uv2[p_dp], uv1[q_dp])`.
   For each preimage segment, also record its `f_here` (the face it lies on, found as a shared element of the two endpoints' relevant `A` sets) and `f_other` (the other sheet's face).
2. Run `sweep_segments(...)` (P5) on the union of all preimage segments in self-mode with `domain` (G3 — close-aware).
3. Filter hits to those where `segs["f_here"][a] == segs["f_here"][b]` (a hit must occur on a single shared face — that's how two preimages meet at the same domain point).
4. Group hits by face-pair interlocking (legacy `find_triple_points` lines 425-441): three hits form a TP when their `f_other` values pair up as `(F3, F2), (F3, F1), (F2, F1)` for some `(F1, F2, F3)`.
5. Verify in 3D: `S(P1) ≈ S(P2) ≈ S(P3)` within `xyz_tol`. Discard if not.
6. Populate `sis_indices` from the SIS index of each of the three preimage segments contributing to the triple.

**Test criterion:**
1. Empty `sis_pairs` → empty TPs.
2. Fig-8 cyl: 0 TPs (single SIC, but more importantly: no three preimage segments meet at one domain point on a shared face).
3. Synthetic three-SIS fixture (handcrafted DPs + sis_pairs, no real surface needed): correct face-pair interlocking detected.
4. **"Immersion with triple point"** surface (to be added to DB): produces ≥1 TP; `sis_indices` reference valid SISs; 3D verification passes.
5. No spurious TPs on standard surfaces (helicoid, torus, paraboloid, Möbius).
6. 3D consistency: for every emitted TP, `‖S(P1) − S(P2)‖, ‖S(P1) − S(P3)‖, ‖S(P2) − S(P3)‖ < xyz_tol`.
**Gotchas:** G3.
**Stop & verify:** `pytest surface_play/test_intersections.py::test_find_triple_points` ⇒ 6 green (test 4 conditionally skipped if the fixture isn't in the DB yet).

---

### C13 — `construction.py`: `build_construction` orchestrator

**Spec:** §"The Construction pipeline" (lines 30-66).
**Files:** create [surface_play/construction.py](surface_play/construction.py), [surface_play/test_construction.py](surface_play/test_construction.py).
**API:**
```python
@dataclass
class ConstructionResult:
    mesh: Mesh
    bcs: list[BoundaryCurve]
    dps: np.ndarray
    sis_pairs: np.ndarray
    sics: list[SelfIntersectingCurve]
    tps: np.ndarray

def build_construction(domain: Domain, surface: SurfaceParams,
                       resolution: int, *, jitter: bool = True,
                       seed: int | None = None) -> ConstructionResult:
    ...
```
**Algorithm:** Sequential calls: `build_mesh` → `build_bcs` → `find_double_points` → `build_sis_pairs` → `build_sics` → `find_triple_points`. Pure data-flow; no side effects.
**Test criterion:**
1. Helicoid rect-no-no res=15: `bcs` non-empty, `dps`/`sics`/`tps` empty.
2. Möbius mo-no res=15: `bcs` 1 closed loop, `dps`/`sics`/`tps` empty.
3. Disk res=20: `bcs` 1 closed loop, others empty.
4. Fig-8 cyl res=20: 0 BCs, 1 SIC with multiple DPs and SISs, 0 TPs.
5. Determinism: same seed → identical `ConstructionResult` (deep-equality on numeric arrays, length-equality on lists).
**Gotchas:** all of G1-G16 are upstream; this step adds none.
**Stop & verify:** `pytest surface_play/test_construction.py::test_build_construction` ⇒ green.

---

### C14 — Probe scripts for construction artifacts

**Files:** create at repo root, alongside existing `diag_*.py`:
- `probe_mesh.py` — render mesh in (u, v) space (matplotlib triplot) for rect, rect-cy, rect-mo, disk, annulus.
- `probe_jitter.py` — overlay pre-jitter and post-jitter mesh; verify boundary reprojection.
- `probe_bcs.py` — render BCs colored by curve index over the mesh.
- `probe_bvh.py` — render edge AABBs and face AABBs in 2D projection; overlay BVH partitions; show candidate pairs as edges between AABB centroids.
- `probe_dps.py` — render mesh + DP scatter (uv space, both `uv1` and `uv2` per DP, color-coded by DP index) on fig-8 cyl.
- `probe_sics.py` — render mesh + DPs + SIC chains (as polylines connecting DPs) on fig-8 cyl.
- `probe_3d.py` — render the surface in 3D with DPs as scatter and SICs as polylines, viewed from a couple of angles.

These give visual sanity for each construction artifact without requiring the Django pipeline. Browser-side debug overlay (live, with rotation) is deferred to Layer W.

---

### Layer C integration test

**Files:** [surface_play/test_layer_c_integration.py](surface_play/test_layer_c_integration.py).
**Goal:** end-to-end on the full fixture set: helicoid, torus-cyl, paraboloid, Möbius band, disk-of-functions, **fig-8 immersion**.
- For each: call `build_construction(...)` and assert ConstructionResult shape (BC count, DP count, SIC count) matches expectation.
- Determinism check: re-run with same seed, expect identical result.
- Caching check: when re-invoked with the same args, the function may be wrapped in `functools.lru_cache` — assert it's a hit (timing or `cache_info()`).

---

## Layer C sign-off gate

Before Layer O is drafted:
1. All 14 step tests + integration test green.
2. User has run at least `probe_dps.py`, `probe_sics.py`, `probe_3d.py` on fig-8 cyl and visually confirmed the SIC traces the expected figure-8 curve.
3. User confirms `Mesh`, `ConstructionResult`, `BoundaryCurve`, `SelfIntersectingCurve` data structures match what the Outline layer should consume.

## Layer O — Outline (17 steps)

Per-viewpoint computation: contour curves and cusps, all six split categories from spec §11, helper curves bridging components, SubCurve assembly, resampling, and visibility. Runs on every rotation/zoom; **not cached**. Consumes a `ConstructionResult` (from Layer C) plus a `Projection` (P3); returns a `OutlineResult` (defined in W1).

Each step has the same shape: **Spec → Files → API → Algorithm → Test criterion → Gotchas → Stop & verify.**

The **fixture viewpoint convention** for tests: each test fixture surface in `test_fixtures.py` ships with one or two canonical `(I, J, eye)` configurations (e.g., `helicoid_ortho_view = ([1,0,0], [0,1,0], None)`, `figure8_persp_view = (..., eye=[0,0,5])`), reproducible across runs.

---

### Section 1 — Contour detection (O1-O4)

#### O1 — `contour.py`: `find_contour_points`

**Spec:** §"Contour points (CP)" (lines 222-228), §"Contour curves" (lines 216-220).
**Files:** create [surface_play/contour.py](surface_play/contour.py), create [surface_play/test_contour.py](surface_play/test_contour.py).
**API:**
```python
cp_dtype = np.dtype([
    ("e", "i4"),               # edge index
    ("s", "f8"),               # parametric (0, 1) on the edge
    ("uv", "2f8"),             # domain coord = e.p + s · e.pq
    ("xyz", "3f8"),            # 3D position
    ("d", "2f8"),              # PROJECTED ker_param: proj_vec(uv, kerdS(uv)) — 2D vector in view plane,
                               # oriented "into the surface" (toward the lit side). Used by VP detection (O4).
    ("ptype", "u1"),           # 0 = interior CP, 4 = on boundary edge (becomes a BCP later)
])

def find_contour_points(mesh: Mesh, projection: Projection,
                        *, use_newton: bool = True) -> np.ndarray:
    ...
```
**Algorithm:**
1. Compute per-vertex `dot_v = SN[v] · viewer_direction(xyz[v])` (vectorized; for ortho this is constant axis, for persp depends on xyz).
2. For each edge, mask `sign_changes(dot_v[p_idx], dot_v[q_idx], flip)` (P6) — true where `dot_p · dot_q · flip < 0`.
3. For surviving edges, batched Newton on `s ↦ projection.viewer_direction(S(e.p + s·e.pq)) · SN(e.p + s·e.pq)`. Newton derivative uses `SN`'s chain-rule via `dSN` (precomputed in P2) and projection's `axis` derivative for persp. **Fall back to `s = 0.5` if Newton escapes (0, 1)** (G9).
4. For each CP, compute `xyz = S(uv)` and `d = projection.proj_vec(uv, projection.kerdS(uv))` — the **projected** ker_param in view-plane 2D. Orient `d` consistently across all CPs (toward the lit side of the surface; e.g., flip to make the dot product with a fixed reference direction positive, or use the convention from the legacy `kerdS` orientation that yields `Z > 0` in 3D before projection).
5. Set `ptype = 4` if the edge is boundary (`g == -1`); else `ptype = 0`.
**Test criterion:**
1. Helicoid rect-no-no, ortho top-down: produces a known CP count consistent with resolution; all `s ∈ (0, 1)`.
2. Sphere (excluded fixture, but a synthetic regular patch): CP locations match analytic equator ± tolerance.
3. Newton convergence: with `use_newton=False` (midpoint only), CP positions agree with `use_newton=True` to within `L/10` (G9 sanity).
4. Boundary CPs: on rect-no-no, CPs with `ptype == 4` match the spec's "endpoints of CCs which are not closed" implementation note.
5. Möbius: `flip` correctly captures sign reversal at seam — sign-change mask doesn't double-detect the seam edge.
**Gotchas:** G9.
**Stop & verify:** `pytest surface_play/test_contour.py::test_find_contour_points` ⇒ 5 green.

---

#### O2 — `contour.py`: `build_contour_segments`

**Spec:** §"Contour segments" (lines 230-231).
**Files:** extend [surface_play/contour.py](surface_play/contour.py).
**API:**
```python
cs_dtype = np.dtype([
    ("p_cp", "i4"), ("q_cp", "i4"),  # CP indices (endpoints)
    ("face", "i4"),                  # the face on which both CPs sit
    ("split1", "i4"), ("split2", "i4"),  # SPT indices (G17); -1 initially
])

def build_contour_segments(cps: np.ndarray, mesh: Mesh) -> np.ndarray:
    ...
```
**Algorithm:** For each face, find CPs whose edge belongs to that face. By construction (sign-change on edges of a single face), at most 2 CPs per face. If exactly 2, emit a CS connecting them; record the face index. Initialize splits to -1.
**Test criterion:**
1. Helicoid: CS count matches expected (paired CPs).
2. Each CS's two CPs are on different edges of the same face.
3. No CS crosses through a vertex (CPs are strictly interior to edges).
4. Splits initialized to -1 (G17).
**Gotchas:** G17.
**Stop & verify:** `pytest surface_play/test_contour.py::test_build_contour_segments` ⇒ 4 green.

---

#### O3 — `contour.py`: `build_contour_curves`

**Spec:** §"Contour curves (CC)" (lines 216-220), §"Curves" make_lines convention (lines 174-181).
**Files:** extend [surface_play/contour.py](surface_play/contour.py).
**API:**
```python
@dataclass
class ContourCurve:
    cs_indices: np.ndarray   # 1D int array; sign-encoded as in make_lines output
    is_closed: bool

def build_contour_curves(css: np.ndarray, cps: np.ndarray) -> list[ContourCurve]:
    ...
```
**Algorithm:** Two CSs share a chain endpoint iff they share a CP index. Build `(K, 2)` array of `(p_cp, q_cp)` from `css`; pass to `make_lines` (P4). Wrap each chain. CCs are closed iff their first and last CP indices match.
**Test criterion:**
1. Empty input → empty list.
2. Helicoid ortho: known CC count; closed/open distribution matches the surface's contour topology at the test viewpoint.
3. Each CC's CPs satisfy adjacency on shared faces.
4. Boundary CPs (ptype=4) appear only as endpoints of open CCs (per spec line 282).
**Gotchas:** none new.
**Stop & verify:** `pytest surface_play/test_contour.py::test_build_contour_curves` ⇒ 4 green.

---

#### O4 — `contour.py`: `find_vps`

**Spec:** §"Cusp points (VPs)" (lines 339-378).
**Files:** extend [surface_play/contour.py](surface_play/contour.py).
**API:**
```python
vp_dtype = np.dtype([
    ("cs", "i4"),                # CS index where VP sits
    ("s", "f8"),                 # parametric position on the CS (0, 1)
    ("uv", "2f8"),               # domain coord
    ("xyz", "3f8"),              # 3D position
    ("vis_change", "i1"),        # ±1 (used by O10 split logic)
])

def find_vps(ccs: list[ContourCurve], css: np.ndarray, cps: np.ndarray,
             surface: SurfaceParams, projection: Projection,
             *, refine: bool = False) -> np.ndarray:
    ...
```
**Algorithm:**
1. **Initial detection (spec lines 343-361):** for each CC chain, take the `d` field at consecutive CPs along the chain. Sign change in `dot(d[i], d[i+1])` flags a VP between them. Initial VP location = midpoint (`s = 0.5` on the CS between the two CPs).
2. **Optional Newton-bisection refinement (spec lines 364-371, behind `refine=True`, maps to settings `NEWTON_CUSP`):** at each candidate, take the midpoint, run Newton on the line orthogonal to the CS through the midpoint to find a true CP nearby. Replace midpoint, splitting the CS into two halves; iterate on the half with a sign change. Terminate when high precision reached.
3. Compute `vis_change` per spec: `vis = +1 if dot(dS(uv, q − p), axis) > 0 else -1` (legacy line 357).
**Test criterion:**
1. No-cusp surface (sphere patch): empty VP array.
2. Surface with known cusps (e.g., Torus, has 4 cusps in the inital viewpoint of the DB): VP count matches expectation; positions on the CS within `L/10`.
3. With `refine=False` and `refine=True`: positions agree to within `L/10`.
4. `vis_change ∈ {-1, +1}` everywhere (G10).
5. Möbius: VPs may appear at orientation-reversal points; verify that `flip`-aware sign-change is used.
**Gotchas:** G9, G10.
**Stop & verify:** `pytest surface_play/test_contour.py::test_find_vps` ⇒ 5 green.

---

### Section 2 — Splitting setup (O5)

#### O5 — `splitting.py`: SP/SPT data structures + slot management

**Spec:** §"Split Points" (lines 236-237), §"Splits" (lines 239-240), §"Splitting of segments" (lines 242-245).
**Files:** create [surface_play/splitting.py](surface_play/splitting.py), [surface_play/test_splitting.py](surface_play/test_splitting.py).
**API:**
```python
sp_dtype = np.dtype([
    ("uv",   "2f8"),            # domain coord
    ("xyz",  "3f8"),            # 3D position
    ("xy",   "2f8"),            # projected 2D position
    ("type", "U3"),             # "cn" | "bcp" | "bdp" | "tp" | "vp" | "cdp" | "ha"
])

spt_dtype = np.dtype([
    ("sp_idx",   "i4"),         # index into SPs array
    ("bary",     "f8"),         # barycentric position on the segment, in [0, 1]
    ("vis_chge", "i1"),         # visibility change crossing in forward direction; ∈ {-1, 0, +1} (G10)
])

class SplitArrays:
    sps:  list   # appended-to; length == n_sps
    spts: list   # appended-to; length == n_spts

    def add_sp(self, uv, xyz, xy, sp_type) -> int:
        """Append SP, return index."""

    def add_spt(self, sp_idx, bary, vis_chge) -> int:
        """Append SPT, return index."""

    def attach_to_segment(self, segment_array, seg_idx, spt_idx):
        """Write spt_idx into the first available split slot (split1 or split2)
        of segment_array[seg_idx]. Asserts split2 == -1 before overwriting (G17)."""
```
**Algorithm:**
- Plain Python lists of dicts/structured tuples internally; convert to structured arrays at finalize time.
- `attach_to_segment` enforces G17: if both `split1` and `split2` are already set, raise `SplitSlotOverflowError` with a diagnostic of the segment, the SP types already attached, and the new SP type.
**Test criterion:**
1. `add_sp` returns sequential indices starting at 0.
2. `add_spt` returns sequential indices starting at 0.
3. `attach_to_segment` writes to `split1` first, then `split2`, then raises on the third try (G17).
4. The error message names the segment kind (BE/CS/SIS), index, and existing SP types.
**Gotchas:** G17.
**Stop & verify:** `pytest surface_play/test_splitting.py::test_sp_spt_primitives` ⇒ 4 green.

---

### Section 3 — Split rules (O6-O11)

Each step in this section adds SPs/SPTs for one category from spec §11. Layer C produced segments with `split1=split2=-1`; these steps populate them via `SplitArrays.attach_to_segment`. The splitting mechanism is to create the SP first, then an SPT for each segment being split at this SP. It follows that **All SPTs sharing the same SP can be identified unambiguously**.

#### O6 — `splitting.py`: corner splits on BCs (rect no-id only)

**Spec:** lines 47-48 (corners as `cn` SPs).
**API:**
```python
def split_bcs_at_corners(mesh: Mesh, bcs: list[BoundaryCurve],
                         splits: SplitArrays) -> None:
    """Mutates `mesh.edges` (writes split1/split2) and `splits` (appends SPs/SPTs)."""
```
**Algorithm:** For each `corner_idx` from `mesh.corner_idx` (only nonempty on rect no-id), create an SP, then find the two BEs incident to it. For each such BE, create one SPT with `bary=0` if the corner is the BE's `p_idx`, else `bary=1`, with `type="cn"`, `xyz = S(uv)`, `xy = projection.XY(xyz)`, `vis_chge=0`. Attach the SPT to the BE.
**Test criterion:**
1. Rect no-id `resolution=10`: 4 SPs created (one per corner); 8 SPTs (each corner attached to 2 BEs).
2. Each corner SP's BE attachments have `bary ∈ {0, 1}`.
3. Rect cy-no, mo-no, etc.: 0 corner SPs (no `corner_idx`).
4. Disk: 0 corner SPs.
**Gotchas:** G10 (`vis_chge=0` is required for corners), G17.
**Stop & verify:** `pytest surface_play/test_splitting.py::test_corner_splits` ⇒ 4 green.

---

#### O7 — `splitting.py`: BCP splits (BC × CC at boundary CP)

**Spec:** §"Splitting of BCs at Contour Points (CPs)" (lines 249-282).
**API:**
```python
def split_bcs_at_bcps(mesh: Mesh, bcs, ccs, css, cps, splits, surface, projection) -> None:
    ...
```
**Algorithm:** For each CP with `ptype == 4` (on boundary edge):

1. Create SP with `type="bcp"`, `xyz = cp.xyz`, `xy = projection.XY(cp.xyz)`.
1. Find the CS containing the CP (only one, since the CP is a CC endpoint per spec line 282).
3. SPT on CS: `bary = 0 or 1` per spec (depending on CC traversal direction); `vis_chge = 0`.
4. SPT on the boundary edge `cp.e`: `bary = cp.s`; `vis_chge = bvis_chge(e, s)` per spec lines 258-280 (G18).
5. Attach SPTs to BE and CS via `splits.attach_to_segment`.
**Test criterion:**
1. Helicoid ortho viewpoint with boundary CPs: SP count = number of CPs with `ptype=4`.
2. Each BCP creates 2 SPTs (one BE, one CS).
3. **`bvis_chge` regression:** synthetic configuration with known `kerdS`, `Tb`, `Np` orientations gives the spec's stipulated `vis ∈ {-1, 0, +1}` (G18).
4. Open CCs end at boundary CPs (consistency with spec line 282).
**Gotchas:** G10, G17, G18.
**Stop & verify:** `pytest surface_play/test_splitting.py::test_bcp_splits` ⇒ 4 green.

---

#### O8 — `splitting.py`: BDP splits (BC × SIC at boundary DP)

**Spec:** §"Splitting of BCs at Double Points (DPs)" (lines 284-296).
**API:**
```python
def split_bcs_at_bdps(mesh, bcs, sics, sis_pairs, dps, splits, surface, projection) -> None:
    ...
```
**Algorithm:** For each DP with `on_boundary == True`:
1. Create SP with `type="bdp"`.
2. Find the BE (or BEs) containing the DP (it is either E1, or E2, or both). Find the SIS containing the DP. There is only one since it is on the boundary.
3. SPT on SIS: `bary = 0 or 1` (per which end of the SIS the DP is); `vis_chge = 0`.
4. SPT on BE: For each of the boundary edges `E` (usually there will only be 1, exceptionally 2) `bary = position of DP on E`; `vis_chge` per spec cases 1 and 2 (lines 293-296). Case 1 = DP also on a second BE (other sheet's E1' is also boundary); case 2 = DP's other sheet is interior.
5. Attach SPTs.
**Test criterion:**
1. Surfaces with no BDPs (helicoid, torus, etc.): no `bdp` SPs created.
2. Surface with known BDP (e.g., fig-8): correct count.
3. `vis_chge ∈ {-1, 0, +1}` (G10).
4. Case 1 vs case 2 distinguished correctly (assert via mock DP configurations).
**Gotchas:** G10, G17.
**Stop & verify:** `pytest surface_play/test_splitting.py::test_bdp_splits` ⇒ 4 green.

---

#### O9 — `splitting.py`: TP splits (3 SISs at triple point)

**Spec:** §"Splitting of SICs at Triple Points (TPs)" (lines 299-336).
**API:**
```python
def split_sics_at_tps(tps, sis_pairs, dps, splits, surface, projection) -> None:
    ...
```
**Algorithm:** For each TP:
1. Create SP with `type="tp"`, `xyz = tp.xyz`, `xy = projection.XY(tp.xyz)`.
2. For each of the three SISs (`sis_indices`):
   - Compute `bary` from the TP's location on the SIS preimage segment in domain.
   - Compute `vis_chge` per spec lines 327-334: `vis = +1 if dot(T_Si, N) > 0 else -1` where `N = SN(P_l)` oriented so `Z(N) > 0`, `T_Si = dS(P_h, dir_along_SIS_h)`. (`P_h` is the preimage point of the SIS on its own face; `P_l` is the preimage on the third face — not the SIS's pair.)
   - Attach SPT to the SIS via `splits.attach_to_segment`.
**Test criterion:**
1. Surfaces with no TPs: no `tp` SPs.
2. "Immersion with triple point" surface (when DB fixture exists): exactly the expected number of TPs; each creates 3 SPTs (one per SIS).
3. `vis_chge ∈ {-1, +1}` (G10 — TP splits never produce 0).
4. Synthetic 3-SIS fixture: SPT directions match handcrafted ground truth.
**Gotchas:** G7, G10, G17.
**Stop & verify:** `pytest surface_play/test_splitting.py::test_tp_splits` ⇒ 4 green.

---

#### O10 — `splitting.py`: VP splits (CC at cusp)

**Spec:** §"Splitting of CCs at cusps" (lines 336-378).
**API:**
```python
def split_ccs_at_vps(vps, css, ccs, splits, surface, projection) -> None:
    ...
```
**Algorithm:** For each VP:
1. Create SP with `type="vp"`, `xyz=vp.xyz`, `xy=projection.XY(vp.xyz)`.
2. SPT on CS: `bary = vp.s`, `vis_chge = vp.vis_change` (already computed in O4).
3. Attach SPT to the CS.
**Test criterion:**
1. No-cusp surface: 0 VP SPs.
2. Cusp-bearing surface: VP SP count matches `len(vps)`.
3. `vis_chge ∈ {-1, +1}`.
4. Each VP creates exactly 1 SPT (only on its CS).
**Gotchas:** G10, G17.
**Stop & verify:** `pytest surface_play/test_splitting.py::test_vp_splits` ⇒ 4 green.

---

#### O11 — `splitting.py`: CDP splits (CS × SIS at domain intersection)

**Spec:** §"Splitting of SISs and CSs at their domain-intersections" (lines 380-393).
**API:**
```python
def split_at_cdps(mesh, css, sis_pairs, dps, splits, domain, surface, projection) -> None:
    ...
```
**Algorithm:**
1. Run `sweep_segments(...)` (P5) **cross mode** between CS preimage segments (uv0/uv1 from CS endpoints' `cps[*].uv`) and SIS preimage segments (one per SIS's two preimages). `domain` set for close-aware (G3).
2. For each hit (`uv`, `t_a`, `t_b`):
   - Create SP with `type="cdp"`, `xyz=S(uv)`, `xy=projection.XY(xyz)`.
   - SPT on CS: `bary = t_a`, `vis_chge` per spec line 388 (uses 3D tangent + other-sheet normal).
   - SPT on SIS: `bary = t_b`, `vis_chge` per spec line 392 (uses domain tangent + CS-normal-toward-front-sheet).
   - Attach SPTs to CS and SIS.
**Test criterion:**
1. No-intersection cases (helicoid, torus): 0 CDP SPs.
2. Fig-8 'immersion sans bord': CDP count matches expected (one per CS that crosses a SIS preimage segment in domain).
3. `vis_chge ∈ {-1, +1}` for both SPTs (G10).
4. **Close-aware regression** (G3): seam-crossing CS×SIS preimage intersections are not double-counted.
**Gotchas:** G3, G10, G17.
**Stop & verify:** `pytest surface_play/test_splitting.py::test_cdp_splits` ⇒ 4 green.

---

### Section 4 — Helpers + SubCurve assembly (O12-O13)

#### O12 — `helpers.py`: `build_helper_curves`

**Spec:** §"Construction of helper curves" (lines 394-424). Spec line 408: HC starts as a `SubCurve(kind="HC")` with **only 2 points** — `p0 = qi_sp_idx`, `pN = qj_sp_idx`, `internal = []`. Subdivision into intermediate samples is the job of O14 (resampling).
**Files:** create [surface_play/helpers.py](surface_play/helpers.py), [surface_play/test_helpers.py](surface_play/test_helpers.py).
**API:**
```python
def build_helper_curves(bcs, ccs, sics, css, sis_pairs, cps, dps, mesh,
                        splits: SplitArrays, projection, surface, domain
                        ) -> list[SubCurve]:
    """
    Returns a list of SubCurves with kind="HC", parent_idx=-1, internal=[].
    Each HC carries its two HA SP endpoints (start, end), is_closed=False,
    and vc_in/vc_out per spec lines 416-424. No subdivision; resampling fills it in O14.
    Mutates `splits` (appends HA SPs and SPTs).
    """
```
**Algorithm:**
1. **Component identification:** union-find over BC/CC/SIC sub-curves via SplitPoint identity (G5). Two curves are in the same component iff they share an SP via SPT. SIC two-sheet preimages count as one node automatically (both sheets encoded in the same chain).
2. **Greedy bridging** (Prim-like, spec lines 401-414):
   - Compute pairwise component distances `d(i, j)` in **domain** space (G19) — minimum over sample-cloud pairwise distances. For SICs, iterate both preimages. Use **original (pre-merge) rectangle coordinates** on identified rects (G19) — no `close()` for the distance metric.
   - Pick the pair with smallest distance; `(qi, qj)` = the achievers.
   - Create two SPs of `type="ha"` with `xyz = S(uv)`, `xy = projection.XY(xyz)`. **Dedup at same target sample (G6):** if a previously-emitted HA on this iteration argmin'd to the same sample of the same target SubCurve, share the SP.
   - Create two SPTs (`vis_chge=0`) — one on the segment of the curve `qi` is attached to, one on the segment of the curve `qj` is attached to. Attach via `splits.attach_to_segment`.
   - Compute `vc_in, vc_out` per spec lines 416-424 (case-by-case: BC → 0; CC → 0 or -1 by `dot(T, N)`; SIC → 0 or -1 by `inner(SN', axis)·inner(T, SN')`).
   - **Emit HC SubCurve:** `SubCurve(kind="HC", parent_idx=-1, is_closed=False, start=qi_sp_idx, end=qj_sp_idx, internal=[], vc_in, vc_out)`. Just the two SP endpoints — no `uv_path`, no intermediate samples (those come in O14).
   - Merge components.
3. Iterate until all components are joined.
**Test criterion:**
1. Single-component starting state: empty HC list.
2. Multi-component case (helicoid with disconnected BC + CC): exactly the right number of HCs to bridge them.
3. HA SPs have `type == "ha"`.
4. **Each HC SubCurve has `internal == []` and `is_closed == False`** (spec line 408 — N=2 at construction).
5. **HA dedup (G6):** synthetic test with two HCs argmin'ing to the same target sample of the same SubCurve produces a shared SP, not two coincident ones.
6. `vc_in, vc_out ∈ {-1, 0}` per spec.
7. **HC endpoints are Python-identity-shared SPs (G5):** an HC's `start`/`end` SP is the same Python object as the SPs the SPTs reference.
**Gotchas:** G5, G6, G17, G19.
**Stop & verify:** `pytest surface_play/test_helpers.py::test_build_helper_curves` ⇒ 7 green.

---

#### O13 — `splitting.py`: SubCurve assembly

**Spec:** §"Construction of Split Curves (SCs)" (lines 426-434).
**API:**
```python
@dataclass
class SubCurve:
    kind: Literal["BC", "CC", "SIC", "HC"]
    is_closed: bool
    start: int            # SP index (even closed curve start/end at SP)
    end:   int            # SP index (even closed curve start/end at SP)
    internal: list[tuple[int, float]]   # point chain , each point specified by [segment_idx, bary_coord] 
    vc_in:  int           # ∈ {-1, 0} per spec line 434
    vc_out: int           # ∈ {-1, 0}

def assemble_subcurves(bcs, ccs, sics, hcs, mesh, css, sis_pairs, splits) -> list[SubCurve]:
    ...
```
**Algorithm:**
1. For each BC/CC/SIC, walk the chain. Each segment's `split1, split2` points to SPTs, which point to SPs. Cut the chain at every SPT, producing one SubCurve per arc between consecutive SPTs.
2. For closed BC/CC/SIC : the curve will always contains an SPT. If it doesn't, an error is raised. The walk starts at the first encountered SPT,and goes on, past the last point of the curve to the beginning, until the starting SPT is encountered again. In this construction, the segment containing the starting SPT is cut in two pieces (except if the SPT's `bary`is `0` or `1`). One of the pieces belongs to the first sub-curve, and the other one to the last.  
3. For open parents: the first and last segment must contain SPTs at the two endpoints (corners or BCPs).
4. Compute `vc_in, vc_out` per spec line 434:
   - `vc_in = min(0, ±1 · spt.vis_chge)` — the sign depends on whether the initial segment is reversed (sign-encoded in chain).
   - `vc_out = min(0, ±1 · spt.vis_chge)` — analogous for last segment.
5. HCs from O12 are already SubCurves (HCs are split-at-construction per spec line 76 with `internal=[]`, `N=2`); just append them to the output list. No re-splitting.
6. For BCs, each point inherits 
**Test criterion:**
1. No-split BC (closed loop with no SPTs): 1 closed SubCurve.
2. BC with corner splits (rect no-id): 4 SubCurves per closed BC (one per corner-to-corner arc).
3. CC with cusps and BCPs: SubCurve count matches `splits + 1` for open CCs, `splits` for closed.
4. **Each SubCurve's start/end SP is the SAME PYTHON OBJECT as a neighbor's start/end** (G5): identity check passes for adjacent SubCurves.
5. `vc_in, vc_out ∈ {-1, 0}`.
**Gotchas:** G5, G10, G17.
**Stop & verify:** `pytest surface_play/test_splitting.py::test_assemble_subcurves` ⇒ 5 green.

---

### Section 5 — Resampling (O14)

#### O14 — `curves.py`: `resample_all`

**Spec:** §"Curve resampling" (lines 78-82, 437-449).
**Files:** extend [surface_play/curves.py](surface_play/curves.py), [surface_play/test_curves_resampling.py](surface_play/test_curves_resampling.py).
**API:**
```python
@dataclass
class ResampledCurve:
    kind:   Literal["BC", "CC", "SIC", "HC"] # kind of source Subcurve
    start:  int                   # SP index (into splits.sps); -1 for closed RCs with no SPs
    end:    int                   # SP index (into splits.sps); -1 for closed RCs with no SPs
    depth:  np.ndarray            # (N, 1) — depth of 3D point
    xy:     np.ndarray            # (N, 2) — projected
    dir:    np.ndarray            # (N, 2) — proj dir of surface sheet(s); None for SIC/HC kinds

def resample_all(subcurves: list[SubCurve], surface, projection,
                 resolution: int, *, project_resampled: bool = False) -> list[ResampledCurve]:
    ...
```
**Algorithm:**
1. Compute `M` = approximate diameter of surface in projected space (`max(‖xyz_i − xyz_j‖)` for vertices in mesh, projected via `projection.XY`).
2. For each SP: enumerate SubCurves originating/ending at it; compute each SubCurve's view-plane length `L_k`. Let `L = min(L_k)`. Compute `d = min(L/10, M / resolution)` (G20).
3. For each SubCurve: split at midpoint; on each half, sample at arclengths `0, d, 2d, ...` until `L/2 − d`. Closed curves: uniform `d` spacing.
4. The sampled points usually lie between two points of the original curve. Then the `(uv)` coordinates of these points are interpolated, after `close`, and the projected image of the result is taken as the `xy` of the interpolated point.
5. **SIC two-sheet tracking (G21):** for SIC SubCurves, one needs to use `(uv)` coordinates of the interpolants in **the same preimage**, which can be done using the `flip` field.
6. **Annular BC reprojection** (spec line 446): if the SubCurve's parent BC is on disk/annulus boundary then there is an additional step which is to snap the interpolated  `(u, v)` to `‖(u, v)‖ = r_min` or `r_max` before computing the projected image. Each interpolated point's `dir` field is taken to be the `dir` field of the boundary edge containing both interpolants.
7. **CC Newton refinement** (spec line 448, behind `project_resampled=True`, maps to settings `PROJECT_RESAMPLED`): If the SubCurve's parent is a CC, there is also an additional step which is to run Newton along the line normal to the CC at the interpolated `(u, v)` to find a nearby `dot(SN, axis) = 0`. The `dir` field is either computed as the `dir`field of CPs, or interpolated from the interpolant's `dir` field.
8. The `dir` field for SIC of HC resampled curves is not used, it can be set to `None`.
9. **Endpoint pinning (G5, G21):** `rc.start = sub.start` and `rc.end = sub.end` (same SP index — SP Python identity is preserved by reference equality through `splits.sps[rc.start]`). The first sample's `xy` exactly equals `projection.XY(splits.sps[rc.start].xyz)`; the last sample's `xy` exactly equals the corresponding for `rc.end`. Same for `depth`. (For closed SubCurves with `sub.start == -1`, the RC is also closed: `rc.start = rc.end = -1`.)
**Test criterion:**
1. SubCurve with one SP at each end: at least 4 sample points per half-curve, total ≥ 8 samples.
2. Closed SubCurve (`sub.start == -1`): uniform spacing; sample count ≈ `round(L / d)`; `rc.start == rc.end == -1`.
3. Annular BC: post-resampling, the projected sample positions correspond to uvs whose `‖(u, v)‖` is within `1e-12` of `r_min` or `r_max`.
4. **SIC two-sheet tracking (G21):** for a SIC SubCurve, the resampling loop does not jump preimages mid-walk (verified by an internal-uv tracker that stays close-aware between consecutive samples).
5. **Endpoint SP indices (G5):** `rc.start == sub.start` and `rc.end == sub.end`. Furthermore `splits.sps[rc.start]` is the same Python object as `splits.sps[sub.start]` (identity preserved through index dereferencing).
6. **Endpoint xy pinning:** `rc.xy[0] == projection.XY(splits.sps[rc.start].xyz)` exactly; same for `rc.xy[-1]` and `rc.end` (when `rc.start, rc.end != -1`).
**Gotchas:** G5, G20, G21.
**Stop & verify:** `pytest surface_play/test_curves_resampling.py::test_resample_all` ⇒ 6 green.

---

### Section 6 — Visibility (O15-O17)

#### O15 — `visibility.py`: `compute_projection_breaks`

**Spec:** §"Projection intersections" (lines 452-467).
**Files:** create [surface_play/visibility.py](surface_play/visibility.py), [surface_play/test_visibility.py](surface_play/test_visibility.py).
**API:**
```python
break_dtype = np.dtype([
    ("rc_idx",      "i4"),       # ResampledCurve being occluded
    ("sample_idx",  "i4"),       # sample in rc.xy preceding the BK
    ("t",           "f8"),       # parametric (0, 1) on the segment between sample_idx and sample_idx+1
    ("xy",          "2f8"),      # 2D projected position of the BK
    ("delta_v",     "i1"),       # ∈ {-2, -1, +1, +2} per spec line 461
])

def compute_projection_breaks(rcs: list[ResampledCurve], surface, projection
                              ) -> np.ndarray:
    ...
```
**Algorithm:**
1. Flatten all sample-to-sample segments of all RCs in projected space; record `(rc_idx, sample_idx)` per segment.
2. Run `sweep_segments(...)` (P5) **self-mode** with `domain=None` (flat 2D in view plane).
3. For each hit, evaluate `Z` depth at the crossing point on both segments. The OCCLUDED side is the one with greater `Z` (farther from viewer).
4. Compute `delta_v` magnitude per spec line 461: `2` if occluder is CC type, `1` if BC, `0` otherwise. Sign per spec lines 463-465: `-` if `dot(T, dir) > 0` where `T` is occluded's tangent and `dir` is the occluder's view-plane normal toward the surface.
5. Do NOT skip same-RC self-crossings; skip HC-as-occluder (HCs don't occlude — they're invisible bridging lines); skip SIC-as-occluder.
**Test criterion:**
1. Non-overlapping projected curves: empty BKs.
2. Single sphere-patch silhouette in front of an inner curve: BKs on the inner curve, magnitude 2.
3. Fig-8 cyl ortho: BK count matches expected; `delta_v` signs distribute as expected.
4. **Same-RC skipped:** RC self-intersection in projection doesn't emit a BK.
5. **HC skipped:** an HC crossing a BC produces no BK.
**Gotchas:** G22.
**Stop & verify:** `pytest surface_play/test_visibility.py::test_compute_projection_breaks` ⇒ 5 green.

---

#### O16 — `visibility.py`: `bfs_visibility`

**Spec:** §"Anchor and BFS propagation" (lines 469-481).
**API:**
```python
def bfs_visibility(rcs, breaks, splits, *, anchor_mode: Literal["leftmost", "extremes"] = "leftmost"
                   ) -> dict[int, np.ndarray]:
    """Returns {id(rc): vis_per_sample_array}."""
```
**Algorithm:**
1. **Build SP-to-RCs index:** for each RC with `rc.start != -1`, register `rc` (and which end) under `splits.sps[rc.start]` — the SP key is `id(splits.sps[rc.start])` so identity (G5) drives the lookup. Same for `rc.end`. Closed RCs (`start == end == -1`) have no SP anchors and are reachable only via projection breaks or are their own islands until O17 LP rescues them.
2. Find anchor sample(s): `anchor_mode="leftmost"` picks the sample with smallest `xy[0]` across all RCs, ties broken by `(rc_id, sample_idx)` (G23). `"extremes"` picks 4: leftmost, rightmost, topmost, bottommost.
3. Initialize `vis = 0` at anchor sample(s) of the anchor RC(s).
4. BFS frontier of `(rc, sample_idx)`:
   - **Within an RC:** propagate sample-by-sample, applying break `delta_v` at the crossings (cumulative).
   - **At an endpoint sample (sample 0 or last):** look up `splits.sps[rc.start]` (or `rc.end`) in the SP-to-RCs index. For each neighbor RC sharing this SP, `vis_neighbor_endpoint = vis_self_endpoint − delta_self + delta_neighbor` where `delta_self, delta_neighbor` come from the SPTs at this SP for each curve, applied in the appropriate direction relative to whether the neighbor RC begins or ends at this SP.
5. Return per-RC visibility arrays.
**Test criterion:**
1. Sphere ortho: anchor `vis = 0`, all samples `vis ≤ 0`.
2. Single closed BC, no breaks: all samples have `vis = 0`.
3. Synthetic break: anchor `vis = 0`; samples after the break have `vis = -2` (CC-type occluder).
4. Multi-curve through SPs: visibility propagates correctly across the SP interface (Python-identity-shared SP).
5. **Determinism (G23):** repeated runs produce identical `vis` arrays.
**Gotchas:** G5, G22, G23.
**Stop & verify:** `pytest surface_play/test_visibility.py::test_bfs_visibility` ⇒ 5 green.

---

#### O17 — `visibility.py`: `lp_refine_visibility`

**Spec:** §"LP pass" (lines 483-486).
**API:**
```python
def lp_refine_visibility(rcs, breaks, splits, vis_bfs,
                         *, anchors: Literal["leftmost", "extremes"] = "leftmost"
                         ) -> dict[int, np.ndarray]:
    ...
```
**Algorithm:** Build a sparse LP per G24:
- Variables: `vis[k]` per sample (≤ 0), `vc_i` per break (free), `s_i ≥ 0` (slack).
- Equality rows:
  - Per-segment propagation: `vis[k+1] − vis[k] − Σ vc_i (over breaks in segment) = 0`.
  - **SP coupling** at each SP shared by two or more RCs (looked up via the same SP-to-RCs index built in O16, keyed on `id(splits.sps[rc.start])` / `id(splits.sps[rc.end])`): `vis_o_endpoint − vis_ref_endpoint = delta_o − delta_ref`, where the SPT `vis_chge` values from each curve incident at that SP set the constants.
  - Anchor: `vis = 0` at anchor sample(s).
- Inequality rows: `vc_i − s_i ≤ vc_i⁰` and `−vc_i − s_i ≤ −vc_i⁰` for L1 slack on the objective.
- Objective: minimize `Σ s_i`.
- Solve via `scipy.optimize.linprog(method="highs")` on a sparse CSR matrix.
- Round `vis` to integers at exit.
- On status=2 (infeasible): raise `LPInfeasibleError` (G7).
**Test criterion:**
1. **No-correction case:** input visibility is consistent → LP returns it unchanged.
2. **Single-bad-break correction:** synthetic 3-sample open RC with a bogus `+3` break → LP corrects to all-zero.
3. **`anchors="extremes"`:** all 4 anchor samples receive `vis = 0`.
4. **Infeasibility (G7):** input with deliberately inconsistent SP δv values raises `LPInfeasibleError`.
5. **Sphere ortho:** LP path matches BFS path (no correction needed).
**Gotchas:** G7, G24.
**Stop & verify:** `pytest surface_play/test_visibility.py::test_lp_refine_visibility` ⇒ 5 green.

---

### Layer O integration test

**Files:** [surface_play/test_layer_o_integration.py](surface_play/test_layer_o_integration.py).
**Goal:** end-to-end on full fixture set with canonical viewpoints. For each fixture:
1. `build_construction(...)` once (Layer C).
2. For each canonical viewpoint:
   - Run all of O1-O17.
   - Assert: SubCurve counts ≥ 0, ResampledCurve sample counts > 0, BK counts plausible, visibility distribution (`vis ≤ 0` everywhere; `vis = 0` at anchor; min `vis` no worse than the maximum surface-sheet count for that surface).
3. Smoke: pass `vis_mode ∈ {"bfs", "lp1", "lp4"}`; all three return.

---

### Layer O probe scripts

**Files:** at repo root:
- `probe_cps.py` — render mesh + CP scatter for each fixture/viewpoint.
- `probe_ccs.py` — render mesh + CC chains.
- `probe_vps.py` — render CCs + VP markers.
- `probe_splits.py` — render all SPs colored by type (`cn, bcp, bdp, tp, vp, cdp, ha`).
- `probe_subcurves.py` — render SubCurves colored by `vc_in`/`vc_out`.
- `probe_resampled.py` — render resampled samples on top of curves.
- `probe_breaks.py` — render projected curves + BK markers.
- `probe_visibility.py` — render projected curves colored by visibility integer.

---

## Layer O sign-off gate

Before Layer W is drafted:
1. All 17 step tests + integration test green.
2. User has run at least `probe_visibility.py` on helicoid, fig-8 cyl, and Möbius and confirmed the visibility distribution is correct visually.
3. User confirms `OutlineResult` (or whatever pipeline.py exposes in W1) carries everything the existing frontend contract needs (`lines_by_visibility`, `si_lines_by_visibility`, `origin`).

## Layer W — Wiring (5 steps)

Plug the modular pipeline into Django, retire `silhouette.py`. **Layer C/O sign-off is a precondition**: this layer assumes `build_construction` and the full Layer O chain are stable. Each step preserves the existing views.py contract (this roadmap, lines 65-69) and the frontend wire format (`I, J, O, eye, debug`) bit-for-bit.

Each step has the same shape: **Spec → Files → API → Algorithm → Test criterion → Gotchas → Stop & verify.**

---

### W1 — `pipeline.py`: orchestrator

**Spec:** §"The Construction pipeline" (lines 30-66), §"The Outline pipeline" (referenced from Layer O), §"Debug panel" (this roadmap line 70).
**Files:** create [surface_play/pipeline.py](surface_play/pipeline.py), create [surface_play/test_pipeline.py](surface_play/test_pipeline.py).
**API:**
```python
@dataclass
class SurfaceInit:
    record_pk: int
    cache_key: str                 # hash of (X, Y, Z, parameter values, domain spec, GRID_RESOLUTION)
    domain: Domain
    surface: SurfaceParams
    construction: ConstructionResult
    threejs: dict                  # ready-to-serialize three.js mesh data (projection-independent)

@dataclass
class OutlineResult:
    lines_by_visibility:    dict[int, list[list[tuple[float, float]]]]   # vis → list of xy polylines
    si_lines_by_visibility: dict[int, list[list[tuple[float, float]]]]
    origin:                 tuple[float, float]                           # XY(world origin) under the view's (I, J, O)

def build_surface_init(record: SurfaceRecord, *,
                       resolution: int | None = None,
                       jitter: bool = True, seed: int | None = None) -> SurfaceInit: ...

def mesh_to_threejs(construction: ConstructionResult) -> dict: ...

def build_outline(init: SurfaceInit, I, J, O, eye, *,
                  newton_cusp:        bool = True,
                  canvas_resolution:  int | None = None,
                  project_resampled:  bool = False,
                  propagation:        Literal["BFS", "LP1", "LP4"] = "BFS",
                  ) -> OutlineResult: ...
```
**Algorithm:**
1. `build_surface_init`:
   - Parse `record` fields (X, Y, Z expressions, parameter values, domain spec) → build `Domain` (P1) and `SurfaceParams` (P2).
   - `resolution = resolution or settings.GRID_RESOLUTION`.
   - `cache_key = sha256(record.X + record.Y + record.Z + str(parameters) + domain_spec + str(resolution))`. Cache the full `SurfaceInit` keyed on `(record.pk, cache_key)` via `functools.lru_cache(maxsize=settings.SURFACE_CACHE_SIZE or 16)`. Cache invalidation is implicit: any change to record fields or `GRID_RESOLUTION` produces a new key, so stale entries are simply unused (LRU evicts).
   - Call `build_construction(domain, surface, resolution, jitter=jitter, seed=seed)` (C13).
   - Pre-compute `threejs` via `mesh_to_threejs(construction)` — projection-independent.
2. `mesh_to_threejs`:
   - Output shape: `{"vertices": [...], "faces": [[i,j,k], ...], "normals": [...], "uvs": [...]}` matching what the frontend currently consumes from `Surface.for_3js()`.
   - Pulls from `construction.mesh.xyz`, `.tris`, `.SN`, `.uv`.
   - **Vertex-class compaction (C4):** under the new C4 contract, `tris` carries pre-id labels and every row of `uv/xyz/SN` is live; paired vertices each have their own evaluation. For three.js, decide per `templates/play.html` whether to send all N pre-id vertices (paired vertices appear as duplicated geometry, harmless overdraw) or compact to one representative per `vertex_class` (rebuild a dense array via the `vertex_class` map). Document the choice in the docstring.
3. `build_outline`:
   - Build `Projection(surface, I, J, O, eye)` (P3). `O` is the world-space view-plane anchor — see W2 for the full semantics. Pass through unmodified.
   - Run Layer O end-to-end on `(init.construction, projection)` with the supplied debug knobs threaded through: `newton_cusp` → O1's `use_newton`; `canvas_resolution` → O14's `resolution` (defaults to `settings.CANVAS_RESOLUTION`); `project_resampled` → O14's `project_resampled`; `propagation` selects O16-only (BFS) vs. O17 LP1/LP4.
   - Assemble `OutlineResult`:
     - `lines_by_visibility[v]`: list of polylines (each polyline = list of `(x, y)` view-plane points) for each ResampledCurve segment whose visibility integer equals `v`. Segment boundaries split at projection breaks (BKs).
     - `si_lines_by_visibility[v]`: same for SIC ResampledCurves.
     - `origin = projection.XY([0, 0, 0])` — the view-plane image of the world origin under the supplied `(I, J, O)`.
   - HCs are excluded from output (G22 — invisible bridges).
**Test criterion:**
1. `build_surface_init` on a fixture record returns a `SurfaceInit` whose `domain.type`, `surface.bbox_diag`, and `construction` match a standalone-call output (no caching weirdness).
2. **Cache hit:** second call with same `(record, resolution)` returns the same Python object as the first (`is`-equality).
3. **Cache miss on resolution change:** different `resolution` produces a different object.
4. `mesh_to_threejs` output passes `json.dumps` round-trip; vertex/face counts match `Mesh`.
5. `build_outline` on helicoid ortho with `I=[1,0,0], J=[0,1,0], O=None, eye=None`: returns non-empty `lines_by_visibility`, no exceptions; all `vis ≤ 0`; `origin == (0.0, 0.0)` (world origin maps to view-plane origin when `O = [0,0,0]`).
6. **Debug-knob plumbing:** `propagation="LP1"` exercises O17; `project_resampled=True` propagates into O14 (verified via spy on the Layer O entry point).
7. **No silhouette import:** running `python -c "import surface_play.pipeline; assert 'silhouette' not in sys.modules"` in a subprocess succeeds.
**Gotchas:** none new (W1 is pure assembly).
**Stop & verify:** `pytest surface_play/test_pipeline.py` ⇒ 7 green.

---

### W2 — `views.py`: rewrite

**Spec:** §"views.py contract (preserved)" (this roadmap, lines 65-68); §"Debug panel" (this roadmap line 70).
**Files:** rewrite [surface_play/views.py](surface_play/views.py); add [surface_play/test_views.py](surface_play/test_views.py).
**API:**
```python
@require_GET
def play(request, pk: int) -> HttpResponse: ...
    # Renders templates/play.html with three.js mesh data injected.

@require_POST
def play_outline(request, pk: int) -> JsonResponse: ...
    # Body:    {"I":[...], "J":[...], "O":..., "eye":[...]|null, "debug":{...}}
    # Returns: {"lines_by_visibility":{...}, "si_lines_by_visibility":{...}, "origin":[...]}
```
**Algorithm:**
1. `play(request, pk)`:
   - `record = get_object_or_404(SurfaceRecord, pk=pk)`.
   - `init = pipeline.build_surface_init(record)`.
   - Render `templates/play.html` with `init.threejs` injected as JSON in the page context (replaces `Surface.for_3js()`).
2. `play_outline(request, pk)`:
   - `record = get_object_or_404(SurfaceRecord, pk=pk)`.
   - Parse JSON body → `(I, J, O, eye, debug)`. All four geometry fields are required; `eye` may be `null` (ortho mode).
   - `init = pipeline.build_surface_init(record)` (cache-hit on subsequent calls per W1).
   - `result = pipeline.build_outline(init, I, J, O, eye, **debug_kwargs(debug))` where `debug_kwargs` validates the `debug` dict and maps keys → kwargs: `PROPAGATION → propagation`, `CANVAS_RESOLUTION → canvas_resolution`, `NEWTON_CUSP → newton_cusp`, `PROJECT_RESAMPLED → project_resampled`. Unknown keys: log a warning, ignore (don't 400 — preserve existing client compat).
   - Return `JsonResponse({"lines_by_visibility": ..., "si_lines_by_visibility": ..., "origin": ...})`.
3. **`O` semantics** (resolved from the frontend at [templates/play.html:330-336](templates/play.html#L330)): `O` is the **world-space point that the frontend wants mapped to view-plane `(0, 0)`**.
   - In ortho mode: `O = unproject(NDC (0, 0, -1))` — the world-space point on the camera's near plane that sits at the screen origin. `I` and `J` are the in-plane direction vectors `unproject(NDC (1,0,-1)) − O` and `unproject(NDC (0,1,-1)) − O`.
   - In persp mode: `O = camera.position = eye`. So `O == eye` is invariant in persp; W2 must validate this and 400 if violated. Persp `XY` is anchored so `XY(eye) = (0, 0)`, consistent with `O = eye`.
   - The contract is consumed by `Projection` in P3 (which now accepts `O` as a constructor argument).
4. No silhouette imports anywhere in the module.
**Test criterion (Django test client):**
1. `GET /play/<pk>/` on a fixture record: 200; body contains expected three.js JSON keys.
2. `POST /play/<pk>/` with `{"I":[1,0,0],"J":[0,1,0],"O":[0,0,0],"eye":null,"debug":{}}`: 200; JSON has all three top-level keys; `lines_by_visibility` non-empty.
3. `POST` with `eye=[0,0,5], O=[0,0,5]`: persp result, no errors, plausible `origin`.
4. **Persp `O ≠ eye` rejected:** `POST` with `eye=[0,0,5], O=[0,0,0]` returns 400 (constructor `ValueError` from P3 surfaced as bad-request).
5. **Ortho `O` translates output:** two `POST`s with the same `(I, J)` but different `O` produce `lines_by_visibility` differing exactly by a constant 2D shift on every polyline (regression for the `−O` shift in P3's `XY`).
6. `POST` with `debug={"PROPAGATION":"LP4"}`: 200; differs from BFS only where the LP correction kicks in.
7. `POST` with `debug={"BOGUS":42}`: 200, warning logged (not a 400).
8. **Caching observable:** second `POST` to same pk faster than first, or `build_surface_init.cache_info().hits ≥ 1` after.
9. **404 on missing pk:** both endpoints return 404.
10. **Static check:** `grep -rn "silhouette" surface_play/views.py` empty.
**Gotchas:** none new (W2 is plumbing).
**Stop & verify:** `pytest surface_play/test_views.py` ⇒ 10 green; manually load `/play/<some_pk>/`, rotate the view, no console errors.

---

### W3 — `thumbnail.py`: migrate to pipeline

**Spec:** preserve current thumbnail file format and dimensions. The frontend wire format is fixed by [templates/play.html:694-708](templates/play.html#L694) — the "save thumbnail" button POSTs `{I, J}` (ortho) or `{I, J, eye}` (persp) to `surface-save-thumbnail`, the server renders a PNG at that viewpoint and stores it on the `SurfaceRecord`.
**Files:** rewrite [surface_play/thumbnail.py](surface_play/thumbnail.py); add [surface_play/test_thumbnail.py](surface_play/test_thumbnail.py); update the `surface-save-thumbnail` view (in W2's views.py rewrite) to call into it.
**API:**
```python
def render_thumbnail(record: SurfaceRecord,
                     I: list[float], J: list[float],
                     eye: list[float] | None = None,
                     *, size: tuple[int, int] = (256, 256)) -> bytes:
    """Render the supplied viewpoint to PNG bytes. Viewpoint is always supplied by the
    caller (the frontend's save-thumbnail payload); there is no default viewpoint."""
```
**Algorithm:**
1. `init = pipeline.build_surface_init(record)` — shares the LRU cache with the views.
2. Compute `O` from the frontend convention (the save-thumbnail payload omits `O`): `O = eye` if `eye is not None`, else `O = [0, 0, 0]`. This matches the implicit anchor the frontend uses for thumbnails.
3. `outline = pipeline.build_outline(init, I, J, O, eye, propagation="BFS")`.
4. Render to a `size[0] × size[1]` PNG: matplotlib `Figure(dpi=...)` with the same line widths, colors, and visibility-integer-to-alpha mapping the current thumbnails use. Plot:
   - For each `vis` group in `lines_by_visibility`: plot polylines with the conventional color/alpha for that visibility level.
   - For each `vis` group in `si_lines_by_visibility`: same with the SIC color.
   - Crop to surface bbox in view plane.
5. Return PNG bytes via `BytesIO`.
**Test criterion:**
1. Returns valid PNG bytes (`PIL.Image.open(BytesIO(...))` succeeds).
2. Dimensions match `size`.
3. Non-trivial pixel content (not all-white): pixel-value standard deviation above threshold.
4. **Cache hit shared with views:** calling `render_thumbnail` after `build_surface_init` from a view (same `record`) is fast; `build_surface_init.cache_info().hits ≥ 1`.
5. PNG is accepted by the storage backend / `ImageField` without error (regression for the existing thumbnail consumer).
6. **Ortho/persp `O` derivation:** for `eye=None`, the projection used has `O=[0,0,0]`; for `eye=[a,b,c]`, the projection has `O=eye` (verified via spy on `Projection.__init__`).
7. **No silhouette import:** `grep -rn "silhouette" surface_play/thumbnail.py` empty.
**Gotchas:** none new.
**Stop & verify:** `pytest surface_play/test_thumbnail.py` ⇒ 7 green; manually click "save thumbnail" in the browser on 3-4 surfaces, confirm the saved PNG resembles the live view (modulo expected fig-8 etc. fixes).

---

### W4 — Smoke + visual regression

**Spec:** §3 "Test strategy" → "Visual regression (final)" (this roadmap line 125).
**Files:** create [surface_play/test_smoke.py](surface_play/test_smoke.py), [surface_play/test_visual_regression.py](surface_play/test_visual_regression.py), and a `regression_baselines/` directory at repo root.
**API:** none — these are test modules.
**Algorithm:**
1. **Smoke** (`test_smoke.py`):
   - For every `SurfaceRecord` in the dev DB: call `pipeline.build_surface_init` then `pipeline.build_outline` at the canonical viewpoint; assert no exceptions and that the response shape is correct (`lines_by_visibility` and `si_lines_by_visibility` are dicts of `int → list`, `origin` is `(2,)`).
   - For every fixture in `test_fixtures.py`: same.
   - Hit the Django views via the test client too (smoke parity with W2's tests, but on the full DB rather than a single record).
2. **Visual regression** (`test_visual_regression.py`):
   - Two modes: `--mode=record` writes baselines; `--mode=check` (default) compares against existing baselines.
   - **Baselines:** run `silhouette.Surface.traitement(...)` on each fixture surface at canonical viewpoints; serialize the resulting `(lines_by_visibility, si_lines_by_visibility, origin)` to `regression_baselines/<fixture>_<viewpoint>.json`.
   - **Check:** run `pipeline.build_outline` on the same inputs; compare to baseline.
   - **Tolerances per channel:**
     - `origin`: `‖new − old‖ < 1e-6`.
     - For each `vis` key: number of curves matches; per-curve sample count matches within ±2; per-sample `xy` matches within `0.01 · M` (M = projected bbox diag).
     - **Known divergences allowed:** the fixtures listed in success criteria 5-6 (Onde radiale, vagues, fig-8 immersion) carry per-fixture override flags that loosen tolerances or skip channels with documented reasons (textual `"reason":` field in the override JSON).
   - **Failure output:** writes a side-by-side matplotlib PNG of `(silhouette, pipeline, diff)` to `test_artifacts/<fixture>_<viewpoint>.png` for human review.
**Test criterion:**
1. Smoke passes on every DB record and every test fixture (no exceptions, valid shapes).
2. Visual-regression baselines exist for every test fixture × canonical viewpoint pair.
3. `--mode=check` passes for all fixtures except those with documented known divergences; for those, the divergence is in the loosened channel only and is explained by the textual override reason.
4. **Determinism of new pipeline:** two runs in a row with the same seed produce byte-identical regression output.
5. CI integration: `pytest surface_play/test_smoke.py surface_play/test_visual_regression.py` runs and exits 0 on the reference dev DB.
**Gotchas:**
- Baselines are committed to the repo and treated as a contract. If a code change legitimately moves a baseline, the change must be reviewed and the baseline regenerated in the same PR.
- `silhouette.py` must still be importable in this step — it's the source of baselines. Retired only in W5.
**Stop & verify:** smoke + visual regression both green; `git status` shows the new baseline files and the (small) set of expected-divergence overrides.

---

### W5 — Retire `silhouette.py`

**Spec:** §0 "Non-goals" line 21 (no retirement until feature parity); §0 "Success criteria" item 1.
**Files:** move `surface_play/silhouette.py` → `old stuff/silhouette_legacy.py` (preserve history via `git mv`); search-and-replace any remaining imports across the project; update [Reorganization_of_the_surface_app.md](Reorganization_of_the_surface_app.md) status note.
**API:** none.
**Algorithm:**
1. `git mv surface_play/silhouette.py "old stuff/silhouette_legacy.py"` — moved out of the import path, not deleted (kept for forensic reference).
2. `grep -rn "from surface_play.silhouette" .` and `grep -rn "import silhouette" .` — both empty (one allowed exception: the visual-regression baseline-recording path, which now imports from the legacy location).
3. Update the visual-regression script to import the legacy `Surface` from its new location (use `importlib` if the hyphen-in-folder breaks plain `import`). Re-run visual regression in `check` mode against existing baselines: must pass.
4. Update [Reorganization_of_the_surface_app.md](Reorganization_of_the_surface_app.md) header to reflect "silhouette.py retired YYYY-MM-DD; legacy at old stuff/silhouette_legacy.py".
5. Run the full sweep: `pytest surface_play/`. All green.
6. Manual smoke: spin up the dev server, browse `/play/<pk>/` for 3-4 records (helicoid, fig-8 immersion, Möbius, paraboloid), rotate, switch to perspective, toggle each `debug` knob. No regressions vs. baseline; fig-8 visibly correct (success criterion 6).
**Test criterion:**
1. `grep` confirms no live import of `silhouette.py` from `surface_play/`.
2. Full pytest pass.
3. Visual regression still passes against the W4 baselines.
4. Manual smoke (4 records × 4 viewpoints each = 16 page loads) shows no visible regressions.
**Gotchas:**
- **Pre-condition: W4 visual regression must be green and reviewed.** Do not retire silhouette.py while any baseline divergence remains unexplained.
- If a record in the production DB depends on a silhouette-only feature not yet covered (e.g., a custom surface type), surface it as a blocker — do not silently fall back to silhouette.
**Stop & verify:** all four test criteria pass; manual smoke documented in the PR description with screenshots.

---

### Layer W sign-off gate (final — project complete)

1. All 5 step tests + smoke + visual regression green.
2. User has visually inspected at least 4 surfaces (helicoid, fig-8 immersion, Möbius band, one complex like Onde radiale or vagues) at multiple viewpoints in the browser, including fig-8 (which should now render correctly per success criterion 6).
3. User confirms no frontend regression on the existing surfaces (success criterion 4).
4. `silhouette.py` retired; only `surface_play/` modules are imported by views/thumbnail (success criterion 1).
5. All Layer P/C/O unit + integration tests still green on the final tree (success criterion 3).
