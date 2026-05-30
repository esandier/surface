# General goal of the app: plot a contour outline of a parametrized surface S.

> **Status (W5, 2026-05-30): `silhouette.py` retired.** The modular rewrite (Layers C/O/W) is complete and is now the sole production path (`surface_play/pipeline.py`, `views.py`, `thumbnail.py`). The legacy monolith was moved to `old stuff/silhouette_legacy.py` (history preserved via `git mv`); it is no longer importable as part of the package. It is loaded by file path only by the visual-regression *record* path (`test_visual_regression.py`) to regenerate baselines — check-mode pytest does not need it.

This file describes an app to compute the outline of the surface. Currently, it is partially and mainly implemented in the file silhouette.py, with the html templates located elsewhere. However several attempts to start from this file to add self-intersection logic have failed.

I would like Claude to refactor the code, only using the current silhouette.py as extra information. The present file remains the true specification of the app, the bible. Claude should not hesitate to diverge from silhouette.py implementation choices if it serves unification of functions, streamlining and performance. However I am quite happy with the frontend. It should not be  modified in a way that would significantly change the user experience. Also, silhouette.py actually works, and is quite good at what it does. Therefore, for logic which does not involve self-intersections, it can be used to clarify points which are not detailed in this document. 

I insist that computations should be unified in functions when they are called several times, and as vectorized as possible for performance. Among unified functions:
+ ``make_lines`` (see below).
+ Mesh generation.
+ 2D line sweep for intersections, which is called in the domain and in the projection. One subtle point here is that there are possible identification of opposite domain rectangle sides. In that case, segments in the domain with very different coordinates can be very close, and thus before testing intersection, one must choose coordinates for the segment ends which are close. This is explained more below.
+ Detection of sign changing on a collection of segments, which is used for detection of cusps and on contour points.

Several data structures are needed to insure the visualization of the surface. Some depend only of the surface and are computed when it is initialized, other ones depend on the viewpoint and type of projection, they are computed each time these change. The former are not re-computed at each change of viewpoint (i.e. rotation or zoom).

## The app from the user point of view

From the user point of view, the app has a home page which is a list of surfaces with thumbnails, from which the user may choose to add an item to the list, edit an existing item, or simply visualize an existing surface. The top banner has buttons:
+ A '+' to add an item.
+ A camera button which, when clicked, makes the current view the inital view of the surface, and updates the thumbnail accordingly
+ A 'controls' button which triggers the display of a panel to control rotation/zoom of the surface.
+ An 'eye' button which toggles between perspective mode and orthographic mode
+ A 'debug' button, which is displayed only on the dev server, and which triggers the display of a debug panel.

When an item is selected for display, a 3D view of the surface is displayed (with a resolution set by CANVAS_RESOLUTION in settings.py). Upon pressing and releasing mouse, two pipelines are triggered: the construction pipeline, then the outline pipeline (see below). These pass to the front-end the outline curves and their visibility data, which are displayed.

The user can then rotate/zoom the surface by mouse control (or the control panel, see above). Upon releasing the mouse, the outline pipeline is triggered, the outline curves are computed anew, and the display is updated. During rotation/zoom, the outline is not displayed, it reappears upon mouse release.

Note that zoom in orthographic projection does not trigger the outline pipeline, since the outline curves are unchanged.

## The Construction pipeline

The steps of the construction pipeline are:
+ Generate a mesh according to the surface's 'domain' data: the shape of the domain (rectangular or disk/annulus), the possible identifications of sides (cylinder type or Moebius type, see below) in the rectangle case. The mesh resolution is determined by the GRID_RESOLUTION variable in settings.py. The mesh is generated, then identification of the vertices and edges takes place.

+ Then the mesh is perturbed, respecting the boundary (boundary vertices are perturbed so as to stay on the boundary), to avoid non-generic situations.

+ Add to the surface parametrized equations a perturbation which should be small (not visible, of the order of .0001 times the bounding box of the surface), random (not really, coefficients are chosen once and for all, like .000342876545), and respect the identification (cyl, Moebius) if there is one. It should not be identically zero on the boundary, though, so as to effectively perturb the surface. 
+ Generate lambdified functions for the surface parametrization $S(u,v)$, the normal to the surface $N(u,v) = \partial_u S\times \partial_v S$, partial derivatives of $S$, second partial derivatives of $S$, partial derivatives of $N$. 

+ Compute the normal at every vertex of the mesh.

+ Mark the mesh edges which are boundaries. Add to each boundary point the direction in (uv) space pointing inside the surface. 

+ For each edge, add a 'flip' flag which is +1 if the inner product of the normals at the edge ends is positive, and '-1' otherwise. A -1 indicates a change of local orientation, which can happen if the case of 'moebius' type identification.

+ In the case of a rectangle domain, store in a collection the indices of corner vertices.

+ Make Boundary Curves (BCs) by chaining boundary edges (BEs), either closed (disk/annulus case or rectangle with at least one identification) or open (rectangle case, no identification) with corners serving as ends. Create a Split Point (SP) for each corner, and record a Split (SPT) on each BE containing a corner, with a `bary`field equal to 0 or 1 according to wether the corner is the origin or the end of the BE. 

+ Compute Double Points (DPs).

+ Make Self Intersecting segments (SISs) joining DPs which belong to the same face on each sheet of the surface.
+ Make Self-Intersecting Curves (SICs) by chaining SISs. They are either closed, or open with ends at DPs belonging to the boundary (BDPs).
+ Record a Split Point (SP) for each BDP, and add a Split (SPT) on the segment containing the BDP (see below for the definition of a Split). This SPT has a `bary` field equal to either 0 or 1, according to wether the BDP is the origin or the end of the SIS.
+ Compute intersections in (uv) space of SISs. Intersection of two SISs happen at Triple points (TPs). 
+ Record splits of SISs at TPs.

### implementation note.
Curves of types CC, BC or SIC are stored in different arrays, or the same array but then each has a 'type' field, in addition to an 'open' boolean field. 

Each curve, in addition to these fields, contains an array of indices of chained segments, with a negative index corresponding to a reversed segment.

The type of the curve indicates to which array of segments  the indices refer to. For BCs, the segments are mesh edges, for CCs, it is Contour Segments (CSs), for SICs, it is Self intersecting Segments (SISs). Segments all have a `split` field which is equal to -1 for no split, and otherwise is an index to an array SPLITS. An element of that array is a triple `(sp_index, bary, vis_chge)`, where `sp_index` is an index into an array of Split Points (SPs), `bary`is the barycentric coordinate of the SP in the segment, and `vis_chge` is an integer equal to minus the change in the number of surface sheets occluding the point when it crosses the SP in the forward direction of the segment. 

We will describe below the computation of `vis_chge` for all the different types of splits. 


## The Outline pipeline
+ Compute the Contour Points (CPs). This involves testing every edge of the mesh. 
+ Make contour segments (CSs), each joining 2 CPs belonging to the same face.
+ Make Contours Curves (CCs) by chaining the CSs. They are either closed, or ending at CPs which are on the boundary (BCPs). 
+ For every BCP, create an SP, and record two SPTs, one for the BE containing the BCP, and one for the CS which originates, or ends, at the BCP.
+ Detect Cusp Points (VPs) on CSs. For each VP, create a corresponding SP, and record an SPT on the CS containing the VP.
+ Compute intersections in (uv) space of SISs with CSs. The intersections are Contour Double Points (CDPs). For each CDP, create a corresponding SP, and record two SPTs: one on the SIS, and one on the CS.
+ Add Helper curves (HCs), if needed, so that all curves belong to the same component (see below for the definition of a component.) In practice, an HC is a line segment in `(uv)` coordinates connecting two previously defined curves. Their endpoints are called Helper Attachments (HAs).
+ For each HA, create a corresponding SP, and record an SPT on the curve the HA is attached to. In particular, compute Visibility Changes on HCs at HAs. Note that the curve to which the HC is attached has no visibility change at an HA.
+ Do the splitting: trasverse every BC, CC, SIC, and build out of them the Splitted Curves, denoted SBCs, SCCs, and SSICs. The HCs are splitted at construction, they do not need resplitting. A splitted curve is `{type, {sp_indexp, sp_indexq}, {vcin,vcout}, [p1,...,pn-1]}`.
+ Resample curves:
-  At every SP, measure the lengths of the curves of which it is an endpoint, let `l` be the smallest of these lengths. Let `d=min(l/10, M/GRID_RESOLUTION)`, where `M` is an approximation of the diameter of the surface in projected space. 
- Then, for each half-curve originating or ending at the SP, resample it so that the segments have length `d`.
- In the case of SCCs, reproject the interpolated points on the CC by running a newton method on the line orthogonal to the SCC through the interpolated point.
- In the case of SBCs, if the boundary is a disk or an annulus, reproject the interpolated points on the boundary by adjusting their norm. 
+ Compute intersections of curves (SBCs, SCCs, SSICs, HCs) in projected space. Each intersection creates a Break Point (BKs) on the curve which is behind from the viewer point-of-view.
+ Compute the visibility change (VC) associated to each BK.
+ Determine anchor points (APs). These are the points with min or max, x or y coordinate in projected space. There are 4 of them.
+ Determine the visibility of portions of curves between BKs or endpoints: this is **minus** the number of sheets of the surface hiding the curve portion from view. It is `0` at anchor points, and then is propagated to the whole curve network either by BFS propagation, or by solving a linear program with the constraint of `0` visibility at either one or four APs. The method used is given by the variable PROPAGATION in settings.py, which can take the value 'BFS', 'LP1' or LP4'
+ The projected curves and their visibility info is forwarded to the front-end for display. 


# Domain of the parametrization

## User specification

The user specifies the type of domain used which can be
- A rectangular domain, specified by the x and y coordinates range. For each pair of opposite edges of the rectangle, the user may specify wether if they are: 
    - not identified,
    - identified keeping orientation (cylinder-like),
    - identified reversing orientation (Moebius band-like).
- An annular domain, specified by the inner and outer radius (inner radius is zero for a disk)

The data specified by the user can be a sympy expression which evaluates to a real number, like $2\pi$ for instance.

## Data stored

The mesh keeps **pre-identification** vertex indices throughout. Vertices have `(u,v)` coordinates in the original rectangular (or polar) domain; identifications do not rewrite indices in `tris`. A separate `vertex_class[i]` array gives the canonical equivalence class of each pre-id vertex (used for 3D evaluation and for traversal across the seam); identification of higher-dimensional cells is then derived geometrically (see below), **not** by collapsing labels in `tris`.

The edge `(i,j)` contains a `p` and `pq` field equal to `p` and `q-p`, where `p` and `q` are the pre-id `(u,v)` coordinates of vertices `i` and `j`. Because both endpoints of any single edge live in the same fundamental domain (a triangle is never split across the seam), no `close()` arithmetic is needed here. Edges also contain the indices of the two adjacent faces; boundary edges have only 1 adjacent face, and the second index is `-1`.

Boundary edges contain a vector `dir`, orthogonal to `pq`, oriented so that the inner product of `r-p` and `dir` is positive, where `r` is the third vertex of the adjacent face. An array contains the indices of boundary edges, together with their `dir` field.

Faces contain the indices of their vertices and of their edges, as well as fields `p`, `pq`, and `pr` equal to `p`, `q-p` and `r-p`, where `p`, `q`, `r` are the pre-id `(u,v)` coordinates of the face's vertices. All three live in the same fundamental domain by construction (a face is one pre-id triangle), so these are direct subtractions.

### Identification of edges and faces — geometric, not by vertex labels

Identifications are interpreted as a gluing of *sides of the rectangle*, not as a relabeling of vertex indices. The rules:

- **Vertices**: paired by the side-identification rule (cy: same coordinate on the other side; mo: mirrored coordinate). The pairing produces equivalence classes, recorded in `vertex_class`. Vertices are otherwise kept distinct with their original pre-id indices.
- **Edges**: two pre-id edges `{a, b}` and `{c, d}` merge into one post-id edge **iff** both edges lie on identified sides of the rectangle (both endpoints of `{a, b}` on one identified side, both endpoints of `{c, d}` on the partner side) **and** their endpoints correspond under that side's pairing (`a ~ c, b ~ d` or `a ~ d, b ~ c`). All other edges keep their own identity, even when their endpoint pairs happen to share canonical labels with another edge.
- **Faces (triangles)**: never merged. Each pre-id triangle is its own 2-cell in the post-id mesh.

This rule is required because under reverse identifications on both axes (mo-mo), opposite-corner pre-id triangles can have all three of their vertices pairwise identified, producing identical canonical vertex labels for geometrically distinct cells. A label-based dedup would silently delete one of those cells and tear the manifold. The post-id complex is best viewed as a Δ-complex (semi-simplicial): multiple cells may share the same vertex set without being the same cell.

## Domain API

The domain object provides a `close(p,q)` function which, in the case where there is identification returns (uv) coordinates for `q` which are close to those of `p`. 

It also provides a function `bary`, of signature `bary(e,s)`, where `e` is an edge,  which returns the `uv` point `e.p + s e.pq`, or `bary(f,s,t)`, where `f` is a face, which returns `f.p + s f.pq + t f.pr`.

It also provides `boundary_tangent(uv, edge_dp)` which returns the **analytic** unit tangent of the smooth boundary curve at `uv`, sign-matched to `edge_dp`'s direction. For a rect domain each side is axis-aligned, so the chord is already analytic — `boundary_tangent` returns the normalized `edge_dp`. For disk/annulus the boundary is a circle, and the analytic tangent at `(u, v)` is `(-v, u)/sqrt(u²+v²)`. This is used in `bvis_chge` and elsewhere where a smooth boundary tangent is needed; the discretized chord can deviate from the analytic tangent by O(1/res) which is enough to flip sign-discriminators in near-edge-on views.


# Parametrization of the surface

## Parametric expression.
The surface S is specified parametrically, by sympy expressions for the coordinates in space of a point. The coordinates can be cylindrical or cartesian, and are given by a sympy expression depending on the coordinates of a point in the domain. These coordinates are cartesian for a rectangular domain, but may be either cartesian or polar for an annular domain. In the case of polar coordinates, these are then translated to cartesian internally for the rest of the pipeline. 

## Perturbation.
In order to avoid degenerate situations, it is useful to add a perturbation to the specified parametrization of the surface which is not visible to the user but resolves non-generic situations.The expression we use for the perturbation is (see silhouette.py file)

``` 
        freq_u = 2 if self.u_identify != "mo" else 1
        freq_v = 2 if self.v_identify != "mo" else 1
        du, dv = u_max - u-min, v_max - v_min


        S_sp = S_pur + (sp.sin((u - self.u_min) * freq_u * sp.pi / d_u) * sp.sin((v - self.v_min) * freq_v * sp.pi / d_v)) * sp.Array(
            [
                0.0005 * dX * .346,
                0.0005 * dY * .632,
                0.0005 * dZ * .693,
            ]
        )
```

## Derivatives, normals.

From the parametric expressions, a number of lambdified functions are generated for later use. These include:
+ First and second derivatives with respect to the parameters.
+ Normal at a point, defined as the cross-product of the two first derivatives. It is denoted SN(p), where p is a point in the domain.

**Implementation note** Use cse() to regroup S, Su, Sv, Suu, etc. in a single evaluation, for speed. Also, use pure vectorization for evaluation of these functions, something like 
```
def _eval_vec(fn, u, v):
   	return np.array(fn(u, v))
```    



# Projection of the surface

The 3D surface S is projected in the vision plane either by orthographic projection, or a perspective projection.

When a viewpoint, zoom factor is selected, the component of the normal to the surface in the direction of the viewer is computed in a vectorized way at each vertex of the triangulation.


# Mesh
A mesh is computed from the domain data. It is important that faces are triangles.

In the rectangle case, the mesh is generated on the *raw* rectangle, then identifications are applied geometrically:

- **Vertex equivalence** is recorded in a `vertex_class` array: two boundary vertices on identified sides are placed in the same class (canonical = the smaller pre-id index). Vertex indices in `tris` are **not** rewritten.
- **Edges** are merged across identified sides: a pre-id edge on the u-min side merges with its partner on the u-max side iff both endpoints pair under the side's identification rule, and similarly for v-sides. Edges that do not lie entirely on identified sides retain their own identity even if their endpoints' canonical labels coincide with another edge's. (See §"Data stored".)
- **Faces (triangles)** are never merged. Each pre-id triangle is its own 2-cell in the post-id mesh.

Each mesh vertex contains as data its coordinates in the domain, the normal to the surface (SN) computed at this point, and the 3D coordinates of the point. Paired vertices each carry their own evaluations; consistency across the seam is handled by the `flip` flag on edges (see below) and by `vertex_class` lookups when needed.

In the rectangular case, the corners are special points; an array contains the indices of the corner vertex classes (after merging via `vertex_class`).

In the rectangular case with identifications which reverse orientation, it may happen that the computed normals to the surface at the ends of an edge have negative scalar product instead of being positive. Thus a `flip` flag is associated to edges, which is `true` if this inner product is negative. For a side-identified edge with two adjacent faces (one in each fundamental-domain copy), `flip` compares the SN at corresponding paired endpoints across the seam.

# Curves

Various types of curves are associated to the surface. Each curve is an array of segment indices. The index has a `-` sign if the segment is reversed in the curve. The curve may be closed, if its endpoints are equal, or open.

Each curve has a field `closed` which is either `true` or `false`, and a `type` field, which will come up below. The `type` field indicates to what array of segments the indices refer to. It can be `bc`, `cc` or `sic`, for boundary curves (BCs), Contour Curves (CCs) or self-intersecting curves (SICs). 

Their construction is unified: a `make_lines(segments)` function exists, where `segments` is a collection of segments (i.e. pairs of point indices). The algorithm for this function is to choose a segment, and extend it in both directions by segments which share an endpoint with it, until a point with degree 1 is reached, or the curve has closed. Once the segment collection is consumed, a collection of curves is returned.

Note: It is left to Claude to decide the precise data structure for the collection, but speed is important.

## Boundary curves (BC)

Each boundary of the surface gives rise to a number of boundary curves, which are  closed in the annular domain case or in the rectangular case. In the rectangle case with both pairs of opposite edges identified, there is no boundary curve. 

The initial boundary curves are generated by calling make_lines on the collection of boundary edges (those which belong to a single face).

Boundary curves are computed in the construction pipeline, they are not updated when the surface is rotated/zoomed.

## Self-intersecting curves (SIC)

SICs are the curves where the 3D surface intersects itself. They are polylines of double points (DP, to be defined below). They may be closed, or open with endpoints being DPs which belong to the boundary of the surface.

They are generated by calling make_lines on the collection of self-intersecting segments (SIS, to be defined below). 

SICs are computed at initialization time, they need not be updated after rotation/zoom.

To be clear: an SIC is stored as one curve only. It does have two preimages in the domain, but these are deduced by looking at the two domain points associated to a double point (DP) of the SIC, and then tracking which DP belongs to which preimage by using the 'flip' flag on SISs.

**Double points** are either the intersection in 3D space of a triangulation edge with a face it does not belong to, or the intersection in 3D space of two triangulation edges. Each DP has associated to it :
- 3D space coordinates
- Two domain coordinates.
- The intersection pair (E1,F2) or (E1,E2)
- Two sets of faces : A1 = faces to which E1 belongs, there are either 1 or 2. A2 = {F2} for an (E1,F2) intersection, or A2 = the set of faces E2 belongs to in the case of an (E1, E2) intersection.

The algorithm to determine the collection of double points consists in consuming (edge,face) pairs of the triangulation, each pair being treated as described in the appendix. This is **the** most time-consuming routine of the app. It must be super-optimized, it is left to Claude's appreciation how to achieve this, but the goal is at most a few seconds for a 300x300 grid. A crucial step is to select candidates from bounding-box data. This will require a tree approach, and possibly numba to optimize for speed.

**Implementation Note** When computing the domain-coordinates of double points, one will interpolate using the barycentric coordinates of the intersection in the image. This requires to use domain coordinates of the implied vertices in the domain, **which are close**, and thus use the ``close(p,q)`` function to avoid problems coming from identification. 

**Self-intersecting segments** Two DPs $P$ and $P'$ form an SIS if and only if, with the obvious notation,  $A_1\cap A_1'\neq \varnothing$ and $A_2\cap A_2'\neq \varnothing$, or $A_1\cap A_2'\neq \varnothing$ and $A_2\cap A_1'\neq \varnothing$. A `flip` flag is associated to each segment. In the first case it is `false` and in the second case it is `true`. This flag is important when chaining segment in the `(uv)` space to form SIC pre-image curves.

**Implementation note (MT vectorization):** The Möller-Trumbore test over all candidate (edge, face) pairs is batched with NumPy before per-hit disambiguation, reducing Python-loop overhead from O(n_candidates) to O(n_hits). Timing: ~0.1s for ~1400 candidates, 56 DPs. 


## Contour curves (CC)

Contour curves are viewpoint dependent, and zoom factor-dependent in the case of a perspective view. They are the curves on the surface where the normal is orthogonal to the viewing direction, which varies in the case of perspective view and is constant in the orthographic case. It is denoted as the ``axis`` (a 3D vector pointing in the viewer's direction).

They are constructed by calling make_lines on the collection of contour segments (CS, to be defined below).

**Contour points (CP)** are computed as follows:
- Detect edges where ``inner(SN, axis)*flip`` change sign on the edge. Here `inner` denotes the inner product of two vectors. 
- For each edge of the resulting collection, run a Newton-algorithm to determine a point where ``inner(SN, axis) = 0``.

This results in a collection of CPs, each having 3D coordinates, and domain coordinates, as well as an edge index, and a barycentric coordinate on that edge. The Newton search should be vectorized for speed, as well as the selection of sign-changing edges.

Note that the newton method on edge`e` is applied to the function `s -> inner(SN(e.p+se.pq), axis)`, with $s\in(0,1)$, where `e.p`and `e.pq` have been previously defined. 

**Contour segments (CS)** Two contour points form a contour segment if and only if the edges they belong to  are adjacent to a common face. 

# Splitting of curves.

Curves are split at specific points. The first case is when two curves intersect in the domain. The second case is that of CCs which are split at cusp points. 

## Split Points
The Split Points (SPs) are  a pair `(u,v)` of domain coordinates, and a type which describes if they are an intersection, and of which type, or a cusp. The possible types are `cdp` (contour double-point), `bcp` (boundary contour point), `bdp` (boudary double-point), `vp` (cusps), `cn` (corners), `ha` (helper curve attachment).  All split points are held in a common array.

## Splits
A Split (SPT), record the effect of the split point on one curve it belongs to. If a SP is the intersection of two curves, it will give rise to two SPTs. An SPT has the structure `[sp_index, bary, vis_chge]`, that is the index of the corresponding split point, the barycentric coordinate of the point on the segment, the visibility change when crossing the split point in the forward direction on the segment. All splits are held in a common array.

## Splitting of segments
Splitting a segment (either a BE, a CS or an SIS) consists in adding to the `split` field the index of the split. 

**Important** Every BE, CS and SIS has a `split1` and `split2` field (not mentionned before), which contains indices to SPTs, or `-1` if the segment is not split.

**Visibility** The visibility of a point of the surface is minus the number of sheets of the surface which hide the point from the viewer. It is a nonpositive number, equal to zero if and only if the point is visible to the user. The  *visibility change*, denoted `vis_chge` is minus the number of additional surface sheets occluding the point from the viewer after crossing the split.

## Splitting of BCs at Contour Points (CPs).
When a contour point (CP) belongs to a boundary edge, a splitting occurs.

+ An SP is created of the type `bcp`.
+ A split (SPT) is created on the boundary edge (see below for the visibility change).
+ A split is created on the CS the point belongs to (there is only one such segment, since we are on the boundary.) The `bary` of the SPT is either 0 or 1, depending on wether the point is the origin or the end of the CS.

**Visibility change** The `vis_chge` field of the second split is `0`. The `vis_chge`of the first split (the one on the BE) is computed as follows:

```
    def bvis_chge(self, e, s):  # e est une arête de bord, s est la coord barycentrique du point sur e.
        p, dp, Nb = e["p"], e["pq"], e["dir"]  # segment dans le domaine

        p0 = p + s * dp  # le point où la visibilité peut changer

        # Tb = ANALYTIC boundary tangent at p0, sign-matched to dp.
        # NOTE: deviates from earlier `dp / norm(dp)` (the chord). The chord
        # can deviate from the smooth tangent by O(1/res); for near-edge-on
        # views, this is enough to flip the sign of the f3 discriminator.
        Tb = domain.boundary_tangent(p0, dp)

        ker = self.kerdS(*p0)[0]
        # oriented toward the upper sheet (axis · dS(ker) > 0)
        if self.Z(self.dS(*p0, *ker)) < 0:
            ker = -ker

        # Np = ∇_uv(axis · SN), the SILHOUETTE-FUNCTION gradient.
        # NOTE: deviates from earlier `SD_jac(*p0, *self.axis) = ∇(axis·S)`
        # (the depth gradient). The depth gradient is WRONG for the silhouette
        # crossing test: it gave the same sign at both BCPs of the same CC
        # (cycle did not close on paraboloid trial 2). The silhouette is
        # `axis·SN = 0`, so its gradient is the right quantity.
        Np = np.array([
            (np.cross(Suu, Sv) + np.cross(Su, Suv)) @ axis,
            (np.cross(Suv, Sv) + np.cross(Su, Svv)) @ axis,
        ])
        if np.inner(Np, ker) < 0:
            Np = -Np  # oriented toward upper sheet

        # Tp = silhouette tangent in uv (perpendicular to Np), oriented inward.
        Tp = np.array([-Np[1], Np[0]])
        if np.inner(Tp, Nb) < 0:
            Tp = -Tp

        # Image projections of dS(Tb), dS(Tp) onto the screen plane.
        Tb_proj = proj_vec(dS(Tb))
        Tp_proj = proj_vec(dS(Tp))

        f1 = np.inner(Tb, Np)         # silhouette-function change along Tb
        f2 = np.inner(Np, ker)        # = +|Np||ker| after orientation
        f3 = np.inner(Tb_proj, Tp_proj)

        # NOTE: this 3-way `f1·f2·f3` gate REPLACES the earlier
        # `np.inner(ker, Nb) > 0` gate, which was numerically brittle in
        # near-edge-on views (ker · Nb fired on FP noise).
        if f1 * f2 * f3 > 0:
            return 0    # no observable visibility change at this BCP
        return 1 if f1 * f2 > 0 else -1
```

**Known degeneracy** for the formula above: when the BCP coincides with a depth-aligned boundary point (i.e., `dS(Tb)` is nearly parallel to `axis`, so `|Tb_proj| ≈ 0`), the `f3` sign is FP-noise. Failure in this configuration is accepted as non-generic.
**Implementation note** The contour points belonging to boundary edge are easily found: they are the endpoints of CCs which are not closed. 

## Splitting of BCs at Double Points (DPs). 
When a double point belongs to a boundary edge, a splitting occurs.

+ An SP is created of the type `bdp`.
+ A split (SPT) is created on the boundary edge (see below for the visibility change).
+ A split is created on the SIS the point belongs to (there is only one such segment, since we are on the boundary.) The `bary` of the SPT is either 0 or 1, depending on wether the point is the origin or the end of the SIS. For this split, `vis_chge = 0`.


**Visibility change** The visibility change on the BE is computed as follows.
+ The DP belongs to two sheets of the surface. On one of the sheets, it belongs to the boundary edge `e`. On the other sheet, it belongs either 1) to a boundary edge `e'`, or  2) to a face or edge which is not on the boundary. In both cases, let `Tb` be the tangent vector to `e` in 3D space and let `SN'` be the normal to the surface's other sheet at the point in 3D space. Also, let `Tb_proj` be the projection of `Tb` in screen space. 

+ In case 1), let in addition `dir'` be the normal to `Tb_proj` in projected space oriented toward the surface. If `inner(SN',axis)*inner(Tb,SN')*inner(Tb_proj,dir)>0`, then the visibility change for `e` at the point is $0$. If `inner(SN',axis)*inner(Tb,SN')*inner(Tb_proj,dir)<0`, then the visibility change is +1 if  `inner(SN',axis)*inner(Tb,SN') > 0`, otherwise it is -1.
+ In case 2),  the visibility change depends is +1 if `inner(SN',axis)*inner(Tb,SN') > 0`, otherwise it is -1.
 

## Splitting of SICs at Triple Points (TPs). 

### Triple points.

A triple point is the intersection in 3D space of three SISs.

**Labelling convention:**
Three SISs  `S1, S2, S3` meet at a triple point. Each SIS has two preimage segments;
the preimage segment NOT belonging to `Si` is called `pSi`. Labelling:
- `pS1` has face-pair `(F3, F2)` — intersection of preimages of `S2` and `S3`
- `pS2` has face-pair `(F3, F1)` — intersection of preimages of `S1` and `S3`
- `pS3` has face-pair `(F2, F1)` — intersection of preimages of `S1` and `S2`

The three domain intersection points `P1 ∈ F1`, `P2 ∈ F2`, `P3 ∈ F3` all map to the
same 3D point (the triple point).

**Algorithm:**
1. Run a 2D domain sweep between all SIS preimage segments to find all pairwise
   domain intersections. Note that when testing for intersection two segments in the domain, one must take into account potential identification. Therefore, domain coordinates of the segment edges are chosen close, using ``close(p,q)``.
2. Group intersections into triples by matching their face-pair data: three intersections
   form a TP when their face-pairs interlock as `(F3,F2)`, `(F3,F1)`, `(F2,F1)`. 

### Splits
For each triple point we create an SP, and an SPT for each of the three SISs involved.

###  Visibility changes.
For `Si` with face-pair `(Fj, Fk)` at the TP (i.e. the face-pair of `Si`'s own preimage
segments at the TP):
- Let `Fl` be the third face (not in `Si`'s face-pair), `Pl ∈ Fl` the corresponding
  domain preimage point.
- `N = SN(Pl)`, oriented so `Z(N) > 0`.
- `T_Si` = 3D tangent to `Si` at the TP (from `dS` along `Si`'s preimage segment).
- `vis_change on Si = +1 if dot(T_Si, N) > 0 else -1`

Same formula applies to `S2` and `S3` with their respective third faces.


## Splitting of CCs at cusps.
Cusps are (apart from contour points), the only generic singularities of maps from the plane to the plane. They belong to CCs. We call them VPs because the shape V reminds of a cusp.

### Cusp points (VPs).
Just as CPs are detected first by a sign change on triangulation edges, and then with a refinement by Newton ; VPs are first detected by a sign change on CSs, and then with a Newton method refinement, which is more complex then for CPs.

**Sign change** The sign change logic is the following.
```
directions = dir_vec[:, l]          # shape (2, len(l)) — 2D "into-surface" direction at each contour pt

# Dot product between each consecutive pair of direction vectors
signes = np.sum(directions[:, 0:-1] * directions[:, 1:], axis=0)  # shape (len(l)-1,)

# A negative dot product means the direction flipped → cusp between pts i and i+1
indices = np.where(signes < 0)[0]
for i in indices:
    p = np.array([u_vec[l[i]], v_vec[l[i]]])               # param coords of left pt
    q = self.close(p, [u_vec[l[i + 1]], v_vec[l[i + 1]]]) # right pt, adjusted for domain wrap
    r = (p + q) / 2                                         # cusp location = midpoint

    # Visibility: is the surface moving toward or away from the viewer along this edge?
    vis = 1 if np.inner(self.dS(*p, *(q - p)), self.axis) > 0 else -1

    cusps.append((l[i], l[i + 1], 1 / 2, r))
    self.breaks["c"][l[i]] = (0.5, vis, len(cusps) - 1)  # register break at this edge
```


**Newton refinement** Then we need to depart from the current version. To find the VP more precisely, we proceed by bisection:
+ On an contour segment (CS) with a sign change, take the midpoint. 
+ Run Newton on the line orthogonal to the CS through the midpoint to find a CP on this line. 
+ This CP replaces the midpoint which was not a true CP. Then the original CS is replaced by two segments having the new CP as one endpoint.
+ Detect which of the 2 segments has a sign change.
+ Iterate until VP is found with high precision. 

This refinement may or may not prove useful, we'll see. It is used, or not, according to the value of the variables "NEWTON_CUSP" of settings.py, which can be true or false.

Once the cusp is determined, it is recorded as a SP, and then a SPT is assigned to the original CS, before bisection. 

### Visibility change at VP.

The visibility change is computed as in the current silhouette.py file.


## Splitting of SISs and CSs at their domain-intersections.

### Computation of domain intersections.
SIS preimages and CSs can intersect in the domain. The intersection points are determined by a 2D sweep (taking into account the possible identifications, thus using the ``close(p,q)`` function). In that case, the SIS itself, as well as the CS, are split at the intersection point. At this point, both the SIS and CS get visibility changes.

### Visibility change of CS at SIS intersection.
The visibility change of the CS at an SIS preimage intersection is computed exactly as in the case of BDPs. Let `T` be the tangent vector to the CS in 3D space and let `SN'` be the normal to the surface's other sheet at the point in 3D space. 

The visibility change is +1 if `inner(SN',axis)*inner(T,SN') > 0`, otherwise it is -1. 

### Visibility change of SIS at CS intersection.

The visibility change of the SIS at the point with preimage P is determined as follows: Let `T` be the tangent to the SIS in `(uv)` space, let `N` be the normal to the CS in `(uv)` space, oriented toward the front sheet. Then the visibility change on the SIS is +1 if `inner(T,N)>0`, and -1 otherwise. 

# Construction of helper curves (HC) bridging components.
Curves are grouped into components, and these components are then connected by *helper curves* (HCs), which are straight lines in the domain.

## Definition of components.

Two curves having a Split (SPT) pointing to the same split Point (SP) belong to the same component. Note that in our terminology, two preimages of an SIC are not different curves, the curve is the SIC, hence the preimages automatically belong to the same component. 

## Algorithm for bridging.

For components $i$, $j$, the distance $d(i,j)$ in the domain between them is computed. Note that this minimum is achieved at points `qi`and `qj` which are endpoints of segments of curves in the components. It is very important, in the case where there is identification, that the domain coordinates used are the coordinates in the original domain rectangle !

**step 1** The components $(i,j)$ for which $d(i,j)$ is minimal are joined by the line segment `[qi,qj]` in the domain. The meaning of this is that SPs `qi` and `qj` are created, a SP at `qi` is created on either segment of the curve it belongs to, with `bary` equal to 0 or 1 according to what end of the segment `qi` is.
The same is done for `qj`. The visibility changes for these SPTs is 0. The SPs created are called Helper Attachments (HAs), and this the type of the SPs `qi` and `qj` which are created

Then a *helper curve* (HC) is created. This is a curve of a new type, a *Split Curve* (SC). An SC is a structure `[open, type, (vcin, vcout), [p0,..., pN]]` where `open` is boolean, `type` is `cc`, `bc`, `sic` or `hc`, where `vcin`and `vcout`are the visibility changes when entering the curve  from `p0` and from  `pN`, respectiveley. Here `p0` and `pN`are SP indices, and  the other `pi`s are point indices, pointing to an array determined by `type`.

For the helper curve we create, we have `open = true`, `type = hc`, `p0` is the index of `qi`, `pN` is the index of `qj`. The number of points is `N=2`. The visibility changes will be described below. 

Then components $(i,j)$ are grouped in a single component $n$, the distance of $n$ to component $k$ is the minimum of $d(i,k)$ and $d(j,k)$.

**step 2** If there is only one component after merging, *stop*, otherwise go to *step 1* to iterate with the new set of componants and distances.

## Computation of `(vcin, vcout)` for the HC `[qi,qj]`

The computation of both is similar. Let `T_proj` be either `qj-qi` (for `vcin`) or  `qi-qj` (for `vcout`). The point `q` is `qi` for the computation of `vcin` and is `qj` for the computation of `vcout`. 

If the curve `q` is attached to is a BC, then the visibility change is 0. 

If the curve `q`is attached to is a CC, we let `N` be the normal to the CS in `(uv)` space, oriented toward the front sheet. If `inner(T,N)>0`, the visibility change is 0, otherwise it is -1.

If the curve `q` is attached to is a SIC, let `SN'` be the normal to the surface's other sheet at the point in 3D space. Let `T` be the image of `T_proj` by the differential of `S` at `q`. The visibility change is 0 if `inner(SN',axis)*inner(T,SN') > 0`, otherwise it is -1. 

# Construction of Split Curves (SCs)

Apart from HCs, which are already of SC type, we need to cut the other curves into pieces, each of which is of SC type. 

The algorithm is as follows, and is very simple: every curve is cut at SPTs, and everything in between two consecutive SPTs make a SC. Several precisions must be made, though.
+ In case the curve is open, it is necessarily the case the there are SPTs at each end of the curve: a SPT in the first segment with `bary=0` (or `bary = 1` if the segment is reversed in the curve), and one in the last segment with `bary=1` (or `bary = 0` if the segment is reversed in the curve).
+ If the curve is closed, it contains necessarily a SPT, then the algorithm starts at this SPT, creates SCs between consecutive SPT until it returns to the startpoint. If only one SC has been created, the result is a closed curve.
+ It is quite clear who `p0,...,pN` are.
+ `vcin` is computed by taking the `vis_chge` field of the starting SPT, multiplying it by -1 if the initial segment of the curve is reversed, and taking the `min` of the result and 0. Thus `vcin` can be either `0` or `-1`. For `vcout`, take the `vis_chge` field of the ending SPT, multiply it by -1 if the last segment of the curve is **not** reversed, and take the `min` of the result and 0. Here too, `vcout` can be either `0` or `-1`. 

# Curve resampling.
## Resampled points.
Curves need to be resampled because we will compute their intersections in projected space, and these intersections will have discretization artefacts if the segments have non uniform lengths. Chord-alignment between curves sharing an SP suppresses near-tangent graze phantoms (same-sign break clusters and `+1/-1` pairs that BFS accumulates as visibility errors).

The resampling rule is **per kind**, with arclength inheritance at certain SP types to align chord boundaries between curves that meet there:

**CC (contour curve).** Samples = the SubCurve's own raw CPs verbatim. The CPs come from `find_contour_points` (Newton-refined), so they lie on the true contour. The CC defines the authoritative arclength ladder at every SP it ends at.

**BC (boundary curve).** Samples = the SubCurve's own raw CPs (mesh-boundary vertices crossed by the chain), **plus** inherited CC arclengths near each endpoint that is a **BCP**. Specifically: at each endpoint SP of kind `bcp`, find the unique CC ending at that BCP whose outgoing image-space tangent has positive dot product with the BC's outgoing tangent (the "tangent-aligned" BC arc at the BCP). On that winning BC arc, add samples at xy-arclengths from the BCP equal to the CC's CP xy-arclengths from the BCP, capped at `L_BC / 2` to avoid overlap with the inheritance from the other endpoint. The losing BC arc (whose tangent is antiparallel to the CC's) has no positive-aligned neighbour at this BCP and keeps only its own CPs.

The rule does **not** apply at BCs' corner endpoints, BDPs, CDPs, or any non-BCP SP.

**HC (helper curve).** At each HC endpoint (an HA), tangent-pick among all non-HC SubCurves incident at the HA — typically the two halves of the CC the HA is anchored to. The winner is the half whose outgoing tangent has the largest positive dot product with the HC's outgoing tangent. Resample the HC's emanating half at xy-arclengths from the HA equal to the winning CC half's CP xy-arclengths from the HA, capped at `L_HC / 2`. If no neighbour aligns (degenerate), fall back on the standard tapered ladder `d_k = min(k·δ, ell)` with `δ = (min incident SC length)/10` and `ell = M / RESOLUTION` on that half.

**SP-less SubCurves.** A SubCurve with no Split Points at all (closed parent whose component was alone — e.g., a flat disk with only its boundary curve, no silhouette and no helper curve) is passed through to the resampled output without subdivision. Its polyline is taken verbatim from the chain.


## Projection of the new points.
If the resampled curve is a BC, and we are in the case of an annular domain, the resampled points must be projected on the boundary. This is done by modifying the norm of the resampled point so that it lies on the inner or outer boundary.

If the resampled curve is a CC, then also the resampled points can be reprojected on the CC. This is optional, determined by the PROJECT_RESAMPLED variable in settings.py. If true, the reprojection is done by applying a Newton method on the line normal to the CC at the resampled point in `(uv)` coordinates to find a point nearby where `inner(SN(u,v),axis) = 0`.

## Per-sample data stored on a ResampledCurve.
Each `ResampledCurve` stores, per sample: `xy` (2D image-space position), `depth` (`axis · S`), `dir` (inward image normal, BC/CC only — lifted from uv via `dS` then projected), and `tan` (the **analytic image-space tangent**, BC/CC only).

`tan` is computed as follows:
- BC sample: `Tb_uv = domain.boundary_tangent(uv, edge_dp)`, then `tan = proj_vec(dS(Tb_uv))`.
- CC sample: `Np_uv = ∇_uv(axis · SN)`, `Tp_uv = (-Np_uv[1], Np_uv[0])`, then `tan = proj_vec(dS(Tp_uv))`.

After endpoint pinning, `tan` is sign-aligned with the local chord in image (`xy[j+1] - xy[j-1]` at interior samples, or the boundary chord at endpoints) so that `tan[j]` points along chain-forward direction. Without this alignment, an rc that traverses a CS in reverse of its native `p_cp → q_cp` direction gets `tan` pointing backwards relative to chain direction, which flips downstream sign-discriminators (e.g. projection-break signs).



# Projection intersections.
Until now, we computed visibility change at the endpoints of curves. Now we'll compute visibility changes at intersection points of curves *in the view plane*, i.e. intersections of the projected curves.

## Computation of intersections.
As before, intersections are computed using a 2D sweep of the collection of segments. Only intersections involving a BC of CC are of interest, the other intersections do not cause visibility changes.

## Visibility changes at intersections. 
A visibility change on curve C occurs if it intersects a BC/CC, and if the BC/CC point in 3D space is closer to the viewer than the curve point.

The visibility change is $\pm 2$ in the case of a CC, it is $\pm 1$ in the case of a BC.

The sign is determined as follows:
+ CC case : the sign is $-$ if ``inner(T,dir)>0``, where T is the tangent to the OCCLUDED curve in the view plane at the crossing, and dir is the normal to the CC in the view plane, pointing in the direction of the surface.
+ BC case : the sign is $-$ if ``inner(T,dir)>0``, where T is the tangent to the OCCLUDED curve in the view plane at the crossing, and dir is the normal to the BC in the view plane, pointing in the direction of the surface.

**Implementation notes**

- `T` is the **analytic image-space tangent** stored in `rc.tan` at each sample (populated by `resample_all`), interpolated to the crossing parameter. The earlier chord-based `T = xy[si+1] - xy[si]` is sample-jitter-noisy at fine resolutions and can flip `T·dir`'s sign when the curve is nearly parallel to the occluder's normal in image space.
- `dir` is computed as the **perpendicular-to-image-tangent component** of `rc.dir` (the lifted-uv-inward direction projected to image, stored on the occluder RC). Concretely:
  ```
  proj_coeff = (rc.dir @ rc.tan) / (rc.tan @ rc.tan)
  dir = rc.dir - proj_coeff * rc.tan
  ```
  When the surface tilts strongly in the depth direction at the boundary (e.g. paraboloid edges in a near-edge-on view), the raw `rc.dir` is NOT perpendicular to the curve's image tangent — the dS-lift picks up a depth contribution that contaminates the in-image direction. Projecting away the along-tangent component restores the geometrically-correct in-image inward normal while preserving the sign convention from `rc.dir`.

Visibility changes on curves resulting from these intersections (called BKs) are stored but curves are not split at these intersection points. 

# Visibility computation.
Once all the visibility changes are computed, the visibility of segments is propagated through the curves, starting from an anchor.

## Anchor and BFS propagation.
Consider the collection of all the points in all BC/CC/SIC curves, in the view plane. The leftmost point is the anchor; the **shared visibility** at its SP (or its visibility if it is not an SP) is zero.

**SP-aware anchor convention**: when the leftmost sample lands AT an endpoint of its RC (sample 0 or sample N-1) and that endpoint is an SP, the anchor's *sample value* is `rc.vc_in` (if at sample 0) or `rc.vc_out` (if at sample N-1), NOT zero. This sets the SP's shared visibility to 0 (the geometric meaning: an SP IS on the silhouette boundary, so its shared visibility is 0 by definition). If the anchor sample is INTERIOR to its RC (not at an SP), then its value is 0 directly. The same convention applies in the LP pass. Without this rule, several `disk_paraboloid_ca` cases were LP-INFEASIBLE because the leftmost sample landed at an SP with non-zero `vc_in`/`vc_out`.

Then the visibility is propagated:
+ Inside a curve, it changes at breaks.
+ Reaching the endpoint of the curve, the visibility of the endpoint is computed by *substracting* the visibility change at the endpoint.
+ When two curves share an endpoint, the visibility is propagated from one to the other.
+ When entering a curve from an endpoint, the visibility change at the endpoint is applied.

Thanks to helper curves, every curve is treated by this propagation (BFS) algorithm.

## LP pass to fix visibility inconsistencies.
Because numerical artifacts can mess things up, visibility can be determined otherwise, by solving a linear program. This involves treating the visibility changes at **breaks only**, not at curve endpoints, as variables, as well as the visibility of segments. The constraints are that visibility at anchor is $0$, and the propagation rule. The quantity to minimize is the sum of $|vc_i - vc_i^0|$, where $vc_i$ is the visibility change at break $i$ (a variable) and the visibility change $vc_i^0$ that was computed before. 

A variant could be to use 4 anchors with visibility $0$: the leftmost, upmost, rightmost, downmost points in the view plane. The choice between BFS, LP with one anchor, LP with 4 anchors is recorded in the PROPAGATION variable of settings.py.  



# Appendix : Edge-face intersections
## Problem
 An edge E1 can intersect a face F2 at its boundary (say edge E2), and then it is numerically ambiguous to determine wether E1 intersects F2, or another face F2' having E2 as a boundary. We ignore intersections at vertices, a situation which is too exceptional.

## Algorithm 
For every E1 and F2 that intersect at P:
- If dist(P, boundary of F2)>= delta, store a DP (E1,F2).
- Else, let E2 be the edge to which P is closest. Let F1, F1' be the faces to which E1 belongs. Let F2' be the other face to which E2 belongs.
+ If E2 intersects neither F1 or F1', store a DP (E1,E2). (note: this should happen only if E1 is a boundary edge). 
+ If E2 intersects, say, F1 at Q, and if dist(Q,E1) >= delta, store a DP (E2, F1).
+ If E2 intersects, say, F1 at Q, and if dist(Q,E1) < delta, store a DP (E1, E2).
- If a DP (E1, E2) was stored, then remove (E1, F2), (E1, F2'), (E2, F1), (E2, F1') from pairs to be tested.
- If a DP (E1, F2) or a DP(E1, F2') was stored, then remove both (E1, F2), (E1, F2') from pairs to be tested.

# TODO

- **Resampling of very small curves produces too many points.** When a
  SubCurve has projected length much smaller than the surface diameter
  `M`, the per-SP density `d = M/RESOLUTION` is much larger than the SC
  itself, but the current arclength sampler still emits an endpoint pair
  plus internal density-based samples, which oversamples short SCs.
  Revisit `curves.py:_sample_arclengths` so that short SCs degrade
  gracefully to a 2-point line.

# Appendix : Debug panel
The debug panel of the front-end allows to redefined the settings.py variables mentionned above. When a change is made, the construction pipeline and the outline pipeline are triggered. The changes survive when another surface is chosen.

The debug panel can be displayed using a button on the app banner, which is displayed only in the dev server, not on the deployed app. There, the variables are set by environment variables with the same name. 




