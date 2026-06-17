"""curves.py — P4: make_lines, P6: sign_changes, C7: build_bcs, O14: resample_all"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, Optional

import numpy as np

if TYPE_CHECKING:
    from surface_play.mesh import Mesh
    from surface_play.projection import Projection
    from surface_play.splitting import SplitArrays, SubCurve
    from surface_play.surface import SurfaceParams


def sign_changes(vals_p: np.ndarray, vals_q: np.ndarray,
                 flip: np.ndarray | None = None) -> np.ndarray:
    """
    Returns boolean mask shape (N,): True where vals_p[i] * vals_q[i] * flip[i] < 0.
    flip defaults to all +1. flip ∈ {-1, +1} per segment (Möbius mesh edges have -1).
    """
    vals_p = np.asarray(vals_p)
    vals_q = np.asarray(vals_q)
    if flip is None:
        return (vals_p * vals_q) < 0
    return (vals_p * vals_q * np.asarray(flip)) < 0


def make_lines(segments: np.ndarray) -> list[np.ndarray]:
    """
    segments: (N, 2) int array — each row is a pair of endpoint indices.
    Returns: list of 1D int arrays. Each array is a chain of segment indices,
             negative values denote reversed traversal. Closed chains end where
             they began (first index == last index in absolute value).
    Only valid when every vertex has degree ≤ 2 (no branch points).
    """
    segments = np.asarray(segments, dtype=np.intp)
    if segments.ndim != 2 or segments.shape[1] != 2:
        if segments.size == 0:
            return []
        raise ValueError("segments must be (N, 2)")
    N = len(segments)
    if N == 0:
        return []

    # Build half-edge table: for each endpoint, store (vertex, signed_seg_id)
    # signed_seg_id = +(i+1)  means "segment i entered from p-side (forward)"
    # signed_seg_id = -(i+1)  means "segment i entered from q-side (reversed)"
    # We use i+1 (1-indexed) so that 0 is not ambiguous with sign.
    p = segments[:, 0]
    q = segments[:, 1]

    # For vertex v, forward half-edge of seg i: arriving at v=p[i], we came
    # from q[i], so the segment sign to record in the chain is -(i+1) (reversed).
    # For vertex v, backward half-edge of seg i: arriving at v=q[i], we came
    # from p[i], so chain sign is +(i+1) (forward).
    #
    # Table columns: [vertex, signed_seg_id]
    # signed_seg_id encodes both which segment and in which direction we
    # *leave* that vertex along this segment.
    # Leaving v=p[i] → forward (+): signed = +(i+1)
    # Leaving v=q[i] → backward (-): signed = -(i+1)

    idx = np.arange(N, dtype=np.intp)
    he_vertex = np.concatenate([p, q])               # (2N,)
    he_signed = np.concatenate([idx + 1, -(idx + 1)])  # (2N,) leaving p fwd, leaving q rev

    # Sort half-edges by vertex
    order = np.argsort(he_vertex, kind="stable")
    he_vertex = he_vertex[order]
    he_signed = he_signed[order]

    # For each vertex collect its half-edges (its adjacency list)
    # We use searchsorted to find groups
    unique_verts, counts = np.unique(he_vertex, return_counts=True)

    # Build adjacency: adj[v] = list of signed seg ids leaving v
    adj: dict[int, list[int]] = {}
    pos = 0
    for v, c in zip(unique_verts, counts):
        adj[int(v)] = he_signed[pos:pos + c].tolist()
        pos += c

    # Degree of each vertex
    degree = {v: len(lst) for v, lst in adj.items()}

    visited = np.zeros(N, dtype=bool)
    chains: list[np.ndarray] = []

    def seg_id(signed: int) -> int:
        return abs(signed) - 1

    def other_end(seg: int, from_v: int) -> int:
        """Return the other endpoint of segment seg given we came from from_v."""
        if segments[seg, 0] == from_v:
            return int(segments[seg, 1])
        return int(segments[seg, 0])

    def traverse_from(start_v: int, start_signed: int) -> np.ndarray | None:
        """Walk a chain starting at start_v along start_signed."""
        chain = []
        prev_v = start_v
        cur_signed = start_signed

        while True:
            s = seg_id(cur_signed)
            if visited[s]:
                break
            visited[s] = True
            chain.append(cur_signed)

            next_v = other_end(s, prev_v)

            # Find the continuation half-edge at next_v (the one that is NOT
            # the reverse of the edge we just traversed)
            # The reverse of cur_signed leaving prev_v is the same edge leaving
            # next_v with opposite sign.
            reverse_signed = -(cur_signed) if cur_signed > 0 else -cur_signed
            # Wait, let me think again.
            # cur_signed leaving prev_v: if cur_signed = +(i+1), that means leaving p[i]
            # The corresponding half-edge arriving at q[i]=next_v is -(i+1)
            # So the "incoming" token at next_v is -(cur_signed) ... no.
            # Token +(i+1) = leaving p[i] forward.
            # Token -(i+1) = leaving q[i] backward.
            # If we used +(i+1) to leave prev_v (= p[i]), then next_v = q[i],
            # and the token -(i+1) represents leaving q[i] (= next_v) backward.
            # That is the reverse, which we want to skip.
            incoming_at_next = -(cur_signed)  # the reverse token at next_v

            neighbors = adj.get(next_v, [])
            # continuation = the other neighbor (not the incoming reverse)
            continuations = [t for t in neighbors if t != incoming_at_next]

            if len(continuations) == 0:
                # degree-1 endpoint: chain ends here
                break
            if len(continuations) == 1:
                next_signed = continuations[0]
                s_next = seg_id(next_signed)
                if visited[s_next]:
                    # closed loop: the chain closes back on itself
                    chain.append(next_signed)
                    break
                prev_v = next_v
                cur_signed = next_signed
            else:
                # branch point (degree ≥ 3) — stop
                break

        return np.array(chain, dtype=np.intp) if chain else None

    # First pass: start from degree-1 vertices (open chain endpoints)
    for v in list(adj.keys()):
        if degree[v] == 1:
            signed_edge = adj[v][0]
            s = seg_id(signed_edge)
            if not visited[s]:
                chain = traverse_from(v, signed_edge)
                if chain is not None and len(chain):
                    chains.append(chain)

    # Second pass: closed loops (all remaining unvisited segments)
    for i in range(N):
        if not visited[i]:
            # Start from p[i], forward
            chain = traverse_from(int(p[i]), i + 1)
            if chain is not None and len(chain):
                chains.append(chain)

    return chains


@dataclass
class BoundaryCurve:
    edge_indices: np.ndarray  # 1D int array; signs encode reversal as in make_lines output
    is_closed: bool


def build_bcs(mesh: Mesh) -> list[BoundaryCurve]:
    """
    Assemble boundary curves from mesh boundary edges. Vertex indices on edges
    are already canonical post-compaction, so chain nodes are taken directly
    from `p_idx`/`q_idx` without a `vertex_class` indirection.
    """
    if len(mesh.boundary_edge_idx) == 0:
        return []

    bnd_edges = mesh.edges[mesh.boundary_edge_idx]
    segments = np.column_stack([
        bnd_edges["p_idx"].astype(np.intp),
        bnd_edges["q_idx"].astype(np.intp),
    ])

    chains = make_lines(segments)

    result = []
    for chain in chains:
        is_closed = bool(chain[0] == chain[-1])
        result.append(BoundaryCurve(
            edge_indices=np.asarray(chain, dtype=np.intp),
            is_closed=is_closed,
        ))
    return result


# ── O14: resample_all ────────────────────────────────────────────────────────

@dataclass
class ResampledCurve:
    """Resampled SubCurve in projected space.

    `start`/`end` mirror the parent SubCurve's SP indices (G5 — Python-identity
    preserved through `splits.sps[idx]`); `-1` for SP-less closed RCs.
    `vc_in`/`vc_out` are carried over from the parent SubCurve for O16 BFS
    (visibility change applied at the start/end SP).
    """

    kind: Literal["BC", "CC", "SIC", "HC"]
    start: int
    end: int
    depth: np.ndarray             # (N,)   z along view axis (anchor-relative)
    xy: np.ndarray                # (N, 2) projected
    dir: Optional[np.ndarray]     # (N, 2) inward NORMAL in image; BC/CC only
    tan: Optional[np.ndarray] = None  # (N, 2) analytic TANGENT in image; BC/CC only
    vc_in: int = 0
    vc_out: int = 0
    uv: Optional[np.ndarray] = None   # (N, 2) canonical domain preimage of each sample


def _seg_uv_at_bary(sub_kind: str, seg_idx: int, bary: float,
                    mesh, css, sis_pairs, cps, dps) -> np.ndarray:
    """uv at bary on a BC/CC chain segment. `bary` is segment-local (0..1).

    SIC segments are NOT handled here: a SIC is one curve of DPs whose two
    preimages are deduced per-SIS via `flip` (spec §SIC, lines 219/231) — see
    `_sic_preAB_at`. `sis_pairs`/`dps` are kept in the signature for caller
    compatibility (probes pass the same positional args for all kinds).
    """
    if sub_kind == "BC":
        edge = mesh.edges[seg_idx]
        return np.asarray(edge["p"], dtype=float) + float(bary) * np.asarray(edge["pq"], dtype=float)
    if sub_kind == "CC":
        cs = css[seg_idx]
        p_uv = np.asarray(cps[int(cs["p_cp"])]["uv"], dtype=float)
        q_uv = np.asarray(cps[int(cs["q_cp"])]["uv"], dtype=float)
        return (1.0 - float(bary)) * p_uv + float(bary) * q_uv
    raise ValueError(f"_seg_uv_at_bary: unsupported kind {sub_kind!r}")


def _seg_uv_at_bary_batch(sub_kind: str, seg_idxs: np.ndarray, barys: np.ndarray,
                          mesh, css, cps) -> np.ndarray:
    """Vectorized `_seg_uv_at_bary` over many (seg_idx, bary) on a BC/CC chain.

    Byte-identical to looping `_seg_uv_at_bary` (same per-element float ops) but
    one array pass instead of thousands of scalar calls — the `_build_polyline`
    hot path. Returns (M, 2).
    """
    barys = np.asarray(barys, dtype=float)
    if sub_kind == "BC":
        edges = mesh.edges[np.asarray(seg_idxs, dtype=np.int64)]
        return (np.asarray(edges["p"], dtype=float)
                + barys[:, None] * np.asarray(edges["pq"], dtype=float))
    if sub_kind == "CC":
        rows = css[np.asarray(seg_idxs, dtype=np.int64)]
        p_uv = np.asarray(cps["uv"][rows["p_cp"]], dtype=float)
        q_uv = np.asarray(cps["uv"][rows["q_cp"]], dtype=float)
        return (1.0 - barys)[:, None] * p_uv + barys[:, None] * q_uv
    raise ValueError(f"_seg_uv_at_bary_batch: unsupported kind {sub_kind!r}")


def _needs_close(domain) -> bool:
    """True iff `domain.close` can actually move a point — i.e. a rect domain
    with at least one identified axis. For unidentified rect (and disk/annulus)
    `close` is a no-op, so callers skip the call entirely (it dominated the
    resample profile at ~100k no-op invocations on non-periodic surfaces)."""
    return domain is not None and bool(getattr(domain, "needs_close", False))


def _close_aware_lerp(a_uv: np.ndarray, b_uv: np.ndarray, t: float, domain) -> np.ndarray:
    if domain is None:
        return np.asarray(a_uv, dtype=float) + float(t) * (np.asarray(b_uv, dtype=float) - np.asarray(a_uv, dtype=float))
    return domain.interpolate(a_uv, b_uv, t)


def _build_polyline(
    sub: "SubCurve", splits: "SplitArrays", mesh: "Mesh",
    css: np.ndarray, sis_pairs: np.ndarray,
    cps: np.ndarray, dps: np.ndarray,
    surface: "SurfaceParams", projection: "Projection",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (uv_poly, xyz_poly, xy_poly) — vertices of the SC's polyline.

    HC: just [start, end] (2 vertices). SP-less closed SC: walk internal only.
    """
    # SIC is one curve of DPs — node skeleton from DP/SP xyz (sheet-independent);
    # preimages are deduced per-SIS via `flip` at resample time, not here.
    if sub.kind == "SIC":
        return _sic_node_polyline(sub, sis_pairs, dps, splits, projection)

    # Internal samples, vectorized (was a per-segment _seg_uv_at_bary loop —
    # the dominant _build_polyline cost). HC has none.
    if sub.kind != "HC" and sub.internal:
        _seg = np.fromiter((int(e[0]) for e in sub.internal), dtype=np.int64,
                           count=len(sub.internal))
        _bar = np.fromiter((float(e[1]) for e in sub.internal), dtype=float,
                           count=len(sub.internal))
        iu = _seg_uv_at_bary_batch(sub.kind, _seg, _bar, mesh, css, cps)
    else:
        iu = np.zeros((0, 2), dtype=float)

    if sub.kind == "HC":
        uv_arr = np.array([splits.sps[sub.start][0], splits.sps[sub.end][0]],
                          dtype=float)
    elif sub.start == -1 and sub.end == -1:
        # SP-less closed SC — internal verbatim; close the loop.
        uv_arr = np.vstack([iu, iu[:1]]) if len(iu) else np.zeros((0, 2), dtype=float)
    else:
        # Normal BC/CC SC: start SP → internal → end SP.
        s0 = np.asarray(splits.sps[sub.start][0], dtype=float).reshape(1, 2)
        s1 = np.asarray(splits.sps[sub.end][0], dtype=float).reshape(1, 2)
        uv_arr = np.vstack([s0, iu, s1])

    # Close-aware adjust consecutive vertices, then lift to xyz/xy.
    # Only for GLOBAL periodicities (rect cy/mo), where S is periodic so the
    # aliased vertex is the same surface point. NOT for antipodal: there close()
    # is a boundary-only identification (S(x) ≠ S(-x) in the interior), so
    # reflecting an interior arc would relocate it onto different surface points
    # (e.g. a CC tail near a cusp jumping to the antipodal preimage). A
    # seam-crossing polyline keeps its true uv; its two boundary endpoints P and
    # -P already lift to one xy via S(P) = S(-P), so no long chord appears.
    domain = getattr(mesh, "domain", None)
    if (_needs_close(domain) and not getattr(domain, "is_antipodal", False)
            and len(uv_arr) > 1):
        for i in range(1, len(uv_arr)):
            uv_arr[i] = domain.interpolate(uv_arr[i - 1], uv_arr[i], 1.0)
    if len(uv_arr):
        S_batch = np.asarray(
            surface.S(uv_arr[:, 0], uv_arr[:, 1]), dtype=float,
        )  # (3, N) from vectorized lambdified call.
        xyz_arr = np.ascontiguousarray(S_batch.T)
        xy_arr = projection.XY(xyz_arr)
    else:
        xyz_arr = np.zeros((0, 3), dtype=float)
        xy_arr = np.zeros((0, 2), dtype=float)
    return uv_arr, xyz_arr, xy_arr


def _arclengths(xy_poly: np.ndarray) -> np.ndarray:
    """Cumulative xy-arclength of a polyline; shape (N,) with [0]==0."""
    if len(xy_poly) < 2:
        return np.zeros(len(xy_poly), dtype=float)
    diffs = np.linalg.norm(np.diff(xy_poly, axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(diffs)])


def _interp_along_polyline(
    uv_poly: np.ndarray, xy_poly: np.ndarray, cum: np.ndarray,
    s_target: float, domain,
) -> tuple[np.ndarray, int, float]:
    """Return (uv, seg_index, alpha) at xy-arclength `s_target` on polyline.

    Uses xy-linear interpolation for the parameter α, then close-aware-lerps uv.
    """
    s_target = float(np.clip(s_target, 0.0, cum[-1]))
    seg = int(np.searchsorted(cum, s_target, side="right") - 1)
    seg = max(0, min(seg, len(cum) - 2))
    span = cum[seg + 1] - cum[seg]
    if span <= 0:
        alpha = 0.0
    else:
        alpha = (s_target - cum[seg]) / span
    uv_a = uv_poly[seg]
    uv_b = uv_poly[seg + 1]
    uv = _close_aware_lerp(uv_a, uv_b, alpha, domain)
    return uv, seg, alpha


def _densify_bc_polyline(uv_p: np.ndarray, surface, projection, domain,
                         n_sub: int):
    """Densify a BC build-polyline for arclength-accurate resampling.

    The BC build-polyline has one vertex per boundary mesh edge — far too
    coarse to track the true projected curve near a *projection fold* (an
    apparent-contour extremity / BCP), where uv→image is near-singular and the
    surface bends sharply away from the straight build-polyline chord. The
    coarse `cum = _arclengths(xy_p)` then mis-measures arclength inside the
    first segment, so the BCP arclength-match (`_bc_s_targets` inheriting the
    CC ladder) reprojects to the wrong positions and the BC/CC image polylines
    weave (spurious projection-break visibility changes). See
    [[bfs-foldtip-weaving-2026-05-30]].

    Each original segment is subdivided into exactly `n_sub` pieces *linear in
    uv* (a BC segment is a straight uv line along the domain boundary), so the
    original segment that owns dense segment `d` is simply `d // n_sub` — this
    preserves the seg→mesh-edge mapping the BC inward-normal `dir` relies on.
    Mirrors the dense-reparam the HC branch already uses (curves.py HC path).

    Returns `(uv_dense, xy_dense, cum_dense)` where `cum_dense` is an accurate
    cumulative xy-arclength. The original vertex k sits at dense index
    `k * n_sub`.
    """
    M = len(uv_p)
    if M < 2 or n_sub < 1:
        xy = projection.XY(np.ascontiguousarray(
            np.asarray(surface.S(uv_p[:, 0], uv_p[:, 1]), dtype=float).T))
        return uv_p.copy(), xy, _arclengths(xy)
    # Build dense uv: start vertex + n_sub interior+end points per segment.
    pieces = [uv_p[:1]]
    for k in range(M - 1):
        a = uv_p[k]; b = uv_p[k + 1]
        t = (np.arange(1, n_sub + 1) / float(n_sub))[:, None]
        pieces.append(a[None, :] + t * (b - a)[None, :])
    uv_dense = np.ascontiguousarray(np.vstack(pieces))
    S = np.asarray(surface.S(uv_dense[:, 0], uv_dense[:, 1]), dtype=float)
    xyz_dense = np.ascontiguousarray(S.T)
    xy_dense = projection.XY(xyz_dense)
    cum_dense = _arclengths(xy_dense)
    return uv_dense, xy_dense, cum_dense


def _snap_annular_bc(uv: np.ndarray, mesh) -> np.ndarray:
    """Snap a uv to ‖uv‖ = r_min or r_max on a disk/annulus boundary.

    Assumes domain.coord_type == 'ca' (cartesian uv with norm = radius).
    """
    domain = mesh.domain
    r_min = float(domain.bounds[0])
    r_max = float(domain.bounds[1])
    r = float(np.linalg.norm(uv))
    if r == 0.0:
        return uv
    target = r_max if abs(r - r_max) <= abs(r - r_min) else r_min
    return uv * (target / r)


def _uv_for_bc_lift(edge, mesh, uv_sample: np.ndarray) -> np.ndarray:
    """Bring `uv_sample` into the canonical p-copy of `edge` for Su/Sv lift.

    Used in resample_all Phase 2 to build `uv_eval`, so the batched BC dir/tan
    share one seam-canonical uv (cf. [[bc_lift_patch_match_2026_05_26]]).
    """
    domain = getattr(mesh, "domain", None)
    if _needs_close(domain):
        p_canonical = mesh.uv[int(edge["p_idx"])]
        return domain.interpolate(p_canonical, uv_sample, 1.0)
    return uv_sample


def _tan_for_cc_samples_batched(
    uv_samples: np.ndarray,    # (N, 2)
    chain_segs: np.ndarray,    # (N,) int — CC chain segment index per sample
    S_all: np.ndarray,         # (N, 3)
    Su_all: np.ndarray,        # (N, 3)
    Sv_all: np.ndarray,        # (N, 3)
    Suu_all: np.ndarray,       # (N, 3)
    Suv_all: np.ndarray,       # (N, 3)
    Svv_all: np.ndarray,       # (N, 3)
    css: np.ndarray, cps: np.ndarray,
    projection,
) -> np.ndarray:
    """Analytic image-space TANGENT of the CC silhouette at each uv, sign-matched
    to the cs chain direction (p_cp → q_cp).

    Math: the silhouette curve in uv is `axis · SN = 0`, so its uv-tangent is
    perpendicular to the gradient
        Np = ∇_uv(axis · SN)
           = ((Suu×Sv + Su×Suv)·axis, (Suv×Sv + Su×Svv)·axis),
    giving Tp_uv = (-Np_y, Np_x). Lifted via dS (Tp_uv[0]·Su + Tp_uv[1]·Sv) and
    projected to the image; the sign is flipped if Tp_uv·(q_uv - p_uv) < 0.
    In persp `axis = viewer_direction(S(uv))`; the extra gradient terms
    −Su·SN / −Sv·SN vanish identically (triple product with a repeated vector),
    so the Np formula is unchanged.

    Implementation: uses the scalar-triple-product identity
    `(a × b) · c = det([a, b, c])` to avoid the per-sample `np.cross + @`
    overhead, and inlines `proj_vec` so every sample is processed in one
    numpy pass.
    """
    N = uv_samples.shape[0]
    if N == 0:
        return np.empty((0, 2), dtype=float)

    axis = np.asarray(projection.viewer_direction(S_all), dtype=float)  # (N, 3)

    def _stp(a, b, c):  # scalar triple product, batched over leading axis
        return (
            a[:, 0] * (b[:, 1] * c[:, 2] - b[:, 2] * c[:, 1])
            + a[:, 1] * (b[:, 2] * c[:, 0] - b[:, 0] * c[:, 2])
            + a[:, 2] * (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])
        )

    Np_x = _stp(Suu_all, Sv_all, axis) + _stp(Su_all, Suv_all, axis)
    Np_y = _stp(Suv_all, Sv_all, axis) + _stp(Su_all, Svv_all, axis)
    Tp_uv = np.empty((N, 2), dtype=float)
    Tp_uv[:, 0] = -Np_y
    Tp_uv[:, 1] = Np_x

    cs_rows = css[chain_segs]
    p_uv = np.asarray(cps[cs_rows["p_cp"]]["uv"], dtype=float)
    q_uv = np.asarray(cps[cs_rows["q_cp"]]["uv"], dtype=float)
    chain_dir = q_uv - p_uv
    dots = Tp_uv[:, 0] * chain_dir[:, 0] + Tp_uv[:, 1] * chain_dir[:, 1]
    flip = dots < 0.0
    if flip.any():
        Tp_uv[flip] = -Tp_uv[flip]

    Tp_3d = Tp_uv[:, 0:1] * Su_all + Tp_uv[:, 1:2] * Sv_all  # (N, 3)

    # Batched equivalent of projection.proj_vec(uv, Tp_3d).
    I, J = projection.I, projection.J
    if projection.mode == "ortho":
        out = np.empty((N, 2), dtype=float)
        out[:, 0] = Tp_3d @ I
        out[:, 1] = Tp_3d @ J
        return out
    eye = projection.eye
    n_axis = projection._axis
    d = S_all - eye
    z = -(d @ n_axis)
    if np.any(z == 0.0):
        bad = int(np.argmax(z == 0.0))
        raise ValueError(
            f"proj_vec undefined: S(uv)={S_all[bad]} lies on the image plane through eye"
        )
    a = d @ I
    b = d @ J
    nv = Tp_3d @ n_axis
    out = np.empty((N, 2), dtype=float)
    out[:, 0] = ((Tp_3d @ I) + (a / z) * nv) / z
    out[:, 1] = ((Tp_3d @ J) + (b / z) * nv) / z
    return out


def _boundary_tangent_batched(domain, uv: np.ndarray, edge_dp: np.ndarray) -> np.ndarray:
    """Batched `domain.boundary_tangent` over (M, 2) uv / edge_dp.

    Byte-identical to the per-sample scalar version (same formulas, numpy
    elementwise): rect → normalized edge_dp; disk/annulus → boundary-circle
    tangent (-v, u)/r sign-matched to edge_dp (r == 0 falls back to edge_dp).
    """
    M = uv.shape[0]
    if domain is None or getattr(domain, "type", None) == "rect":
        n = np.linalg.norm(edge_dp, axis=1, keepdims=True)
        return np.where(n > 0.0, edge_dp / np.where(n > 0.0, n, 1.0), 0.0)
    u = uv[:, 0]
    v = uv[:, 1]
    r = np.sqrt(u * u + v * v)
    Tan = np.empty((M, 2), dtype=float)
    nz = r != 0.0
    Tan[nz, 0] = -v[nz] / r[nz]
    Tan[nz, 1] = u[nz] / r[nz]
    if (~nz).any():
        edp = edge_dp[~nz]
        n = np.linalg.norm(edp, axis=1, keepdims=True)
        Tan[~nz] = np.where(n > 0.0, edp / np.where(n > 0.0, n, 1.0), 0.0)
    dots = Tan[:, 0] * edge_dp[:, 0] + Tan[:, 1] * edge_dp[:, 1]
    flip = dots < 0.0
    Tan[flip] = -Tan[flip]
    return Tan


def _proj_vec2_batched(projection, S: np.ndarray, v1: np.ndarray, v2: np.ndarray):
    """Batched `projection.proj_vec` for two vector fields v1, v2 (M, 3) at the
    same surface points S (M, 3). Returns (proj_v1, proj_v2), each (M, 2).

    `S` is the precomputed surface position (== surface.S(uv), which is
    _eval_all(uv)[0]); reusing it keeps the persp depth byte-identical to the
    per-sample scalar proj_vec without a second eval. Inlines the same formula
    as `_tan_for_cc_samples_batched`'s proj step.
    """
    I, J = projection.I, projection.J
    M = S.shape[0]
    if projection.mode == "ortho":
        o1 = np.empty((M, 2), dtype=float)
        o2 = np.empty((M, 2), dtype=float)
        o1[:, 0] = v1 @ I; o1[:, 1] = v1 @ J
        o2[:, 0] = v2 @ I; o2[:, 1] = v2 @ J
        return o1, o2
    eye = projection.eye
    n_axis = projection._axis
    d = S - eye
    z = -(d @ n_axis)
    if np.any(z == 0.0):
        bad = int(np.argmax(z == 0.0))
        raise ValueError(
            f"proj_vec undefined: S(uv)={S[bad]} lies on the image plane through eye"
        )
    a = d @ I
    b = d @ J

    def _pv(v):
        nv = v @ n_axis
        o = np.empty((M, 2), dtype=float)
        o[:, 0] = ((v @ I) + (a / z) * nv) / z
        o[:, 1] = ((v @ J) + (b / z) * nv) / z
        return o

    return _pv(v1), _pv(v2)


def _newton_cc_refine(uv: np.ndarray, surface, projection, *, max_iter: int = 50
                      ) -> np.ndarray:
    """Run Newton in uv to land the point on the silhouette `axis·SN = 0`.

    Step = uv-gradient descent of `f(uv) = axis·SN`, using ANALYTIC derivatives
    `dSN/du = Suu×Sv + Su×Suv`, `dSN/dv = Suv×Sv + Su×Svv` (same as
    `_newton_orthogonal_cp`). In perspective the `d(axis)/duv` term vanishes
    because `SN ⊥ Su, Sv`, so the gradient is just `axis · dSN/d·`.

    Iterates to convergence — until the Newton step no longer moves `uv` at
    double precision (or the residual is negligible) — not a fixed count. The
    analytic gradient lets it reach machine precision in a few iterations; a
    finite-difference gradient would plateau at the fd-step floor. `max_iter`
    is a generous safety cap against a non-converging / oscillating point.

    Persp (2026-05-27): `axis = viewer_direction(S(uv))` recomputed each
    iteration as `uv` moves.
    """
    uv_cur = uv.copy()
    for _ in range(max_iter):
        u, v = float(uv_cur[0]), float(uv_cur[1])
        vals = surface._eval_all(u, v)
        S_p = np.asarray(vals[0], dtype=float).reshape(3)
        Su  = np.asarray(vals[1], dtype=float).reshape(3)
        Sv  = np.asarray(vals[2], dtype=float).reshape(3)
        Suu = np.asarray(vals[3], dtype=float).reshape(3)
        Suv = np.asarray(vals[4], dtype=float).reshape(3)
        Svv = np.asarray(vals[5], dtype=float).reshape(3)
        SN  = np.asarray(vals[6], dtype=float).reshape(3)
        axis = projection.viewer_direction(S_p).reshape(3)
        f = float(axis @ SN)
        if not np.isfinite(f) or abs(f) < 1e-13:
            break
        dSN_du = np.cross(Suu, Sv) + np.cross(Su, Suv)
        dSN_dv = np.cross(Suv, Sv) + np.cross(Su, Svv)
        gu = float(axis @ dSN_du)
        gv = float(axis @ dSN_dv)
        g2 = gu * gu + gv * gv
        if g2 < 1e-20:
            break
        du = f * gu / g2
        dv = f * gv / g2
        uv_cur[0] -= du
        uv_cur[1] -= dv
        # Converged: the step no longer moves uv at double precision.
        if du * du + dv * dv < 1e-26:
            break
    return uv_cur


def _sample_arclengths(L_total: float, ell: float,
                       delta_start: float, delta_end: float,
                       is_closed: bool) -> np.ndarray:
    """Pick xy-arclength positions for resampling (spec change 2026-05-27).

    Per-half growing-spacing rule:
      - Each sub is split into a start-half and end-half at L_total/2.
      - In the START half, walking inward from the start SP, the k-th
        segment (k=1,2,...) has length `d_k = min(k * delta_start, ell)`.
      - Symmetrically in the END half from the end SP, with `delta_end`.
      - `delta_X = L_min_incident_at_X / 10` (per-SP).
      - `ell = M / resolution` (global).

    Effect: segments grow linearly from delta near the SP until saturating
    at ell. Catches near-SP events densely and transitions smoothly to the
    coarse interior. No "blind zone" between dense and coarse regions.

    Closed curves: uniform `ell` spacing.
    """
    if L_total <= 0:
        return np.zeros(0, dtype=float)

    if is_closed:
        coarse = ell if ell > 0 else L_total / 30.0
        n = max(int(round(L_total / coarse)), 4)
        return np.linspace(0.0, L_total, n, endpoint=False)

    coarse = ell if ell > 0 else max(delta_start, delta_end, L_total / 30.0)
    d_s = delta_start if delta_start > 0 else coarse
    d_e = delta_end   if delta_end   > 0 else coarse
    half = L_total / 2.0

    # Forward from start: 0, d_1, d_1+d_2, ... until reaching `half`.
    start: list[float] = [0.0]
    s = 0.0
    k = 1
    while s < half:
        dk = min(k * d_s, coarse)
        if dk <= 0:
            break
        s += dk
        if s > half:
            break
        start.append(s)
        k += 1

    # Backward from end: L_total, L_total - d_1, L_total - (d_1+d_2), ...
    end_rev: list[float] = [L_total]
    s = 0.0
    k = 1
    while s < half:
        dk = min(k * d_e, coarse)
        if dk <= 0:
            break
        s += dk
        if s > half:
            break
        end_rev.append(L_total - s)
        k += 1
    end = list(reversed(end_rev))

    return np.array(sorted(set(start + end)))


def _sub_outgoing_xy_tangent(sub, sp_idx: int, xy_p: np.ndarray):
    """Image-space outgoing tangent of a non-HC sub at SP `sp_idx`, normalised
    to point AWAY from SP into the curve. Returns None if the polyline has no
    non-zero-length step on the relevant side.
    """
    if len(xy_p) < 2:
        return None
    if sp_idx == sub.start:
        anchor = xy_p[0]
        for k in range(1, len(xy_p)):
            d = xy_p[k] - anchor
            n = float(np.linalg.norm(d))
            if n > 0:
                return d / n
    elif sp_idx == sub.end:
        anchor = xy_p[-1]
        for k in range(len(xy_p) - 2, -1, -1):
            d = xy_p[k] - anchor
            n = float(np.linalg.norm(d))
            if n > 0:
                return d / n
    return None


def _pick_neighbour_arclengths(
    sp_idx: int, T_self: np.ndarray, subcurves, polys, half_L: float,
    exclude_kinds=("HC",), self_idx: int | None = None,
    return_length: bool = False,
):
    """Tangent-pick: among subs incident at sp_idx whose kind is NOT in
    `exclude_kinds` and whose index is not `self_idx`, pick the one with the
    largest T_sub · T_self (must be > 0). Return its near-SP CP arclengths in
    (0, half_L]. Returns None when no positive-aligned neighbour exists
    (degenerate; caller falls back to its own ladder).

    If `return_length` is True, return `(arcs, L_neighbour)` instead of `arcs`
    (and `(None, None)` on the degenerate paths) so the caller can cap the
    matched region at the neighbour's own half (`_bc_s_targets`).
    """
    if T_self is None:
        return (None, None) if return_length else None
    exclude_set = set(exclude_kinds)
    candidates = []
    for i, sub in enumerate(subcurves):
        if sub.kind in exclude_set:
            continue
        if self_idx is not None and i == self_idx:
            continue
        if sub.start != sp_idx and sub.end != sp_idx:
            continue
        xy_p = polys[i][2]
        T = _sub_outgoing_xy_tangent(sub, sp_idx, xy_p)
        if T is None:
            continue
        dotp = float(np.dot(T, T_self))
        candidates.append((dotp, i, sub, xy_p))
    if not candidates:
        return (None, None) if return_length else None
    candidates.sort(key=lambda x: (-x[0], x[1]))
    best_dot, best_i, best_sub, xy_p = candidates[0]
    if best_dot <= 0:
        return (None, None) if return_length else None
    diffs = xy_p[1:] - xy_p[:-1]
    seg_len = np.sqrt(np.einsum("ij,ij->i", diffs, diffs))
    cum = np.concatenate(([0.0], np.cumsum(seg_len)))
    L_sub = float(cum[-1])
    if best_sub.start == sp_idx:
        arc_from_sp = cum
    else:
        arc_from_sp = L_sub - cum
    arcs = arc_from_sp[(arc_from_sp > 0) & (arc_from_sp <= half_L)]
    arcs_sorted = np.sort(np.unique(arcs))
    return (arcs_sorted, L_sub) if return_length else arcs_sorted


def _hc_s_targets(
    sub_hc, i_hc: int, subcurves, polys,
    surface, projection,
    L_xy_total: float, ell: float, delta_s: float, delta_e: float,
) -> np.ndarray:
    """HC sample arclengths using the tangent-pick / arclength-match rule.

    At each of the HC's two endpoint SPs, find the non-HC sub incident there
    whose outgoing image tangent is closest to T_HC (largest dot product).
    Inherit its near-SP CP arclengths in (0, L_xy_total/2]. If no such
    neighbour exists (no positive-dot alignment), fall back to the standard
    `_sample_arclengths` ladder on that half. See [[resume-hc-match-cc]].
    """
    uv_q0 = polys[i_hc][0][0]
    uv_q1 = polys[i_hc][0][1]
    duv = uv_q1 - uv_q0

    def _T_HC_at(uv_sp: np.ndarray, duv_signed: np.ndarray):
        Su = np.asarray(surface.Su(uv_sp[0], uv_sp[1]), dtype=float).ravel()
        Sv = np.asarray(surface.Sv(uv_sp[0], uv_sp[1]), dtype=float).ravel()
        T_3d = duv_signed[0] * Su + duv_signed[1] * Sv
        T_xy = projection.proj_vec(uv_sp, T_3d)
        n = float(np.linalg.norm(T_xy))
        return None if n == 0 else T_xy / n

    T_start = _T_HC_at(uv_q0, +duv)
    T_end   = _T_HC_at(uv_q1, -duv)
    half_L = L_xy_total / 2.0

    start_arcs = _pick_neighbour_arclengths(
        sub_hc.start, T_start, subcurves, polys, half_L,
        exclude_kinds=("HC",), self_idx=i_hc,
    )
    end_arcs = _pick_neighbour_arclengths(
        sub_hc.end, T_end, subcurves, polys, half_L,
        exclude_kinds=("HC",), self_idx=i_hc,
    )

    s_set = {0.0, float(L_xy_total)}
    std_ladder = None
    if start_arcs is None:
        std_ladder = _sample_arclengths(L_xy_total, ell, delta_s, delta_e, False)
        for s in std_ladder:
            if 0.0 < float(s) <= half_L:
                s_set.add(float(s))
    else:
        for a in start_arcs:
            s_set.add(float(a))
    if end_arcs is None:
        if std_ladder is None:
            std_ladder = _sample_arclengths(L_xy_total, ell, delta_s, delta_e, False)
        for s in std_ladder:
            if half_L <= float(s) < L_xy_total:
                s_set.add(float(s))
    else:
        for a in end_arcs:
            s_set.add(float(L_xy_total - a))

    return np.array(sorted(s_set), dtype=float)


def _bc_s_targets(
    sub_bc, i_bc: int, subcurves, polys, splits,
    L_total: float, own_cum: np.ndarray,
) -> np.ndarray:
    """BC sample arclengths via tangent-pick at BCP endpoints only.

    At each BC endpoint that is a **BCP** (BC × CC tangent boundary contact),
    tangent-pick the CC ending at the BCP whose outgoing image tangent aligns
    with the BC's outgoing tangent, and adopt that CC's vertex arclengths-from-
    the-SP as the BC's samples (REPLACING — not unioning — the BC's own raw
    vertices). Beyond that matched extent (the un-matched middle) the BC reverts
    to its own raw vertices (`own_cum`).

    Matched extent:
    - **Both ends match** a CC: each end takes its own half, capped at
      `min(L_BC/2, L_CC/2)`, so the two ladders meet at the midpoint without
      colliding. (For a lune — both endpoints share the same CC, near-equal
      lengths — the halves meet and the BC becomes a vertex-exact copy of the
      CC.)
    - **Only one end matches** (the other is a corner / BDP / CDP / non-BCP, or
      has no tangent-aligned CC): there is no competing rule at the far end, so
      the single match extends over the WHOLE BC, capped only by the CC's own
      length `min(L_BC, L_CC)`. Beyond `L_CC` (CC shorter than BC) the BC
      reverts to its own vertices.

    Rationale: at a BCP the BC and the tangent CC are near-coincident in the
    image; sampling them at *identical* arclengths-from-the-shared-origin makes
    their polylines coincide, so they cannot weave (the old UNION kept the BC's
    own vertices in the matched region, which bulged off the CC's chord and
    produced spurious grazing crossings — the Onde fold-tip "lune" ±4 defect).

    The rule does NOT apply at corners, BDPs, CDPs, or other non-BCP SPs.
    See [[resume-hc-match-cc]].
    """
    if sub_bc.start < 0 or sub_bc.end < 0:
        return own_cum.copy()
    xy_p = polys[i_bc][2]
    half_L = L_total / 2.0

    def _sp_is_bcp(sp_idx: int) -> bool:
        try:
            t = str(splits.sps[sp_idx][3])
        except (IndexError, TypeError):
            return False
        return t == "bcp"

    # Probe both endpoints for a tangent-aligned CC. Cap the picker at the full
    # BC length (L_total) so a single-end match can extend over the whole BC;
    # the per-end extent is applied below once we know whether the OTHER end
    # also has a rule.
    start_arcs = end_arcs = None
    L_cc_start = L_cc_end = 0.0
    if _sp_is_bcp(sub_bc.start):
        T_start = _sub_outgoing_xy_tangent(sub_bc, sub_bc.start, xy_p)
        start_arcs, L_cc_start = _pick_neighbour_arclengths(
            sub_bc.start, T_start, subcurves, polys, L_total,
            exclude_kinds=("BC", "HC", "SIC"), self_idx=i_bc,
            return_length=True,
        )
    if _sp_is_bcp(sub_bc.end):
        T_end = _sub_outgoing_xy_tangent(sub_bc, sub_bc.end, xy_p)
        end_arcs, L_cc_end = _pick_neighbour_arclengths(
            sub_bc.end, T_end, subcurves, polys, L_total,
            exclude_kinds=("BC", "HC", "SIC"), self_idx=i_bc,
            return_length=True,
        )

    has_start = start_arcs is not None
    has_end = end_arcs is not None

    # Decide the matched extent from each end. Both ends matching is the only
    # "conflict": split the BC at the midpoint (each half-capped at its CC's
    # own half). A lone match owns the whole BC, capped by its CC's length.
    ext_start = ext_end = 0.0
    if has_start and has_end:
        ext_start = min(half_L, L_cc_start / 2.0)
        ext_end = min(half_L, L_cc_end / 2.0)
    elif has_start:
        ext_start = min(L_total, L_cc_start)
    elif has_end:
        ext_end = min(L_total, L_cc_end)

    matched: list[float] = []  # arclengths-from-the-BC-start of matched samples
    if has_start:
        matched.extend(float(a) for a in start_arcs if a <= ext_start)
    if has_end:
        matched.extend(float(L_total - a) for a in end_arcs if a <= ext_end)

    s_set = {0.0, float(L_total)}
    s_set.update(matched)
    # Un-matched middle: BC's own raw vertices strictly outside both matched
    # regions [0, ext_start] and [L_total - ext_end, L_total].
    lo, hi = ext_start, L_total - ext_end
    s_set.update(float(c) for c in own_cum if lo < float(c) < hi)
    return np.array(sorted(s_set), dtype=float)


# ── SIC resampling (single curve of DPs; preimages deduced per-SIS via flip) ──
#
# Spec §SIC (lines 219, 231): an SIC is ONE curve; its two domain preimages are
# deduced per-DP by `flip`, NOT extracted as two resampled polylines. We mirror
# BC/CC resampling: the arclength skeleton is the DP/SP `xyz` polyline (a DP's
# `xyz` is sheet-independent, so it is single-valued and continuous regardless
# of how `flip` alternates). Each resampled sample is itself a DP: within a SIS
# we interpolate BOTH preimages (close-aware) using `flip` to pair the two DP
# endpoints' sheets, then lift `xyz = S(preimage-A)`. No sheet is ever "picked"
# by proximity, and no continuous preimage curve is materialised.

def _dp_of_internal(entry, sis_pairs) -> int:
    """DP index at the chain-forward end of an SIC `internal` segment entry.

    `entry = (sis_idx, end_bary)`; `end_bary` is 1.0 if the segment is forward
    in the chain (its q-end is the join vertex) else 0.0 (its p-end).
    """
    s = int(entry[0])
    row = sis_pairs[s]
    return int(row["q_dp"]) if float(entry[1]) > 0.5 else int(row["p_dp"])


_SIC_BARY_EPS = 1e-9


def _sic_preAB_at(s, bary, sis_pairs, dps, domain):
    """(preimage-A, preimage-B) uv of SIS `s` at native bary `bary`.

    Preimage-A is the sheet through `dps[p].uv1`; preimage-B through
    `dps[p].uv2`; `flip` selects the q-end sheet. Both preimage segments are
    interpolated close-aware at `bary` (bary 0 → p-end, 1 → q-end).
    """
    row = sis_pairs[s]
    p = int(row["p_dp"]); q = int(row["q_dp"]); f = int(row["flip"])
    uv1p = np.asarray(dps[p]["uv1"], dtype=float)
    uv2p = np.asarray(dps[p]["uv2"], dtype=float)
    qA = np.asarray(dps[q]["uv1" if f == 1 else "uv2"], dtype=float)
    qB = np.asarray(dps[q]["uv2" if f == 1 else "uv1"], dtype=float)
    A = _close_aware_lerp(uv1p, qA, float(bary), domain)
    B = _close_aware_lerp(uv2p, qB, float(bary), domain)
    return A, B


def _sic_owning_sis(vk, vk1, sis_pairs, dp_pair_to_sis, sp_to_sis):
    """Return `(sis_idx, bary_k, bary_k1)` for the polyline segment `vk → vk1`.

    `bary_*` are native barys on the OWNING SIS (0 = p-end, 1 = q-end). An SP
    sitting at a DP (its SPT bary is 0 or 1) is resolved to that DP and the
    segment's owning SIS is the DP-pair edge; an SP at interior bary keeps its
    own SIS. This is structural — no proximity.
    """
    def _dp_bary(dp, s):
        row = sis_pairs[s]
        if int(row["p_dp"]) == dp:
            return 0.0
        if int(row["q_dp"]) == dp:
            return 1.0
        return None

    def _as_dp(vert):
        """DP index if `vert` is a DP, or an SP located at a DP (bary 0/1)."""
        if vert[0] == "DP":
            return int(vert[1])
        for s, t in sp_to_sis.get(int(vert[1]), ()):
            if t <= _SIC_BARY_EPS:
                return int(sis_pairs[s]["p_dp"])
            if t >= 1.0 - _SIC_BARY_EPS:
                return int(sis_pairs[s]["q_dp"])
        return None

    def _isp(vert):
        """Interior-SP candidates [(sis, t)] for an SP not sitting at a DP."""
        if vert[0] == "DP":
            return []
        return [(s, t) for s, t in sp_to_sis.get(int(vert[1]), ())
                if _SIC_BARY_EPS < t < 1.0 - _SIC_BARY_EPS]

    dk, dk1 = _as_dp(vk), _as_dp(vk1)
    if dk is not None and dk1 is not None:
        s = dp_pair_to_sis[(min(dk, dk1), max(dk, dk1))]
        return s, _dp_bary(dk, s), _dp_bary(dk1, s)
    if dk is None and dk1 is not None:
        for s, t in _isp(vk):
            b = _dp_bary(dk1, s)
            if b is not None:
                return s, t, b
        raise RuntimeError(f"SIC: interior SP {vk[1]} shares no SIS with DP {dk1}")
    if dk is not None and dk1 is None:
        for s, t in _isp(vk1):
            b = _dp_bary(dk, s)
            if b is not None:
                return s, b, t
        raise RuntimeError(f"SIC: interior SP {vk1[1]} shares no SIS with DP {dk}")
    map_k1 = {s: t for s, t in _isp(vk1)}
    for s, t in _isp(vk):
        if s in map_k1:
            return s, t, map_k1[s]
    raise RuntimeError(
        f"SIC: interior SPs {vk[1]} and {vk1[1]} share no SIS")


def _sic_vertices(sub, sis_pairs):
    """Ordered vertices of an SIC SubCurve: list of `("SP", idx)` / `("DP", idx)`.

    SP-less closed subs are pure DP loops; otherwise SP-start, internal DPs,
    SP-end.
    """
    if sub.start == -1 and sub.end == -1:
        return [("DP", _dp_of_internal(e, sis_pairs)) for e in sub.internal]
    verts = [("SP", int(sub.start))]
    verts.extend(("DP", _dp_of_internal(e, sis_pairs)) for e in sub.internal)
    verts.append(("SP", int(sub.end)))
    return verts


def _sic_vert_xyz(vert, splits, dps) -> np.ndarray:
    if vert[0] == "DP":
        return np.asarray(dps[int(vert[1])]["xyz"], dtype=float)
    return np.asarray(splits.sps[int(vert[1])][1], dtype=float)


def _sic_build(sub, sis_pairs, dps, splits, dp_pair_to_sis, sp_to_sis, domain):
    """Build the SIC arclength skeleton + per-segment preimage endpoints.

    Returns `(node_xyz, seg_A, seg_B)` where `node_xyz` is the (M+1, 3) DP/SP
    xyz polyline and `seg_A[k] = (A0, A1)` / `seg_B[k] = (B0, B1)` are the two
    preimage-A / preimage-B endpoints of polyline segment k (already brought
    close-aware-consistent). A sample at fraction α in segment k has preimage-A
    `A0 + α·(A1 - A0)` (and likewise B), so the sample is a DP.
    """
    verts = _sic_vertices(sub, sis_pairs)
    closed = (sub.start == -1 and sub.end == -1)
    n = len(verts)
    seg_A: list = []
    seg_B: list = []
    node_xyz: list = []
    npair = n if closed else n - 1
    for k in range(npair):
        vk = verts[k]
        vk1 = verts[(k + 1) % n]
        s, b_k, b_k1 = _sic_owning_sis(vk, vk1, sis_pairs, dp_pair_to_sis, sp_to_sis)
        A0, B0 = _sic_preAB_at(s, b_k, sis_pairs, dps, domain)
        A1, B1 = _sic_preAB_at(s, b_k1, sis_pairs, dps, domain)
        if _needs_close(domain):
            # localize (NOT interpolate): a seam-spanning SIS segment must keep
            # its endpoint as the inside→outside σ-rep, so the lerp renders the
            # short on-surface arc (S(σ(A1)) = S(A1)); interpolate's map-back
            # would return the far A1 and re-create a disk-spanning spike.
            A1 = domain.localize(A0, A1)
            B1 = domain.localize(B0, B1)
        seg_A.append((A0, A1))
        seg_B.append((B0, B1))
        node_xyz.append(_sic_vert_xyz(vk, splits, dps))
    node_xyz.append(_sic_vert_xyz(verts[0] if closed else verts[-1], splits, dps))
    return np.asarray(node_xyz, dtype=float), seg_A, seg_B


def _sic_node_polyline(sub, sis_pairs, dps, splits, projection):
    """Node-skeleton (uv, xyz, xy) for an SIC SubCurve — used by `_build_polyline`
    for the arclength of `L_per_sub` and the SP-less verbatim pass-through.

    `uv` is the canonical (sheet-1 / SP-uv) node coords; it is metadata only —
    SIC visibility/rendering consume `xy`/`depth`, never `uv`.
    """
    verts = _sic_vertices(sub, sis_pairs)
    if sub.start == -1 and sub.end == -1 and verts:
        verts = verts + [verts[0]]   # close the loop
    uvs: list = []
    xyzs: list = []
    for v in verts:
        if v[0] == "DP":
            uvs.append(np.asarray(dps[int(v[1])]["uv1"], dtype=float))
            xyzs.append(np.asarray(dps[int(v[1])]["xyz"], dtype=float))
        else:
            sp = splits.sps[int(v[1])]
            uvs.append(np.asarray(sp[0], dtype=float))
            xyzs.append(np.asarray(sp[1], dtype=float))
    uv_arr = np.asarray(uvs, dtype=float) if uvs else np.zeros((0, 2), dtype=float)
    xyz_arr = np.asarray(xyzs, dtype=float) if xyzs else np.zeros((0, 3), dtype=float)
    xy_arr = projection.XY(xyz_arr) if len(xyz_arr) else np.zeros((0, 2), dtype=float)
    return uv_arr, xyz_arr, xy_arr


def _invert_distance_to_arclength(
    poly, which_end, targets, vp_xy, surface, projection, domain,
):
    """Arclength positions on a CC half-curve at given Euclidean distances to a VP.

    Walks the branch's CP polyline from its VP end; for each target distance
    `d` finds the point at Euclidean (image) distance `d` from `vp_xy` by
    bisecting on the uv-interpolation parameter and REPROJECTING through the
    surface (`S(uv) → XY`) — never a straight image-space lerp. Returns the
    matching arclengths (in the polyline's start→end `cum` frame) so the caller
    can feed them straight into the CC resampler's `s_targets` ladder.

    Returns None if a target can't be bracketed monotonically from the VP
    (the branch curves back within the master's reach) → caller falls back to
    native sampling for that VP.
    """
    uv_p, xyz_p, xy_p = poly
    N = len(xy_p)
    if N < 2:
        return None
    tgt = np.asarray(targets, dtype=float)
    if not len(tgt):
        return []
    cum = _arclengths(xy_p)
    dist_v = np.linalg.norm(xy_p - vp_xy, axis=1)
    order = np.arange(N) if which_end == "start" else np.arange(N - 1, -1, -1)
    d_ord = dist_v[order]                          # distances from the VP, in VP order

    # First bracketing segment (in VP order) per target — the segment whose
    # endpoint distances straddle the target. Vectorized over targets; no
    # surface evals here. Bail to native sampling if any target is unbracketable
    # (matches the old per-target scan's `return None`).
    seg_lo = np.minimum(d_ord[:-1], d_ord[1:])
    seg_hi = np.maximum(d_ord[:-1], d_ord[1:])
    inb = (seg_lo[None, :] - 1e-12 <= tgt[:, None]) & (tgt[:, None] <= seg_hi[None, :] + 1e-12)
    if not inb.any(axis=1).all():
        return None
    oi = inb.argmax(axis=1)                         # first bracketing segment
    va = order[oi]
    vb = order[oi + 1]
    lo = np.minimum(va, vb)
    hi = np.maximum(va, vb)
    uva = uv_p[lo]
    uvb = uv_p[hi]

    def dist_at(alpha):                             # alpha (M,) → distances (M,)
        if domain is None:
            uv = uva + alpha[:, None] * (uvb - uva)
        else:
            uv = np.asarray(domain.interpolate(uva, uvb, alpha), dtype=float)
        S = np.asarray(surface.S(uv[:, 0], uv[:, 1]), dtype=float)
        xy = projection.XY(np.ascontiguousarray(S.T))
        return np.linalg.norm(xy - vp_xy, axis=1)

    # Vectorized bisection within each bracket (distance is monotone in alpha
    # over the small segment, either direction).
    a0 = np.zeros(len(tgt)); a1 = np.ones(len(tgt))
    inc = dist_at(a1) >= dist_at(a0)
    for _ in range(30):
        am = 0.5 * (a0 + a1)
        below = (dist_at(am) < tgt) == inc
        a0 = np.where(below, am, a0)
        a1 = np.where(below, a1, am)
    alpha = 0.5 * (a0 + a1)
    s = cum[lo] + alpha * (cum[hi] - cum[lo])
    return s.tolist()


def _cc_vp_match_equidist(subcurves, polys, L_per_sub, splits, surface,
                          projection, domain, vp_trim):
    """Equal-distance VP matching + cusp trim (spec 2026-06-16).

    For each VP shared by exactly two CC half-curves, the SHORTER one (by
    projected length) is the *master*; the longer branch is resampled so its
    points sit at the SAME Euclidean distances to the VP (in image space) as
    the master's, up to the master's far end (beyond that the longer branch
    keeps its native CP arclengths). The two near-cusp branches are then
    mirror-matched.

    The innermost `vp_trim` points are then dropped from BOTH branches: at a
    cusp the projected contour velocity → 0, so the near-tip points are
    hypersensitive and straddle the VP (a phantom inter-branch occlusion
    crossing — the Klein bottle bug). Dropping them removes the straddle while
    keeping the two branches matched (they drop the same distance band).

    Returns `{sub_index: s_targets ndarray}` overriding the CC `s_targets = cum`
    default. Each VP end is handled independently, so a branch that is master
    at one end and longer at the other is modified near both ends with its
    native vertices kept in between.
    """
    # VP-SP → CC half-curves touching it.
    vp_eps: dict[int, list[tuple[int, str]]] = {}
    for i, sub in enumerate(subcurves):
        if sub.kind != "CC":
            continue
        for end, sp in (("start", sub.start), ("end", sub.end)):
            if sp >= 0 and splits.sps[int(sp)][3] == "vp":
                vp_eps.setdefault(int(sp), []).append((i, end))

    # Per sub, accumulate modifications (role at each shared VP end).
    mods: dict[int, list[dict]] = {}
    for sp, eps in vp_eps.items():
        if len(eps) != 2:
            continue
        (ia, ea), (ib, eb) = eps
        if ia == ib:
            continue  # closed CC through a single VP — skip
        if L_per_sub[ia] <= L_per_sub[ib]:
            im, em, io, eo = ia, ea, ib, eb
        else:
            im, em, io, eo = ib, eb, ia, ea
        vp_xy = np.asarray(splits.sps[int(sp)][2], dtype=float)
        # Master target distances: its CP vertices from the VP outward.
        m_xy = polys[im][2]
        m_d = np.linalg.norm(m_xy - vp_xy, axis=1)
        if em == "end":
            m_d = m_d[::-1]
        targets = []
        for d in m_d[1:]:
            if d > 1e-12 and (not targets or d > targets[-1] + 1e-12):
                targets.append(float(d))
        if not targets:
            continue
        reach = targets[-1]
        s_list = _invert_distance_to_arclength(
            polys[io], eo, targets, vp_xy, surface, projection, domain)
        if s_list is None:
            continue
        mods.setdefault(im, []).append(dict(role="master", end=em, vp_xy=vp_xy, reach=reach))
        mods.setdefault(io, []).append(dict(role="longer", end=eo, vp_xy=vp_xy,
                                             reach=reach, s_list=s_list))

    out: dict[int, np.ndarray] = {}
    for i, items in mods.items():
        uv_p, xyz_p, xy_p = polys[i]
        N = len(xy_p)
        cum = _arclengths(xy_p)
        L = float(cum[-1]) if N else 0.0
        drop: set[int] = set()
        add_s: list[float] = []
        for m in items:
            if m["role"] == "master":
                # Drop the vp_trim native vertices nearest the VP (the VP
                # endpoint vertex itself is kept / pinned).
                order = range(1, N) if m["end"] == "start" else range(N - 2, -1, -1)
                for cnt, v in enumerate(order):
                    if cnt >= vp_trim:
                        break
                    drop.add(int(v))
            else:  # longer
                # Replace the CONTIGUOUS near-VP run (from the VP outward, in
                # arclength, until the distance first exceeds the master's
                # reach) with the matched points. Walking the run — rather than
                # dropping every vertex within Euclidean `reach` — is essential:
                # a long contour can loop back near the VP at a distant
                # arclength, and dropping those would leave a gap (a spurious
                # long skip segment). Drop the innermost vp_trim matched points.
                order = range(N) if m["end"] == "start" else range(N - 1, -1, -1)
                for v in order:
                    if float(np.linalg.norm(xy_p[v] - m["vp_xy"])) <= m["reach"]:
                        drop.add(int(v))
                    else:
                        break
                add_s.extend(m["s_list"][vp_trim:])
        drop.discard(0)
        drop.discard(N - 1)  # never drop the SP endpoints
        s_all = [0.0, L]
        s_all.extend(float(cum[v]) for v in range(N) if v not in drop)
        s_all.extend(add_s)
        out[i] = np.unique(np.round(np.clip(s_all, 0.0, L), 9))
    return out


def _cc_vp_match_asym(subcurves, polys, splits, vp_trim):
    """Asymmetric cusp trim at VPs (spec 2026-06-17, experimental).

    At a cusp (VP) the two CC half-curves are tangent and overlap in the image
    near the tip; sampled independently, their near-tip polylines straddle the
    VP and inject a phantom inter-branch occlusion break (the Klein-bottle bug).
    Like `_cc_vp_match_equidist` this kills the straddle, but by DROPPING
    near-cusp vertices only — no resampling, no reprojection:

      * `reference` branch — the one whose `vp_trim`-th vertex from the VP is
        FARTHER from the VP in the image — drops exactly `vp_trim` vertices
        nearest the VP. Call the image distance of its first surviving vertex
        `d_ref`.
      * `matched` branch drops as many near-VP vertices as makes its first
        surviving vertex's image distance to the VP closest to `d_ref`, so both
        branches leave the cusp at the same radius (symmetric — no straddle).
        It always drops at least `vp_trim` (choosing the farther branch as the
        reference guarantees the matched branch only ever drops MORE, never
        fewer, to reach `d_ref`).

    The VP endpoint vertex (index 0 from the VP) and the far SP endpoint are
    always kept. Returns `{sub_index: s_targets ndarray}` overriding the CC
    `s_targets = cum` default. Each VP end is handled independently, so a branch
    sharing a VP at each end is trimmed near both, keeping its mid vertices.
    """
    # VP-SP → CC half-curves touching it.
    vp_eps: dict[int, list[tuple[int, str]]] = {}
    for i, sub in enumerate(subcurves):
        if sub.kind != "CC":
            continue
        for end, sp in (("start", sub.start), ("end", sub.end)):
            if sp >= 0 and splits.sps[int(sp)][3] == "vp":
                vp_eps.setdefault(int(sp), []).append((i, end))

    def _vp_order(i, end):
        """Original vertex indices ordered from the VP outward, + their
        image distances to the VP (index 0 = the VP endpoint)."""
        xy = polys[i][2]
        N = len(xy)
        order = list(range(N)) if end == "start" else list(range(N - 1, -1, -1))
        vp_xy = xy[order[0]]
        d = np.linalg.norm(xy[order] - vp_xy, axis=1)
        return order, d

    # Per sub, accumulate the near-VP vertex drops from each shared VP end.
    mods: dict[int, set[int]] = {}
    for sp, eps in vp_eps.items():
        if len(eps) != 2:
            continue
        (ia, ea), (ib, eb) = eps
        if ia == ib:
            continue  # closed CC through a single VP — skip
        oa, da = _vp_order(ia, ea)
        ob, db = _vp_order(ib, eb)
        if min(len(da), len(db)) < 2:
            continue
        # First-surviving-vertex distance when each branch trims exactly
        # vp_trim (clamped to keep the far endpoint).
        ka = min(vp_trim, len(da) - 1)
        kb = min(vp_trim, len(db) - 1)
        # Reference = farther first-survivor (fixed trim); matched grows toward it.
        if da[ka] >= db[kb]:
            ref_i, ref_o, ref_k, d_ref = ia, oa, ka, float(da[ka])
            mat_i, mat_o, d_m = ib, ob, db
        else:
            ref_i, ref_o, ref_k, d_ref = ib, ob, kb, float(db[kb])
            mat_i, mat_o, d_m = ia, oa, da
        # Matched branch: grow its trim from vp_trim until its first survivor's
        # distance is closest to d_ref, stopping once it passes d_ref along the
        # increasing near-cusp run (so a far loop-back can't masquerade as the
        # match).
        k0 = min(vp_trim, len(d_m) - 1)
        mat_k, best_err = k0, abs(float(d_m[k0]) - d_ref)
        for k in range(k0 + 1, len(d_m)):
            err = abs(float(d_m[k]) - d_ref)
            if err < best_err:
                mat_k, best_err = k, err
            if d_m[k] >= d_ref:
                break
        # Drop the innermost trim vertices on each branch (keep the VP endpoint).
        mods.setdefault(ref_i, set()).update(int(v) for v in ref_o[1:ref_k + 1])
        mods.setdefault(mat_i, set()).update(int(v) for v in mat_o[1:mat_k + 1])

    out: dict[int, np.ndarray] = {}
    for i, drop in mods.items():
        xy_p = polys[i][2]
        N = len(xy_p)
        cum = _arclengths(xy_p)
        L = float(cum[-1]) if N else 0.0
        drop = set(drop)
        drop.discard(0)
        drop.discard(N - 1)  # never drop the SP endpoints
        s_all = [0.0, L]
        s_all.extend(float(cum[v]) for v in range(N) if v not in drop)
        out[i] = np.unique(np.round(np.clip(s_all, 0.0, L), 9))
    return out


def _cc_vp_match_targets(subcurves, polys, L_per_sub, splits, surface,
                         projection, domain, vp_trim):
    """Dispatch VP cusp-matching by `settings.VP_MATCH_MODE`.

    "match" (default, committed) → `_cc_vp_match_equidist` (equal-distance
    resampling of both branches). "trim" (experimental) → `_cc_vp_match_asym`
    (asymmetric near-cusp trim, no resampling). Selectable per-request from the
    debug panel ("VP match mode") for A/B comparison on the fixture zoo.
    """
    from surface_play import settings as _settings
    mode = getattr(_settings, "VP_MATCH_MODE", "match")
    if mode == "trim":
        return _cc_vp_match_asym(subcurves, polys, splits, vp_trim)
    if mode != "match":
        raise ValueError(
            f"VP_MATCH_MODE must be 'match' or 'trim', got {mode!r}")
    return _cc_vp_match_equidist(
        subcurves, polys, L_per_sub, splits, surface, projection, domain, vp_trim)


def resample_all(
    subcurves: "list[SubCurve]",
    surface: "SurfaceParams",
    projection: "Projection",
    splits: "SplitArrays",
    mesh: "Mesh",
    css: np.ndarray,
    sis_pairs: np.ndarray,
    cps: np.ndarray,
    dps: np.ndarray,
    *,
    resolution: int | None = None,
    project_resampled: bool | None = None,
) -> list[ResampledCurve]:
    """Resample each SubCurve in projected space.

    Spec: §"Curve resampling" (lines 452-465).
    `resolution` / `project_resampled` default to `surface_play.settings`.
    """
    from surface_play import settings as _settings
    if resolution is None:
        resolution = _settings.RESOLUTION
    if project_resampled is None:
        project_resampled = _settings.PROJECT_RESAMPLED
    domain = getattr(mesh, "domain", None)

    # Step A — M (mesh-xy bbox diagonal).
    mesh_xy = projection.XY(mesh.xyz)
    bbox = mesh_xy.max(axis=0) - mesh_xy.min(axis=0) if len(mesh_xy) else np.array([1.0, 1.0])
    M = float(np.hypot(bbox[0], bbox[1])) or 1.0

    # Step B — polylines and per-SC arclength.
    polys = []
    L_per_sub = []
    for sub in subcurves:
        uv_p, xyz_p, xy_p = _build_polyline(
            sub, splits, mesh, css, sis_pairs, cps, dps, surface, projection,
        )
        polys.append((uv_p, xyz_p, xy_p))
        L_per_sub.append(float(_arclengths(xy_p)[-1]) if len(xy_p) > 1 else 0.0)

    # Step C — resampling scales (spec change 2026-05-27, per-SP variant):
    #   ell        = M / resolution                  (global coarse spacing)
    #   delta[sp]  = (min incident SubCurve length at sp) / 10
    # `_sample_arclengths` takes (L_total, ell, delta_start, delta_end,
    # is_closed) and produces 5 dense segments at each endpoint, ell-spaced
    # interior. Local-per-SP delta avoids a globally tiny sub from
    # contaminating distant subs' near-SP resolution.
    ell = M / float(resolution)
    # Shared BC/HC densification oversampling factor (dense spacing = ell/this).
    densify_subdiv = float(_settings.DENSIFY_SUBDIV)
    # Degenerate-SubCurve handling (2026-06-17). A SubCurve whose xy polyline is
    # shorter than `1e-4 * ell` is collapsed to essentially a single image point
    # — e.g. a tiny CC loop at a bump tip under some jitter+view, or a
    # non-generic axis-aligned view (the helicoid CC at u=0 under the Z-axis
    # view collapses to (0, 0)). History: originally a hard raise (aborted the
    # whole outline → intermittent HTTP 400 on bumpy surfaces), then a verbatim
    # pass-through. Now such subs are DELETED from the pipeline entirely (no RC
    # emitted): they carry no visible line, but their near-coincident points
    # destabilise the client cubic-Bézier fit (`_genBezier`'s least-squares
    # det → 0 → blown-up control arms → image spikes). Dropping them also keeps
    # their tiny L out of the per-SP delta so `_sample_arclengths` can't spend
    # 10⁷+ iterations climbing back to `ell` for a neighbour sharing the SP.
    # SP-less tiny closed loops are included (they spike the same way).
    _L_FLOOR = 1e-4 * ell
    degenerate = {
        i for i in range(len(subcurves)) if 0 < L_per_sub[i] < _L_FLOOR
    }
    L_per_sp: dict[int, float] = {}
    for i, sub in enumerate(subcurves):
        L = L_per_sub[i]
        if L <= 0 or i in degenerate:
            continue
        for sp in (sub.start, sub.end):
            if sp >= 0:
                L_per_sp[sp] = min(L_per_sp.get(sp, float("inf")), L)
    delta_per_sp = {sp: L / 10.0 for sp, L in L_per_sp.items()}

    # VP cusp matching: for each VP shared by two CC half-curves, override the
    # CC arclength ladder near the cusp so the two overlapping near-cusp
    # polylines can't spuriously cross (which injects a phantom occlusion
    # break). Algorithm selected by settings.VP_MATCH_MODE ("match" equal-
    # distance resample, default / "trim" asymmetric trim). See
    # _cc_vp_match_targets dispatcher.
    cc_match = _cc_vp_match_targets(
        subcurves, polys, L_per_sub, splits, surface, projection, domain,
        int(_settings.VP_TRIM))

    # SIC interpolation maps (built once): DP-pair → SIS, and SP → [(SIS, bary)].
    # Used to resolve each SIC polyline segment's owning SIS (and SP barys)
    # structurally — no proximity heuristic.
    dp_pair_to_sis: dict[tuple[int, int], int] = {}
    sp_to_sis: dict[int, list[tuple[int, float]]] = {}
    if len(sis_pairs):
        for _s in range(len(sis_pairs)):
            _row = sis_pairs[_s]
            _p = int(_row["p_dp"]); _q = int(_row["q_dp"])
            dp_pair_to_sis[(min(_p, _q), max(_p, _q))] = _s
            for _slot_name in ("split1", "split2"):
                _slot = int(_row[_slot_name])
                if _slot >= 0:
                    _spt = splits.spts[_slot]
                    sp_to_sis.setdefault(int(_spt[0]), []).append((_s, float(_spt[1])))

    # Step E — sample each SC.
    out: list[ResampledCurve] = []
    for i, sub in enumerate(subcurves):
        uv_p, xyz_p, xy_p = polys[i]
        L_total = L_per_sub[i]

        # Degenerate SubCurve (collapsed to ~a point) → DELETE from the pipeline:
        # emit no RC at all. It carries no visible line, and its near-coincident
        # points blow up the client Bézier fit into image spikes. Checked before
        # the SP-less branch so tiny closed loops are dropped too. Its SPs, if
        # shared, stay anchored by the real curves that use them.
        if i in degenerate:
            continue

        # SP-less closed SC → verbatim copy.
        if sub.start == -1 and sub.end == -1:
            depth = (np.asarray(projection.Z(xyz_p), dtype=float)
                     if len(xyz_p) else np.zeros(0, dtype=float))
            out.append(ResampledCurve(
                kind=sub.kind, start=-1, end=-1,
                depth=depth, xy=xy_p.copy(), dir=None,
                vc_in=int(sub.vc_in), vc_out=int(sub.vc_out),
                uv=uv_p.copy() if len(uv_p) else uv_p,
            ))
            continue

        # SIC → arclength resample on the DP/SP xyz skeleton; each sample is a
        # DP whose preimage-A is interpolated per-SIS (flip + close-aware). No
        # sheet picked by proximity; `dir/tan = None`.
        if sub.kind == "SIC":
            node_xyz, seg_A, _seg_B = _sic_build(
                sub, sis_pairs, dps, splits, dp_pair_to_sis, sp_to_sis, domain)
            node_xy = (projection.XY(node_xyz) if len(node_xyz)
                       else np.zeros((0, 2), dtype=float))
            cum = _arclengths(node_xy)
            L_xy = float(cum[-1]) if len(cum) else 0.0
            delta_s = delta_per_sp.get(sub.start, ell)
            delta_e = delta_per_sp.get(sub.end, ell)
            s_targets = (_sample_arclengths(L_xy, ell, delta_s, delta_e, sub.is_closed)
                         if L_xy > 0 and len(seg_A) else np.array([0.0, L_xy]))
            N = len(s_targets)
            sample_uv = np.zeros((N, 2), dtype=float)
            for j, s in enumerate(s_targets):
                st = float(np.clip(s, 0.0, cum[-1])) if len(cum) else 0.0
                seg = int(np.searchsorted(cum, st, side="right") - 1)
                seg = max(0, min(seg, len(seg_A) - 1))
                span = cum[seg + 1] - cum[seg]
                alpha = (st - cum[seg]) / span if span > 0 else 0.0
                A0, A1 = seg_A[seg]
                sample_uv[j] = A0 + alpha * (A1 - A0)
            if N > 0:
                S_all = np.asarray(surface.S(sample_uv[:, 0], sample_uv[:, 1]),
                                   dtype=float)
                sample_xyz = np.ascontiguousarray(S_all.T)
                sample_xy = projection.XY(sample_xyz)
            else:
                sample_xyz = np.zeros((0, 3), dtype=float)
                sample_xy = np.zeros((0, 2), dtype=float)
            # Pin endpoints to their SP positions (image-exact).
            if N >= 1 and sub.start >= 0:
                sp0 = splits.sps[sub.start]
                sample_xyz[0] = np.asarray(sp0[1], dtype=float)
                sample_xy[0] = np.asarray(sp0[2], dtype=float)
            if N >= 1 and sub.end >= 0:
                sp1 = splits.sps[sub.end]
                sample_xyz[-1] = np.asarray(sp1[1], dtype=float)
                sample_xy[-1] = np.asarray(sp1[2], dtype=float)
            depth = (np.asarray(projection.Z(sample_xyz), dtype=float)
                     if N else np.zeros(0, dtype=float))
            out.append(ResampledCurve(
                kind="SIC", start=sub.start, end=sub.end,
                depth=depth, xy=sample_xy, dir=None, tan=None,
                vc_in=int(sub.vc_in), vc_out=int(sub.vc_out),
                uv=sample_uv.copy(),
            ))
            continue

        # HC → resample the uv straight line using the same arclength
        # logic as BC/CC, so the HC has comparable per-segment xy density.
        # Step 1: pre-sample uniformly in uv (dense) → curved xy polyline +
        # cumulative xy arclength. Step 2: use `_sample_arclengths` to pick
        # target arclengths just like BC. Step 3: interpolate back to uv.
        if sub.kind == "HC":
            uv_q0 = uv_p[0]; uv_q1 = uv_p[1]
            # Unified densification (spec 2026-06-03): N points =
            # resolution·DENSIFY_SUBDIV·L/M — samples per unit image-arclength,
            # unit = M (mesh-xy bbox diagonal), so spacing ≈ ell/DENSIFY_SUBDIV
            # and BC/HC share one density knob. The HC is a straight uv line but
            # curved in the image, so bootstrap L from a coarse pass before N.
            t_boot = np.linspace(0.0, 1.0, 33)
            uv_boot = uv_q0[None, :] + t_boot[:, None] * (uv_q1 - uv_q0)[None, :]
            xy_boot = projection.XY(np.ascontiguousarray(
                np.asarray(surface.S(uv_boot[:, 0], uv_boot[:, 1]), dtype=float).T))
            L_boot = float(_arclengths(xy_boot)[-1])
            N_dense = max(2, int(round(resolution * densify_subdiv * L_boot / M)))
            t_dense = np.linspace(0.0, 1.0, N_dense)
            uv_dense = uv_q0[None, :] + t_dense[:, None] * (uv_q1 - uv_q0)[None, :]
            # Batched dense xyz/xy: one `surface.S` call for all N_dense uv pairs.
            S_dense = np.asarray(
                surface.S(uv_dense[:, 0], uv_dense[:, 1]), dtype=float,
            )  # (3, N_dense)
            xyz_dense = np.ascontiguousarray(S_dense.T)
            xy_dense = projection.XY(xyz_dense)
            cum_dense = _arclengths(xy_dense)
            L_xy_total = float(cum_dense[-1])

            delta_s = delta_per_sp.get(sub.start, ell)
            delta_e = delta_per_sp.get(sub.end, ell)
            s_targets = _hc_s_targets(
                sub, i, subcurves, polys,
                surface, projection,
                L_xy_total, ell, delta_s, delta_e,
            )

            sample_uv = np.empty((len(s_targets), 2), dtype=float)
            sample_xyz = np.empty((len(s_targets), 3), dtype=float)
            sample_xy = np.empty((len(s_targets), 2), dtype=float)
            for j, s in enumerate(s_targets):
                if s <= 0.0:
                    t = 0.0
                elif s >= L_xy_total:
                    t = 1.0
                else:
                    idx = int(np.searchsorted(cum_dense, s) - 1)
                    idx = max(0, min(idx, N_dense - 2))
                    denom = cum_dense[idx + 1] - cum_dense[idx]
                    frac = (s - cum_dense[idx]) / denom if denom > 0 else 0.0
                    t = t_dense[idx] + frac * (t_dense[idx + 1] - t_dense[idx])
                uv_s = uv_q0 + t * (uv_q1 - uv_q0)
                sample_uv[j] = uv_s
            # Batched sample xyz/xy for all selected samples.
            sample_tan = None
            if len(sample_uv):
                S_samples = np.asarray(
                    surface.S(sample_uv[:, 0], sample_uv[:, 1]), dtype=float,
                )
                sample_xyz[:] = S_samples.T
                sample_xy[:] = projection.XY(sample_xyz)
                # Analytic image-space tangent: lift the constant uv-line
                # direction via S_u·Δu + S_v·Δv at each sample, then project.
                # Replaces the chord fallback in compute_projection_breaks
                # (xy[k+1]-xy[k]), which shares one coarse tangent across all
                # crossings inside the same segment and breaks dv-sign
                # alternation when multiple breaks land in one segment.
                duv = uv_q1 - uv_q0
                Su_samples = np.asarray(
                    surface.Su(sample_uv[:, 0], sample_uv[:, 1]),
                    dtype=float,
                ).T  # (N, 3)
                Sv_samples = np.asarray(
                    surface.Sv(sample_uv[:, 0], sample_uv[:, 1]),
                    dtype=float,
                ).T  # (N, 3)
                T_3d = duv[0] * Su_samples + duv[1] * Sv_samples
                sample_tan = np.empty((len(sample_uv), 2), dtype=float)
                for j in range(len(sample_uv)):
                    sample_tan[j] = projection.proj_vec(sample_uv[j], T_3d[j])
            depth = (np.asarray(projection.Z(sample_xyz), dtype=float)
                     if len(sample_xyz) else np.zeros(0, dtype=float))
            out.append(ResampledCurve(
                kind="HC", start=sub.start, end=sub.end,
                depth=depth, xy=sample_xy, dir=None, tan=sample_tan,
                vc_in=int(sub.vc_in), vc_out=int(sub.vc_out),
                uv=sample_uv.copy(),
            ))
            continue

        # Pick sample arclengths (two-phase tapered spacing).
        delta_s = delta_per_sp.get(sub.start, ell)
        delta_e = delta_per_sp.get(sub.end, ell)
        cum = _arclengths(xy_p)
        # Per-kind interpolation polyline. BC densifies for arclength accuracy
        # near projection folds (see _densify_bc_polyline); CC/SIC interpolate
        # against their own (already-fine) build-polyline. `interp_nsub` maps a
        # dense segment index back to the original segment via `// interp_nsub`,
        # preserving the seg→mesh-edge mapping the BC `dir` lookup needs.
        interp_uv, interp_xy, interp_cum = uv_p, xy_p, cum
        interp_nsub = 1
        if sub.kind == "CC":
            # CC samples = polyline vertices verbatim. The CPs are already on
            # the true contour (find_contour_points uses Newton). The CC
            # defines the authoritative arclength ladder at each shared SP;
            # BCs and HCs inherit it via tangent-pick at BCPs / HAs.
            # See [[resume-hc-match-cc]].
            # Exception: a longer branch at a shared VP uses the equal-distance
            # matched ladder (Phase-1 interpolates uv at these arclengths and
            # Phase-2 reprojects, so the matched points lie on the contour).
            s_targets = cc_match[i] if i in cc_match else cum.copy()
        elif sub.kind == "BC":
            # Unified densification (spec 2026-06-03): total dense points ≈
            # resolution·DENSIFY_SUBDIV·L/M (samples per unit image-arclength,
            # unit = M), matching HC. `_densify_bc_polyline` subdivides each of
            # the n_seg mesh-edge segments uniformly, so map the target total to
            # a per-segment subdivision count (≥ 1).
            n_seg = max(1, len(uv_p) - 1)
            interp_nsub = max(1, int(round(resolution * densify_subdiv * L_total / M / n_seg)))
            interp_uv, interp_xy, interp_cum = _densify_bc_polyline(
                uv_p, surface, projection, domain, interp_nsub)
            # Native-vertex arclengths in the ACCURATE metric (dense vertex k
            # lives at index k*interp_nsub) — so own CPs and inherited CC arcs
            # share one consistent arclength scale.
            native_arc = interp_cum[np.arange(len(uv_p)) * interp_nsub]
            # BC samples = own CPs + inherited CC arclengths near BCPs.
            # See `_bc_s_targets` for the tangent-pick rule.
            s_targets = _bc_s_targets(sub, i, subcurves, polys, splits,
                                       float(interp_cum[-1]), native_arc)
        else:
            s_targets = _sample_arclengths(L_total, ell, delta_s, delta_e, sub.is_closed)
        N = len(s_targets)
        sample_uv = np.zeros((N, 2), dtype=float)
        sample_xyz = np.zeros((N, 3), dtype=float)
        sample_xy = np.zeros((N, 2), dtype=float)
        sample_dir = np.zeros((N, 2), dtype=float) if sub.kind in ("BC", "CC") else None
        sample_tan = np.zeros((N, 2), dtype=float) if sub.kind in ("BC", "CC") else None
        # Track which internal seg we hit (for BC dir, annular snap, CC dir/newton).
        internal_seg_of_poly_seg = []
        # Polyline index → (sub.internal entry idx) for non-endpoint segments.
        # poly vertex 0 = start SP; poly vertex k for k in 1..len(internal) = internal[k-1];
        # poly vertex len(internal)+1 = end SP. The polyline SEGMENT k (between vertex k and k+1)
        # corresponds to sub.internal[k] for k in 0..len(internal)-1, and the last segment
        # corresponds to sub.internal[-1] (same).
        # Simpler: per poly-segment k, the "owning" sub.internal entry that holds its END is
        # sub.internal[k] (if k < len(sub.internal)); the last poly-segment ends at end SP,
        # which sits on sub.internal[-1]'s segment (if any) or directly the start segment.

        # Phase 1 (vectorized) — resolve all sample uv in one batched polyline
        # walk. The former per-sample `_interp_along_polyline` (clip /
        # searchsorted / close-aware lerp) is array-friendly, and
        # `domain.interpolate` already accepts (N,2)/(N,) inputs, so the whole
        # `s_targets` ladder resolves in a handful of numpy ops instead of N
        # Python-level calls. No surface evals here so Phase 2 can do a single
        # batched _eval_all.
        chain_segs = np.full(N, -1, dtype=np.int64)
        if N > 0 and len(interp_cum) >= 2:
            st = np.clip(np.asarray(s_targets, dtype=float), 0.0, interp_cum[-1])
            dseg = np.searchsorted(interp_cum, st, side="right") - 1
            dseg = np.clip(dseg, 0, len(interp_cum) - 2)
            span = interp_cum[dseg + 1] - interp_cum[dseg]
            alphas = np.where(span > 0, (st - interp_cum[dseg]) / span, 0.0)
            uv_a = interp_uv[dseg]
            uv_b = interp_uv[dseg + 1]
            if domain is None:
                sample_uv = uv_a + alphas[:, None] * (uv_b - uv_a)
            else:
                sample_uv = np.asarray(
                    domain.interpolate(uv_a, uv_b, alphas), dtype=float)
            seg_ps = (dseg // interp_nsub).astype(np.int64)
        else:
            alphas = np.zeros(N, dtype=float)
            seg_ps = np.zeros(N, dtype=np.int64)
        # Per-sample geometric corrections (rare paths; absent for plain rect /
        # non-projected resampling — kept looped to preserve exact behaviour).
        if (sub.kind == "BC" and domain is not None
                and getattr(domain, "type", None) in ("disk", "annulus")):
            for j in range(N):
                sample_uv[j] = _snap_annular_bc(sample_uv[j], mesh)
        if sub.kind == "CC" and project_resampled:
            for j in range(N):
                sample_uv[j] = _newton_cc_refine(sample_uv[j], surface, projection)
        # Map each sample's original build-polyline segment to its mesh/CS index.
        if sample_dir is not None and sub.internal:
            internal0 = np.array([int(e[0]) for e in sub.internal], dtype=np.int64)
            chain_segs = internal0[np.clip(seg_ps, 0, len(sub.internal) - 1)]

        # Phase 2 — single batched `_eval_all` over all samples. The cse'd
        # lambdified callable computes S, Su, Sv, Suu, Suv, Svv, SN for
        # every uv in one numpy pass — ~5 ms total vs ~150 µs × N
        # scalar calls. For BC samples, uv must first be brought into
        # canonical p's identification copy so Su, Sv match edge["dir"]
        # / edge["pq"]'s frame at the mo seam.
        uv_eval = sample_uv.copy()
        if sub.kind == "BC":
            for j in range(N):
                cs = int(chain_segs[j])
                if cs >= 0:
                    edge = mesh.edges[cs]
                    uv_eval[j] = _uv_for_bc_lift(edge, mesh, sample_uv[j])
        if N > 0:
            S_all, Su_all, Sv_all, Suu_all, Suv_all, Svv_all, _SN_all = \
                surface._eval_all(uv_eval[:, 0], uv_eval[:, 1])
            # Convert (3, N) → (N, 3) once.
            S_all = np.ascontiguousarray(np.asarray(S_all, dtype=float).T)
            Su_all = np.ascontiguousarray(np.asarray(Su_all, dtype=float).T)
            Sv_all = np.ascontiguousarray(np.asarray(Sv_all, dtype=float).T)
            Suu_all = np.ascontiguousarray(np.asarray(Suu_all, dtype=float).T)
            Suv_all = np.ascontiguousarray(np.asarray(Suv_all, dtype=float).T)
            Svv_all = np.ascontiguousarray(np.asarray(Svv_all, dtype=float).T)
            sample_xyz[:] = S_all
            sample_xy[:] = projection.XY(S_all)

        # Phase 3 — dir/tan. CC tangents AND CC dirs are batched; BC still
        # loops per-sample (per-edge lookups + domain.interpolate / boundary
        # tangent that aren't vectorized here).
        if sample_dir is not None and sub.kind == "CC" and N > 0:
            cc_mask = chain_segs >= 0
            if cc_mask.any():
                idx = np.flatnonzero(cc_mask)
                sample_tan[idx] = _tan_for_cc_samples_batched(
                    sample_uv[idx], chain_segs[idx],
                    S_all[idx], Su_all[idx], Sv_all[idx],
                    Suu_all[idx], Suv_all[idx], Svv_all[idx],
                    css, cps, projection,
                )
                # CC dir = CP d-field linearly interpolated along the cs:
                # (1-α)·d_p + α·d_q, where d_p/d_q are the CPs' d vectors and α
                # is the Phase-1 along-segment fraction. Batched (same float ops
                # per element as the former per-sample loop).
                cs_rows = css[chain_segs[idx]]
                d_p = np.asarray(cps[cs_rows["p_cp"]]["d"], dtype=float)
                d_q = np.asarray(cps[cs_rows["q_cp"]]["d"], dtype=float)
                a = np.asarray(alphas[idx], dtype=float).reshape(-1, 1)
                sample_dir[idx] = (1.0 - a) * d_p + a * d_q
        if sample_dir is not None and sub.kind == "BC" and N > 0:
            # BC dir = projected inward 2D normal: edge["dir"] lifted via Su, Sv
            #          then proj_vec.
            # BC tan = projected analytic boundary tangent: domain.boundary_tangent
            #          (edge["pq"]-signed) lifted via Su, Sv then proj_vec. The
            #          analytic tangent replaces the sample chord, which is
            #          jitter-noisy at fine res and can flip the projection-break
            #          discriminator's sign.
            # Both edge["dir"] / edge["pq"] live in the canonical p's
            # identification copy; `uv_eval` (Phase 2) already holds the
            # seam-canonical lift uv (_uv_for_bc_lift), so Su/Sv/S there match
            # that frame (e.g. Möbius Su(v+2π) = -Su(v) would otherwise corrupt
            # the lift). See [[bc_lift_patch_match_2026_05_26]].
            bc_mask = chain_segs >= 0
            if bc_mask.any():
                idx = np.flatnonzero(bc_mask)
                edge_rows = mesh.edges[chain_segs[idx]]
                dir_uv = np.asarray(edge_rows["dir"], dtype=float).reshape(-1, 2)
                edge_pq = np.asarray(edge_rows["pq"], dtype=float).reshape(-1, 2)
                Su_i, Sv_i = Su_all[idx], Sv_all[idx]
                uv_lift = uv_eval[idx]
                inward_3d = dir_uv[:, 0:1] * Su_i + dir_uv[:, 1:2] * Sv_i
                Tb_uv = _boundary_tangent_batched(domain, uv_lift, edge_pq)
                Tb_3d = Tb_uv[:, 0:1] * Su_i + Tb_uv[:, 1:2] * Sv_i
                sample_dir[idx], sample_tan[idx] = _proj_vec2_batched(
                    projection, S_all[idx], inward_3d, Tb_3d)

        # Endpoint pinning (G5 / G21) — first/last sample anchored to SP positions.
        if N >= 1 and sub.start >= 0:
            sp_start = splits.sps[sub.start]
            sample_xyz[0] = np.asarray(sp_start[1], dtype=float)
            sample_xy[0] = np.asarray(sp_start[2], dtype=float)
        if N >= 1 and sub.end >= 0:
            # Always pin the end to its SP, even for closed subs (start == end):
            # otherwise the last sample drifts from the SP by up to one
            # spacing-step, leaving a visible gap that looks like a spurious
            # segment in tilted views. `sub.end >= 0` already excludes SP-less
            # closed subs (which use start = end = -1).
            sp_end = splits.sps[sub.end]
            sample_xyz[-1] = np.asarray(sp_end[1], dtype=float)
            sample_xy[-1] = np.asarray(sp_end[2], dtype=float)

        # Sign-align analytic tangents with local chain-forward direction in
        # image (chord between neighboring samples). This is robust because
        # only the SIGN matters: the chord direction is dominated by the
        # macroscopic chain direction. Without this, the analytic tangent
        # may point against chain direction when the rc traverses a CS in
        # reverse of its native (p_cp → q_cp) order — that flipped half the
        # CC break signs on torus trial 7.
        if sample_tan is not None and N >= 2:
            for j in range(N):
                if j == 0:
                    chord = sample_xy[1] - sample_xy[0]
                elif j == N - 1:
                    chord = sample_xy[N - 1] - sample_xy[N - 2]
                else:
                    chord = sample_xy[j + 1] - sample_xy[j - 1]
                if float(sample_tan[j] @ chord) < 0.0:
                    sample_tan[j] = -sample_tan[j]

        depth = (np.asarray(projection.Z(sample_xyz), dtype=float)
                 if N else np.zeros(0, dtype=float))
        out.append(ResampledCurve(
            kind=sub.kind, start=sub.start, end=sub.end,
            depth=depth, xy=sample_xy, dir=sample_dir, tan=sample_tan,
            vc_in=int(sub.vc_in), vc_out=int(sub.vc_out),
            uv=sample_uv.copy(),
        ))

    return out
