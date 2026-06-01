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

    uvs: list[np.ndarray] = []

    if sub.kind == "HC":
        uvs.append(np.asarray(splits.sps[sub.start][0], dtype=float))
        uvs.append(np.asarray(splits.sps[sub.end][0], dtype=float))
    elif sub.start == -1 and sub.end == -1:
        # SP-less closed SC — walk internal verbatim; close the loop.
        for seg_idx, bary in sub.internal:
            uvs.append(_seg_uv_at_bary(sub.kind, int(seg_idx), float(bary),
                                       mesh, css, sis_pairs, cps, dps))
        if uvs:
            uvs.append(uvs[0].copy())
    else:
        # Normal BC/CC SC: start SP → internal → end SP.
        uvs.append(np.asarray(splits.sps[sub.start][0], dtype=float))
        for seg_idx, bary in sub.internal:
            uvs.append(_seg_uv_at_bary(sub.kind, int(seg_idx), float(bary),
                                       mesh, css, sis_pairs, cps, dps))
        uvs.append(np.asarray(splits.sps[sub.end][0], dtype=float))

    # Close-aware adjust consecutive vertices, then lift to xyz/xy.
    # Only for GLOBAL periodicities (rect cy/mo), where S is periodic so the
    # aliased vertex is the same surface point. NOT for antipodal: there close()
    # is a boundary-only identification (S(x) ≠ S(-x) in the interior), so
    # reflecting an interior arc would relocate it onto different surface points
    # (e.g. a CC tail near a cusp jumping to the antipodal preimage). A
    # seam-crossing polyline keeps its true uv; its two boundary endpoints P and
    # -P already lift to one xy via S(P) = S(-P), so no long chord appears.
    domain = getattr(mesh, "domain", None)
    if _needs_close(domain) and not getattr(domain, "is_antipodal", False):
        for i in range(1, len(uvs)):
            uvs[i] = domain.interpolate(uvs[i - 1], uvs[i], 1.0)

    uv_arr = np.asarray(uvs, dtype=float) if uvs else np.zeros((0, 2), dtype=float)
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

    Extracted so `_dir_for_bc_sample` and `_tan_for_bc_sample` share the
    same seam-canonical uv (cf. [[bc_lift_patch_match_2026_05_26]]).
    """
    domain = getattr(mesh, "domain", None)
    if _needs_close(domain):
        p_canonical = mesh.uv[int(edge["p_idx"])]
        return domain.interpolate(p_canonical, uv_sample, 1.0)
    return uv_sample


def _dir_for_bc_sample(seg_idx: int, mesh, surface, projection,
                       uv_sample: np.ndarray, *,
                       Su=None, Sv=None) -> np.ndarray:
    """Projected inward 2D normal of a BC sample (from edge['dir'] at uv_sample).

    `edge["dir"]` is expressed in the canonical p's identification copy (fixed
    at mesh-build). `uv_sample` may sit in a different copy if the chain's
    polyline was close()-extended across a seam. Lifting via Su, Sv at the
    extended uv yields the wrong 3D vector (e.g. for Möbius, Su(v+2π) = -Su(v)).
    Bring uv_sample into canonical p's copy before evaluating Su, Sv so the
    lift agrees with edge["dir"]'s frame.

    `Su`, `Sv` may be supplied pre-computed at the seam-canonical lift uv
    (one `_eval_all` per sample, shared with `_tan_for_bc_sample` and the
    outer sample_xyz fill).
    """
    edge = mesh.edges[int(seg_idx)]
    dir_uv = np.asarray(edge["dir"], dtype=float)
    uv_for_lift = _uv_for_bc_lift(edge, mesh, uv_sample)
    if Su is None or Sv is None:
        u, v = float(uv_for_lift[0]), float(uv_for_lift[1])
        Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
        Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    inward_3d = dir_uv[0] * Su + dir_uv[1] * Sv
    return projection.proj_vec(uv_for_lift, inward_3d)


def _dir_for_cc_sample(seg_idx_poly: int, alpha: float,
                       sub_internal_seg_idx: int, css: np.ndarray, cps: np.ndarray
                       ) -> np.ndarray:
    """Interpolated CP d-field at a CC sample."""
    cs = css[int(sub_internal_seg_idx)]
    d_p = np.asarray(cps[int(cs["p_cp"])]["d"], dtype=float)
    d_q = np.asarray(cps[int(cs["q_cp"])]["d"], dtype=float)
    return (1.0 - float(alpha)) * d_p + float(alpha) * d_q


def _tan_for_bc_sample(seg_idx: int, mesh, surface, projection, domain,
                       uv_sample: np.ndarray, *,
                       Su=None, Sv=None) -> np.ndarray:
    """Analytic 2D image-space TANGENT of BC at uv, sign-matched to edge chord.

    Replaces the chord `xy[si+1] - xy[si]` used previously: that chord direction
    is sample-jitter-noisy at fine resolutions, which can flip the sign of the
    projection-break discriminator. The analytic tangent is resolution-stable.

    `edge["pq"]` is in canonical p's identification copy; bring `uv_sample`
    into the same copy before evaluating Su, Sv so the lift is frame-consistent
    (matters at the mo seam).

    `Su`, `Sv` may be supplied pre-computed at the seam-canonical lift uv
    (one `_eval_all` per sample, shared with `_dir_for_bc_sample`).
    """
    edge = mesh.edges[int(seg_idx)]
    edge_dp = np.asarray(edge["pq"], dtype=float)
    uv_for_lift = _uv_for_bc_lift(edge, mesh, uv_sample)
    Tb_uv = domain.boundary_tangent(uv_for_lift, edge_dp)
    if Su is None or Sv is None:
        u, v = float(uv_for_lift[0]), float(uv_for_lift[1])
        Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
        Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    Tb_3d = Tb_uv[0] * Su + Tb_uv[1] * Sv
    return projection.proj_vec(uv_for_lift, Tb_3d)


def _tan_for_cc_sample(uv_sample: np.ndarray, css: np.ndarray, cps: np.ndarray,
                       chain_seg: int, surface, projection, *,
                       S_p=None, Su=None, Sv=None,
                       Suu=None, Suv=None, Svv=None) -> np.ndarray:
    """Analytic 2D image-space TANGENT of CC silhouette at uv, sign-matched to
    the chain direction (from p_cp to q_cp of the cs).

    Silhouette curve in uv: `axis · SN = 0`. Its tangent is perpendicular to
    `Np = ∇_uv(axis · SN) = ((Suu×Sv + Su×Suv)·axis, (Suv×Sv + Su×Svv)·axis)`.
    Lifted via dS and projected to image.

    Persp (2026-05-27): `axis = viewer_direction(S(uv_sample))`. The
    extra gradient terms `−Su·SN` / `−Sv·SN` vanish identically (triple
    product with repeated vector), so the Np formula is unchanged.

    `S_p`, `Su`, `Sv`, `Suu`, `Suv`, `Svv` may be supplied pre-computed
    (one `_eval_all` per sample, shared with the outer sample_xyz fill).
    """
    if S_p is None:
        u, v = float(uv_sample[0]), float(uv_sample[1])
        S_p, Su, Sv, Suu, Suv, Svv, _SN = surface._eval_all(u, v)
        S_p = np.asarray(S_p, dtype=float).reshape(3)
        Su = np.asarray(Su, dtype=float).reshape(3)
        Sv = np.asarray(Sv, dtype=float).reshape(3)
        Suu = np.asarray(Suu, dtype=float).reshape(3)
        Suv = np.asarray(Suv, dtype=float).reshape(3)
        Svv = np.asarray(Svv, dtype=float).reshape(3)
    axis = projection.viewer_direction(S_p).reshape(3)
    Np = np.array([
        float(np.cross(Suu, Sv) @ axis + np.cross(Su, Suv) @ axis),
        float(np.cross(Suv, Sv) @ axis + np.cross(Su, Svv) @ axis),
    ])
    Tp_uv = np.array([-Np[1], Np[0]])
    # Sign-match against the cs chain direction in uv (p_cp → q_cp).
    cs = css[int(chain_seg)]
    p_uv = np.asarray(cps[int(cs["p_cp"])]["uv"], dtype=float)
    q_uv = np.asarray(cps[int(cs["q_cp"])]["uv"], dtype=float)
    chain_dir = q_uv - p_uv
    if float(Tp_uv @ chain_dir) < 0.0:
        Tp_uv = -Tp_uv
    Tp_3d = Tp_uv[0] * Su + Tp_uv[1] * Sv
    return projection.proj_vec(uv_sample, Tp_3d)


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
    """Batched _tan_for_cc_sample. See its docstring for the math.

    Uses the scalar-triple-product identity `(a × b) · c = det([a, b, c])` to
    avoid the per-sample `np.cross + @` overhead, and inlines `proj_vec` so
    every sample is processed in one numpy pass.
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
    vertices) out to `min(L_BC/2, L_CC/2)`, whichever half comes first. Beyond
    that matched extent (the un-matched middle) the BC reverts to its own raw
    vertices (`own_cum`).

    Rationale: at a BCP the BC and the tangent CC are near-coincident in the
    image; sampling them at *identical* arclengths-from-the-shared-origin makes
    their polylines coincide, so they cannot weave (the old UNION kept the BC's
    own vertices in the matched region, which bulged off the CC's chord and
    produced spurious grazing crossings — the Onde fold-tip "lune" ±4 defect).
    For a lune (both endpoints share the same CC, near-equal lengths) the two
    matched halves meet → the BC becomes a vertex-exact copy of the CC.

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

    # Matched region from each BCP end: CC vertex arclengths-from-SP, capped at
    # min(BC half, CC half). `ext_*` is the matched extent (0 = no match).
    ext_start = 0.0
    ext_end = 0.0
    matched: list[float] = []  # arclengths-from-the-BC-start of matched samples

    if _sp_is_bcp(sub_bc.start):
        T_start = _sub_outgoing_xy_tangent(sub_bc, sub_bc.start, xy_p)
        start_arcs, L_cc = _pick_neighbour_arclengths(
            sub_bc.start, T_start, subcurves, polys, half_L,
            exclude_kinds=("BC", "HC", "SIC"), self_idx=i_bc,
            return_length=True,
        )
        if start_arcs is not None:
            ext_start = min(half_L, L_cc / 2.0)
            matched.extend(float(a) for a in start_arcs if a <= ext_start)
    if _sp_is_bcp(sub_bc.end):
        T_end = _sub_outgoing_xy_tangent(sub_bc, sub_bc.end, xy_p)
        end_arcs, L_cc = _pick_neighbour_arclengths(
            sub_bc.end, T_end, subcurves, polys, half_L,
            exclude_kinds=("BC", "HC", "SIC"), self_idx=i_bc,
            return_length=True,
        )
        if end_arcs is not None:
            ext_end = min(half_L, L_cc / 2.0)
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
    # Collapsed-SubCurve guard (2026-05-27). A SubCurve whose xy polyline
    # is shorter than `1e-4 * ell` has both SPs at essentially the same
    # image-space point — typically a non-generic axis-aligned view that
    # projects an entire curve to a single point (e.g., the helicoid CC
    # at u=0 under the Z-axis view collapses to (0, 0)). Without this
    # guard, delta_per_sp inherits the tiny L → `_sample_arclengths`
    # spends 10⁷+ iterations climbing back to `ell` for every other sub
    # sharing the SP. Surface the degenerate input loudly rather than
    # hanging silently.
    _L_FLOOR = 1e-4 * ell
    for i, sub in enumerate(subcurves):
        L = L_per_sub[i]
        if 0 < L < _L_FLOOR and not (sub.start == -1 and sub.end == -1):
            raise ValueError(
                f"resample_all: SubCurve {i} (kind={sub.kind}, "
                f"start={sub.start}, end={sub.end}) has xy-length "
                f"{L:.3e} < 1e-4·ell={_L_FLOOR:.3e}. Both SPs project to "
                f"essentially the same image point — likely a non-generic "
                f"axis-aligned view collapsing a curve to a single point. "
                f"Use a generic random view, or skip this fixture."
            )
    L_per_sp: dict[int, float] = {}
    for i, sub in enumerate(subcurves):
        L = L_per_sub[i]
        if L <= 0:
            continue
        for sp in (sub.start, sub.end):
            if sp >= 0:
                L_per_sp[sp] = min(L_per_sp.get(sp, float("inf")), L)
    delta_per_sp = {sp: L / 10.0 for sp, L in L_per_sp.items()}

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
            N_dense = _settings.HC_DENSIFY_N
            t_dense = np.linspace(0.0, 1.0, N_dense)
            uv_dense = uv_q0[None, :] + t_dense[:, None] * (uv_q1 - uv_q0)[None, :]
            # Batched dense xyz/xy: one `surface.S` call for all 200 uv pairs.
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
            s_targets = cum.copy()
        elif sub.kind == "BC":
            interp_nsub = _settings.BC_DENSIFY_NSUB
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

        # Phase 1 — per-sample uv resolution (polyline walk + snap + Newton).
        # No surface evals here so Phase 2 can do a single batched _eval_all.
        chain_segs = np.full(N, -1, dtype=np.int64)
        seg_ps = np.zeros(N, dtype=np.int64)
        alphas = np.zeros(N, dtype=float)
        for j, s in enumerate(s_targets):
            uv_s, dseg, alpha = _interp_along_polyline(
                interp_uv, interp_xy, interp_cum, s, domain)
            # Map dense segment back to the original build-polyline segment.
            seg_p = dseg // interp_nsub
            if (sub.kind == "BC" and domain is not None
                    and getattr(domain, "type", None) in ("disk", "annulus")):
                uv_s = _snap_annular_bc(uv_s, mesh)
            if sub.kind == "CC" and project_resampled:
                uv_s = _newton_cc_refine(uv_s, surface, projection)
            sample_uv[j] = uv_s
            seg_ps[j] = seg_p
            alphas[j] = alpha
            if sample_dir is not None:
                int_idx = seg_p
                if int_idx >= len(sub.internal):
                    int_idx = len(sub.internal) - 1
                if int_idx < 0:
                    int_idx = 0
                chain_segs[j] = (int(sub.internal[int_idx][0])
                                  if sub.internal else -1)

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

        # Phase 3 — per-sample dir/tan helpers (using precomputed derivs).
        # CC tangents are batched (see _tan_for_cc_samples_batched); BC and
        # CC dirs still loop per-sample (the BC formulas need per-sample
        # edge lookups; CC dir is a simple per-segment interp).
        if sample_dir is not None and sub.kind == "CC" and N > 0:
            cc_mask = chain_segs >= 0
            if cc_mask.any():
                idx = np.flatnonzero(cc_mask)
                cc_tans = _tan_for_cc_samples_batched(
                    sample_uv[idx], chain_segs[idx],
                    S_all[idx], Su_all[idx], Sv_all[idx],
                    Suu_all[idx], Suv_all[idx], Svv_all[idx],
                    css, cps, projection,
                )
                sample_tan[idx] = cc_tans
        if sample_dir is not None:
            for j in range(N):
                cs = int(chain_segs[j])
                if cs < 0:
                    continue
                if sub.kind == "BC":
                    sample_dir[j] = _dir_for_bc_sample(
                        cs, mesh, surface, projection, sample_uv[j],
                        Su=Su_all[j], Sv=Sv_all[j],
                    )
                    sample_tan[j] = _tan_for_bc_sample(
                        cs, mesh, surface, projection, domain, sample_uv[j],
                        Su=Su_all[j], Sv=Sv_all[j],
                    )
                elif sub.kind == "CC":
                    sample_dir[j] = _dir_for_cc_sample(
                        int(seg_ps[j]), float(alphas[j]), cs, css, cps,
                    )
                    # sample_tan[j] computed via batched call above.

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
