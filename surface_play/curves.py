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


def _seg_uv_at_bary(sub_kind: str, seg_idx: int, bary: float,
                    mesh, css, sis_pairs, cps, dps,
                    sis_preimage: int | None = None) -> np.ndarray:
    """uv at bary on a chain segment. `bary` is segment-local (0..1)."""
    if sub_kind == "BC":
        edge = mesh.edges[seg_idx]
        return np.asarray(edge["p"], dtype=float) + float(bary) * np.asarray(edge["pq"], dtype=float)
    if sub_kind == "CC":
        cs = css[seg_idx]
        p_uv = np.asarray(cps[int(cs["p_cp"])]["uv"], dtype=float)
        q_uv = np.asarray(cps[int(cs["q_cp"])]["uv"], dtype=float)
        return (1.0 - float(bary)) * p_uv + float(bary) * q_uv
    if sub_kind == "SIC":
        row = sis_pairs[seg_idx]
        p = int(row["p_dp"]); q = int(row["q_dp"]); f = int(row["flip"])
        # Preimage 1: uv1[p] → (uv1 if f=+1 else uv2)[q].
        # Preimage 2: uv2[p] → (uv2 if f=+1 else uv1)[q].
        if sis_preimage == 2:
            p_uv = np.asarray(dps[p]["uv2"], dtype=float)
            q_uv = np.asarray(dps[q]["uv2" if f == 1 else "uv1"], dtype=float)
        else:
            p_uv = np.asarray(dps[p]["uv1"], dtype=float)
            q_uv = np.asarray(dps[q]["uv1" if f == 1 else "uv2"], dtype=float)
        return (1.0 - float(bary)) * p_uv + float(bary) * q_uv
    raise ValueError(f"_seg_uv_at_bary: unsupported kind {sub_kind!r}")


def _close_aware_lerp(a_uv: np.ndarray, b_uv: np.ndarray, t: float, domain) -> np.ndarray:
    if domain is not None and getattr(domain, "type", None) == "rect":
        b_loc = domain.close(a_uv, b_uv)
    else:
        b_loc = b_uv
    return a_uv + float(t) * (b_loc - a_uv)


def _build_polyline(
    sub: "SubCurve", splits: "SplitArrays", mesh: "Mesh",
    css: np.ndarray, sis_pairs: np.ndarray,
    cps: np.ndarray, dps: np.ndarray,
    surface: "SurfaceParams", projection: "Projection",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (uv_poly, xyz_poly, xy_poly) — vertices of the SC's polyline.

    HC: just [start, end] (2 vertices). SP-less closed SC: walk internal only.
    """
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
        # Normal SC: start SP → internal → end SP.
        # For SIC subs we pick the preimage of each endpoint that lies on the
        # same sheet as the chain walk (using SP.uv_alt = other preimage uv).
        # Start SP: defer the decision until the first internal sample is
        # known, then choose the start preimage closer to it.
        sp_start_uv = np.asarray(splits.sps[sub.start][0], dtype=float)
        sp_start_uv_alt = (
            np.asarray(splits.sps[sub.start][4], dtype=float)
            if sub.kind == "SIC" else None
        )
        uvs.append(sp_start_uv)
        domain_local = getattr(mesh, "domain", None)
        first_internal_for_start_pick = None
        for seg_idx, bary in sub.internal:
            if sub.kind == "SIC":
                # Pick whichever preimage's vertex is closer (close-aware) to
                # the previous polyline vertex — preserves G21 sheet identity.
                prev = uvs[-1]
                cand1 = _seg_uv_at_bary("SIC", int(seg_idx), float(bary),
                                        mesh, css, sis_pairs, cps, dps,
                                        sis_preimage=1)
                cand2 = _seg_uv_at_bary("SIC", int(seg_idx), float(bary),
                                        mesh, css, sis_pairs, cps, dps,
                                        sis_preimage=2)
                if domain_local is not None and getattr(domain_local, "type", None) == "rect":
                    cand1_close = domain_local.close(prev, cand1)
                    cand2_close = domain_local.close(prev, cand2)
                else:
                    cand1_close, cand2_close = cand1, cand2
                d1 = float(np.linalg.norm(cand1_close - prev))
                d2 = float(np.linalg.norm(cand2_close - prev))
                picked = cand1 if d1 <= d2 else cand2
                uvs.append(picked)
                if first_internal_for_start_pick is None:
                    first_internal_for_start_pick = picked
            else:
                uvs.append(_seg_uv_at_bary(sub.kind, int(seg_idx), float(bary),
                                           mesh, css, sis_pairs, cps, dps))
        # Retroactively fix the SIC start SP if uv_alt is closer to the first
        # internal sample (start uv was just a placeholder above).
        if (sub.kind == "SIC"
                and sp_start_uv_alt is not None
                and not np.any(np.isnan(sp_start_uv_alt))
                and first_internal_for_start_pick is not None):
            target = first_internal_for_start_pick
            if domain_local is not None and getattr(domain_local, "type", None) == "rect":
                a_close = domain_local.close(target, sp_start_uv)
                b_close = domain_local.close(target, sp_start_uv_alt)
            else:
                a_close, b_close = sp_start_uv, sp_start_uv_alt
            if (np.linalg.norm(a_close - target)
                    > np.linalg.norm(b_close - target)):
                uvs[0] = sp_start_uv_alt
        # End SP: same picking against the last internal sample (or start SP
        # if there is no internal sample).
        sp_end_uv = np.asarray(splits.sps[sub.end][0], dtype=float)
        if sub.kind == "SIC":
            sp_end_uv_alt = np.asarray(splits.sps[sub.end][4], dtype=float)
            if not np.any(np.isnan(sp_end_uv_alt)):
                prev = uvs[-1]
                if domain_local is not None and getattr(domain_local, "type", None) == "rect":
                    a_close = domain_local.close(prev, sp_end_uv)
                    b_close = domain_local.close(prev, sp_end_uv_alt)
                else:
                    a_close, b_close = sp_end_uv, sp_end_uv_alt
                if (np.linalg.norm(a_close - prev)
                        > np.linalg.norm(b_close - prev)):
                    sp_end_uv = sp_end_uv_alt
        uvs.append(sp_end_uv)

    # Close-aware adjust consecutive vertices, then lift to xyz/xy.
    domain = getattr(mesh, "domain", None)
    if domain is not None and getattr(domain, "type", None) == "rect":
        for i in range(1, len(uvs)):
            uvs[i] = domain.close(uvs[i - 1], uvs[i])

    uv_arr = np.asarray(uvs, dtype=float) if uvs else np.zeros((0, 2), dtype=float)
    xyz_arr = np.zeros((len(uv_arr), 3), dtype=float)
    for i in range(len(uv_arr)):
        xyz_arr[i] = surface.S(float(uv_arr[i, 0]), float(uv_arr[i, 1]))
    xy_arr = projection.XY(xyz_arr) if len(xyz_arr) else np.zeros((0, 2), dtype=float)
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


def _dir_for_bc_sample(seg_idx: int, mesh, surface, projection,
                       uv_sample: np.ndarray) -> np.ndarray:
    """Projected inward 2D normal of a BC sample (from edge['dir'] at uv_sample).

    `edge["dir"]` is expressed in the canonical p's identification copy (fixed
    at mesh-build). `uv_sample` may sit in a different copy if the chain's
    polyline was close()-extended across a seam. Lifting via Su, Sv at the
    extended uv yields the wrong 3D vector (e.g. for Möbius, Su(v+2π) = -Su(v)).
    Bring uv_sample into canonical p's copy before evaluating Su, Sv so the
    lift agrees with edge["dir"]'s frame.
    """
    edge = mesh.edges[int(seg_idx)]
    dir_uv = np.asarray(edge["dir"], dtype=float)
    domain = getattr(mesh, "domain", None)
    if domain is not None and getattr(domain, "type", None) == "rect":
        p_canonical = mesh.uv[int(edge["p_idx"])]
        uv_for_lift = domain.close(p_canonical, uv_sample)
    else:
        uv_for_lift = uv_sample
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
                       uv_sample: np.ndarray) -> np.ndarray:
    """Analytic 2D image-space TANGENT of BC at uv, sign-matched to edge chord.

    Replaces the chord `xy[si+1] - xy[si]` used previously: that chord direction
    is sample-jitter-noisy at fine resolutions, which can flip the sign of the
    projection-break discriminator. The analytic tangent is resolution-stable.

    `edge["pq"]` is in canonical p's identification copy; bring `uv_sample`
    into the same copy before evaluating Su, Sv so the lift is frame-consistent
    (matters at the mo seam).
    """
    edge = mesh.edges[int(seg_idx)]
    edge_dp = np.asarray(edge["pq"], dtype=float)
    if domain is not None and getattr(domain, "type", None) == "rect":
        p_canonical = mesh.uv[int(edge["p_idx"])]
        uv_for_lift = domain.close(p_canonical, uv_sample)
    else:
        uv_for_lift = uv_sample
    Tb_uv = domain.boundary_tangent(uv_for_lift, edge_dp)
    u, v = float(uv_for_lift[0]), float(uv_for_lift[1])
    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    Tb_3d = Tb_uv[0] * Su + Tb_uv[1] * Sv
    return projection.proj_vec(uv_for_lift, Tb_3d)


def _tan_for_cc_sample(uv_sample: np.ndarray, css: np.ndarray, cps: np.ndarray,
                       chain_seg: int, surface, projection) -> np.ndarray:
    """Analytic 2D image-space TANGENT of CC silhouette at uv, sign-matched to
    the chain direction (from p_cp to q_cp of the cs).

    Silhouette curve in uv: `axis · SN = 0`. Its tangent is perpendicular to
    `Np = ∇_uv(axis · SN) = ((Suu×Sv + Su×Suv)·axis, (Suv×Sv + Su×Svv)·axis)`.
    Lifted via dS and projected to image.

    Persp (2026-05-27): `axis = viewer_direction(S(uv_sample))`. The
    extra gradient terms `−Su·SN` / `−Sv·SN` vanish identically (triple
    product with repeated vector), so the Np formula is unchanged.
    """
    u, v = float(uv_sample[0]), float(uv_sample[1])
    S_p = np.asarray(surface.S(u, v), dtype=float).reshape(3)
    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    Suu = np.asarray(surface.Suu(u, v), dtype=float).reshape(3)
    Suv = np.asarray(surface.Suv(u, v), dtype=float).reshape(3)
    Svv = np.asarray(surface.Svv(u, v), dtype=float).reshape(3)
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


def _newton_cc_refine(uv: np.ndarray, surface, projection, *, n_iter: int = 5
                      ) -> np.ndarray:
    """Run Newton along chord-normal to find nearby axis·SN = 0.

    Direction = uv-grad of `axis·SN`. This is rough but sufficient for the
    G7 / PROJECT_RESAMPLED test (the test only asserts convergence to a sign
    change, not a tight tolerance).

    Persp (2026-05-27): `axis = viewer_direction(S(uv))` recomputed each
    iteration as `uv` moves.
    """
    uv_cur = uv.copy()
    for _ in range(n_iter):
        u, v = float(uv_cur[0]), float(uv_cur[1])
        S_p = np.asarray(surface.S(u, v), dtype=float).reshape(3)
        axis = projection.viewer_direction(S_p).reshape(3)
        SN = np.asarray(surface.SN(u, v), dtype=float).reshape(3)
        f = float(axis @ SN)
        if abs(f) < 1e-10:
            break
        # Numerical gradient in uv (small fd step). In persp the axis
        # itself depends on S(u,v), so the finite-difference SN_u/SN_v
        # are evaluated with the same axis at uv_cur (consistent with
        # Newton's Jacobian approximation).
        h = 1e-5
        SN_u = np.asarray(surface.SN(u + h, v), dtype=float).reshape(3)
        SN_v = np.asarray(surface.SN(u, v + h), dtype=float).reshape(3)
        gu = (float(axis @ SN_u) - f) / h
        gv = (float(axis @ SN_v) - f) / h
        g2 = gu * gu + gv * gv
        if g2 < 1e-20:
            break
        uv_cur[0] -= f * gu / g2
        uv_cur[1] -= f * gv / g2
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

    # Step E — sample each SC.
    out: list[ResampledCurve] = []
    for i, sub in enumerate(subcurves):
        uv_p, xyz_p, xy_p = polys[i]
        L_total = L_per_sub[i]

        # SP-less closed SC → verbatim copy.
        if sub.start == -1 and sub.end == -1:
            depth = np.array([projection.Z(p) for p in xyz_p], dtype=float)
            out.append(ResampledCurve(
                kind=sub.kind, start=-1, end=-1,
                depth=depth, xy=xy_p.copy(), dir=None,
                vc_in=int(sub.vc_in), vc_out=int(sub.vc_out),
            ))
            continue

        # HC → resample the uv straight line using the same arclength
        # logic as BC/CC, so the HC has comparable per-segment xy density.
        # Step 1: pre-sample uniformly in uv (dense) → curved xy polyline +
        # cumulative xy arclength. Step 2: use `_sample_arclengths` to pick
        # target arclengths just like BC. Step 3: interpolate back to uv.
        if sub.kind == "HC":
            uv_q0 = uv_p[0]; uv_q1 = uv_p[1]
            N_dense = 200
            t_dense = np.linspace(0.0, 1.0, N_dense)
            uv_dense = uv_q0[None, :] + t_dense[:, None] * (uv_q1 - uv_q0)[None, :]
            xyz_dense = np.empty((N_dense, 3), dtype=float)
            for j, uv_s in enumerate(uv_dense):
                xyz_dense[j] = np.asarray(
                    surface.S(float(uv_s[0]), float(uv_s[1])), dtype=float
                ).reshape(3)
            xy_dense = np.array([projection.XY(p) for p in xyz_dense], dtype=float)
            cum_dense = _arclengths(xy_dense)
            L_xy_total = float(cum_dense[-1])

            delta_s = delta_per_sp.get(sub.start, ell)
            delta_e = delta_per_sp.get(sub.end, ell)
            s_targets = _sample_arclengths(L_xy_total, ell, delta_s, delta_e, False)

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
                sample_xyz[j] = np.asarray(
                    surface.S(float(uv_s[0]), float(uv_s[1])), dtype=float
                ).reshape(3)
                sample_xy[j] = projection.XY(sample_xyz[j])
            depth = np.array([projection.Z(p) for p in sample_xyz], dtype=float)
            out.append(ResampledCurve(
                kind="HC", start=sub.start, end=sub.end,
                depth=depth, xy=sample_xy, dir=None,
                vc_in=int(sub.vc_in), vc_out=int(sub.vc_out),
            ))
            continue

        # Pick sample arclengths (two-phase tapered spacing).
        delta_s = delta_per_sp.get(sub.start, ell)
        delta_e = delta_per_sp.get(sub.end, ell)
        s_targets = _sample_arclengths(L_total, ell, delta_s, delta_e, sub.is_closed)

        cum = _arclengths(xy_p)
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

        for j, s in enumerate(s_targets):
            uv_s, seg_p, alpha = _interp_along_polyline(uv_p, xy_p, cum, s, domain)

            # Annular BC snap.
            if (sub.kind == "BC" and domain is not None
                    and getattr(domain, "type", None) in ("disk", "annulus")):
                uv_s = _snap_annular_bc(uv_s, mesh)

            # CC Newton refinement.
            if sub.kind == "CC" and project_resampled:
                uv_s = _newton_cc_refine(uv_s, surface, projection)

            sample_uv[j] = uv_s
            sample_xyz[j] = np.asarray(
                surface.S(float(uv_s[0]), float(uv_s[1])), dtype=float
            ).reshape(3)
            sample_xy[j] = projection.XY(sample_xyz[j])

            if sample_dir is not None:
                # Resolve owning seg_idx (BC: the chain segment containing poly-segment seg_p;
                # CC: the css idx). poly-segment seg_p goes from poly-vertex seg_p (SP or internal)
                # to seg_p+1.
                int_idx = seg_p  # internal[seg_p] holds the END vertex of poly-segment seg_p,
                                  # which is the chain segment THIS sample sits on.
                if int_idx >= len(sub.internal):
                    int_idx = len(sub.internal) - 1
                if int_idx < 0:
                    int_idx = 0
                chain_seg = int(sub.internal[int_idx][0]) if sub.internal else -1
                if sub.kind == "BC" and chain_seg >= 0:
                    sample_dir[j] = _dir_for_bc_sample(
                        chain_seg, mesh, surface, projection, uv_s,
                    )
                    sample_tan[j] = _tan_for_bc_sample(
                        chain_seg, mesh, surface, projection, domain, uv_s,
                    )
                elif sub.kind == "CC" and chain_seg >= 0:
                    sample_dir[j] = _dir_for_cc_sample(
                        seg_p, alpha, chain_seg, css, cps,
                    )
                    sample_tan[j] = _tan_for_cc_sample(
                        uv_s, css, cps, chain_seg, surface, projection,
                    )

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

        depth = np.array([projection.Z(p) for p in sample_xyz], dtype=float)
        out.append(ResampledCurve(
            kind=sub.kind, start=sub.start, end=sub.end,
            depth=depth, xy=sample_xy, dir=sample_dir, tan=sample_tan,
            vc_in=int(sub.vc_in), vc_out=int(sub.vc_out),
        ))

    return out
