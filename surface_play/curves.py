"""curves.py — P4: make_lines, P6: sign_changes, C7: build_bcs"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from surface_play.mesh import Mesh


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
    Assemble boundary curves from mesh boundary edges.
    Vertex indices are mapped through vertex_class so that mo-identified
    corners are recognized as the same chain node.
    """
    if len(mesh.boundary_edge_idx) == 0:
        return []

    bnd_edges = mesh.edges[mesh.boundary_edge_idx]
    p_class = mesh.vertex_class[bnd_edges["p_idx"].astype(np.intp)]
    q_class = mesh.vertex_class[bnd_edges["q_idx"].astype(np.intp)]
    segments = np.column_stack([p_class, q_class]).astype(np.intp)

    chains = make_lines(segments)

    result = []
    for chain in chains:
        is_closed = bool(chain[0] == chain[-1])
        result.append(BoundaryCurve(
            edge_indices=np.asarray(chain, dtype=np.intp),
            is_closed=is_closed,
        ))
    return result
