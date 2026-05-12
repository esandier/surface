"""Contour curves for a surface at a given viewpoint.

O1: find_contour_points  — edges → contour points (CP)
O2: build_contour_segments — CPs → contour segments (CS)
O3: build_contour_curves — CSs → ContourCurve list
O4: find_vps — CCs → cusp points (VP)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from surface_play.curves import sign_changes, make_lines
from surface_play.mesh import Mesh
from surface_play.projection import Projection


# ── dtype definitions ─────────────────────────────────────────────────────────

cp_dtype = np.dtype([
    ("e",     "i4"),    # index into mesh.edges
    ("s",     "f8"),    # parametric position on the edge, in (0, 1)
    ("uv",    "2f8"),   # domain coord = e.p + s·e.pq
    ("xyz",   "3f8"),   # 3D position
    ("d",     "2f8"),   # projected ker_param (view-plane 2D), oriented toward viewer
    ("ptype", "u1"),    # 0 = interior CP; 4 = on boundary edge
])

cs_dtype = np.dtype([
    ("p_cp",   "i4"),   # CP index (endpoint)
    ("q_cp",   "i4"),   # CP index (endpoint)
    ("face",   "i4"),   # face on which both CPs sit
    ("split1", "i4"),   # SPT index (G17); -1 initially
    ("split2", "i4"),   # SPT index (G17); -1 initially
])


# ── O1 ────────────────────────────────────────────────────────────────────────

def find_contour_points(
    mesh: Mesh, projection: Projection, *, use_newton: bool = True
) -> np.ndarray:
    """Find contour points: positions on edges where SN·viewer_direction = 0.

    Returns a cp_dtype structured array, one entry per contour point.
    G9: Newton steps that escape (0, 1) fall back to s = 0.5.
    """
    edges = mesh.edges

    # 1. Per-vertex SN·viewer_direction
    dot_v = projection.per_vertex_viewer_dot(mesh)  # (K,)

    # 2. Sign-change detection, respecting edge flip (Möbius seam edges have flip=-1)
    flip = edges["flip"].astype(float)
    mask = sign_changes(dot_v[edges["p_idx"]], dot_v[edges["q_idx"]], flip)

    cand_idx = np.nonzero(mask)[0]
    if len(cand_idx) == 0:
        return np.zeros(0, dtype=cp_dtype)

    cand = edges[cand_idx]
    pq_u = cand["pq"][:, 0]  # (N,)
    pq_v = cand["pq"][:, 1]  # (N,)
    s = np.full(len(cand_idx), 0.5)

    # 3. Vectorised Newton iterations over candidate edges
    if use_newton:
        surface = mesh.surface
        for _ in range(10):
            u = cand["p"][:, 0] + s * pq_u
            v = cand["p"][:, 1] + s * pq_v

            vals = surface._eval_all(u, v)
            # each component is (3, N) from _eval_all
            SN_arr  = vals[6].T   # (N, 3)
            Su_arr  = vals[1].T
            Sv_arr  = vals[2].T
            Suu_arr = vals[3].T
            Suv_arr = vals[4].T
            Svv_arr = vals[5].T

            # dSN/ds = (Suu×Sv + Su×Suv)·pq_u + (Suv×Sv + Su×Svv)·pq_v
            dSN_ds = (
                (np.cross(Suu_arr, Sv_arr) + np.cross(Su_arr, Suv_arr)) * pq_u[:, None]
                + (np.cross(Suv_arr, Sv_arr) + np.cross(Su_arr, Svv_arr)) * pq_v[:, None]
            )  # (N, 3)

            if projection.mode == "ortho":
                axis = projection._axis
                f_s  = SN_arr  @ axis   # (N,)
                fp_s = dSN_ds  @ axis   # (N,)
            else:
                # f(s)  = SN · (S - eye);  f'(s) = dSN/ds · (S - eye)
                # (SN · dS/ds = 0 because SN ⊥ Su, Sv)
                xyz_arr = vals[0].T  # (N, 3)
                vd = xyz_arr - projection.eye  # (N, 3)
                f_s  = np.einsum("ij,ij->i", SN_arr, vd)
                fp_s = np.einsum("ij,ij->i", dSN_ds, vd)

            safe = np.abs(fp_s) > 1e-15
            s = np.where(safe, s - f_s / np.where(safe, fp_s, 1.0), s)

        # G9: reset any s that escaped (0, 1) to midpoint
        s = np.where((s > 0.0) & (s < 1.0), s, 0.5)

    # 4. Final domain / 3D positions
    u_f = cand["p"][:, 0] + s * pq_u
    v_f = cand["p"][:, 1] + s * pq_v
    uv_arr  = np.stack([u_f, v_f], axis=-1)   # (N, 2)
    xyz_arr = mesh.surface.S(u_f, v_f).T       # (N, 3)

    # 5. d field: project kerdS, orient toward viewer (positive depth component)
    d_arr = np.zeros((len(cand_idx), 2))
    for i in range(len(cand_idx)):
        uv_i     = uv_arr[i]
        kerdS_3d = projection.kerdS(uv_i)
        if projection._axis @ kerdS_3d < 0.0:
            kerdS_3d = -kerdS_3d
        d_arr[i] = projection.proj_vec(uv_i, kerdS_3d)

    # 6. ptype: 4 on boundary edge (g == -1), else 0
    ptype_arr = np.where(cand["g"] == -1, np.uint8(4), np.uint8(0))

    result = np.zeros(len(cand_idx), dtype=cp_dtype)
    result["e"]     = cand_idx
    result["s"]     = s
    result["uv"]    = uv_arr
    result["xyz"]   = xyz_arr
    result["d"]     = d_arr
    result["ptype"] = ptype_arr
    return result


# ── O2 ────────────────────────────────────────────────────────────────────────

def build_contour_segments(cps: np.ndarray, mesh: Mesh) -> np.ndarray:
    """Pair contour points within each face to form contour segments (CS).

    For each face: collect CPs whose edge belongs to that face; if exactly 2,
    emit one CS.  By the sign-change parity argument a face always has 0 or 2
    such CPs, never 1 or 3.  Splits initialised to -1 (G17).
    """
    if len(cps) == 0:
        return np.zeros(0, dtype=cs_dtype)

    # edge index → CP index; -1 where no CP (each edge carries at most one CP)
    edge_to_cp = np.full(len(mesh.edges), -1, dtype=np.int32)
    edge_to_cp[cps["e"]] = np.arange(len(cps), dtype=np.int32)

    # For every face look up its 3 edges' CP indices: shape (M, 3)
    face_edges  = mesh.faces["edges"]          # (M, 3) — edge indices per face
    cp_in_face  = edge_to_cp[face_edges]       # (M, 3) — CP index or -1

    # Faces with exactly 2 CPs
    cp_count = (cp_in_face >= 0).sum(axis=1)  # (M,)
    paired   = np.nonzero(cp_count == 2)[0]   # face indices

    if len(paired) == 0:
        return np.zeros(0, dtype=cs_dtype)

    # Extract the two valid CP indices per paired face.
    # np.sort puts the -1 entry first (smallest), so columns 1 and 2 are valid.
    cp_vals    = cp_in_face[paired]              # (n, 3)
    sorted_cp  = np.sort(cp_vals, axis=1)        # (n, 3): [-1, p, q]

    out = np.zeros(len(paired), dtype=cs_dtype)
    out["p_cp"]   = sorted_cp[:, 1]
    out["q_cp"]   = sorted_cp[:, 2]
    out["face"]   = paired
    out["split1"] = -1
    out["split2"] = -1
    return out


# ── O3 ────────────────────────────────────────────────────────────────────────

@dataclass
class ContourCurve:
    cs_indices: np.ndarray  # 1D int array; sign-encoded as in make_lines output
    is_closed: bool


def build_contour_curves(css: np.ndarray, cps: np.ndarray) -> list[ContourCurve]:
    """Chain contour segments into contour curves via make_lines.

    Two CSs share an endpoint iff they share a CP index (p_cp / q_cp).
    CCs are closed iff make_lines returns a chain whose first and last
    signed indices have the same absolute value.
    """
    if len(css) == 0:
        return []

    segments = np.column_stack([
        css["p_cp"].astype(np.intp),
        css["q_cp"].astype(np.intp),
    ])
    chains = make_lines(segments)

    result = []
    for chain in chains:
        ch = np.asarray(chain, dtype=np.intp)
        is_closed = len(ch) > 1 and abs(int(ch[0])) == abs(int(ch[-1]))
        result.append(ContourCurve(cs_indices=ch, is_closed=is_closed))
    return result
