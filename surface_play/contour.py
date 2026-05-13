"""Contour curves for a surface at a given viewpoint.

O1: find_contour_points  â€” edges â†’ contour points (CP)
O2: build_contour_segments â€” CPs â†’ contour segments (CS)
O3: build_contour_curves â€” CSs â†’ ContourCurve list
O4: find_vps â€” CCs â†’ cusp points (VP)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from surface_play.curves import sign_changes, make_lines
from surface_play.mesh import Mesh
from surface_play.projection import Projection
from surface_play.surface import SurfaceParams


# â”€â”€ dtype definitions â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

cp_dtype = np.dtype([
    ("e",     "i4"),    # index into mesh.edges
    ("s",     "f8"),    # parametric position on the edge, in (0, 1)
    ("uv",    "2f8"),   # domain coord = e.p + sÂ·e.pq
    ("xyz",   "3f8"),   # 3D position
    ("d",     "2f8"),   # unit 2D curvature direction of the projected contour at the
                        # CP (see _compute_d_at): non-zero at non-cusp CPs and reverses
                        # across cusps.  Drives O4 sign-change detection.
    ("ptype", "u1"),    # 0 = interior CP; 4 = on boundary edge
])

cs_dtype = np.dtype([
    ("p_cp",   "i4"),   # CP index (endpoint)
    ("q_cp",   "i4"),   # CP index (endpoint)
    ("face",   "i4"),   # face on which both CPs sit
    ("split1", "i4"),   # SPT index (G17); -1 initially
    ("split2", "i4"),   # SPT index (G17); -1 initially
])


# â”€â”€ O1 helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def _compute_d_batch(
    uv_arr: np.ndarray,
    vals: tuple[np.ndarray, ...],
    projection: Projection,
) -> np.ndarray:
    """Vectorized 2D curvature direction of the projected contour at N uv points.

    `vals` is the 7-tuple returned by `SurfaceParams._eval_all(u, v)` evaluated
    at the same N points â€” passed in so callers can share one evaluation
    between xyz, the d-field, and any other consumer (Layer-O steps run once
    per viewpoint; recomputing Su/Sv/Suu/Suv/Svv per-CP in a python loop is
    the previous hot path we are eliminating).

    Algorithm (legacy silhouette.py "dir_vec", lines 2025-2034, vectorized):
      ker  = right-singular vector of the 2Ã—2 view-plane Jacobian
             [IÂ·Su  IÂ·Sv ; JÂ·Su  JÂ·Sv]  (with persp rank-1 correction)
             for the smallest singular value â€” at every CP at once.
      dÂ²S  = SuuÂ·k0Â² + 2Â·SuvÂ·k0Â·k1 + SvvÂ·k1Â²       (3D, quadratic so the
                                                    sign ambiguity of ker
                                                    is irrelevant)
      tan  = SuÂ·(-k1) + SvÂ·k0                        (image-plane contour
                                                    tangent in 3D)
      diff2, im = proj_vec(dÂ²S), proj_vec(tan)       (both 2D)
      diff2 â† diff2 âˆ’ (diff2Â·im / |im|Â²) Â· im        (strip tangent)
      d    = diff2 / |diff2|                          (unit vector)
    """
    N = len(uv_arr)
    if N == 0:
        return np.zeros((0, 2), dtype=float)

    S_arr   = vals[0].T   # (N, 3)
    Su_arr  = vals[1].T
    Sv_arr  = vals[2].T
    Suu_arr = vals[3].T
    Suv_arr = vals[4].T
    Svv_arr = vals[5].T

    I = projection.I
    J = projection.J
    axis = projection._axis

    # 2Ã—2 view-plane Jacobian per CP.
    ISu = Su_arr @ I;  ISv = Sv_arr @ I
    JSu = Su_arr @ J;  JSv = Sv_arr @ J

    if projection.mode == "ortho":
        J_arr = np.empty((N, 2, 2), dtype=float)
        J_arr[:, 0, 0] = ISu;  J_arr[:, 0, 1] = ISv
        J_arr[:, 1, 0] = JSu;  J_arr[:, 1, 1] = JSv
    else:
        # Persp: J[r, c] = (basis_r Â· Sc) âˆ’ (basis_r Â· d / z) Â· (axis Â· Sc)
        # with d = S âˆ’ eye, z = axis Â· d.  The uniform 1/z prefactor is dropped
        # (irrelevant for the SVD right-singular vector â€” same direction).
        d3 = S_arr - projection.eye                   # (N, 3)
        z  = d3 @ axis                                # (N,)
        if (z == 0.0).any():
            raise ValueError(
                "ker_param undefined: a CP lies on the image plane through eye"
            )
        a_over_z = (d3 @ I) / z                       # (N,)
        b_over_z = (d3 @ J) / z
        nSu = Su_arr @ axis                           # (N,)
        nSv = Sv_arr @ axis
        J_arr = np.empty((N, 2, 2), dtype=float)
        J_arr[:, 0, 0] = ISu - a_over_z * nSu;  J_arr[:, 0, 1] = ISv - a_over_z * nSv
        J_arr[:, 1, 0] = JSu - b_over_z * nSu;  J_arr[:, 1, 1] = JSv - b_over_z * nSv

    _U, _S, Vh = np.linalg.svd(J_arr)                  # batched 2Ã—2 SVD
    ker = Vh[:, -1, :]                                 # (N, 2)

    k0 = ker[:, 0:1]   # (N, 1) â€” broadcasts over the 3D vectors
    k1 = ker[:, 1:2]

    # dÂ²S in kernel direction, and the contour tangent in 3D.
    d2S = Suu_arr * (k0 ** 2) + 2.0 * Suv_arr * (k0 * k1) + Svv_arr * (k1 ** 2)
    tan = Su_arr  * (-k1)     + Sv_arr  * k0

    # Project both to the view plane (vectorized â€” N at a time).
    if projection.mode == "ortho":
        diff2 = np.stack([d2S @ I, d2S @ J], axis=-1)    # (N, 2)
        im    = np.stack([tan @ I, tan @ J], axis=-1)
    else:
        # Persp proj_vec(uv, v3) = (1/z) Â· ((IÂ·v3) âˆ’ (a/z)(axisÂ·v3),
        #                                   (JÂ·v3) âˆ’ (b/z)(axisÂ·v3))
        a_over_z = (d3 @ I) / z   # noqa: F811  (recomputed for clarity)
        b_over_z = (d3 @ J) / z
        z_inv = 1.0 / z
        nd2S = d2S @ axis
        ntan = tan @ axis
        diff2 = np.stack([
            ((d2S @ I) - a_over_z * nd2S) * z_inv,
            ((d2S @ J) - b_over_z * nd2S) * z_inv,
        ], axis=-1)
        im = np.stack([
            ((tan @ I) - a_over_z * ntan) * z_inv,
            ((tan @ J) - b_over_z * ntan) * z_inv,
        ], axis=-1)

    # Strip the tangent component (Gram-Schmidt), then normalize.
    im_sq = (im * im).sum(axis=1)                      # (N,)
    diff2_dot_im = (diff2 * im).sum(axis=1)
    safe_im = im_sq > 0.0
    coef = np.where(safe_im, diff2_dot_im / np.where(safe_im, im_sq, 1.0), 0.0)
    diff2 = diff2 - coef[:, None] * im

    n = np.linalg.norm(diff2, axis=1)                  # (N,)
    safe_n = n > 0.0
    d_arr = np.zeros_like(diff2)
    d_arr[safe_n] = diff2[safe_n] / n[safe_n, None]
    return d_arr


def _compute_d_at(uv: np.ndarray, surface: SurfaceParams,
                  projection: Projection) -> np.ndarray:
    """Scalar fallback of `_compute_d_batch` for one uv point.

    Used by O4 refinement's per-Newton-iteration evaluations.  Delegates to
    the batched helper with N = 1 so the algorithm lives in one place.
    """
    uv = np.asarray(uv, dtype=float).reshape(2)
    u_ = np.array([uv[0]])
    v_ = np.array([uv[1]])
    vals = surface._eval_all(u_, v_)
    return _compute_d_batch(uv.reshape(1, 2), vals, projection)[0]


# â”€â”€ O1 â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def find_contour_points(
    mesh: Mesh, projection: Projection, *, use_newton: bool = True
) -> np.ndarray:
    """Find contour points: positions on edges where SNÂ·viewer_direction = 0.

    Returns a cp_dtype structured array, one entry per contour point.
    G9: Newton steps that escape (0, 1) fall back to s = 0.5.
    """
    edges = mesh.edges

    # 1. Per-vertex SNÂ·viewer_direction
    dot_v = projection.per_vertex_viewer_dot(mesh)  # (K,)

    # 2. Sign-change detection, respecting edge flip (MÃ¶bius seam edges have flip=-1)
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

            # dSN/ds = (SuuÃ—Sv + SuÃ—Suv)Â·pq_u + (SuvÃ—Sv + SuÃ—Svv)Â·pq_v
            dSN_ds = (
                (np.cross(Suu_arr, Sv_arr) + np.cross(Su_arr, Suv_arr)) * pq_u[:, None]
                + (np.cross(Suv_arr, Sv_arr) + np.cross(Su_arr, Svv_arr)) * pq_v[:, None]
            )  # (N, 3)

            if projection.mode == "ortho":
                axis = projection._axis
                f_s  = SN_arr  @ axis   # (N,)
                fp_s = dSN_ds  @ axis   # (N,)
            else:
                # f(s)  = SN Â· (S - eye);  f'(s) = dSN/ds Â· (S - eye)
                # (SN Â· dS/ds = 0 because SN âŠ¥ Su, Sv)
                xyz_arr = vals[0].T  # (N, 3)
                vd = xyz_arr - projection.eye  # (N, 3)
                f_s  = np.einsum("ij,ij->i", SN_arr, vd)
                fp_s = np.einsum("ij,ij->i", dSN_ds, vd)

            safe = np.abs(fp_s) > 1e-15
            s = np.where(safe, s - f_s / np.where(safe, fp_s, 1.0), s)

        # G9: reset any s that escaped (0, 1) to midpoint
        s = np.where((s > 0.0) & (s < 1.0), s, 0.5)

    # 4. Final domain coords + one shared vectorized surface evaluation feeding
    #    both xyz_arr and the d-field (avoids 5Â·N redundant scalar derivative
    #    calls per viewpoint that the previous per-CP `_compute_d_at` loop did).
    u_f = cand["p"][:, 0] + s * pq_u
    v_f = cand["p"][:, 1] + s * pq_v
    uv_arr = np.stack([u_f, v_f], axis=-1)              # (N, 2)
    vals_final = mesh.surface._eval_all(u_f, v_f)       # 7-tuple, each (3, N)
    xyz_arr = vals_final[0].T                            # (N, 3)

    # 5. d field: 2D curvature direction of the projected contour (vectorized).
    d_arr = _compute_d_batch(uv_arr, vals_final, projection)

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


# â”€â”€ O2 â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def build_contour_segments(cps: np.ndarray, mesh: Mesh) -> np.ndarray:
    """Pair contour points within each face to form contour segments (CS).

    For each face: collect CPs whose edge belongs to that face; if exactly 2,
    emit one CS.  By the sign-change parity argument a face always has 0 or 2
    such CPs, never 1 or 3.  Splits initialised to -1 (G17).
    """
    if len(cps) == 0:
        return np.zeros(0, dtype=cs_dtype)

    # edge index â†’ CP index; -1 where no CP (each edge carries at most one CP)
    edge_to_cp = np.full(len(mesh.edges), -1, dtype=np.int32)
    edge_to_cp[cps["e"]] = np.arange(len(cps), dtype=np.int32)

    # For every face look up its 3 edges' CP indices: shape (M, 3)
    face_edges  = mesh.faces["edges"]          # (M, 3) â€” edge indices per face
    cp_in_face  = edge_to_cp[face_edges]       # (M, 3) â€” CP index or -1

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


# â”€â”€ O3 â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

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


