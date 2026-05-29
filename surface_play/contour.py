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
    ("e",        "i4"),    # index into mesh.edges
    ("s",        "f8"),    # parametric position on the edge, in (0, 1)
    ("uv",       "2f8"),   # domain coord = e.p + sÂ·e.pq
    ("xyz",      "3f8"),   # 3D position
    ("d",        "2f8"),   # unit 2D curvature direction of the projected contour at the
                           # CP (see _compute_d_at): non-zero at non-cusp CPs and reverses
                           # across cusps.  Drives O4 sign-change detection.
    ("front_uv", "2f8"),   # 2D uv-space vector pointing toward the FRONT sheet of the
                           # surface. Computed as oriented Np = grad_uv(axis.SN) per the
                           # 4-step recipe (compute ker, orient so axis.dS(ker)>0, compute
                           # Np, orient so Np.ker>0). Used by O4 (cusp detection between
                           # CS endpoints), O7 BCP vis_chge, and O11 CDP CS-side vis_chge.
    ("ptype",    "u1"),    # 0 = interior CP; 4 = on boundary edge
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
) -> tuple[np.ndarray, np.ndarray]:
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
        # Persp: J[r, c] = (basis_r · Sc) + (basis_r · d / z) · (axis · Sc)
        # with d = S − eye, z = axis · (eye − S) > 0 (under "axis toward viewer").
        # The uniform 1/z prefactor is dropped (SVD right-singular vector is
        # sign/scale-blind).
        d3 = S_arr - projection.eye                   # (N, 3); into-scene
        z  = -(d3 @ axis)                             # (N,); = axis·(eye−S) > 0 in-scene
        if (z == 0.0).any():
            raise ValueError(
                "ker_param undefined: a CP lies on the image plane through eye"
            )
        a_over_z = (d3 @ I) / z                       # (N,)
        b_over_z = (d3 @ J) / z
        nSu = Su_arr @ axis                           # (N,)
        nSv = Sv_arr @ axis
        J_arr = np.empty((N, 2, 2), dtype=float)
        J_arr[:, 0, 0] = ISu + a_over_z * nSu;  J_arr[:, 0, 1] = ISv + a_over_z * nSv
        J_arr[:, 1, 0] = JSu + b_over_z * nSu;  J_arr[:, 1, 1] = JSv + b_over_z * nSv

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
        # Persp proj_vec(uv, v3) = (1/z) · ((I·v3) + (a/z)(axis·v3),
        #                                   (J·v3) + (b/z)(axis·v3))
        # with z = axis·(eye − S) > 0 (in-scene); d = S − eye, a = I·d, b = J·d.
        a_over_z = (d3 @ I) / z   # noqa: F811  (recomputed for clarity)
        b_over_z = (d3 @ J) / z
        z_inv = 1.0 / z
        nd2S = d2S @ axis
        ntan = tan @ axis
        diff2 = np.stack([
            ((d2S @ I) + a_over_z * nd2S) * z_inv,
            ((d2S @ J) + b_over_z * nd2S) * z_inv,
        ], axis=-1)
        im = np.stack([
            ((tan @ I) + a_over_z * ntan) * z_inv,
            ((tan @ J) + b_over_z * ntan) * z_inv,
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

    # ── front_uv: oriented Np = grad_uv(axis . SN) ──────────────────────────
    # `axis` is per-point in persp (viewer_direction(S)), constant in ortho.
    # Matches the convention used by `_bvis_chge` in splitting.py and the
    # other vis_chge sites. The Np simplification still holds: any extra
    # gradient term `d(viewer_direction)/duv · SN = -dS/duv · SN` vanishes
    # via the triple product (Su · SN = Sv · SN = 0).
    if projection.mode == "ortho":
        axis_per = np.broadcast_to(axis, S_arr.shape)   # (N, 3)
    else:
        axis_per = projection.eye - S_arr               # (N, 3)

    # Step A. Orient `ker` so axis . dS(ker) > 0.
    dS_ker = k0 * Su_arr + k1 * Sv_arr                  # (N, 3)
    axis_dot_dS_ker = np.einsum("ij,ij->i", dS_ker, axis_per)
    sign_ker = np.where(axis_dot_dS_ker >= 0.0, 1.0, -1.0)
    ker_oriented = ker * sign_ker[:, None]              # (N, 2)

    # Step B. Compute Np = grad_uv(axis . SN) in uv space (per _bvis_chge).
    cross_a = np.cross(Suu_arr, Sv_arr) + np.cross(Su_arr, Suv_arr)  # (N, 3)
    cross_b = np.cross(Suv_arr, Sv_arr) + np.cross(Su_arr, Svv_arr)  # (N, 3)
    Np_u = np.einsum("ij,ij->i", cross_a, axis_per)     # (N,)
    Np_v = np.einsum("ij,ij->i", cross_b, axis_per)     # (N,)
    Np = np.stack([Np_u, Np_v], axis=-1)                # (N, 2)

    # Step C. Orient Np so Np . ker_oriented > 0.
    Np_dot_ker = (Np * ker_oriented).sum(axis=-1)
    sign_Np = np.where(Np_dot_ker >= 0.0, 1.0, -1.0)
    front_uv = Np * sign_Np[:, None]                    # (N, 2)

    return d_arr, front_uv


def _compute_d_at(uv: np.ndarray, surface: SurfaceParams,
                  projection: Projection) -> np.ndarray:
    """Scalar fallback of `_compute_d_batch` for one uv point (d only).

    Used by O4 refinement's per-Newton-iteration evaluations.  Delegates to
    the batched helper with N = 1 so the algorithm lives in one place.
    Returns just the `d` field; if you need `front_uv` too, call
    `_compute_d_batch` directly.
    """
    uv = np.asarray(uv, dtype=float).reshape(2)
    u_ = np.array([uv[0]])
    v_ = np.array([uv[1]])
    vals = surface._eval_all(u_, v_)
    d_arr, _front = _compute_d_batch(uv.reshape(1, 2), vals, projection)
    return d_arr[0]


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
                # f(s)  = SN · viewer_dir;  f'(s) = dSN/ds · viewer_dir
                # viewer_dir = eye - S (toward viewer); SN · dS/ds = 0 because
                # SN ⊥ Su, Sv. Sign-blind: Newton finds the same zero either way,
                # but kept consistent with per_vertex_viewer_dot.
                xyz_arr = vals[0].T  # (N, 3)
                vd = projection.eye - xyz_arr  # (N, 3)
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

    # 5. d field + front_uv: 2D curvature direction of the projected contour
    #    and the uv-space front-sheet vector (both vectorized in one pass that
    #    shares the SVD-derived `ker` between them).
    d_arr, front_uv_arr = _compute_d_batch(uv_arr, vals_final, projection)

    # 6. ptype: 4 on boundary edge (g == -1), else 0
    ptype_arr = np.where(cand["g"] == -1, np.uint8(4), np.uint8(0))

    result = np.zeros(len(cand_idx), dtype=cp_dtype)
    result["e"]        = cand_idx
    result["s"]        = s
    result["uv"]       = uv_arr
    result["xyz"]      = xyz_arr
    result["d"]        = d_arr
    result["front_uv"] = front_uv_arr
    result["ptype"]    = ptype_arr
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


# ── O4 ────────────────────────────────────────────────────────────────────────

vp_dtype = np.dtype([
    ("cs",         "i4"),   # CS index where the VP sits
    ("s",          "f8"),   # parametric position on the CS — kept at 0.5 (the
                            # chord midpoint between the two d-sign-change CPs);
                            # refinement adjusts uv off the chord but s is the
                            # initial bracket midpoint that O10 consumes.
    ("uv",         "2f8"),  # domain coord of the VP (refined if refine=True)
    ("xyz",        "3f8"),  # 3D position S(uv)
    ("vis_change", "i1"),   # ±1 — G10
])


def _cp_sequence(cc: ContourCurve, css: np.ndarray) -> list[int]:
    """Ordered CP indices along a CC chain.

    Open chain of N CSs → N+1 entries; closed → N+1 with seq[0] == seq[-1].
    Chain position j connects ``seq[j] → seq[j+1]`` via ``cs_indices[j]``
    (sign of cs_indices[j] tells which endpoint of the CS comes first).
    """
    seq: list[int] = []
    for j, si in enumerate(cc.cs_indices):
        cs = css[abs(int(si)) - 1]
        start_cp = int(cs["p_cp"]) if si > 0 else int(cs["q_cp"])
        end_cp   = int(cs["q_cp"]) if si > 0 else int(cs["p_cp"])
        if j == 0:
            seq.append(start_cp)
        seq.append(end_cp)
    return seq


def _newton_orthogonal_cp(
    uv_mid: np.ndarray, n_hat: np.ndarray,
    surface: SurfaceParams, projection: Projection, *,
    n_iter: int = 10,
) -> np.ndarray:
    """Run Newton from ``uv_mid`` along ``n_hat`` to find a CP (SN·viewer_dir = 0).

    f(t)  = SN(uv_mid + t·n_hat) · viewer_dir
    f'(t) = (dSN/du · n_hat[0] + dSN/dv · n_hat[1]) · viewer_dir
        with dSN/du = Suu×Sv + Su×Suv, dSN/dv = Suv×Sv + Su×Svv.
    For perspective viewer_dir = eye − S(uv) (toward viewer); SN · dS/dt = 0
    because SN ⊥ Su, Sv, so the extra term vanishes.

    No (0,1) clamping — the caller's bisection brackets the cusp.
    """
    t = 0.0
    for _ in range(n_iter):
        uv_t = uv_mid + t * n_hat
        u_t, v_t = float(uv_t[0]), float(uv_t[1])
        vals = surface._eval_all(u_t, v_t)
        Su   = np.asarray(vals[1], dtype=float).reshape(3)
        Sv   = np.asarray(vals[2], dtype=float).reshape(3)
        Suu  = np.asarray(vals[3], dtype=float).reshape(3)
        Suv  = np.asarray(vals[4], dtype=float).reshape(3)
        Svv  = np.asarray(vals[5], dtype=float).reshape(3)
        SN_  = np.asarray(vals[6], dtype=float).reshape(3)

        dSN_du = np.cross(Suu, Sv) + np.cross(Su, Suv)
        dSN_dv = np.cross(Suv, Sv) + np.cross(Su, Svv)
        dSN_dt = dSN_du * n_hat[0] + dSN_dv * n_hat[1]

        if projection.mode == "ortho":
            axis = projection._axis
            f_t  = float(SN_ @ axis)
            fp_t = float(dSN_dt @ axis)
        else:
            S_t = np.asarray(vals[0], dtype=float).reshape(3)
            vd  = projection.eye - S_t  # toward viewer; sign-blind for Newton
            f_t  = float(SN_ @ vd)
            fp_t = float(dSN_dt @ vd)

        if abs(fp_t) < 1e-15:
            break
        step = f_t / fp_t
        t -= step
        if abs(step) < 1e-14:
            break
    return uv_mid + t * n_hat


def _refine_cusp(
    p_uv: np.ndarray, q_uv: np.ndarray,
    surface: SurfaceParams, projection: Projection, *,
    max_iter: int = 20, tol: float = 1e-10,
) -> np.ndarray:
    """Newton-bisection refinement of a cusp bracketed by two CPs.

    Pre: d(p_uv)·d(q_uv) < 0 (caller verified the d-field sign change).

    Each iteration takes the midpoint M of [p, q], runs Newton orthogonally
    to (q-p) at M to land on a true CP M', then keeps whichever of [p, M']
    or [M', q] still brackets the d sign change.  Terminates on
    ``|q - p| < tol`` or after ``max_iter`` halvings.
    """
    p = p_uv.copy()
    q = q_uv.copy()
    d_p = _compute_d_at(p, surface, projection)
    d_q = _compute_d_at(q, surface, projection)

    for _ in range(max_iter):
        diff = q - p
        if float(np.linalg.norm(diff)) < tol:
            break
        mid = 0.5 * (p + q)
        perp = np.array([-diff[1], diff[0]])  # 90° rotation of (q-p) in 2D
        n_hat = perp / float(np.linalg.norm(perp))

        new_uv = _newton_orthogonal_cp(mid, n_hat, surface, projection)
        d_new  = _compute_d_at(new_uv, surface, projection)

        if float(d_p @ d_new) < 0.0:
            q, d_q = new_uv, d_new
        else:
            p, d_p = new_uv, d_new

    return 0.5 * (p + q)


def find_vps(
    ccs: list[ContourCurve],
    css: np.ndarray,
    cps: np.ndarray,
    surface: SurfaceParams,
    projection: Projection,
    *,
    refine: bool = False,
) -> np.ndarray:
    """Find cusp points (VPs) on contour curves.

    Walks each CC chain in order, flags a sign change of the d field between
    consecutive CPs (``d_i · d_j < 0``); each such CS hosts one VP.

    Note (2026-05-27): the `front_uv` field on CPs is also available and
    detects cusps via the analogous criterion (`front_i · front_j < 0`),
    but `d` (image-curvature) is more numerically robust in degenerate
    configurations such as the torus side view where the silhouette factors
    into perpendicular contour-line families in uv (the `front_uv` arrows at
    the two arms are then orthogonal, dot ≈ 0, noisy sign). `front_uv` is
    consumed by O7 BCP and O11 CDP instead.

    refine=False : VP location = midpoint of the CS in uv.
    refine=True  : Newton-bisection refinement orthogonal to the chord
                   (spec §"Cusp points (VPs)" — Newton refinement bullet).

    Visibility change (G10): ``+1 if dS(p, q-p) · axis > 0 else -1``, with
    dS evaluated at the *left* CP p (Reorganization §"Cusp points" line 373).
    axis = `projection.viewer_direction(xyz_vp)` — per-point viewer direction
    at the cusp's 3D location. In ortho this returns the constant principal
    axis (matching the original spec); in persp it returns `eye - xyz_vp`,
    which differs from the camera axis whenever the cusp is off the optical
    axis. Using the constant axis in persp flipped vis_change on cusps far
    from the principal direction, producing phantom BFS updates that LP had
    to repair.

    The dot product is computed in the CS's NATIVE direction (p_cp → q_cp),
    not the CC's chain direction, because downstream `_emit_subcurves` applies
    its own chain-sign correction (`sign_in * spt.vis_chge`) when emitting
    SubCurves — matching the convention used by BCP / BDP / CDP SPTs. Storing
    `vis_change` in chain direction here would double-correct on CSs that
    appear with a negative sign in the CC chain.
    """
    rows: list[tuple] = []

    for cc in ccs:
        seq = _cp_sequence(cc, css)
        # For a closed chain seq[-1] == seq[0]; drop the wrap-around pair.
        n_pairs = len(cc.cs_indices) - (1 if cc.is_closed else 0)

        for j in range(n_pairs):
            i_cp = seq[j]
            k_cp = seq[j + 1]

            d_i = cps[i_cp]["d"]
            d_k = cps[k_cp]["d"]
            if float(d_i @ d_k) >= 0.0:
                continue

            cs_idx = abs(int(cc.cs_indices[j])) - 1

            p_uv = cps[i_cp]["uv"].copy()
            q_uv = surface.domain.close(p_uv, cps[k_cp]["uv"])  # wrap-adjusted

            if refine:
                uv_vp = _refine_cusp(p_uv, q_uv, surface, projection)
            else:
                uv_vp = 0.5 * (p_uv + q_uv)

            u_, v_ = float(uv_vp[0]), float(uv_vp[1])
            xyz_vp = np.asarray(surface.S(u_, v_), dtype=float).reshape(3)

            # vis_change in NATIVE CS direction (p_cp → q_cp).
            cs = css[cs_idx]
            p_nat = int(cs["p_cp"])
            q_nat = int(cs["q_cp"])
            p_uv_nat = cps[p_nat]["uv"].copy()
            q_uv_nat = surface.domain.close(p_uv_nat, cps[q_nat]["uv"])
            dir_uv_nat = q_uv_nat - p_uv_nat
            up_, vp_ = float(p_uv_nat[0]), float(p_uv_nat[1])
            Su_p = np.asarray(surface.Su(up_, vp_), dtype=float).reshape(3)
            Sv_p = np.asarray(surface.Sv(up_, vp_), dtype=float).reshape(3)
            dS_3d = Su_p * dir_uv_nat[0] + Sv_p * dir_uv_nat[1]
            axis_at_vp = projection.viewer_direction(xyz_vp).reshape(3)
            vis = np.int8(1) if float(dS_3d @ axis_at_vp) > 0.0 else np.int8(-1)

            rows.append((cs_idx, 0.5, uv_vp, xyz_vp, vis))

    if not rows:
        return np.zeros(0, dtype=vp_dtype)

    out = np.zeros(len(rows), dtype=vp_dtype)
    for k, (cs_idx, s_vp, uv_vp, xyz_vp, vis) in enumerate(rows):
        out[k]["cs"]         = cs_idx
        out[k]["s"]          = s_vp
        out[k]["uv"]         = uv_vp
        out[k]["xyz"]        = xyz_vp
        out[k]["vis_change"] = vis
    return out

