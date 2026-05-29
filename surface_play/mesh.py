import math
from dataclasses import dataclass

import numpy as np
import triangle

from surface_play.domain import Domain
from surface_play.surface import SurfaceParams


edge_dtype = np.dtype([
    ("p_idx", "i4"), ("q_idx", "i4"),
    ("p", "2f8"),
    ("pq", "2f8"),
    ("f", "i4"), ("g", "i4"),
    ("dir", "2f8"),
    ("flip", "i1"),
    ("split1", "i4"), ("split2", "i4"),
])

face_dtype = np.dtype([
    ("verts", "3i4"),
    ("edges", "3i4"),
    ("p", "2f8"),
    ("pq", "2f8"),
    ("pr", "2f8"),
])


def _boundary_edge_count(tris: np.ndarray) -> int:
    """Count edges that appear in exactly one triangle (boundary edges)."""
    all_edges = np.sort(
        np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1
    )
    _, counts = np.unique(all_edges, axis=0, return_counts=True)
    return int((counts == 1).sum())


def _generate_rect_mesh(domain: Domain, resolution: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns:
      uv:   (N, 2) float64 — vertex (u, v) coordinates.
      tris: (M, 3) int32   — triangle vertex indices.
    No identifications applied yet (paired vertices have distinct indices). No jitter.
    """
    u_min, u_max, v_min, v_max = domain.bounds
    du = (u_max - u_min) / resolution
    dv = (v_max - v_min) / resolution
    n_bv = 4 * resolution

    bverts = np.empty((n_bv, 2), dtype=np.float64)
    # Bottom: u from u_min..u_max-du, v=v_min
    for i in range(resolution):
        bverts[i] = [u_min + i * du, v_min]
    # Right: u=u_max, v from v_min..v_max-dv
    for i in range(resolution):
        bverts[resolution + i] = [u_max, v_min + i * dv]
    # Top: u from u_max..u_min+du, v=v_max
    for i in range(resolution):
        bverts[2 * resolution + i] = [u_max - i * du, v_max]
    # Left: u=u_min, v from v_max..v_min+dv
    for i in range(resolution):
        bverts[3 * resolution + i] = [u_min, v_max - i * dv]

    idx = np.arange(n_bv, dtype=np.int32)
    bsegs = np.column_stack([idx, np.roll(idx, -1)])

    rect_area = (u_max - u_min) * (v_max - v_min)
    area = rect_area / resolution ** 2

    pslg = {"vertices": bverts, "segments": bsegs}
    # Use fixed-point format (NOT .17g): for small areas, .17g switches to
    # scientific notation (e.g. "7.85e-05") and the C triangle library's
    # option parser treats 'e' as the -e flag, truncating the area number
    # and silently dropping the area constraint. Bites disk domains with
    # r_max=1 at resolution=200 (area=7.85e-5 < 1e-4 threshold).
    result = triangle.triangulate(pslg, opts=f"pYq30a{area:.20f}")

    uv = np.asarray(result["vertices"], dtype=np.float64)
    tris = np.asarray(result["triangles"], dtype=np.int32)
    return uv, tris


def _build_edges_faces(
    uv: np.ndarray,
    tris: np.ndarray,
    uv_pre: np.ndarray,
    tris_pre: np.ndarray,
    SN_per_vertex: np.ndarray,
    on_u_seam: np.ndarray,
    on_v_seam: np.ndarray,
    domain: Domain,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build edges and faces from compacted `uv`/`tris` (one row per equivalence class
    in `uv`, canonical indices in `tris`).  Edges are keyed by their canonical
    vertex pair, so seam-paired edges (whose two physical sides relabel to the
    same canonical pair) merge automatically into a single edge with two
    adjacent faces.

    `p`, `q_idx`, etc. use compacted indices.  `pq` / `pr` are pre-identification
    subtractions: per roadmap line 552, `pq = uv_pre[bp] − uv_pre[ap]` where
    `(ap, bp)` are the pre-id labels of the canonical pair from one incident
    face's row of `tris_pre`.  This keeps the vector representing the *physical*
    edge (short way) even when the canonical endpoints sit on opposite sides of
    a cy/mo seam — without that, `p + s·pq` would interpolate across the
    "long way" and surface evaluation at a periodic alias of the true edge
    point breaks every downstream consumer that cares about which period the
    point lives in (notably O1 contour points on multi-branch silhouettes).

    Flip detection: an edge has `flip = -1` iff its canonical-uv chord
    crosses a Möbius seam (chord component on a `mo` axis exceeds half the
    period). All other edges have `flip = +1`.
    """
    M = len(tris)

    u_is_mo = domain.type == "rect" and domain.u_identify == "mo"
    v_is_mo = domain.type == "rect" and domain.v_identify == "mo"
    if domain.type == "rect":
        _u_min_d, _u_max_d, _v_min_d, _v_max_d = domain.bounds
        _u_range_half = 0.5 * (_u_max_d - _u_min_d)
        _v_range_half = 0.5 * (_v_max_d - _v_min_d)
    else:
        _u_range_half = _v_range_half = float("inf")

    # records[key] : list of (f_idx, third_compacted, a_compacted, b_compacted, pre_pair)
    # where pre_pair = (pre_p, pre_q) ordered to match key = (min, max) on canonical labels.
    records: dict[tuple[int, int], list[tuple[int, int, int, int, tuple[int, int]]]] = {}
    face_edge_idx = np.empty((M, 3), dtype=np.int32)

    for f_idx in range(M):
        i  = int(tris[f_idx, 0]);  ip = int(tris_pre[f_idx, 0])
        j  = int(tris[f_idx, 1]);  jp = int(tris_pre[f_idx, 1])
        k  = int(tris[f_idx, 2]);  kp = int(tris_pre[f_idx, 2])
        local_edges = (
            (i, j, k, ip, jp),
            (j, k, i, jp, kp),
            (k, i, j, kp, ip),
        )
        for a, b, third, ap, bp in local_edges:
            if a < b:
                key = (a, b)
                pre_pair = (ap, bp)
            else:
                key = (b, a)
                pre_pair = (bp, ap)
            rec_list = records.get(key)
            if rec_list is None:
                rec_list = []
                records[key] = rec_list
            rec_list.append((f_idx, third, a, b, pre_pair))

    key_to_idx = {key: idx for idx, key in enumerate(records.keys())}

    for f_idx in range(M):
        i = int(tris[f_idx, 0])
        j = int(tris[f_idx, 1])
        k = int(tris[f_idx, 2])
        local_edges = ((i, j, k), (j, k, i), (k, i, j))
        for local_pos, (a, b, _third) in enumerate(local_edges):
            key = (a, b) if a < b else (b, a)
            face_edge_idx[f_idx, local_pos] = key_to_idx[key]

    n_edges = len(records)
    edges = np.zeros(n_edges, dtype=edge_dtype)

    for key, recs in records.items():
        e_idx = key_to_idx[key]
        if len(recs) > 2:
            raise ValueError(
                f"Edge {key} has {len(recs)} face references (expected 1 or 2). "
                f"Records: {recs}"
            )

        f_idx, third, a, b, pre_pair = recs[0]
        p_idx, q_idx = key  # already sorted (compacted)
        pre_p, pre_q = pre_pair
        p  = uv[p_idx]
        pq = uv_pre[pre_q] - uv_pre[pre_p]   # pre-id subtraction (G15-compliant)

        # mo-seam axis-flip correction for `p + s·pq` interpolation.
        # When canonical p and pre-id pre_p live in different identification
        # copies (canonical chose the OTHER side of the seam as representative),
        # the pre-id chord pq is expressed in pre_p's copy. To keep
        # `p + s·pq` evaluating to the geometrically correct 3D point in
        # canonical p's copy, we must flip the OTHER axis under mo:
        #   - mo on u (identification (0, v) ~ (u_max, -v)): u-seam crossing
        #     reverses v across copies → flip pq[1].
        #   - mo on v: analogous with axes swapped → flip pq[0].
        # cy identifications don't reverse any axis, so no flip.
        diff_p = p - uv_pre[pre_p]
        if u_is_mo and abs(float(diff_p[0])) > _u_range_half:
            pq[1] = -pq[1]
        elif v_is_mo and abs(float(diff_p[1])) > _v_range_half:
            pq[0] = -pq[0]

        edges[e_idx]["p_idx"] = p_idx
        edges[e_idx]["q_idx"] = q_idx
        edges[e_idx]["p"] = p
        edges[e_idx]["pq"] = pq

        if len(recs) == 1:
            edges[e_idx]["f"] = f_idx
            edges[e_idx]["g"] = -1
            d = np.array([-pq[1], pq[0]], dtype=np.float64)
            offset = uv[third] - p  # boundary triangle: never spans a seam
            if float(np.dot(d, offset)) < 0.0:
                d = -d
            nrm = float(np.linalg.norm(d))
            if nrm > 0.0:
                d = d / nrm
            edges[e_idx]["dir"] = d
        else:
            edges[e_idx]["f"] = f_idx
            edges[e_idx]["g"] = recs[1][0]
            edges[e_idx]["dir"] = (0.0, 0.0)

        # flip = -1 iff the canonical-uv chord straddles a Möbius seam
        # (i.e. its periodic-axis component is longer than half the range).
        # The `on_u_seam[p] and on_u_seam[q]` test was too narrow: across
        # the mo-u seam, the triangle on one side has its third vertex in
        # the interior, and the edge between that interior vertex (canonical
        # close to u=0) and an opposite-side seam vertex (canonical at
        # u=u_max) is a genuine mo-flip edge but only one endpoint is on
        # the seam. A wider `SN·SN<0` heuristic (pre-9f02eab) spuriously
        # fired on highly-curved smooth surfaces (Vagues at res=200, fig8
        # near a fold) — the chord-jump test is the precise replacement.
        if u_is_mo and abs(float(uv[q_idx, 0] - uv[p_idx, 0])) > _u_range_half:
            crosses_mo = True
        elif v_is_mo and abs(float(uv[q_idx, 1] - uv[p_idx, 1])) > _v_range_half:
            crosses_mo = True
        else:
            crosses_mo = False
        edges[e_idx]["flip"] = -1 if crosses_mo else 1

        edges[e_idx]["split1"] = -1
        edges[e_idx]["split2"] = -1

    faces = np.zeros(M, dtype=face_dtype)
    for f_idx in range(M):
        i  = int(tris[f_idx, 0]);  ip = int(tris_pre[f_idx, 0])
        j  = int(tris[f_idx, 1]);  jp = int(tris_pre[f_idx, 1])
        k  = int(tris[f_idx, 2]);  kp = int(tris_pre[f_idx, 2])
        faces[f_idx]["verts"] = (i, j, k)
        faces[f_idx]["edges"] = face_edge_idx[f_idx]
        faces[f_idx]["p"]  = uv[i]
        faces[f_idx]["pq"] = uv_pre[jp] - uv_pre[ip]   # pre-id subtraction
        faces[f_idx]["pr"] = uv_pre[kp] - uv_pre[ip]   # pre-id subtraction

    return edges, faces


def _jitter(
    uv: np.ndarray,
    tris: np.ndarray,
    domain: Domain,
    *,
    seed: int | None = None,
) -> np.ndarray:
    """
    Returns a jittered copy of uv. Independent jitter per vertex by ±0.1% of typical
    edge length. Reprojects true boundary vertices onto the domain boundary. Run BEFORE
    `_apply_identifications` (C4): paired vertices receive independent jitter — the
    discarded one is dropped at identification, so coherence is unnecessary.
    RNG seeded by `seed` (default = id(uv) for repro).
    """
    if seed is None:
        seed = id(uv)

    e0 = np.linalg.norm(uv[tris[:, 1]] - uv[tris[:, 0]], axis=1)
    e1 = np.linalg.norm(uv[tris[:, 2]] - uv[tris[:, 1]], axis=1)
    e2 = np.linalg.norm(uv[tris[:, 0]] - uv[tris[:, 2]], axis=1)
    L = float(np.median(np.concatenate([e0, e1, e2])))

    rng = np.random.default_rng(seed)
    delta = rng.uniform(-1.0, 1.0, uv.shape) * 0.001 * L
    uv_new = uv + delta

    tol = 1e-12
    if domain.type == "rect":
        u_min, u_max, v_min, v_max = domain.bounds
        if domain.u_identify == "no":
            uv_new[np.abs(uv[:, 0] - u_min) < tol, 0] = u_min
            uv_new[np.abs(uv[:, 0] - u_max) < tol, 0] = u_max
        if domain.v_identify == "no":
            uv_new[np.abs(uv[:, 1] - v_min) < tol, 1] = v_min
            uv_new[np.abs(uv[:, 1] - v_max) < tol, 1] = v_max
    elif domain.type in ("disk", "annulus"):
        r_min, r_max = domain.bounds[0], domain.bounds[1]
        r = np.sqrt(uv[:, 0] ** 2 + uv[:, 1] ** 2)
        ring_radii = [r_max] + ([r_min] if r_min > 0 else [])
        for target_r in ring_radii:
            mask = np.abs(r - target_r) < tol
            r_new = np.sqrt(uv_new[mask, 0] ** 2 + uv_new[mask, 1] ** 2)
            uv_new[mask] = uv_new[mask] / r_new[:, None] * target_r

    return uv_new


def _apply_identifications(
    uv_jittered: np.ndarray,
    tris: np.ndarray,
    domain: Domain,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Identify seam-paired vertices, compact `uv` to one row per equivalence class,
    and rewrite `tris` to those compacted indices. Non-canonical vertex copies are
    removed outright.

    (mo, mo) identification is rejected: a smooth immersion cannot be parametrized
    by a square with both Möbius identifications. RP² and similar non-orientables
    are handled later via disk-based parametrizations.

    Returns (uv_compacted, tris_canonical, on_u_seam, on_v_seam, corner_idx):
      - uv_compacted: (K, 2) — one row per equivalence class, in canonical order.
      - tris_canonical: (M, 3) — original tris rewritten to compacted indices.
      - on_u_seam: (K,) bool — True if the canonical class includes any vertex on
        an identified u-side (u_min or u_max). Used by C5 for mo-seam flip detection.
      - on_v_seam: (K,) bool — analogous for v.
      - corner_idx: (≤4,) int32 — compacted indices of the original four rect corners,
        deduplicated via identification. Empty for non-rect.
    """
    N = len(uv_jittered)

    if domain.type != "rect":
        return (
            uv_jittered.copy(),
            tris.astype(np.int32, copy=True),
            np.zeros(N, dtype=bool),
            np.zeros(N, dtype=bool),
            np.array([], dtype=np.int32),
        )

    n = _boundary_edge_count(tris) // 4
    pre_corners = np.array([0, n, 2 * n, 3 * n], dtype=np.int32) if n > 0 else \
                  np.array([], dtype=np.int32)

    if domain.u_identify == "no" and domain.v_identify == "no":
        return (
            uv_jittered.copy(),
            tris.astype(np.int32, copy=True),
            np.zeros(N, dtype=bool),
            np.zeros(N, dtype=bool),
            pre_corners,
        )

    if domain.u_identify == "mo" and domain.v_identify == "mo":
        raise ValueError(
            "(mo, mo) identification is not supported: a smooth immersion cannot "
            "be parametrized this way. Use a disk-based parametrization for RP²."
        )

    pairs: list[tuple[int, int]] = []

    if domain.u_identify == "cy":
        pairs.append((0, n))
        for i in range(n):
            pairs.append((3 * n + i, 2 * n - i))
    elif domain.u_identify == "mo":
        pairs.append((0, 2 * n))
        for i in range(n):
            pairs.append((3 * n + i, n + i))

    if domain.v_identify == "cy":
        pairs.append((n, 2 * n))
        for i in range(n):
            pairs.append((i, 3 * n - i))
    elif domain.v_identify == "mo":
        pairs.append((n, 3 * n))
        for i in range(n):
            pairs.append((i, 2 * n + i))

    parent = np.arange(N, dtype=np.int32)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            if ra > rb:
                ra, rb = rb, ra
            parent[rb] = ra

    canonical_old = np.array([find(i) for i in range(N)], dtype=np.int32)

    # Per-vertex IDENTIFIED-seam membership flags in PRE-ID indexing.
    # Only set for axes whose identification is cy or mo; unidentified-axis
    # sides are plain boundary, not seams.
    pre_on_u = np.zeros(N, dtype=bool)
    pre_on_v = np.zeros(N, dtype=bool)
    u_seam = domain.u_identify in ("cy", "mo")
    v_seam = domain.v_identify in ("cy", "mo")
    if n > 0:
        if v_seam:
            pre_on_v[0:n] = True            # v_min interior side
            pre_on_v[2 * n:3 * n] = True    # v_max interior side
        if u_seam:
            pre_on_u[n:2 * n] = True        # u_max interior side
            pre_on_u[3 * n:4 * n] = True    # u_min interior side
        # Corners belong to both adjacent sides — set per-axis as identified.
        for corner in (0, n, 2 * n, 3 * n):
            if u_seam:
                pre_on_u[corner] = True
            if v_seam:
                pre_on_v[corner] = True

    # Map each canonical to a new compacted index.
    canonical_indices = np.unique(canonical_old)
    K = len(canonical_indices)
    old_to_new = np.full(N, -1, dtype=np.int32)
    old_to_new[canonical_indices] = np.arange(K, dtype=np.int32)
    # Members of a class inherit the new index of their canonical.
    new_idx = old_to_new[canonical_old]

    uv_compacted = uv_jittered[canonical_indices]
    tris_canonical = new_idx[tris].astype(np.int32)

    # Seam membership of each compacted class: True if ANY class member was on
    # that seam. Computed by OR-reducing per pre-id membership into classes.
    on_u_seam = np.zeros(K, dtype=bool)
    on_v_seam = np.zeros(K, dtype=bool)
    np.logical_or.at(on_u_seam, new_idx, pre_on_u)
    np.logical_or.at(on_v_seam, new_idx, pre_on_v)

    corner_idx = np.unique(new_idx[pre_corners]).astype(np.int32)

    return uv_compacted, tris_canonical, on_u_seam, on_v_seam, corner_idx


@dataclass
class Mesh:
    domain: Domain
    surface: SurfaceParams
    uv: np.ndarray           # (K, 2) — compacted to one row per equivalence class
    tris: np.ndarray         # (M, 3) — compacted indices into uv
    edges: np.ndarray        # structured array, edge_dtype
    faces: np.ndarray        # structured array, face_dtype
    SN: np.ndarray           # (K, 3)
    xyz: np.ndarray          # (K, 3)
    boundary_edge_idx: np.ndarray  # indices into edges where g == -1
    corner_idx: np.ndarray         # corner indices in compacted uv (rect only)


def build_mesh(
    domain: Domain,
    surface: SurfaceParams,
    resolution: int,
    *,
    jitter: bool = True,
    seed: int | None = None,
) -> Mesh:
    if domain.type == "rect":
        uv_raw, tris_raw = _generate_rect_mesh(domain, resolution)
    else:
        uv_raw, tris_raw = _generate_disk_mesh(domain, resolution)

    uv_jittered = _jitter(uv_raw, tris_raw, domain, seed=seed) if jitter else uv_raw

    uv, tris, on_u_seam, on_v_seam, corner_idx = _apply_identifications(
        uv_jittered, tris_raw, domain
    )

    xyz = surface.S(uv[:, 0], uv[:, 1]).T
    SN = surface.SN(uv[:, 0], uv[:, 1]).T

    edges, faces = _build_edges_faces(
        uv, tris, uv_jittered, tris_raw, SN, on_u_seam, on_v_seam, domain
    )

    boundary_edge_idx = np.nonzero(edges["g"] == -1)[0]

    return Mesh(
        domain=domain,
        surface=surface,
        uv=uv,
        tris=tris,
        edges=edges,
        faces=faces,
        SN=SN,
        xyz=xyz,
        boundary_edge_idx=boundary_edge_idx,
        corner_idx=corner_idx,
    )


def _generate_disk_mesh(domain: Domain, resolution: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Disk (r_min == 0) or annulus (r_min > 0). Returns (uv, tris) in cartesian (u, v).
    """
    r_min, r_max = domain.bounds[0], domain.bounds[1]
    bbox_diag = domain.bbox_diag()  # = 2 * r_max

    n_outer = max(3, round(2 * math.pi * r_max * resolution / bbox_diag))
    theta_outer = np.linspace(0, 2 * math.pi, n_outer, endpoint=False)
    outer_verts = np.column_stack(
        [r_max * np.cos(theta_outer), r_max * np.sin(theta_outer)]
    )
    idx_o = np.arange(n_outer, dtype=np.int32)
    outer_segs = np.column_stack([idx_o, np.roll(idx_o, -1)])

    if r_min > 0:
        n_inner = max(3, round(2 * math.pi * r_min * resolution / bbox_diag))
        theta_inner = np.linspace(0, 2 * math.pi, n_inner, endpoint=False)
        inner_verts = np.column_stack(
            [r_min * np.cos(theta_inner), r_min * np.sin(theta_inner)]
        )
        idx_i = n_outer + np.arange(n_inner, dtype=np.int32)
        inner_segs = np.column_stack([idx_i, np.roll(idx_i, -1)])

        all_verts = np.vstack([outer_verts, inner_verts])
        all_segs = np.vstack([outer_segs, inner_segs])
        domain_area = math.pi * (r_max ** 2 - r_min ** 2)
        pslg = {
            "vertices": all_verts,
            "segments": all_segs,
            "holes": np.array([[0.0, 0.0]]),
        }
    else:
        all_verts = outer_verts
        all_segs = outer_segs
        domain_area = math.pi * r_max ** 2
        pslg = {"vertices": all_verts, "segments": all_segs}

    area = domain_area / resolution ** 2
    result = triangle.triangulate(pslg, opts=f"pYq30a{area:.20f}")

    uv = np.asarray(result["vertices"], dtype=np.float64)
    tris = np.asarray(result["triangles"], dtype=np.int32)
    return uv, tris
