import math

import numpy as np
import triangle

from surface_play.domain import Domain


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
    result = triangle.triangulate(pslg, opts=f"pYq30a{area:.17g}")

    uv = np.asarray(result["vertices"], dtype=np.float64)
    tris = np.asarray(result["triangles"], dtype=np.int32)
    return uv, tris


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
) -> tuple[np.ndarray, np.ndarray]:
    """
    For rect with cy/mo on either axis: for each pair of identified boundary vertices,
    pick one as canonical (the u_min/v_min-side vertex) and replace occurrences of the
    other in tris. The disposed vertex's row in uv is left in place (unused) — caller
    may compact later. Returns (uv_unchanged, tris_with_merged_indices).
    Vertex pairing matches by the geometric criterion at PRE-jitter positions: cy uses
    same v (or u) coord; mo uses mirrored. (Pairing is determined from C1's deterministic
    boundary placement, not from current jittered uvs.)
    """
    if domain.type != "rect" or (domain.u_identify == "no" and domain.v_identify == "no"):
        return uv_jittered, tris.copy()

    # Infer n = resolution from boundary edge count (= 4*n for a rect PSLG with -Y).
    n = _boundary_edge_count(tris) // 4
    N = len(uv_jittered)

    # C1 boundary layout (indices):
    #   bottom: 0..n-1  at (u_min + k*du, v_min)
    #   right:  n..2n-1 at (u_max, v_min + k*dv)
    #   top:   2n..3n-1 at (u_max - k*du, v_max)
    #   left:  3n..4n-1 at (u_min, v_max - k*dv)
    pairs: list[tuple[int, int]] = []

    if domain.u_identify == "cy":
        # (u_min, v) pairs with (u_max, same-v): (0,n) and (3n+i, 2n-i)
        pairs.append((0, n))
        for i in range(n):
            pairs.append((3 * n + i, 2 * n - i))
    elif domain.u_identify == "mo":
        # (u_min, v) pairs with (u_max, mirrored-v): (0,2n) and (3n+i, n+i)
        pairs.append((0, 2 * n))
        for i in range(n):
            pairs.append((3 * n + i, n + i))

    if domain.v_identify == "cy":
        # (u, v_min) pairs with (same-u, v_max): (n,2n) and (i, 3n-i)
        pairs.append((n, 2 * n))
        for i in range(n):
            pairs.append((i, 3 * n - i))
    elif domain.v_identify == "mo":
        # (u, v_min) pairs with (mirrored-u, v_max): (n,3n) and (i, 2n+i)
        pairs.append((n, 3 * n))
        for i in range(n):
            pairs.append((i, 2 * n + i))

    # Union-find with path halving; smaller index wins as canonical.
    # Sequential assignment would create cycles for mo-mo (e.g. n->3n->n),
    # so union-find is required to resolve all equivalence classes correctly.
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

    remap = np.array([find(i) for i in range(N)], dtype=np.int32)
    new_tris = remap[tris].astype(np.int32)

    # Remove degenerate faces (two vertices collapsed to the same canonical) and
    # duplicate faces (two pre-id triangles that map to the same triple of vertices,
    # which happens at corners under mo-mo identification).
    valid = np.array([len(set(t)) == 3 for t in new_tris])
    new_tris = new_tris[valid]
    _, first_occurrence = np.unique(np.sort(new_tris, axis=1), axis=0, return_index=True)
    new_tris = new_tris[np.sort(first_occurrence)]

    return uv_jittered, new_tris


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
    result = triangle.triangulate(pslg, opts=f"pYq30a{area:.17g}")

    uv = np.asarray(result["vertices"], dtype=np.float64)
    tris = np.asarray(result["triangles"], dtype=np.int32)
    return uv, tris
