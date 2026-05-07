import math

import numpy as np
import triangle

from surface_play.domain import Domain


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
