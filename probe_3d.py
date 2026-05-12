"""C14 — 3D render of the surface, DPs (scatter), and SICs (polylines).

Two viewing angles, side by side. Uses matplotlib 3D (no Three.js).
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: F401, E402

from surface_play.intersections import (  # noqa: E402
    build_sics, build_sis_pairs, find_double_points,
)
from surface_play.mesh import build_mesh  # noqa: E402
from surface_play.test_fixtures import fig8  # noqa: E402

OUT = Path(__file__).resolve().parent / "probes_out"


def _sic_polyline_3d(sic, sis_pairs, dps):
    """3D polyline of a SIC: chain of DP `xyz` values in traversal order."""
    pts = []
    for tok in sic.sis_indices:
        sis_idx = abs(int(tok)) - 1
        p = int(sis_pairs[sis_idx]["p_dp"])
        q = int(sis_pairs[sis_idx]["q_dp"])
        if int(tok) < 0:
            p, q = q, p
        if not pts:
            pts.append(dps["xyz"][p])
        pts.append(dps["xyz"][q])
    return np.asarray(pts)


def _render_axes(ax, mesh, dps, sics, sis_pairs, view):
    xyz = mesh.xyz
    tris = mesh.tris
    tri_verts = xyz[tris]   # (M, 3, 3)
    coll = ax.add_collection3d(
        plt.matplotlib.collections.PolyCollection([])
    )
    coll.remove()

    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    pc = Poly3DCollection(tri_verts, alpha=0.15, facecolor="C0",
                          edgecolor="lightgray", linewidth=0.1)
    ax.add_collection3d(pc)

    if len(dps):
        ax.scatter(dps["xyz"][:, 0], dps["xyz"][:, 1], dps["xyz"][:, 2],
                   color="red", s=10, depthshade=False, label=f"{len(dps)} DPs")
    for i, sic in enumerate(sics):
        pl = _sic_polyline_3d(sic, sis_pairs, dps)
        ax.plot(pl[:, 0], pl[:, 1], pl[:, 2],
                color=f"C{(i + 1) % 10}", lw=2, label=f"SIC{i}")

    mins = xyz.min(axis=0); maxs = xyz.max(axis=0)
    ax.set_xlim(mins[0], maxs[0])
    ax.set_ylim(mins[1], maxs[1])
    ax.set_zlim(mins[2], maxs[2])
    ax.view_init(elev=view[0], azim=view[1])
    ax.set_title(f"elev={view[0]}, azim={view[1]}")
    ax.legend(fontsize=7)


def main():
    surface = fig8()
    mesh = build_mesh(surface.domain, surface, resolution=20)
    dps = find_double_points(mesh, surface)
    sis_pairs = build_sis_pairs(dps)
    sics = build_sics(sis_pairs)

    fig = plt.figure(figsize=(14, 7))
    for k, view in enumerate([(30, -60), (60, 30)]):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        _render_axes(ax, mesh, dps, sics, sis_pairs, view)

    OUT.mkdir(exist_ok=True)
    path = OUT / "probe_3d_fig8.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}  ({len(sics)} SIC, {len(dps)} DPs)")


if __name__ == "__main__":
    main()
