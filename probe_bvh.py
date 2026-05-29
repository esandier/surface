"""C14 — edge/face broad phase visualization in 2D (XY projection).

Renders edge AABBs (blue), face AABBs (gray), and candidate non-adjacent
overlapping pairs as segments connecting the AABB centroids (red). Uses the
fused-grid broad phase (`edge_face_candidates`), so shared-vertex adjacency
pairs are already excluded — this probe is for visual sanity on the
*broad-phase output*, not any internal acceleration-structure layout.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from surface_play.intersections import edge_face_candidates
from surface_play.mesh import build_mesh
from surface_play.test_fixtures import fig8

OUT = Path(__file__).resolve().parent / "probes_out"


def main():
    surface = fig8()
    mesh = build_mesh(surface.domain, surface, resolution=15)

    xyz = mesh.xyz
    edges = mesh.edges
    tris = mesh.tris

    e_pts_p = xyz[edges["p_idx"]]
    e_pts_q = xyz[edges["q_idx"]]
    e_bbox = np.hstack([np.minimum(e_pts_p, e_pts_q),
                        np.maximum(e_pts_p, e_pts_q)])
    f_pts = xyz[tris]
    f_bbox = np.hstack([f_pts.min(axis=1), f_pts.max(axis=1)])

    ei, fi = edge_face_candidates(
        e_bbox, f_bbox, edges["p_idx"], edges["q_idx"], tris
    )

    fig, ax = plt.subplots(figsize=(9, 9))
    for box in f_bbox:
        ax.add_patch(mpatches.Rectangle(
            (box[0], box[1]), box[3] - box[0], box[4] - box[1],
            fill=False, edgecolor="lightgray", lw=0.3,
        ))
    for box in e_bbox:
        ax.add_patch(mpatches.Rectangle(
            (box[0], box[1]), box[3] - box[0], box[4] - box[1],
            fill=False, edgecolor="C0", lw=0.3, alpha=0.5,
        ))

    if len(ei):
        e_cen = 0.5 * (e_bbox[:, :3] + e_bbox[:, 3:])
        f_cen = 0.5 * (f_bbox[:, :3] + f_bbox[:, 3:])
        for k in range(min(len(ei), 2000)):
            p = e_cen[ei[k]]
            q = f_cen[fi[k]]
            ax.plot([p[0], q[0]], [p[1], q[1]], color="red", lw=0.4, alpha=0.4)

    ax.set_title(
        f"fig-8 res=15  AABBs (E={len(edges)}, F={len(tris)})  "
        f"candidate pairs = {len(ei)}"
    )
    ax.set_xlabel("X"); ax.set_ylabel("Y")
    ax.set_aspect("equal")

    OUT.mkdir(exist_ok=True)
    path = OUT / "probe_bvh_fig8_xy.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}  ({len(ei)} candidate pairs)")


if __name__ == "__main__":
    main()
