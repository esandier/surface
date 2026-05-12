"""C14 — boundary curves over the (u, v) mesh, color-coded by curve index.

Walks each BC's signed edge tokens, lifts (p, q) endpoints by `domain.close`
to avoid wrap-induced visual jumps on identified rectangles.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from surface_play.curves import build_bcs
from surface_play.mesh import build_mesh
from surface_play.test_fixtures import (
    disk_paraboloid_po, helicoid, mobius_u,
)

OUT = Path(__file__).resolve().parent / "probes_out"


def _polyline(mesh, bc):
    """Reconstruct the BC polyline by walking edges with sign-encoded reversal.

    `bc.edge_indices` is signed 1-indexed into `mesh.boundary_edge_idx`
    (the subset passed to `make_lines`), NOT into `mesh.edges` directly.
    """
    pts = []
    domain = mesh.domain
    for s in bc.edge_indices:
        bnd_pos = abs(int(s)) - 1
        idx = int(mesh.boundary_edge_idx[bnd_pos])
        e = mesh.edges[idx]
        p_uv = mesh.uv[int(e["p_idx"])]
        q_uv = mesh.uv[int(e["q_idx"])]
        q_local = domain.close(p_uv, q_uv) if hasattr(domain, "close") else q_uv
        if int(s) > 0:
            pts.append(p_uv)
            pts.append(q_local)
        else:
            pts.append(q_local)
            pts.append(p_uv)
    return np.asarray(pts)


def _render(name, surface, resolution=15):
    mesh = build_mesh(surface.domain, surface, resolution=resolution, jitter=False)
    bcs = build_bcs(mesh)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.triplot(mesh.uv[:, 0], mesh.uv[:, 1], mesh.tris,
               color="lightgray", lw=0.3)
    for i, bc in enumerate(bcs):
        pl = _polyline(mesh, bc)
        ax.plot(pl[:, 0], pl[:, 1], color=f"C{i % 10}", lw=2,
                label=f"BC{i} ({'closed' if bc.is_closed else 'open'})")
    ax.set_title(f"{name}  BCs: {len(bcs)}")
    ax.set_aspect("equal")
    if bcs:
        ax.legend(fontsize=8)

    OUT.mkdir(exist_ok=True)
    path = OUT / f"probe_bcs_{name}.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}  ({len(bcs)} BCs)")


def main():
    _render("rect_no",   helicoid())
    _render("rect_mo_u", mobius_u())
    _render("disk_po",   disk_paraboloid_po())


if __name__ == "__main__":
    main()
