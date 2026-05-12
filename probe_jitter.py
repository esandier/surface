"""C14 — visual sanity for the `_jitter` step: side-by-side pre/post.

Boundary reprojection (true rect outer sides, disk r=r_min/r=r_max ring)
should keep boundary vertices on the boundary after jitter. Inspect the
edges of each pair of plots for that.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from surface_play.mesh import build_mesh
from surface_play.test_fixtures import (
    disk_paraboloid_po, helicoid, mobius_u,
)

OUT = Path(__file__).resolve().parent / "probes_out"


def _render(name, surface, resolution=12, seed=42):
    nj = build_mesh(surface.domain, surface, resolution=resolution, jitter=False)
    j = build_mesh(surface.domain, surface, resolution=resolution,
                   jitter=True, seed=seed)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 6))
    a1.triplot(nj.uv[:, 0], nj.uv[:, 1], nj.tris, color="C0", lw=0.5)
    a1.set_title(f"{name}  no jitter")
    a1.set_aspect("equal")
    a2.triplot(j.uv[:, 0], j.uv[:, 1], j.tris, color="C1", lw=0.5)
    a2.set_title(f"{name}  jitter seed={seed}")
    a2.set_aspect("equal")

    OUT.mkdir(exist_ok=True)
    path = OUT / f"probe_jitter_{name}.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


def main():
    _render("rect_no",   helicoid())
    _render("rect_mo_u", mobius_u())
    _render("disk_po",   disk_paraboloid_po())


if __name__ == "__main__":
    main()
