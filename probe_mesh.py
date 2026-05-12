"""C14 — visual sanity for `build_mesh` in (u, v) space.

Renders the post-identification mesh for representative domain types and
saves a PNG per fixture under probes_out/. No Django.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from surface_play.mesh import build_mesh
from surface_play.test_fixtures import (
    cylinder_cy, disk_paraboloid_po, fig8, helicoid, mobius_u,
)

OUT = Path(__file__).resolve().parent / "probes_out"


def _render(name, surface, resolution=15):
    mesh = build_mesh(surface.domain, surface, resolution=resolution, jitter=False)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.triplot(mesh.uv[:, 0], mesh.uv[:, 1], mesh.tris, color="C0", lw=0.5)
    ax.set_title(f"{name}  res={resolution}  V={len(mesh.uv)}  F={len(mesh.tris)}")
    ax.set_aspect("equal")
    OUT.mkdir(exist_ok=True)
    path = OUT / f"probe_mesh_{name}.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


def main():
    _render("rect_no",   helicoid())
    _render("rect_cy_v", cylinder_cy())
    _render("rect_mo_u", mobius_u())
    _render("rect_cy_cy", fig8())
    _render("disk_po",   disk_paraboloid_po())


if __name__ == "__main__":
    main()
