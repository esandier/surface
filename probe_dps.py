"""C14 — DP preimages in (u, v) space on the fig-8 immersion.

Each DP has two domain preimages `uv1, uv2`. They are plotted with circles
and squares respectively, colored by DP index, with a faint connector
between sibling preimages.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from surface_play.intersections import find_double_points
from surface_play.mesh import build_mesh
from surface_play.test_fixtures import fig8

OUT = Path(__file__).resolve().parent / "probes_out"


def main():
    surface = fig8()
    mesh = build_mesh(surface.domain, surface, resolution=20)
    dps = find_double_points(mesh, surface)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.triplot(mesh.uv[:, 0], mesh.uv[:, 1], mesh.tris,
               color="lightgray", lw=0.3)

    n = max(1, len(dps))
    colors = plt.cm.viridis(np.linspace(0, 1, n))
    for i, dp in enumerate(dps):
        c = colors[i]
        ax.plot([dp["uv1"][0]], [dp["uv1"][1]], "o",
                color=c, ms=5, mec="black", mew=0.3)
        ax.plot([dp["uv2"][0]], [dp["uv2"][1]], "s",
                color=c, ms=5, mec="black", mew=0.3)
        ax.plot([dp["uv1"][0], dp["uv2"][0]],
                [dp["uv1"][1], dp["uv2"][1]],
                "-", color=c, lw=0.4, alpha=0.4)

    ax.set_title(
        f"fig-8 res=20: {len(dps)} DPs   ○=uv1  □=uv2  (color = DP index)"
    )
    ax.set_aspect("equal")

    OUT.mkdir(exist_ok=True)
    path = OUT / "probe_dps_fig8.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}  ({len(dps)} DPs)")


if __name__ == "__main__":
    main()
