"""C14 — SIC chains over the (u, v) mesh on the fig-8 immersion.

Two subplots, one per sheet:
- Sheet A: walk each SIC in chain order, plotting `uv1[p_dp] → uv1/uv2[q_dp]`
  depending on the SIS's `flip`.
- Sheet B: the complementary preimage.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from surface_play.intersections import (
    build_sics, build_sis_pairs, find_double_points,
)
from surface_play.mesh import build_mesh
from surface_play.test_fixtures import fig8

OUT = Path(__file__).resolve().parent / "probes_out"


def _sheet_polylines(sics, sis_pairs, dps, start_sheet):
    """Trace each SIC as a *single* continuous preimage polyline.

    Sheet identity is local to each SIS (uv1 vs uv2 of a DP), and `flip=-1`
    swaps sheets across an SIS. To produce a single continuous preimage, we
    track the current sheet as state, flipping it at every `flip=-1` SIS. The
    second preimage is obtained by starting from the opposite sheet.

    `start_sheet`: 0 → start on uv1 at the first DP of the chain; 1 → uv2.
    """
    polys = []
    for sic in sics:
        if len(sic.sis_indices) == 0:
            polys.append(np.empty((0, 2)))
            continue

        first_tok = int(sic.sis_indices[0])
        first_rec = sis_pairs[abs(first_tok) - 1]
        start_dp = int(first_rec["p_dp"] if first_tok > 0 else first_rec["q_dp"])

        cur_sheet = start_sheet
        key = "uv1" if cur_sheet == 0 else "uv2"
        pts = [dps[key][start_dp]]

        for tok in sic.sis_indices:
            rec = sis_pairs[abs(int(tok)) - 1]
            p = int(rec["p_dp"])
            q = int(rec["q_dp"])
            flip = int(rec["flip"])
            end_dp = q if int(tok) > 0 else p
            if flip == -1:
                cur_sheet = 1 - cur_sheet
            key = "uv1" if cur_sheet == 0 else "uv2"
            pts.append(dps[key][end_dp])

        polys.append(np.asarray(pts))
    return polys


def main():
    surface = fig8()
    mesh = build_mesh(surface.domain, surface, resolution=20)
    dps = find_double_points(mesh, surface)
    sis_pairs = build_sis_pairs(dps)
    sics = build_sics(sis_pairs)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for start_sheet, (ax, label) in enumerate(
        zip(axes, ("preimage starting on uv1", "preimage starting on uv2"))
    ):
        ax.triplot(mesh.uv[:, 0], mesh.uv[:, 1], mesh.tris,
                   color="lightgray", lw=0.3)
        polys = _sheet_polylines(sics, sis_pairs, dps, start_sheet)
        for i, pl in enumerate(polys):
            ax.plot(pl[:, 0], pl[:, 1], "-o",
                    color=f"C{i % 10}", lw=1.5, ms=2.5, label=f"SIC{i}")
        ax.set_title(f"fig-8 res=20: {len(sics)} SIC(s)  —  {label}")
        ax.set_aspect("equal")
        if polys:
            ax.legend(fontsize=8)

    OUT.mkdir(exist_ok=True)
    path = OUT / "probe_sics_fig8.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}  ({len(sics)} SIC, {len(dps)} DPs, {len(sis_pairs)} SIS)")


if __name__ == "__main__":
    main()
