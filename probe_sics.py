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


def _sheet_polylines(sics, sis_pairs, dps, sheet):
    """Return list of polyline arrays, one per SIC, on the requested sheet (0=A, 1=B).

    For each SIC chain token (signed 1-indexed SIS id), append the preimage
    points on the requested sheet, respecting the SIS `flip` so that
    sheet-A endpoints map to the A-sheet of both endpoint DPs.
    """
    polys = []
    for sic in sics:
        pts = []
        for tok in sic.sis_indices:
            sis_idx = abs(int(tok)) - 1
            rec = sis_pairs[sis_idx]
            p = int(rec["p_dp"])
            q = int(rec["q_dp"])
            flip = int(rec["flip"])
            if int(tok) < 0:
                p, q = q, p
            # On sheet A: A-end at p is uv1[p]; A-end at q is uv1[q] (flip=+1)
            # or uv2[q] (flip=-1). Sheet B swaps 1↔2 throughout.
            if sheet == 0:
                p_uv = dps["uv1"][p]
                q_uv = dps["uv1"][q] if flip == 1 else dps["uv2"][q]
            else:
                p_uv = dps["uv2"][p]
                q_uv = dps["uv2"][q] if flip == 1 else dps["uv1"][q]
            if not pts:
                pts.append(p_uv)
            pts.append(q_uv)
        polys.append(np.asarray(pts))
    return polys


def main():
    surface = fig8()
    mesh = build_mesh(surface.domain, surface, resolution=20)
    dps = find_double_points(mesh, surface)
    sis_pairs = build_sis_pairs(dps)
    sics = build_sics(sis_pairs)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for sheet, (ax, label) in enumerate(zip(axes, ("sheet A (uv1)", "sheet B (uv2)"))):
        ax.triplot(mesh.uv[:, 0], mesh.uv[:, 1], mesh.tris,
                   color="lightgray", lw=0.3)
        polys = _sheet_polylines(sics, sis_pairs, dps, sheet)
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
