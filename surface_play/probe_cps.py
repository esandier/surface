"""Visual probe for Contour Points (O1: find_contour_points).

Usage:
    python -m surface_play.probe_cps [FIXTURE] [TRIAL] [RES]
"""

import sys

import numpy as np
import matplotlib.pyplot as plt

from surface_play._probe_common import (
    build_context, draw_domain_box, draw_outline_uv, draw_outline_xy,
    draw_surface_3d, parse_cli, project_xyz, save_two_panels, uv_extent,
)


PTYPE_COLOR = {0: "tab:blue", 4: "tab:red"}
PTYPE_LABEL = {0: "interior", 4: "boundary"}


def probe(fixture_name, trial, resolution):
    R = build_context(fixture_name, trial, resolution)
    cps = R['cps']
    domain = R['surf'].domain

    n_int = int(np.sum(cps['ptype'] == 0))
    n_bnd = int(np.sum(cps['ptype'] == 4))
    print(f"\n=== probe_cps: {fixture_name} view={R['view_label']} res={resolution} ===")
    print(f"axis = {R['axis_v'].round(4).tolist()}")
    print(f"n_cps={len(cps)}  (interior={n_int}, boundary={n_bnd})")
    for i, cp in enumerate(cps):
        print(f"  cp[{i:>3}] ptype={int(cp['ptype'])}  e={int(cp['e'])}  "
              f"s={float(cp['s']):.3f}  "
              f"uv=({cp['uv'][0]:+.3f},{cp['uv'][1]:+.3f})  "
              f"d=({cp['d'][0]:+.3f},{cp['d'][1]:+.3f})  "
              f"front_uv=({cp['front_uv'][0]:+.3f},{cp['front_uv'][1]:+.3f})")

    # ── context: UV + 3D + summary ─────────────────────────────────────────
    fig = plt.figure(figsize=(18, 9))

    ax_uv = fig.add_subplot(1, 3, 1)
    ax_uv.set_title("UV  —  CPs; d (blue) and front_uv (green) arrows")
    draw_domain_box(ax_uv, domain)
    draw_outline_uv(ax_uv, R)
    if len(cps):
        scale = 0.04 * uv_extent(domain, R)
        for cp in cps:
            uv = cp['uv']; d = cp['d']; fu = cp['front_uv']
            ax_uv.arrow(uv[0], uv[1], d[0]*scale, d[1]*scale,
                        head_width=scale*0.25, color="tab:blue",
                        alpha=0.7, length_includes_head=True)
            ax_uv.arrow(uv[0], uv[1], fu[0]*scale, fu[1]*scale,
                        head_width=scale*0.25, color="tab:green",
                        alpha=0.7, length_includes_head=True)
        for pt in (0, 4):
            mask = cps['ptype'] == pt
            if np.any(mask):
                ax_uv.plot(cps['uv'][mask, 0], cps['uv'][mask, 1], "o",
                           c=PTYPE_COLOR[pt], ms=5,
                           label=f"ptype={pt} ({PTYPE_LABEL[pt]})")
        ax_uv.legend(loc="upper right", fontsize=8)
    ax_uv.set_xlabel("u"); ax_uv.set_ylabel("v")
    ax_uv.grid(True, alpha=0.3)

    ax_3d = fig.add_subplot(1, 3, 2, projection="3d")
    ax_3d.set_title("3D camera-frame  —  CPs by ptype")
    draw_surface_3d(ax_3d, R)
    if len(cps):
        xyzC = project_xyz(cps['xyz'], R)
        for pt in (0, 4):
            mask = cps['ptype'] == pt
            if np.any(mask):
                ax_3d.scatter(xyzC[mask, 0], xyzC[mask, 1], xyzC[mask, 2],
                              c=PTYPE_COLOR[pt], s=20, label=f"ptype={pt}")
        ax_3d.legend(loc="upper right", fontsize=8)
    ax_3d.view_init(elev=85, azim=-90)
    ax_3d.set_xlabel("I·S"); ax_3d.set_ylabel("J·S"); ax_3d.set_zlabel("axis·S")

    ax_t = fig.add_subplot(1, 3, 3); ax_t.axis("off")
    ax_t.set_title("Summary")
    lines = [f"fixture: {fixture_name}",
             f"view:    {R['view_label']}",
             f"res:     {resolution}",
             f"axis:    {R['axis_v'].round(3).tolist()}",
             "",
             f"n_cps   = {len(cps)}",
             f"  interior (ptype=0): {n_int}",
             f"  boundary (ptype=4): {n_bnd}",
             "",
             "First 25 CPs:"]
    for i, cp in enumerate(cps[:25]):
        lines.append(
            f" [{i:>3}] p={int(cp['ptype'])} e={int(cp['e']):>4} "
            f"uv=({cp['uv'][0]:+.2f},{cp['uv'][1]:+.2f})")
    ax_t.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=8,
              va="top", transform=ax_t.transAxes)

    # ── XY wide ────────────────────────────────────────────────────────────
    fig_xy = plt.figure(figsize=(18, 7))
    ax_xy = fig_xy.add_subplot(1, 1, 1)
    ax_xy.set_title(f"XY projection — CPs by ptype\n"
                    f"fixture={fixture_name}  view={R['view_label']}  "
                    f"res={resolution}  n_cps={len(cps)} "
                    f"(int={n_int}, bnd={n_bnd})")
    draw_outline_xy(ax_xy, R)
    if len(cps):
        xyzC = project_xyz(cps['xyz'], R)
        for pt in (0, 4):
            mask = cps['ptype'] == pt
            if np.any(mask):
                ax_xy.plot(xyzC[mask, 0], xyzC[mask, 1], "o",
                           c=PTYPE_COLOR[pt], ms=6,
                           label=f"ptype={pt} ({PTYPE_LABEL[pt]})  "
                                 f"n={int(mask.sum())}")
        ax_xy.legend(loc="upper left", fontsize=9)
    ax_xy.set_xlabel("x"); ax_xy.set_ylabel("y")
    ax_xy.set_aspect("equal")
    ax_xy.grid(True, alpha=0.3)

    save_two_panels(fig, fig_xy, R, "cps")


if __name__ == "__main__":
    probe(*parse_cli(sys.argv))
