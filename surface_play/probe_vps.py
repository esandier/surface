"""Visual probe for cusp points / VPs (O4: find_vps).

Each VP: position (uv, xyz), `vis_change` (±1), host CS index.

Usage:
    python -m surface_play.probe_vps [FIXTURE] [TRIAL] [RES]
"""

import sys

import numpy as np
import matplotlib.pyplot as plt

from surface_play._probe_common import (
    build_context, draw_domain_box, draw_outline_uv, draw_outline_xy,
    draw_surface_3d, parse_cli, project_xyz, save_two_panels,
)


VIS_COLOR = {-1: "tab:red", +1: "tab:blue"}


def probe(fixture_name, trial, resolution):
    R = build_context(fixture_name, trial, resolution)
    vps = R['vps']; css = R['css']
    domain = R['surf'].domain

    n_neg = int(np.sum(vps['vis_change'] == -1)) if len(vps) else 0
    n_pos = int(np.sum(vps['vis_change'] == +1)) if len(vps) else 0
    print(f"\n=== probe_vps: {fixture_name} view={R['view_label']} res={resolution} ===")
    print(f"n_vps={len(vps)}  (vc=-1: {n_neg}, vc=+1: {n_pos})")
    for i, vp in enumerate(vps):
        print(f"  vp[{i}] cs={int(vp['cs'])}  s={float(vp['s']):.3f}  "
              f"uv=({vp['uv'][0]:+.3f},{vp['uv'][1]:+.3f})  "
              f"xyz=({vp['xyz'][0]:+.3f},{vp['xyz'][1]:+.3f},{vp['xyz'][2]:+.3f})  "
              f"vis_change={int(vp['vis_change']):+d}")

    # ── context ────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 9))

    ax_uv = fig.add_subplot(1, 3, 1)
    ax_uv.set_title("UV  —  VPs colored by vis_change")
    draw_domain_box(ax_uv, domain)
    draw_outline_uv(ax_uv, R, alpha=0.20)
    if len(vps):
        for sgn in (-1, +1):
            mask = vps['vis_change'] == sgn
            if np.any(mask):
                ax_uv.plot(vps['uv'][mask, 0], vps['uv'][mask, 1], "*",
                           c=VIS_COLOR[sgn], ms=14,
                           label=f"vc={sgn:+d}  n={int(mask.sum())}")
        ax_uv.legend(loc="upper right", fontsize=8)
    ax_uv.set_xlabel("u"); ax_uv.set_ylabel("v")
    ax_uv.grid(True, alpha=0.3)

    ax_3d = fig.add_subplot(1, 3, 2, projection="3d")
    ax_3d.set_title("3D camera-frame  —  VPs")
    draw_surface_3d(ax_3d, R)
    if len(vps):
        xyzC = project_xyz(vps['xyz'], R)
        for sgn in (-1, +1):
            mask = vps['vis_change'] == sgn
            if np.any(mask):
                ax_3d.scatter(xyzC[mask, 0], xyzC[mask, 1], xyzC[mask, 2],
                              c=VIS_COLOR[sgn], s=80, marker="*",
                              label=f"vc={sgn:+d}")
        ax_3d.legend(loc="upper right", fontsize=8)
    ax_3d.view_init(elev=85, azim=-90)
    ax_3d.set_xlabel("I·S"); ax_3d.set_ylabel("J·S"); ax_3d.set_zlabel("axis·S")

    ax_t = fig.add_subplot(1, 3, 3); ax_t.axis("off")
    ax_t.set_title("Summary")
    lines = [f"fixture: {fixture_name}",
             f"view:    {R['view_label']}",
             f"res:     {resolution}",
             "",
             f"n_vps   = {len(vps)}",
             f"  vc=-1: {n_neg}",
             f"  vc=+1: {n_pos}",
             f"n_css   = {len(css)}",
             "",
             "VPs (idx | host CS | s | vc):"]
    for i, vp in enumerate(vps):
        lines.append(f" [{i:>2}] cs={int(vp['cs']):>4} s={float(vp['s']):.2f} "
                     f"vc={int(vp['vis_change']):+d}")
    ax_t.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=8,
              va="top", transform=ax_t.transAxes)

    # ── XY wide ────────────────────────────────────────────────────────────
    fig_xy = plt.figure(figsize=(18, 7))
    ax_xy = fig_xy.add_subplot(1, 1, 1)
    ax_xy.set_title(f"XY projection — {len(vps)} VPs  "
                    f"(vc=-1: {n_neg}, vc=+1: {n_pos})\n"
                    f"fixture={fixture_name}  view={R['view_label']}  "
                    f"res={resolution}")
    draw_outline_xy(ax_xy, R, alpha=0.20)
    if len(vps):
        xyzC = project_xyz(vps['xyz'], R)
        for sgn in (-1, +1):
            mask = vps['vis_change'] == sgn
            if np.any(mask):
                ax_xy.plot(xyzC[mask, 0], xyzC[mask, 1], "*",
                           c=VIS_COLOR[sgn], ms=18,
                           label=f"vc={sgn:+d}  n={int(mask.sum())}")
        for i, vp in enumerate(vps):
            ax_xy.text(xyzC[i, 0], xyzC[i, 1], f" {i}",
                       color="black", fontsize=8, alpha=0.7)
        ax_xy.legend(loc="upper left", fontsize=9)
    ax_xy.set_xlabel("x"); ax_xy.set_ylabel("y")
    ax_xy.set_aspect("equal")
    ax_xy.grid(True, alpha=0.3)

    save_two_panels(fig, fig_xy, R, "vps")


if __name__ == "__main__":
    probe(*parse_cli(sys.argv))
