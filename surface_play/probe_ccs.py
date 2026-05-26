"""Visual probe for ContourCurves (O3: build_contour_curves).

Each CC plotted as a colored chain of CS segments; label by CC index;
closed CCs drawn solid, open CCs dashed.

Usage:
    python -m surface_play.probe_ccs [FIXTURE] [TRIAL] [RES]
"""

import sys

import numpy as np
import matplotlib.pyplot as plt

from surface_play._probe_common import (
    build_context, draw_domain_box, draw_outline_uv, draw_outline_xy,
    draw_surface_3d, parse_cli, project_xyz, save_two_panels,
)


def _cc_cp_uv_chain(cc, css, cps):
    """uv-space polyline of CP positions along a CC chain."""
    pts = []
    for j, si in enumerate(cc.cs_indices):
        cs = css[abs(int(si)) - 1]
        p_cp = int(cs["p_cp"]) if si > 0 else int(cs["q_cp"])
        q_cp = int(cs["q_cp"]) if si > 0 else int(cs["p_cp"])
        if j == 0:
            pts.append(np.asarray(cps[p_cp]["uv"], dtype=float))
        pts.append(np.asarray(cps[q_cp]["uv"], dtype=float))
    return np.asarray(pts)


def _cc_cp_xyz_chain(cc, css, cps):
    pts = []
    for j, si in enumerate(cc.cs_indices):
        cs = css[abs(int(si)) - 1]
        p_cp = int(cs["p_cp"]) if si > 0 else int(cs["q_cp"])
        q_cp = int(cs["q_cp"]) if si > 0 else int(cs["p_cp"])
        if j == 0:
            pts.append(np.asarray(cps[p_cp]["xyz"], dtype=float))
        pts.append(np.asarray(cps[q_cp]["xyz"], dtype=float))
    return np.asarray(pts)


def probe(fixture_name, trial, resolution):
    R = build_context(fixture_name, trial, resolution)
    ccs = R['ccs']; css = R['css']; cps = R['cps']
    domain = R['surf'].domain

    n_closed = sum(1 for cc in ccs if cc.is_closed)
    n_open = len(ccs) - n_closed
    print(f"\n=== probe_ccs: {fixture_name} view={R['view_label']} res={resolution} ===")
    print(f"n_ccs={len(ccs)}  (closed={n_closed}, open={n_open})  "
          f"n_css={len(css)}  n_cps={len(cps)}")
    for i, cc in enumerate(ccs):
        flag = "closed" if cc.is_closed else "open"
        print(f"  cc[{i}] {flag}  n_cs={len(cc.cs_indices)}  "
              f"cs_indices={cc.cs_indices.tolist()}")

    cmap = plt.get_cmap("tab20", max(len(ccs), 1))

    # ── context ────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 9))

    ax_uv = fig.add_subplot(1, 3, 1)
    ax_uv.set_title(f"UV  —  {len(ccs)} ContourCurves "
                    f"(solid=closed, dashed=open)")
    draw_domain_box(ax_uv, domain)
    draw_outline_uv(ax_uv, R, alpha=0.12)
    for i, cc in enumerate(ccs):
        col = cmap(i)
        uv = _cc_cp_uv_chain(cc, css, cps)
        ls = "-" if cc.is_closed else "--"
        # Split at seam discontinuities (pu/pv jumps in identified domains)
        pu = domain.period_u if not np.isnan(domain.period_u) else 0.0
        pv = domain.period_v if not np.isnan(domain.period_v) else 0.0
        breaks = []
        for k in range(len(uv) - 1):
            du = abs(uv[k+1, 0] - uv[k, 0])
            dv = abs(uv[k+1, 1] - uv[k, 1])
            if (pu > 0 and du > 0.5 * pu) or (pv > 0 and dv > 0.5 * pv):
                breaks.append(k)
        seg_starts = [0] + [b + 1 for b in breaks] + [len(uv)]
        for s, e in zip(seg_starts[:-1], seg_starts[1:]):
            if e - s >= 2:
                ax_uv.plot(uv[s:e, 0], uv[s:e, 1], ls, c=col, lw=1.6)
        if len(uv):
            mid = uv[len(uv) // 2]
            ax_uv.text(mid[0], mid[1], f"{i}", color=col, fontsize=9,
                       weight="bold")
    ax_uv.set_xlabel("u"); ax_uv.set_ylabel("v")
    ax_uv.grid(True, alpha=0.3)

    ax_3d = fig.add_subplot(1, 3, 2, projection="3d")
    ax_3d.set_title("3D camera-frame  —  CCs colored by index")
    draw_surface_3d(ax_3d, R)
    for i, cc in enumerate(ccs):
        col = cmap(i)
        xyz = _cc_cp_xyz_chain(cc, css, cps)
        xyzC = project_xyz(xyz, R)
        ax_3d.plot(xyzC[:, 0], xyzC[:, 1], xyzC[:, 2], "-", c=col, lw=1.8)
    ax_3d.view_init(elev=85, azim=-90)
    ax_3d.set_xlabel("I·S"); ax_3d.set_ylabel("J·S"); ax_3d.set_zlabel("axis·S")

    ax_t = fig.add_subplot(1, 3, 3); ax_t.axis("off")
    ax_t.set_title("Summary")
    lines = [f"fixture: {fixture_name}",
             f"view:    {R['view_label']}",
             f"res:     {resolution}",
             "",
             f"n_ccs    = {len(ccs)}",
             f"  closed = {n_closed}",
             f"  open   = {n_open}",
             f"n_css    = {len(css)}",
             f"n_cps    = {len(cps)}",
             "",
             "CCs (idx | closed | n_cs):"]
    for i, cc in enumerate(ccs):
        flag = "C" if cc.is_closed else "O"
        lines.append(f" [{i:>2}] {flag}  n_cs={len(cc.cs_indices):>3}")
    ax_t.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=8,
              va="top", transform=ax_t.transAxes)

    # ── XY wide ────────────────────────────────────────────────────────────
    fig_xy = plt.figure(figsize=(18, 7))
    ax_xy = fig_xy.add_subplot(1, 1, 1)
    ax_xy.set_title(f"XY projection — {len(ccs)} ContourCurves "
                    f"(solid=closed, dashed=open)\n"
                    f"fixture={fixture_name}  view={R['view_label']}  "
                    f"res={resolution}")
    draw_outline_xy(ax_xy, R, alpha=0.12)
    for i, cc in enumerate(ccs):
        col = cmap(i)
        xyz = _cc_cp_xyz_chain(cc, css, cps)
        xyzC = project_xyz(xyz, R)
        ls = "-" if cc.is_closed else "--"
        ax_xy.plot(xyzC[:, 0], xyzC[:, 1], ls, c=col, lw=2.2,
                   label=f"cc[{i}] {'closed' if cc.is_closed else 'open'} "
                         f"(n_cs={len(cc.cs_indices)})")
        if len(xyzC):
            mid = xyzC[len(xyzC) // 2]
            ax_xy.text(mid[0], mid[1], f"{i}", color=col, fontsize=11,
                       weight="bold")
    if len(ccs) <= 20:
        ax_xy.legend(loc="upper left", fontsize=8)
    ax_xy.set_xlabel("x"); ax_xy.set_ylabel("y")
    ax_xy.set_aspect("equal")
    ax_xy.grid(True, alpha=0.3)

    save_two_panels(fig, fig_xy, R, "ccs")


if __name__ == "__main__":
    probe(*parse_cli(sys.argv))
