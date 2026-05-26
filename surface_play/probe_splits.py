"""Visual probe for SplitPoints / SPTs.

SPs by type (cn/bcp/bdp/cdp/tp/vp/ha) and SPTs on segments (BC/CC/SIC)
with their `vis_chge` values.

Usage:
    python -m surface_play.probe_splits [FIXTURE] [TRIAL] [RES]
"""

import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from surface_play._probe_common import (
    build_context, draw_domain_box, draw_outline_uv, draw_outline_xy,
    draw_surface_3d, parse_cli, save_two_panels,
)
from surface_play.curves import _seg_uv_at_bary


# Per-type marker / color for SPs
SP_STYLE = {
    "cn":  ("s", "tab:gray"),    # corner
    "bcp": ("o", "tab:blue"),    # boundary contact point
    "bdp": ("D", "tab:cyan"),    # boundary double point
    "cdp": ("v", "tab:purple"),  # contour double point
    "tp":  ("^", "tab:olive"),   # triple point
    "vp":  ("*", "tab:red"),     # cusp / VP
    "ha":  ("P", "tab:green"),   # helper anchor
}


def _list_spts_by_segment(R):
    """Return list of (kind, seg_idx, bary, sp_idx, vis_chge) per SPT host."""
    spts = R['splits'].spts_array()
    mesh = R['mesh']; css = R['css']; sis_pairs = R['sis_pairs']
    out = []

    def _scan(seg_array, kind, idx_list):
        if seg_array is None or len(seg_array) == 0:
            return
        for j in idx_list if idx_list is not None else range(len(seg_array)):
            seg = seg_array[j]
            for slot in ("split1", "split2"):
                spt_idx = int(seg[slot])
                if spt_idx < 0:
                    continue
                spt = spts[spt_idx]
                out.append((kind, int(j), float(spt['bary']),
                            int(spt['sp_idx']), int(spt['vis_chge']), spt_idx))

    # BCs: only BE edges (boundary edges)
    be_idx = mesh.boundary_edge_idx
    _scan(mesh.edges, "BC", be_idx.tolist())
    _scan(css, "CC", None)
    _scan(sis_pairs, "SIC", None)
    return out


def _spt_uv(kind, seg_idx, bary, R):
    """uv position of an SPT on its host segment."""
    return _seg_uv_at_bary(kind, seg_idx, bary,
                           R['mesh'], R['css'], R['sis_pairs'],
                           R['cps'], R['dps'])


KIND_COLOR = {"BC": "tab:blue", "CC": "tab:orange", "SIC": "tab:red"}


def probe(fixture_name, trial, resolution):
    R = build_context(fixture_name, trial, resolution)
    sps = R['sps']
    spt_rows = _list_spts_by_segment(R)
    domain = R['surf'].domain

    # type counts
    type_counts = {t: 0 for t in SP_STYLE}
    for sp in sps:
        t = str(sp['type'])
        if t in type_counts: type_counts[t] += 1
    kind_counts = {"BC": 0, "CC": 0, "SIC": 0}
    for row in spt_rows:
        kind_counts[row[0]] += 1

    print(f"\n=== probe_splits: {fixture_name} view={R['view_label']} res={resolution} ===")
    print(f"n_sps={len(sps)}  by type: " +
          ", ".join(f"{t}={c}" for t, c in type_counts.items() if c))
    print(f"n_spts={len(spt_rows)}  by host: " +
          ", ".join(f"{k}={c}" for k, c in kind_counts.items() if c))
    print(f"\nSPs:")
    for i, sp in enumerate(sps):
        print(f"  sp[{i:>3}] type={str(sp['type']):<3}  "
              f"uv=({sp['uv'][0]:+.3f},{sp['uv'][1]:+.3f})  "
              f"xy=({sp['xy'][0]:+.3f},{sp['xy'][1]:+.3f})")
    print(f"\nSPTs:")
    for kind, seg_idx, bary, sp_idx, vis_chge, spt_idx in spt_rows:
        print(f"  spt[{spt_idx:>3}] host={kind} seg={seg_idx:>4} "
              f"bary={bary:.3f}  sp={sp_idx:>3}  vis_chge={vis_chge:+d}")

    # ── context ────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 9))

    ax_uv = fig.add_subplot(1, 3, 1)
    ax_uv.set_title("UV  —  SPs by type, SPT bary positions marked")
    draw_domain_box(ax_uv, domain)
    draw_outline_uv(ax_uv, R, alpha=0.18)
    for t, (marker, col) in SP_STYLE.items():
        mask = sps['type'] == t
        if np.any(mask):
            ax_uv.plot(sps['uv'][mask, 0], sps['uv'][mask, 1], marker,
                       c=col, ms=8, mew=0.8, mec="black",
                       label=f"{t}  ({int(mask.sum())})")
    # SPT bary markers
    for kind, seg_idx, bary, sp_idx, vis_chge, _ in spt_rows:
        try:
            uv = _spt_uv(kind, seg_idx, bary, R)
        except Exception:
            continue
        ax_uv.plot(uv[0], uv[1], "x", c=KIND_COLOR[kind], ms=7, mew=1.2,
                   alpha=0.7)
    ax_uv.legend(loc="upper right", fontsize=7)
    ax_uv.set_xlabel("u"); ax_uv.set_ylabel("v")
    ax_uv.grid(True, alpha=0.3)

    ax_3d = fig.add_subplot(1, 3, 2, projection="3d")
    ax_3d.set_title("3D camera-frame  —  SPs by type")
    draw_surface_3d(ax_3d, R)
    I_v = R['I_v']; J_v = R['J_v']; axis_v = R['axis_v']
    for t, (marker, col) in SP_STYLE.items():
        mask = sps['type'] == t
        if np.any(mask):
            xyz = sps['xyz'][mask]
            ax_3d.scatter(xyz @ I_v, xyz @ J_v, xyz @ axis_v,
                          c=col, marker=marker, s=40,
                          edgecolors="black", linewidths=0.5,
                          label=f"{t} ({int(mask.sum())})")
    ax_3d.legend(loc="upper right", fontsize=7)
    ax_3d.view_init(elev=85, azim=-90)
    ax_3d.set_xlabel("I·S"); ax_3d.set_ylabel("J·S"); ax_3d.set_zlabel("axis·S")

    ax_t = fig.add_subplot(1, 3, 3); ax_t.axis("off")
    ax_t.set_title("Summary")
    lines = [f"fixture: {fixture_name}",
             f"view:    {R['view_label']}",
             f"res:     {resolution}",
             "",
             f"n_sps   = {len(sps)}"]
    for t, c in type_counts.items():
        if c: lines.append(f"  {t:<4}: {c}")
    lines.append("")
    lines.append(f"n_spts  = {len(spt_rows)}")
    for k, c in kind_counts.items():
        if c: lines.append(f"  on {k:<3}: {c}")
    lines.append("")
    lines.append("SPTs (first 30):")
    for row in spt_rows[:30]:
        kind, seg_idx, bary, sp_idx, vis_chge, spt_idx = row
        lines.append(f" [{spt_idx:>3}] {kind} seg={seg_idx:>4} "
                     f"b={bary:.2f} sp={sp_idx:>3} v={vis_chge:+d}")
    ax_t.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=7,
              va="top", transform=ax_t.transAxes)

    # ── XY wide ────────────────────────────────────────────────────────────
    fig_xy = plt.figure(figsize=(18, 7))
    ax_xy = fig_xy.add_subplot(1, 1, 1)
    ax_xy.set_title(f"XY projection — SPs by type + SPT vis_chge labels\n"
                    f"fixture={fixture_name}  view={R['view_label']}  "
                    f"res={resolution}  n_sps={len(sps)}  n_spts={len(spt_rows)}")
    draw_outline_xy(ax_xy, R, alpha=0.20)
    for t, (marker, col) in SP_STYLE.items():
        mask = sps['type'] == t
        if np.any(mask):
            ax_xy.plot(sps['xy'][mask, 0], sps['xy'][mask, 1], marker,
                       c=col, ms=11, mew=0.8, mec="black",
                       label=f"{t}  ({int(mask.sum())})")
    # SPT labels: vis_chge near sp xy with a tiny offset by host kind
    for kind, seg_idx, bary, sp_idx, vis_chge, spt_idx in spt_rows:
        sp = sps[sp_idx]
        dx = {"BC": -0.02, "CC": 0.02, "SIC": 0.0}[kind]
        dy = {"BC": 0.0,   "CC": 0.0,  "SIC": -0.025}[kind]
        ax_xy.text(sp['xy'][0] + dx, sp['xy'][1] + dy,
                   f"{vis_chge:+d}", color=KIND_COLOR[kind], fontsize=7,
                   alpha=0.9, weight="bold")
    ax_xy.legend(loc="upper left", fontsize=8)
    ax_xy.set_xlabel("x"); ax_xy.set_ylabel("y")
    ax_xy.set_aspect("equal")
    ax_xy.grid(True, alpha=0.3)

    save_two_panels(fig, fig_xy, R, "splits")


if __name__ == "__main__":
    probe(*parse_cli(sys.argv))
