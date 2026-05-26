"""Visual probe for SubCurves (O12 + O13: assemble_subcurves + HCs).

Each SubCurve by kind (BC/CC/SIC/HC); start/end SPs; internal sample points
(seg_idx/bary); vc_in / vc_out.

Usage:
    python -m surface_play.probe_subcurves [FIXTURE] [TRIAL] [RES]
"""

import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from surface_play._probe_common import (
    build_context, draw_domain_box, draw_outline_uv, draw_outline_xy,
    draw_surface_3d, parse_cli, project_xyz, save_two_panels,
)
from surface_play.probe_visibility import _rc_uv, _detect_seam_crossings


KIND_COLOR = {"BC": "tab:blue", "CC": "tab:orange",
              "HC": "tab:green", "SIC": "tab:red"}


def _vc_pair_marker(vc_in, vc_out):
    """Marker style by (vc_in, vc_out) — both -1 / both 0 / mixed."""
    if vc_in == 0 and vc_out == 0:
        return "o", "lightgray"
    if vc_in == -1 and vc_out == -1:
        return "s", "black"
    return "D", "magenta"  # mixed (rare)


def probe(fixture_name, trial, resolution):
    R = build_context(fixture_name, trial, resolution)
    subs = R['subs']; splits = R['splits']
    sps = R['sps']
    mesh = R['mesh']; css = R['css']; sis_pairs = R['sis_pairs']
    cps = R['cps']; dps = R['dps']
    domain = R['surf'].domain

    kind_counts = {"BC": 0, "CC": 0, "SIC": 0, "HC": 0}
    closed_count = 0
    for s in subs:
        kind_counts[s.kind] += 1
        if s.is_closed: closed_count += 1

    print(f"\n=== probe_subcurves: {fixture_name} view={R['view_label']} res={resolution} ===")
    print(f"n_subs={len(subs)}  closed={closed_count}  "
          f"by kind: " + ", ".join(f"{k}={c}" for k, c in kind_counts.items() if c))
    for i, sub in enumerate(subs):
        start_t = str(sps[sub.start]['type']) if sub.start >= 0 else "-"
        end_t = str(sps[sub.end]['type']) if sub.end >= 0 else "-"
        print(f"  sub[{i:>3}] {sub.kind:<3} {'C' if sub.is_closed else 'O'} "
              f"start={sub.start:>3}({start_t}) end={sub.end:>3}({end_t})  "
              f"n_int={len(sub.internal):>3}  "
              f"vc=({sub.vc_in:+d},{sub.vc_out:+d})  "
              f"parent={sub.parent_idx}")

    # ── context ────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 9))

    ax_uv = fig.add_subplot(1, 3, 1)
    ax_uv.set_title(f"UV  —  {len(subs)} SubCurves by kind")
    draw_domain_box(ax_uv, domain)
    for i, sub in enumerate(subs):
        col = KIND_COLOR[sub.kind]
        uv_xy = _rc_uv(sub, splits, mesh, css, sis_pairs, cps, dps)
        if uv_xy is None or len(uv_xy) < 2:
            continue
        crossings = _detect_seam_crossings(uv_xy, domain)
        seg_starts = [0] + [c + 1 for c in crossings] + [len(uv_xy)]
        lw = {"SIC": 2.4, "CC": 1.6, "BC": 1.2, "HC": 1.0}[sub.kind]
        ls = "--" if sub.kind == "HC" else "-"
        for s, e in zip(seg_starts[:-1], seg_starts[1:]):
            if e - s >= 2:
                ax_uv.plot(uv_xy[s:e, 0], uv_xy[s:e, 1], ls, c=col, lw=lw)
        # endpoint SPs
        if sub.start >= 0:
            ax_uv.plot(sps[sub.start]['uv'][0], sps[sub.start]['uv'][1],
                       "ko", ms=3, alpha=0.6)
        if sub.end >= 0:
            ax_uv.plot(sps[sub.end]['uv'][0], sps[sub.end]['uv'][1],
                       "ko", ms=3, alpha=0.6)
    legend = [Line2D([0], [0], color=KIND_COLOR[k], lw=2,
                     ls=("--" if k == "HC" else "-"),
                     label=f"{k} ({kind_counts[k]})")
              for k in ("BC", "CC", "SIC", "HC")]
    ax_uv.legend(handles=legend, loc="upper right", fontsize=8)
    ax_uv.set_xlabel("u"); ax_uv.set_ylabel("v")
    ax_uv.grid(True, alpha=0.3)

    ax_3d = fig.add_subplot(1, 3, 2, projection="3d")
    ax_3d.set_title("3D camera-frame  —  SubCurves by kind")
    draw_surface_3d(ax_3d, R)
    surf = R['surf']
    for sub in subs:
        col = KIND_COLOR[sub.kind]
        uv_xy = _rc_uv(sub, splits, mesh, css, sis_pairs, cps, dps)
        if uv_xy is None or len(uv_xy) < 2: continue
        xyz = np.array([surf.S(float(p[0]), float(p[1])) for p in uv_xy])
        xyzC = project_xyz(xyz, R)
        lw = {"SIC": 2.4, "CC": 1.6, "BC": 1.2, "HC": 1.0}[sub.kind]
        ls = "--" if sub.kind == "HC" else "-"
        ax_3d.plot(xyzC[:, 0], xyzC[:, 1], xyzC[:, 2], ls, c=col, lw=lw)
    ax_3d.view_init(elev=85, azim=-90)
    ax_3d.set_xlabel("I·S"); ax_3d.set_ylabel("J·S"); ax_3d.set_zlabel("axis·S")

    ax_t = fig.add_subplot(1, 3, 3); ax_t.axis("off")
    ax_t.set_title("Summary")
    lines = [f"fixture: {fixture_name}",
             f"view:    {R['view_label']}",
             f"res:     {resolution}",
             "",
             f"n_subs   = {len(subs)}",
             f"  closed = {closed_count}"]
    for k, c in kind_counts.items():
        lines.append(f"  {k:<3}: {c}")
    lines.append("")
    lines.append("SubCurves (first 30 — kind|C/O|n_int|vc):")
    for i, sub in enumerate(subs[:30]):
        co = "C" if sub.is_closed else "O"
        lines.append(f" [{i:>3}] {sub.kind:<3} {co} n_int={len(sub.internal):>3} "
                     f"vc=({sub.vc_in:+d},{sub.vc_out:+d})")
    ax_t.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=7,
              va="top", transform=ax_t.transAxes)

    # ── XY wide ────────────────────────────────────────────────────────────
    fig_xy = plt.figure(figsize=(18, 7))
    ax_xy = fig_xy.add_subplot(1, 1, 1)
    ax_xy.set_title(f"XY projection — SubCurves by kind, vc_in/vc_out at ends\n"
                    f"fixture={fixture_name}  view={R['view_label']}  "
                    f"res={resolution}  n_subs={len(subs)}")
    draw_outline_xy(ax_xy, R, alpha=0.10)
    for i, sub in enumerate(subs):
        col = KIND_COLOR[sub.kind]
        uv_xy = _rc_uv(sub, splits, mesh, css, sis_pairs, cps, dps)
        if uv_xy is None or len(uv_xy) < 2: continue
        xyz = np.array([surf.S(float(p[0]), float(p[1])) for p in uv_xy])
        xyzC = project_xyz(xyz, R)
        lw = {"SIC": 3.0, "CC": 2.0, "BC": 1.5, "HC": 1.0}[sub.kind]
        ls = "--" if sub.kind == "HC" else "-"
        ax_xy.plot(xyzC[:, 0], xyzC[:, 1], ls, c=col, lw=lw)
    # Endpoint markers with vc styling
    for sub in subs:
        for sp_idx, vc in ((sub.start, sub.vc_in), (sub.end, sub.vc_out)):
            if sp_idx < 0: continue
            marker, mcol = _vc_pair_marker(vc, vc)
            sp = sps[sp_idx]
            ax_xy.plot(sp['xy'][0], sp['xy'][1], marker, c=mcol, ms=4,
                       alpha=0.7)
    legend = [Line2D([0], [0], color=KIND_COLOR[k], lw=2,
                     ls=("--" if k == "HC" else "-"),
                     label=f"{k} ({kind_counts[k]})")
              for k in ("BC", "CC", "SIC", "HC")]
    legend += [
        Line2D([0], [0], marker="o", color="lightgray", lw=0, ms=8,
               label="vc=(0,0)"),
        Line2D([0], [0], marker="s", color="black", lw=0, ms=8,
               label="vc=(-1,-1)"),
        Line2D([0], [0], marker="D", color="magenta", lw=0, ms=8,
               label="vc mixed"),
    ]
    ax_xy.legend(handles=legend, loc="upper left", fontsize=8)
    ax_xy.set_xlabel("x"); ax_xy.set_ylabel("y")
    ax_xy.set_aspect("equal")
    ax_xy.grid(True, alpha=0.3)

    save_two_panels(fig, fig_xy, R, "subcurves")


if __name__ == "__main__":
    probe(*parse_cli(sys.argv))
