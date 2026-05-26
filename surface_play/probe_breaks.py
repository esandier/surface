"""Visual probe for projection breaks (O15: compute_projection_breaks).

Each break: xy position, delta_v (±1 if BC occluder, ±2 if CC occluder),
occluded RC index and sample index. The occluder is whichever RC visibly
crosses through the break point — |delta_v| reveals its kind, sign reveals
which side is in front.

Usage:
    python -m surface_play.probe_breaks [FIXTURE] [TRIAL] [RES]
"""

import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from surface_play._probe_common import (
    build_context, draw_outline_xy, draw_outline_uv, draw_domain_box,
    parse_cli, save_two_panels,
)


def _dv_style(dv):
    """Marker / color for a delta_v value."""
    mag = abs(int(dv))
    # magnitude → color: 1 (BC) green, 2 (CC) red, other gray
    col = {1: "tab:green", 2: "tab:red"}.get(mag, "tab:gray")
    # sign → marker
    marker = "v" if dv < 0 else "^"
    return marker, col


def probe(fixture_name, trial, resolution):
    R = build_context(fixture_name, trial, resolution)
    breaks = R['breaks']; rcs = R['rcs']
    domain = R['surf'].domain

    n_total = len(breaks)
    by_mag = {1: 0, 2: 0}
    by_sign = {-1: 0, +1: 0}
    if n_total:
        for b in breaks:
            m = abs(int(b['delta_v']))
            by_mag[m] = by_mag.get(m, 0) + 1
            by_sign[int(np.sign(int(b['delta_v'])))] = (
                by_sign.get(int(np.sign(int(b['delta_v']))), 0) + 1)

    print(f"\n=== probe_breaks: {fixture_name} view={R['view_label']} res={resolution} ===")
    print(f"n_breaks={n_total}  by |dv|: " +
          ", ".join(f"{m}={c}" for m, c in by_mag.items() if c) +
          "  by sign: " +
          ", ".join(f"{('+' if s>0 else '-')}={c}" for s, c in by_sign.items() if c))
    print()
    for bi, b in enumerate(breaks):
        ri = int(b['rc_idx'])
        rc = rcs[ri]
        print(f"  bk[{bi:>3}] occluded=rc[{ri}]({rc.kind}) "
              f"si={int(b['sample_idx']):>4} t={float(b['t']):.3f}  "
              f"xy=({b['xy'][0]:+.3f},{b['xy'][1]:+.3f})  "
              f"delta_v={int(b['delta_v']):+d}  "
              f"(|dv|={abs(int(b['delta_v']))} -> occluder kind = "
              f"{ {1:'BC', 2:'CC'}.get(abs(int(b['delta_v'])), '?') })")

    # ── context ────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 9))

    ax_uv = fig.add_subplot(1, 3, 1)
    ax_uv.set_title("UV — break positions don't have a direct uv (image-space "
                    "only). Outline shown for orientation.")
    draw_domain_box(ax_uv, domain)
    draw_outline_uv(ax_uv, R, alpha=0.25)
    ax_uv.set_xlabel("u"); ax_uv.set_ylabel("v")
    ax_uv.grid(True, alpha=0.3)

    # 3D panel: place breaks at their (xy, occluded-depth) location
    ax_3d = fig.add_subplot(1, 3, 2, projection="3d")
    ax_3d.set_title("3D camera-frame — breaks at occluded depth")
    for rc in rcs:
        ax_3d.plot(rc.xy[:, 0], rc.xy[:, 1], rc.depth, "-", c="lightgray",
                   lw=0.6, alpha=0.5)
    for bi, b in enumerate(breaks):
        ri = int(b['rc_idx']); si = int(b['sample_idx'])
        t = float(b['t'])
        rc = rcs[ri]
        if 0 <= si < len(rc.depth) - 1:
            d = (1.0 - t) * float(rc.depth[si]) + t * float(rc.depth[si + 1])
        else:
            d = float(rc.depth[min(si, len(rc.depth) - 1)])
        marker, col = _dv_style(int(b['delta_v']))
        ax_3d.scatter([b['xy'][0]], [b['xy'][1]], [d],
                      c=col, marker=marker, s=50)
    ax_3d.view_init(elev=85, azim=-90)
    ax_3d.set_xlabel("x = I·S"); ax_3d.set_ylabel("y = J·S")
    ax_3d.set_zlabel("depth = axis·S")

    ax_t = fig.add_subplot(1, 3, 3); ax_t.axis("off")
    ax_t.set_title("Summary")
    lines = [f"fixture: {fixture_name}",
             f"view:    {R['view_label']}",
             f"res:     {resolution}",
             "",
             f"n_breaks = {n_total}",
             ""]
    for m, c in sorted(by_mag.items()):
        if c:
            kind = {1: "BC", 2: "CC"}.get(m, "?")
            lines.append(f"  |dv|={m} ({kind} occluder): {c}")
    lines.append("")
    for s, c in sorted(by_sign.items(), reverse=True):
        if c: lines.append(f"  sign={'+' if s>0 else '-'}: {c}")
    lines.append("")
    lines.append("Breaks (first 30):")
    for bi, b in enumerate(breaks[:30]):
        lines.append(f" [{bi:>3}] rc={int(b['rc_idx']):>3} "
                     f"si={int(b['sample_idx']):>4} "
                     f"t={float(b['t']):.2f} "
                     f"dv={int(b['delta_v']):+d}")
    ax_t.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=7,
              va="top", transform=ax_t.transAxes)

    # ── XY wide ────────────────────────────────────────────────────────────
    fig_xy = plt.figure(figsize=(18, 7))
    ax_xy = fig_xy.add_subplot(1, 1, 1)
    ax_xy.set_title(f"XY projection — {n_total} projection breaks; "
                    f"marker shape = sign, color = occluder kind\n"
                    f"fixture={fixture_name}  view={R['view_label']}  "
                    f"res={resolution}")
    draw_outline_xy(ax_xy, R, alpha=0.30)
    # Highlight occluded RCs in slightly darker line
    occluded_ris = set(int(b['rc_idx']) for b in breaks)
    for ri in occluded_ris:
        rc = rcs[ri]
        ax_xy.plot(rc.xy[:, 0], rc.xy[:, 1], "-", c="dimgray", lw=1.4,
                   alpha=0.55)
    # Break markers
    seen_styles = set()
    for bi, b in enumerate(breaks):
        marker, col = _dv_style(int(b['delta_v']))
        label = None
        key = (marker, col)
        if key not in seen_styles:
            seen_styles.add(key)
            sign = "+" if int(b['delta_v']) > 0 else "-"
            mag = abs(int(b['delta_v']))
            kind = {1: "BC", 2: "CC"}.get(mag, "?")
            label = f"{sign}{mag}  ({kind} occluder)"
        ax_xy.plot(b['xy'][0], b['xy'][1], marker, c=col, ms=10, mew=0.7,
                   mec="black", label=label, alpha=0.9)
        ax_xy.annotate(f"{int(b['delta_v']):+d}",
                       (b['xy'][0], b['xy'][1]),
                       textcoords="offset points", xytext=(6, 4),
                       fontsize=7, color=col, weight="bold")
    if seen_styles:
        ax_xy.legend(loc="upper left", fontsize=8)
    ax_xy.set_xlabel("x"); ax_xy.set_ylabel("y")
    ax_xy.set_aspect("equal")
    ax_xy.grid(True, alpha=0.3)

    save_two_panels(fig, fig_xy, R, "breaks")


if __name__ == "__main__":
    probe(*parse_cli(sys.argv))
