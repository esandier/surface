"""Visual probe for ResampledCurves (O14: resample_all).

Samples colored by depth (axis·S, anchor-relative); per-RC vc_in/vc_out
labels; per-sample count N reported.

Usage:
    python -m surface_play.probe_resampled [FIXTURE] [TRIAL] [RES]
"""

import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib.lines import Line2D

from surface_play._probe_common import (
    build_context, draw_domain_box, draw_outline_xy, draw_surface_3d,
    parse_cli, save_two_panels,
)


def probe(fixture_name, trial, resolution):
    R = build_context(fixture_name, trial, resolution)
    rcs = R['rcs']
    domain = R['surf'].domain

    if not rcs:
        print(f"\n=== probe_resampled: {fixture_name} — no RCs ===")
        return

    # depth range across all samples for shared colormap
    all_depths = np.concatenate([rc.depth for rc in rcs])
    d_lo, d_hi = float(all_depths.min()), float(all_depths.max())
    if d_hi - d_lo < 1e-12:
        d_lo -= 0.5; d_hi += 0.5
    norm = mcolors.Normalize(vmin=d_lo, vmax=d_hi)
    cmap = plt.get_cmap("viridis")

    n_total = sum(len(rc.xy) for rc in rcs)
    kind_counts = {"BC": 0, "CC": 0, "SIC": 0, "HC": 0}
    for rc in rcs:
        kind_counts[rc.kind] += 1

    print(f"\n=== probe_resampled: {fixture_name} view={R['view_label']} res={resolution} ===")
    print(f"n_rcs={len(rcs)}  n_samples_total={n_total}  "
          f"depth=[{d_lo:+.3f},{d_hi:+.3f}]")
    print(f"by kind: " + ", ".join(f"{k}={c}" for k, c in kind_counts.items() if c))
    for i, rc in enumerate(rcs):
        print(f"  rc[{i:>3}] {rc.kind:<3} N={len(rc.xy):>4} "
              f"start={rc.start:>3} end={rc.end:>3}  "
              f"vc=({rc.vc_in:+d},{rc.vc_out:+d})  "
              f"depth=[{float(rc.depth.min()):+.3f},{float(rc.depth.max()):+.3f}]")

    # ── context ────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 9))

    ax_uv = fig.add_subplot(1, 3, 1)
    ax_uv.set_title("UV-space — RCs do not have native uv (skipped); see xy panel")
    draw_domain_box(ax_uv, domain)
    ax_uv.text(0.5, 0.5, "ResampledCurves live in image space.\n"
               "See XY panel for per-sample depth coloring.",
               transform=ax_uv.transAxes, ha="center", va="center",
               fontsize=11, color="gray")
    ax_uv.set_xlabel("u"); ax_uv.set_ylabel("v")
    ax_uv.grid(True, alpha=0.3)

    ax_3d = fig.add_subplot(1, 3, 2, projection="3d")
    ax_3d.set_title("3D camera-frame — RC samples colored by depth")
    draw_surface_3d(ax_3d, R)
    for rc in rcs:
        # We don't have native 3D; rebuild via xy + depth gives (x, y, z)
        ax_3d.scatter(rc.xy[:, 0], rc.xy[:, 1], rc.depth,
                      c=rc.depth, cmap=cmap, norm=norm, s=8, alpha=0.85)
    ax_3d.view_init(elev=85, azim=-90)
    ax_3d.set_xlabel("I·S"); ax_3d.set_ylabel("J·S"); ax_3d.set_zlabel("axis·S")

    ax_t = fig.add_subplot(1, 3, 3); ax_t.axis("off")
    ax_t.set_title("Summary")
    lines = [f"fixture: {fixture_name}",
             f"view:    {R['view_label']}",
             f"res:     {resolution}",
             "",
             f"n_rcs            = {len(rcs)}",
             f"n_samples_total  = {n_total}",
             f"depth range      = [{d_lo:+.3f}, {d_hi:+.3f}]",
             ""]
    for k, c in kind_counts.items():
        lines.append(f"  {k:<3}: {c}")
    lines.append("")
    lines.append("RCs (idx | kind | N | vc):")
    for i, rc in enumerate(rcs[:30]):
        lines.append(f" [{i:>3}] {rc.kind:<3} N={len(rc.xy):>4} "
                     f"vc=({rc.vc_in:+d},{rc.vc_out:+d})")
    ax_t.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=7,
              va="top", transform=ax_t.transAxes)

    # ── XY wide ────────────────────────────────────────────────────────────
    fig_xy = plt.figure(figsize=(18, 7))
    ax_xy = fig_xy.add_subplot(1, 1, 1)
    ax_xy.set_title(f"XY projection — RC samples colored by depth (axis·S)\n"
                    f"fixture={fixture_name}  view={R['view_label']}  "
                    f"res={resolution}  n_rcs={len(rcs)}  n_samples={n_total}  "
                    f"depth=[{d_lo:+.3f},{d_hi:+.3f}]")
    draw_outline_xy(ax_xy, R, alpha=0.10)

    # Plot per-segment lines colored by per-vertex depth (use mean depth)
    for rc in rcs:
        sx = rc.xy[:, 0]; sy = rc.xy[:, 1]; sd = rc.depth
        lw = {"SIC": 3.0, "CC": 2.0, "BC": 1.6, "HC": 1.0}[rc.kind]
        for k in range(len(sx) - 1):
            d_mid = 0.5 * (sd[k] + sd[k+1])
            ax_xy.plot(sx[k:k+2], sy[k:k+2], "-",
                       c=cmap(norm(d_mid)), lw=lw,
                       solid_capstyle="round")
        # endpoint vc labels
        if len(sx):
            ax_xy.annotate(f"{rc.vc_in:+d}", (sx[0], sy[0]),
                           fontsize=7, color="black", alpha=0.6)
            ax_xy.annotate(f"{rc.vc_out:+d}", (sx[-1], sy[-1]),
                           fontsize=7, color="black", alpha=0.6)
        # tiny sample dots
        ax_xy.plot(sx, sy, ".", c="k", ms=1.5, alpha=0.35)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cbar = fig_xy.colorbar(sm, ax=ax_xy, shrink=0.7, pad=0.02)
    cbar.set_label("depth (axis·S) — larger = closer to viewer")
    ax_xy.set_xlabel("x"); ax_xy.set_ylabel("y")
    ax_xy.set_aspect("equal")
    ax_xy.grid(True, alpha=0.3)

    save_two_panels(fig, fig_xy, R, "resampled")


if __name__ == "__main__":
    probe(*parse_cli(sys.argv))
