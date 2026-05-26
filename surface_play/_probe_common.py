"""Shared helpers for the Layer O probe family.

Each `probe_*.py` module:
- imports `parse_cli`, `build_context`, `draw_outline_uv`, `draw_outline_xy`,
  `draw_domain_box`, `save_two_panels` from here;
- builds its own UV / 3D / XY artists for the specific artifact;
- prints a one-page console summary.

`build_context` runs `_build_outline` from probe_visibility once and returns
the result dict + view metadata.
"""

from __future__ import annotations

import sys
from typing import Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from surface_play.probe_visibility import (
    FIXTURE_BUILDERS, _build_outline, _detect_seam_crossings, _kind_color,
    _random_axis, _rc_uv,
)


def parse_cli(argv, default_fixture="mobius_v", default_trial="0",
              default_res=100) -> Tuple[str, str, int]:
    fixture = argv[1] if len(argv) > 1 else default_fixture
    trial = argv[2] if len(argv) > 2 else default_trial
    res = int(argv[3]) if len(argv) > 3 else default_res
    return fixture, trial, res


def build_context(fixture_name, trial, resolution) -> dict:
    if fixture_name not in FIXTURE_BUILDERS:
        raise SystemExit(f"unknown fixture {fixture_name!r}; "
                         f"choose from {list(FIXTURE_BUILDERS)}")
    builder, canonical = FIXTURE_BUILDERS[fixture_name]
    surf = builder(perturb=False)

    if trial == "canonical":
        if canonical is None:
            raise SystemExit(f"no canonical view for {fixture_name}")
        I_in, J_in, _eye = canonical
        I = np.asarray(I_in, dtype=float); J = np.asarray(J_in, dtype=float)
        a = np.cross(I, J)
        view_label = "canonical"
    else:
        t = int(trial)
        a, I, J = _random_axis(t)
        view_label = f"trial {t}"

    R = _build_outline(surf, I, J, resolution)
    R['view_label'] = view_label
    R['axis_v'] = np.asarray(a)
    R['I_v'] = np.asarray(I)
    R['J_v'] = np.asarray(J)
    R['fixture_name'] = fixture_name
    R['trial'] = trial
    R['resolution'] = resolution
    return R


def draw_domain_box(ax, domain):
    if domain.type == "rect":
        u_min, u_max, v_min, v_max = domain.bounds
        ax.add_patch(plt.Rectangle((u_min, v_min), u_max-u_min, v_max-v_min,
                                   fill=False, ec="gray", lw=0.5))


def draw_outline_uv(ax, R, alpha=0.25):
    """Draw all SubCurves in UV at low alpha as background context."""
    splits = R['splits']; mesh = R['mesh']; css = R['css']
    sis_pairs = R['sis_pairs']; cps = R['cps']; dps = R['dps']
    domain = R['surf'].domain
    for sub in R['subs']:
        col = _kind_color(sub.kind)
        uv_xy = _rc_uv(sub, splits, mesh, css, sis_pairs, cps, dps)
        if uv_xy is None or len(uv_xy) < 2:
            continue
        crossings = _detect_seam_crossings(uv_xy, domain)
        seg_starts = [0] + [c + 1 for c in crossings] + [len(uv_xy)]
        for s, e in zip(seg_starts[:-1], seg_starts[1:]):
            if e - s >= 2:
                ax.plot(uv_xy[s:e, 0], uv_xy[s:e, 1], "-", c=col,
                        lw=0.8, alpha=alpha)


def draw_outline_xy(ax, R, alpha=0.25):
    """Draw all ResampledCurves in XY at low alpha as background context."""
    for rc in R['rcs']:
        col = _kind_color(rc.kind)
        ax.plot(rc.xy[:, 0], rc.xy[:, 1], "-", c=col, lw=0.8, alpha=alpha)


def draw_surface_3d(ax, R, n_g=30, alpha=0.12):
    """Faint surface backdrop in the 3D camera-frame panel."""
    surf = R['surf']
    I_v = R['I_v']; J_v = R['J_v']; axis_v = R['axis_v']
    domain = surf.domain
    if domain.type != "rect":
        return  # disk: skip surface backdrop for now
    u_min, u_max, v_min, v_max = domain.bounds
    u_g = np.linspace(u_min, u_max, n_g); v_g = np.linspace(v_min, v_max, n_g)
    UG, VG = np.meshgrid(u_g, v_g, indexing="ij")
    S_grid = np.zeros((n_g, n_g, 3))
    for ii in range(n_g):
        for jj in range(n_g):
            S_grid[ii, jj] = surf.S(UG[ii, jj], VG[ii, jj])
    Xc = (S_grid * I_v).sum(axis=-1)
    Yc = (S_grid * J_v).sum(axis=-1)
    Zc = (S_grid * axis_v).sum(axis=-1)
    ax.plot_surface(Xc, Yc, Zc, alpha=alpha, color="gray", edgecolor="none")


def save_two_panels(fig_ctx, fig_xy, R, probe_label):
    base = (f"probe_{probe_label}_{R['fixture_name']}_"
            f"{R['trial']}_res{R['resolution']}")
    ctx_path = f"{base}_context.png"
    xy_path = f"{base}_xy.png"
    fig_ctx.suptitle(
        f"probe_{probe_label} — {R['fixture_name']}  "
        f"{R['view_label']}  res={R['resolution']}")
    fig_ctx.tight_layout(); fig_xy.tight_layout()
    fig_ctx.savefig(ctx_path, dpi=110)
    fig_xy.savefig(xy_path, dpi=110)
    print(f"\nSaved {ctx_path}  and  {xy_path}")
    plt.close(fig_ctx); plt.close(fig_xy)


def project_xyz(xyz, R):
    """Project Nx3 xyz into NxC camera frame (x, y, depth)."""
    I_v = R['I_v']; J_v = R['J_v']; axis_v = R['axis_v']
    return np.column_stack([xyz @ I_v, xyz @ J_v, xyz @ axis_v])


def uv_extent(domain, R):
    """Characteristic uv length for arrow/marker scales."""
    if domain.type == "rect":
        u_min, u_max, v_min, v_max = domain.bounds
        return max(u_max - u_min, v_max - v_min)
    # disk/annulus: use mesh uv extent
    uv = R['mesh'].uv
    return float(max(np.ptp(uv[:, 0]), np.ptp(uv[:, 1])))
