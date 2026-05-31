"""Visual probe for Layer O visibility pipeline.

Usage:
    python -m surface_play.probe_visibility [FIXTURE] [TRIAL] [RES]

    FIXTURE ∈ {helicoid, fig8, paraboloid, torus, mobius_u, mobius_v,
               disk_paraboloid_po, disk_paraboloid_ca}
    TRIAL   = integer seed for a deterministic random axis (default 0)
    RES     = mesh resolution (default 100)

If TRIAL == "canonical" the canonical viewpoint from test_fixtures is used.

For each run:
- builds the full Layer-O pipeline (mesh → contours → splits → helpers →
  SubCurves → resample → breaks → BFS → LP).
- prints a one-page summary of RCs, SPs, breaks, BFS/LP ranges.
- saves a 4-panel PNG (UV domain, XY projection, camera-frame 3D, text).

Per the Layer O sign-off gate (roadmap line 1513-1518), user reviews this
probe's XY-projection panel on helicoid, fig-8 cyl, and Möbius and
confirms the visibility distribution looks correct.
"""

import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from surface_play.contour import (build_contour_curves, build_contour_segments,
    find_contour_points, find_vps)
from surface_play.curves import (
    build_bcs, resample_all, _seg_uv_at_bary, _sic_vertices,
)
from surface_play.helpers import build_helper_curves
from surface_play.intersections import (build_sics, build_sis_pairs,
    find_double_points, find_triple_points, tp_dtype)
from surface_play.mesh import build_mesh
from surface_play.projection import Projection
from surface_play.splitting import (SplitArrays, assemble_subcurves, split_at_cdps,
    split_bcs_at_bcps, split_bcs_at_bdps, split_bcs_at_corners,
    split_ccs_at_vps, split_sics_at_tps)
from surface_play import test_fixtures as tf
from surface_play.visibility import (bfs_visibility, compute_projection_breaks,
    lp_refine_visibility, LPInfeasibleError, _pick_anchors)


FIXTURE_BUILDERS = {
    "helicoid":             (tf.helicoid,             None),
    "paraboloid":           (tf.paraboloid,           tf.paraboloid_side_view),
    "torus":                (tf.torus,                tf.torus_ortho_view),
    "mobius_u":             (tf.mobius_u,             None),
    "mobius_v":             (tf.mobius_v,             None),
    "fig8":                 (tf.fig8,                 tf.fig8_top_view),
    "disk_paraboloid_po":   (tf.disk_paraboloid_po,   None),
    "disk_paraboloid_ca":   (tf.disk_paraboloid_ca,   None),
}


def _frame_from_axis(a):
    a = np.asarray(a, dtype=float); a /= np.linalg.norm(a)
    helper = np.array([1.0, 0.0, 0.0])
    if abs(a @ helper) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    I = helper - (helper @ a) * a; I /= np.linalg.norm(I)
    J = np.cross(a, I)
    return a, I, J


def _random_axis(trial: int):
    rng = np.random.default_rng(2026 + trial)
    a = rng.standard_normal(3)
    return _frame_from_axis(a)


def _build_outline(surf, I, J, resolution, seed=42):
    mesh = build_mesh(surf.domain, surf, resolution=resolution, jitter=True, seed=seed)
    proj = Projection(surf, I=I.tolist() if isinstance(I, np.ndarray) else I,
                            J=J.tolist() if isinstance(J, np.ndarray) else J)
    bcs = build_bcs(mesh)
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    vps = find_vps(ccs, css, cps, surf, proj)
    dps = find_double_points(mesh, surf)
    sis_pairs = build_sis_pairs(dps)
    sics = build_sics(sis_pairs) if len(sis_pairs) else []
    tps = (find_triple_points(sis_pairs, dps, mesh, surf)
           if len(sis_pairs) >= 3 else np.empty(0, dtype=tp_dtype))
    splits = SplitArrays()
    split_bcs_at_corners(mesh, bcs, splits, proj)
    split_bcs_at_bcps(mesh, bcs, ccs, css, cps, splits, surf, proj)
    split_bcs_at_bdps(mesh, bcs, sics, sis_pairs, dps, splits, surf, proj)
    split_sics_at_tps(tps, sis_pairs, dps, splits, surf, proj)
    split_ccs_at_vps(vps, css, ccs, splits, surf, proj)
    split_at_cdps(mesh, css, cps, sis_pairs, dps, splits, surf, proj)
    hcs = build_helper_curves(bcs, ccs, sics, css, sis_pairs, cps, dps,
                              mesh, splits, proj, surf, surf.domain)
    subs = assemble_subcurves(bcs, ccs, sics, hcs, mesh, css, sis_pairs, splits)
    rcs = resample_all(subs, surf, proj, splits, mesh, css, sis_pairs, cps, dps)
    breaks = compute_projection_breaks(rcs, surf, proj)
    bfs = bfs_visibility(rcs, breaks, splits)
    try:
        lp = lp_refine_visibility(rcs, breaks, splits, vis_bfs=bfs)
        lp_status = "OK"
    except LPInfeasibleError:
        lp = None; lp_status = "INFEASIBLE"
    return dict(mesh=mesh, proj=proj, subs=subs, rcs=rcs, breaks=breaks,
                splits=splits, bfs=bfs, lp=lp, lp_status=lp_status,
                sps=splits.sps_array(), css=css, sis_pairs=sis_pairs,
                cps=cps, dps=dps, ccs=ccs, vps=vps, hcs=hcs, sics=sics,
                tps=tps, bcs=bcs, surf=surf)


def _rc_uv(sub, splits, mesh, css, sis_pairs, cps, dps):
    """uv polyline for a sub: [start SP] + internal samples + [end SP].

    HC has no internal samples and only 2 SPs. SIC is one curve of DPs: its
    node polyline is drawn in the canonical (uv1 / SP-uv) preimage, matching
    `resample_all`'s skeleton (spec §SIC). The canonical preimage may switch
    sheets at a DP where `flip` alternates — that is faithful, not a bug.
    """
    if sub.kind == "SIC":
        verts = _sic_vertices(sub, sis_pairs)
        if sub.is_closed and verts:
            verts = verts + [verts[0]]
        uvs = [
            np.asarray(dps[int(v[1])]["uv1"], dtype=float) if v[0] == "DP"
            else np.asarray(splits.sps[int(v[1])][0], dtype=float)
            for v in verts
        ]
        return np.array(uvs) if uvs else None

    uvs = []
    if sub.start >= 0:
        uvs.append(np.asarray(splits.sps[sub.start][0], dtype=float))
    for entry in sub.internal:
        if len(entry) >= 2:
            try:
                uvs.append(_seg_uv_at_bary(sub.kind, int(entry[0]), float(entry[1]),
                                           mesh, css, sis_pairs, cps, dps))
            except Exception:
                pass
    if sub.end >= 0 and not sub.is_closed:
        uvs.append(np.asarray(splits.sps[sub.end][0], dtype=float))
    elif sub.is_closed and uvs:
        uvs.append(uvs[0].copy())
    return np.array(uvs) if uvs else None


def _detect_seam_crossings(uv_xy, domain):
    if uv_xy is None or len(uv_xy) < 2:
        return []
    crossings = []
    pu = domain.period_u if not np.isnan(domain.period_u) else 0.0
    pv = domain.period_v if not np.isnan(domain.period_v) else 0.0
    for k in range(len(uv_xy) - 1):
        du = abs(uv_xy[k+1, 0] - uv_xy[k, 0])
        dv = abs(uv_xy[k+1, 1] - uv_xy[k, 1])
        if (pu > 0 and du > 0.5 * pu) or (pv > 0 and dv > 0.5 * pv):
            crossings.append(k)
    return crossings


CMAP_VIS = {-4: "#1a0000", -3: "#440000", -2: "#aa0000", -1: "#ff5500",
            0: "#0a8a0a", 1: "#0066ff", 2: "#00aaff", 3: "#88aaff", 4: "#bbddff"}


def _kind_color(kind):
    return {"BC": "tab:blue", "CC": "tab:orange",
            "HC": "tab:green", "SIC": "tab:red"}.get(kind, "gray")


def probe(fixture_name, trial, resolution, out_png=None):
    if fixture_name not in FIXTURE_BUILDERS:
        raise SystemExit(f"unknown fixture {fixture_name!r}; "
                         f"choose from {list(FIXTURE_BUILDERS)}")
    builder, canonical = FIXTURE_BUILDERS[fixture_name]
    surf = builder(perturb=False)

    if trial == "canonical":
        if canonical is None:
            raise SystemExit(f"no canonical view registered for {fixture_name}; "
                             "pass an integer trial instead")
        I_in, J_in, eye = canonical
        I = np.asarray(I_in, dtype=float); J = np.asarray(J_in, dtype=float)
        a = np.cross(I, J)  # axis = I × J, points toward viewer
        view_label = "canonical"
    else:
        t = int(trial)
        a, I, J = _random_axis(t)
        view_label = f"trial {t} (seed=2026+{t})"

    R = _build_outline(surf, I, J, resolution)

    # ── Console summary ─────────────────────────────────────────────────────
    rcs, subs, splits = R['rcs'], R['subs'], R['splits']
    bfs = R['bfs']; lp = R['lp']; breaks = R['breaks']
    sps = R['sps']

    bfs_flat = np.concatenate([bfs[id(rc)] for rc in rcs]) if rcs else np.array([0])
    bfs_lo, bfs_hi = int(bfs_flat.min()), int(bfs_flat.max())
    if lp is not None:
        lp_flat = np.concatenate([lp[id(rc)] for rc in rcs])
        lp_lo, lp_hi = int(lp_flat.min()), int(lp_flat.max())
    else:
        lp_lo = lp_hi = None
    n_diff = (sum(1 for rc in rcs if np.any(bfs[id(rc)] != lp[id(rc)]))
              if lp is not None else None)

    print(f"\n=== probe_visibility: {fixture_name}  view={view_label}  res={resolution} ===")
    print(f"axis = {a.round(4).tolist()}")
    print(f"n_subs={len(subs)} n_rcs={len(rcs)} n_breaks={len(breaks)} "
          f"n_sps={len(sps)} n_cps={len(R['cps'])} n_dps={len(R['dps'])}")
    print(f"BFS range = [{bfs_lo:+d},{bfs_hi:+d}]  LP={R['lp_status']}", end="")
    if lp is not None:
        print(f"  LP range = [{lp_lo:+d},{lp_hi:+d}]  RCs with LP!=BFS = {n_diff}")
    else:
        print()

    print(f"\nRCs:")
    for ri, rc in enumerate(rcs):
        v = bfs[id(rc)]
        lp_str = ""
        if lp is not None:
            lv = lp[id(rc)]
            lp_str = f"  lp[{int(lv.min()):+d},{int(lv.max()):+d}]"
        print(f"  rc[{ri}] {rc.kind} N={len(rc.xy):>4} {rc.start}->{rc.end}  "
              f"vc=({int(rc.vc_in):+d},{int(rc.vc_out):+d})  "
              f"bfs[{int(v.min()):+d},{int(v.max()):+d}]{lp_str}")

    if len(breaks):
        print(f"\nBreaks ({len(breaks)}):")
        for bi, b in enumerate(breaks):
            print(f"  bk[{bi}] rc[{int(b['rc_idx'])}].si={int(b['sample_idx'])} "
                  f"dv={int(b['delta_v']):+d}  xy=({b['xy'][0]:+.3f},{b['xy'][1]:+.3f})")

    # ── Plot ────────────────────────────────────────────────────────────────
    mesh = R['mesh']; proj = R['proj']
    css = R['css']; sis_pairs = R['sis_pairs']
    cps = R['cps']; dps = R['dps']
    axis_v = np.asarray(a); I_v = np.asarray(I); J_v = np.asarray(J)
    domain = surf.domain
    rcs_arr = rcs

    # 3-panel context figure (UV, 3D, info), plus standalone wide XY figure.
    fig = plt.figure(figsize=(18, 9))

    # Panel 1: UV domain
    ax_uv = fig.add_subplot(1, 3, 1)
    ax_uv.set_title(f"UV domain ({domain.type}"
                    f", u_id={domain.u_identify}, v_id={domain.v_identify})")
    if domain.type == "rect":
        u_min, u_max, v_min, v_max = domain.bounds
        ax_uv.add_patch(plt.Rectangle((u_min, v_min), u_max-u_min, v_max-v_min,
                                      fill=False, ec="gray", lw=0.5))
    for ri, sub in enumerate(subs):
        col = _kind_color(sub.kind)
        uv_xy = _rc_uv(sub, splits, mesh, css, sis_pairs, cps, dps)
        if uv_xy is None or len(uv_xy) < 2: continue
        crossings = _detect_seam_crossings(uv_xy, domain)
        seg_starts = [0] + [c + 1 for c in crossings] + [len(uv_xy)]
        for s, e in zip(seg_starts[:-1], seg_starts[1:]):
            if e - s >= 2:
                ax_uv.plot(uv_xy[s:e, 0], uv_xy[s:e, 1], "-", c=col, lw=1.2)
        for c in crossings:
            ax_uv.plot(uv_xy[c, 0], uv_xy[c, 1], "rx", ms=6, mew=1.5)
    for sp_i in range(len(sps)):
        uv = sps[sp_i]["uv"]
        ax_uv.plot(uv[0], uv[1], "ko", ms=4)
    ax_uv.set_xlabel("u"); ax_uv.set_ylabel("v")
    ax_uv.grid(True, alpha=0.3)

    # Panel 2: 3D camera-frame
    ax_3d = fig.add_subplot(1, 3, 2, projection="3d")
    ax_3d.set_title("3D camera-frame")
    if domain.type == "rect":
        u_min, u_max, v_min, v_max = domain.bounds
        n_g = 30
        u_g = np.linspace(u_min, u_max, n_g); v_g = np.linspace(v_min, v_max, n_g)
        UG, VG = np.meshgrid(u_g, v_g, indexing="ij")
        S_grid = np.zeros((n_g, n_g, 3))
        for ii in range(n_g):
            for jj in range(n_g):
                S_grid[ii, jj] = surf.S(UG[ii, jj], VG[ii, jj])
        Xc = (S_grid * I_v).sum(axis=-1)
        Yc = (S_grid * J_v).sum(axis=-1)
        Zc = (S_grid * axis_v).sum(axis=-1)
        ax_3d.plot_surface(Xc, Yc, Zc, alpha=0.15, color="gray", edgecolor="none")
    for ri, sub in enumerate(subs):
        col = _kind_color(sub.kind)
        uv_xy = _rc_uv(sub, splits, mesh, css, sis_pairs, cps, dps)
        if uv_xy is None: continue
        S3 = np.array([surf.S(float(p[0]), float(p[1])) for p in uv_xy])
        Xl = S3 @ I_v; Yl = S3 @ J_v; Zl = S3 @ axis_v
        ax_3d.plot(Xl, Yl, Zl, "-", c=col, lw=2.0 if sub.kind == "CC" else 1.2)
    ax_3d.view_init(elev=85, azim=-90)
    ax_3d.set_xlabel("I·S"); ax_3d.set_ylabel("J·S"); ax_3d.set_zlabel("axis·S")

    # Panel 3: text summary
    ax_t = fig.add_subplot(1, 3, 3); ax_t.axis("off")
    ax_t.set_title("Summary")
    lines = [f"fixture: {fixture_name}",
             f"view:    {view_label}",
             f"res:     {resolution}",
             f"axis:    {axis_v.round(3).tolist()}",
             "",
             f"BFS = [{bfs_lo:+d},{bfs_hi:+d}]"]
    if lp is not None:
        lines.append(f"LP  = [{lp_lo:+d},{lp_hi:+d}]   (status: {R['lp_status']})")
        lines.append(f"RCs with LP != BFS: {n_diff}")
    else:
        lines.append(f"LP  = {R['lp_status']}")
    lines.append("")
    lines.append(f"n_subs   = {len(subs)}")
    lines.append(f"n_rcs    = {len(rcs)}")
    lines.append(f"n_sps    = {len(sps)}")
    lines.append(f"n_breaks = {len(breaks)}")
    lines.append(f"n_cps    = {len(cps)}")
    lines.append(f"n_dps    = {len(dps)}")
    lines.append("")
    lines.append("RCs (kind, N, vc, bfs):")
    for ri, rc in enumerate(rcs):
        v = bfs[id(rc)]
        lp_str = ""
        if lp is not None:
            lv = lp[id(rc)]
            lp_str = f" lp[{int(lv.min()):+d},{int(lv.max()):+d}]"
        lines.append(f" [{ri:>2}] {rc.kind} N={len(rc.xy):>3} "
                     f"({int(rc.vc_in):+d},{int(rc.vc_out):+d}) "
                     f"bfs[{int(v.min()):+d},{int(v.max()):+d}]{lp_str}")
    ax_t.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=8,
              va="top", transform=ax_t.transAxes)

    # ── XY projection — wide standalone figure ──────────────────────────────
    fig_xy = plt.figure(figsize=(18, 7))
    ax_xy = fig_xy.add_subplot(1, 1, 1)
    title_xy = (f"XY projection — visibility colored "
                f"(thicker = CC, thin = BC/HC; SIC red)\n"
                f"fixture={fixture_name}  view={view_label}  res={resolution}  "
                f"axis={axis_v.round(3).tolist()}  "
                f"BFS=[{bfs_lo:+d},{bfs_hi:+d}]  LP={R['lp_status']}")
    ax_xy.set_title(title_xy)

    vis_source = lp if lp is not None else bfs

    # Pre-build per-RC, per-segment break list: (rc_idx, sample_idx) -> sorted
    # list of (t, delta_v). A break with sample_idx=k applies to segment k → k+1.
    breaks_by_seg: dict[tuple[int, int], list[tuple[float, int]]] = {}
    for b in breaks:
        ri = int(b['rc_idx']); si = int(b['sample_idx'])
        t = float(b['t']) if 't' in b.dtype.names else 0.5
        dv = int(b['delta_v'])
        breaks_by_seg.setdefault((ri, si), []).append((t, dv))
    for key in breaks_by_seg:
        breaks_by_seg[key].sort(key=lambda td: td[0])

    for ri, rc in enumerate(rcs_arr):
        v = vis_source[id(rc)]
        sx = rc.xy[:, 0]; sy = rc.xy[:, 1]
        if rc.kind == "SIC":
            lw = 3.5
        elif rc.kind == "CC":
            lw = 3.0
        elif rc.kind == "BC":
            lw = 1.8
        else:
            lw = 1.0
        for k in range(len(v) - 1):
            # Visibility at sample k is the segment's "pre-break" value.
            # If breaks land inside segment k, split the line at each break.
            if rc.kind == "SIC":
                # SIC RCs are not occluded (drawn solid red), no break splits needed.
                ax_xy.plot(sx[k:k+2], sy[k:k+2], "-", c="tab:red", lw=lw,
                           solid_capstyle="round", solid_joinstyle="round")
                continue
            segment_breaks = breaks_by_seg.get((ri, k), [])
            # Walk fractional positions along the segment, plotting each
            # sub-segment with the visibility valid in that sub-interval.
            curr_vis = int(v[k])
            prev_t = 0.0
            for t, dv in segment_breaks:
                if t <= prev_t:
                    curr_vis += dv  # multiple breaks at the same t — accumulate
                    continue
                col = CMAP_VIS.get(curr_vis, "magenta")
                x0 = sx[k] + prev_t * (sx[k+1] - sx[k])
                y0 = sy[k] + prev_t * (sy[k+1] - sy[k])
                x1 = sx[k] + t * (sx[k+1] - sx[k])
                y1 = sy[k] + t * (sy[k+1] - sy[k])
                ax_xy.plot([x0, x1], [y0, y1], "-", c=col, lw=lw,
                           solid_capstyle="round", solid_joinstyle="round")
                curr_vis += dv
                prev_t = t
            # Final stretch from prev_t to t=1.
            col = CMAP_VIS.get(curr_vis, "magenta")
            x0 = sx[k] + prev_t * (sx[k+1] - sx[k])
            y0 = sy[k] + prev_t * (sy[k+1] - sy[k])
            ax_xy.plot([x0, sx[k+1]], [y0, sy[k+1]], "-", c=col, lw=lw,
                       solid_capstyle="round", solid_joinstyle="round")

    # Break markers — below curve so they don't break visual continuity
    for bi, b in enumerate(breaks):
        rcb = rcs_arr[int(b['rc_idx'])]
        si = int(b['sample_idx'])
        if 0 <= si < len(rcb.xy):
            ax_xy.plot(rcb.xy[si, 0], rcb.xy[si, 1], "m^", ms=5,
                       zorder=2, alpha=0.7)

    # SP markers — small, semi-transparent, below curves so endpoints visibly meet
    for sp_i in range(len(sps)):
        xy = sps[sp_i]["xy"]
        ax_xy.plot(xy[0], xy[1], "ks", ms=3, zorder=2, alpha=0.6)

    # Anchor marker — kept on top so it's findable
    ri_a, si_a = _pick_anchors(rcs_arr, "leftmost")[0]
    xy_a = rcs_arr[ri_a].xy[si_a]
    ax_xy.plot(xy_a[0], xy_a[1], "g*", ms=14, zorder=12, alpha=0.85)

    # Legend
    vis_keys = sorted(set(int(v) for rc in rcs_arr
                          for v in vis_source[id(rc)]))
    vis_legend = [Line2D([0], [0], color=CMAP_VIS.get(v, "magenta"),
                          lw=3, label=f"vis={v:+d}") for v in vis_keys]
    vis_legend += [
        Line2D([0], [0], color="tab:red", lw=3, label="SIC"),
        Line2D([0], [0], marker="s", color="k", lw=0, ms=6, label="SP"),
        Line2D([0], [0], marker="^", color="m", lw=0, ms=8, label="break"),
        Line2D([0], [0], marker="*", color="g", lw=0, ms=12, label="anchor"),
    ]
    ax_xy.legend(handles=vis_legend, loc="upper left", fontsize=8)
    ax_xy.set_xlabel("x"); ax_xy.set_ylabel("y")
    ax_xy.set_aspect("equal")
    ax_xy.grid(True, alpha=0.3)

    fig.suptitle(f"probe_visibility — {fixture_name}  {view_label}  res={resolution}")
    fig.tight_layout()
    fig_xy.tight_layout()

    base = out_png or f"probe_visibility_{fixture_name}_{trial}_res{resolution}"
    ctx_path = f"{base}_context.png"
    xy_path = f"{base}_xy.png"
    fig.savefig(ctx_path, dpi=110)
    fig_xy.savefig(xy_path, dpi=110)
    print(f"\nSaved {ctx_path}  and  {xy_path}")
    plt.close(fig); plt.close(fig_xy)


if __name__ == "__main__":
    fixture = sys.argv[1] if len(sys.argv) > 1 else "mobius_v"
    trial = sys.argv[2] if len(sys.argv) > 2 else "0"
    res = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    probe(fixture, trial, res)
