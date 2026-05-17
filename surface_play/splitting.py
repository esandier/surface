"""SP/SPT split-point primitives (Layer O, step O5).

This module is pure plumbing — no geometry. It defines:

- `sp_dtype` / `spt_dtype`: the structured dtypes for split points and splits.
- `SplitArrays`: an append-only buffer that hands out indices and writes them
  into a segment's `split1` / `split2` slot.
- `SplitSlotOverflowError`: raised by `attach_to_segment` when a third SPT is
  attached to a segment that already has both slots filled (G17).

Layer-O steps O6-O11 are the consumers; they decide *what* to split and *why*.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from surface_play.contour import ContourCurve
    from surface_play.curves import BoundaryCurve
    from surface_play.mesh import Mesh
    from surface_play.projection import Projection
    from surface_play.surface import SurfaceParams


sp_dtype = np.dtype([
    ("uv",   "2f8"),
    ("xyz",  "3f8"),
    ("xy",   "2f8"),
    ("type", "U3"),
])

spt_dtype = np.dtype([
    ("sp_idx",   "i4"),
    ("bary",     "f8"),
    ("vis_chge", "i1"),
])


class SplitSlotOverflowError(RuntimeError):
    """G17 violation: a third SPT was attached to a segment with both slots filled."""
    pass


class SplitArrays:
    """Append-only buffer of SPs and SPTs."""

    def __init__(self) -> None:
        self.sps: list = []
        self.spts: list = []

    def add_sp(self, uv, xyz, xy, sp_type) -> int:
        self.sps.append((tuple(uv), tuple(xyz), tuple(xy), sp_type))
        return len(self.sps) - 1

    def add_spt(self, sp_idx, bary, vis_chge) -> int:
        self.spts.append((int(sp_idx), float(bary), int(vis_chge)))
        return len(self.spts) - 1

    def attach_to_segment(self, segment_array, seg_idx, spt_idx,
                          segment_label: str | None = None) -> None:
        seg = segment_array[seg_idx]
        s1 = int(seg["split1"])
        s2 = int(seg["split2"])

        if s1 == -1:
            segment_array[seg_idx]["split1"] = spt_idx
            return
        if s2 == -1:
            segment_array[seg_idx]["split2"] = spt_idx
            return

        # G17 — both slots full.
        kind = segment_label if segment_label is not None else str(segment_array.dtype.names)
        existing_types = (
            self.sps[self.spts[s1][0]][3],
            self.sps[self.spts[s2][0]][3],
        )
        new_type = self.sps[self.spts[spt_idx][0]][3]
        raise SplitSlotOverflowError(
            f"G17: segment kind={kind} seg_idx={seg_idx} already has "
            f"split1/split2 filled with SP types {existing_types}; "
            f"refused to attach a third SPT of type {new_type!r}."
        )

    def sps_array(self) -> np.ndarray:
        return np.array(self.sps, dtype=sp_dtype)

    def spts_array(self) -> np.ndarray:
        return np.array(self.spts, dtype=spt_dtype)


def split_bcs_at_corners(
    mesh: "Mesh",
    bcs: "list[BoundaryCurve]",
    splits: SplitArrays,
    projection: "Projection",
) -> None:
    """Create corner SPs and attach SPTs to incident boundary edges.

    Rect no-id only: in any other domain (disk, annulus, or any identification),
    corners are either absent or are interior points of closed boundary loops,
    and no split is needed (spec line 48).

    Mutates `mesh.edges` (writes split1/split2) and `splits` (appends SPs/SPTs).
    `bcs` is unused — corners are located via `mesh.corner_idx`; we keep the
    parameter to match the roadmap API and for future consistency with O7/O8.
    """
    del bcs  # not needed: corners come from mesh.corner_idx, BEs from mesh.edges.

    domain = mesh.domain
    if (
        domain.type != "rect"
        or domain.u_identify != "no"
        or domain.v_identify != "no"
    ):
        return

    if len(mesh.corner_idx) == 0:
        return

    bnd_idx = mesh.boundary_edge_idx
    bnd_edges = mesh.edges[bnd_idx]
    p_idx_b = bnd_edges["p_idx"]
    q_idx_b = bnd_edges["q_idx"]

    for c in mesh.corner_idx:
        c_int = int(c)
        uv = mesh.uv[c_int]
        xyz = mesh.xyz[c_int]
        xy = projection.XY(xyz)

        sp_idx = splits.add_sp(uv=uv, xyz=xyz, xy=xy, sp_type="cn")

        # Incident boundary edges: those whose p_idx or q_idx equals this corner.
        p_mask = p_idx_b == c_int
        q_mask = q_idx_b == c_int
        incident = np.nonzero(p_mask | q_mask)[0]

        for k in incident:
            e_idx = int(bnd_idx[int(k)])
            bary = 0.0 if bool(p_mask[k]) else 1.0
            spt_idx = splits.add_spt(sp_idx=sp_idx, bary=bary, vis_chge=0)
            splits.attach_to_segment(mesh.edges, e_idx, spt_idx, segment_label="BE")


# ── O7 ────────────────────────────────────────────────────────────────────────

def _bvis_chge(
    surface: "SurfaceParams",
    projection: "Projection",
    edge: np.void,
    s: float,
) -> int:
    """Visibility change on a boundary edge at parametric position s ∈ (0, 1).

    Translation of the spec pseudocode (Reorganization §"Splitting of BCs at
    Contour Points (CPs)", lines 275-296) onto the current refactor:

      - ker_param(p0)     ↔ spec's `kerdS(*p0)[0]` (2D kernel direction in uv).
      - k0·Su + k1·Sv     ↔ `dS(*p0, *ker)` (3D image of the kernel direction).
      - axis · v3         ↔ `Z(v3)` for a tangent vector (depth along view axis;
                            see legacy silhouette.py:434 — the projection.Z
                            method here subtracts an anchor which is wrong for
                            tangents, so we inline the inner product).
      - (axis·Su, axis·Sv) ↔ `SD_jac(*p0, *axis)` — 2D Jacobian in uv of (axis·S).

    Nb sign convention: spec uses `Nb = e["dir"]` directly, and the current
    mesh stores `dir` as the **inward** boundary normal (see _build_edges_faces
    in mesh.py — `d` is flipped to point toward the boundary triangle's apex).
    The legacy mesh (old stuff/silhouette.py:569) also stored an `inward`
    vector under the same name, so the spec formula's `ker · Nb > 0` test is
    the inward-normal version. No flip needed.

    Perspective is not yet supported (the `axis · v3` shortcut for depth is
    ortho-specific; perspective needs a per-point view direction `S - eye`).
    """
    if projection.mode != "ortho":
        raise NotImplementedError("bvis_chge: perspective mode not yet supported")

    p = np.asarray(edge["p"], dtype=float).reshape(2)
    dp = np.asarray(edge["pq"], dtype=float).reshape(2)
    Nb = np.asarray(edge["dir"], dtype=float).reshape(2)

    p0 = p + float(s) * dp
    u, v = float(p0[0]), float(p0[1])

    ker = projection.ker_param(p0)

    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    axis = projection._axis

    dS_ker = ker[0] * Su + ker[1] * Sv
    if float(axis @ dS_ker) < 0.0:
        ker = -ker

    if float(np.dot(ker, Nb)) <= 0.0:
        return 0

    norm_dp = float(np.linalg.norm(dp))
    if norm_dp == 0.0:
        return 0
    Tb = dp / norm_dp

    Np = np.array([float(axis @ Su), float(axis @ Sv)])
    if float(np.dot(Np, ker)) < 0.0:
        Np = -Np

    return 1 if float(np.dot(Tb, Np)) > 0.0 else -1


def split_bcs_at_bcps(
    mesh: "Mesh",
    bcs: "list[BoundaryCurve]",
    ccs: "list[ContourCurve]",
    css: np.ndarray,
    cps: np.ndarray,
    splits: SplitArrays,
    surface: "SurfaceParams",
    projection: "Projection",
) -> None:
    """Split BCs and CSs at every CP with ptype==4 (CP on a boundary edge).

    For each such CP: create one SP of type 'bcp', attach one SPT on the
    unique CS containing the CP (per spec line 270, on the boundary there is
    exactly one), and one SPT on the boundary edge `cp.e`. The CS-side SPT
    has `vis_chge = 0`; the BE-side SPT carries `bvis_chge(e, s)` (G18).

    Mutates `mesh.edges` (split1/split2 on BEs), `css` (split1/split2 on the
    CSs), and `splits`. `bcs` and `ccs` are unused by the algorithm but kept
    in the signature for consistency with O6/O8 and so callers can verify
    sub-assert 4 (every boundary CP is an open-CC endpoint) downstream.
    """
    del bcs, ccs  # API parity only; algorithm uses cps[i]["e"] + css scan.

    if len(cps) == 0:
        return

    bcp_idx = np.nonzero(cps["ptype"] == 4)[0]
    if len(bcp_idx) == 0:
        return

    cs_p = css["p_cp"]
    cs_q = css["q_cp"]

    for i in bcp_idx:
        cp = cps[i]
        uv = np.asarray(cp["uv"], dtype=float)
        xyz = np.asarray(cp["xyz"], dtype=float)
        xy = projection.XY(xyz)
        sp_idx = splits.add_sp(uv=uv, xyz=xyz, xy=xy, sp_type="bcp")

        on_p = cs_p == int(i)
        on_q = cs_q == int(i)
        hits = np.nonzero(on_p | on_q)[0]
        if len(hits) != 1:
            raise RuntimeError(
                f"O7: expected exactly one CS containing boundary CP {int(i)}, "
                f"found {len(hits)} (spec line 270). hits={hits.tolist()}"
            )
        cs_idx = int(hits[0])
        bary_cs = 0.0 if bool(on_p[cs_idx]) else 1.0
        spt_cs = splits.add_spt(sp_idx=sp_idx, bary=bary_cs, vis_chge=0)
        splits.attach_to_segment(css, cs_idx, spt_cs, segment_label="CS")

        e_idx = int(cp["e"])
        bary_be = float(cp["s"])
        vis = _bvis_chge(surface, projection, mesh.edges[e_idx], bary_be)
        spt_be = splits.add_spt(sp_idx=sp_idx, bary=bary_be, vis_chge=vis)
        splits.attach_to_segment(mesh.edges, e_idx, spt_be, segment_label="BE")


# ── O8 ────────────────────────────────────────────────────────────────────────

def _bdp_vis_chge_case2(
    surface: "SurfaceParams",
    projection: "Projection",
    edge: np.void,
    s: float,
    other_uv: np.ndarray,
) -> int:
    """Case 2: BDP's other sheet is interior (face or non-boundary edge).

    Spec line 312: vis_chge = +1 if `inner(SN', axis) * inner(Tb, SN') > 0`
    else -1. SN' is the OTHER sheet's 3D normal at the DP (evaluated at
    `other_uv`), Tb is the 3D tangent to this BE in the p→q direction.
    """
    if projection.mode != "ortho":
        raise NotImplementedError("_bdp_vis_chge: perspective not yet supported")

    p = np.asarray(edge["p"], dtype=float).reshape(2)
    pq = np.asarray(edge["pq"], dtype=float).reshape(2)
    uv = p + float(s) * pq
    u, v = float(uv[0]), float(uv[1])
    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    Tb = pq[0] * Su + pq[1] * Sv

    ou, ov = float(other_uv[0]), float(other_uv[1])
    SN_other = np.asarray(surface.SN(ou, ov), dtype=float).reshape(3)

    axis = projection._axis
    factor = float(axis @ SN_other) * float(Tb @ SN_other)
    return 1 if factor > 0.0 else -1


def _bdp_vis_chge_case1(
    surface: "SurfaceParams",
    projection: "Projection",
    edge_self: np.void,
    s_self: float,
    edge_other: np.void,
    s_other: float,
    other_uv: np.ndarray,
) -> int:
    """Case 1: BDP's other sheet is also on a boundary edge.

    Spec lines 311 (with `dir'` interpreted per legacy splitting.py:401 — the
    perpendicular to the OTHER edge's projected tangent, oriented toward the
    surface via the OTHER edge's inward 2D normal projected to screen).

    Lifted: a = `inner(SN', axis) * inner(Tb, SN')`,
            b = `inner(Tb_proj, dir')`.
    If `a * b > 0`: vis_chge = 0.
    Else:           vis_chge = +1 if a > 0 else -1.
    """
    if projection.mode != "ortho":
        raise NotImplementedError("_bdp_vis_chge: perspective not yet supported")

    # OUR edge geometry.
    p = np.asarray(edge_self["p"], dtype=float).reshape(2)
    pq = np.asarray(edge_self["pq"], dtype=float).reshape(2)
    uv = p + float(s_self) * pq
    u, v = float(uv[0]), float(uv[1])
    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    Tb = pq[0] * Su + pq[1] * Sv
    Tb_proj = projection.proj_vec(uv, Tb)

    # OTHER sheet's normal at the DP.
    ou, ov = float(other_uv[0]), float(other_uv[1])
    SN_other = np.asarray(surface.SN(ou, ov), dtype=float).reshape(3)

    # OTHER edge's projected tangent and inward normal at the DP.
    op = np.asarray(edge_other["p"], dtype=float).reshape(2)
    opq = np.asarray(edge_other["pq"], dtype=float).reshape(2)
    odir = np.asarray(edge_other["dir"], dtype=float).reshape(2)
    o_uv = op + float(s_other) * opq
    ou2, ov2 = float(o_uv[0]), float(o_uv[1])
    Su_o = np.asarray(surface.Su(ou2, ov2), dtype=float).reshape(3)
    Sv_o = np.asarray(surface.Sv(ou2, ov2), dtype=float).reshape(3)
    T_other_3d = opq[0] * Su_o + opq[1] * Sv_o
    T_other_proj = projection.proj_vec(o_uv, T_other_3d)
    inward_3d = odir[0] * Su_o + odir[1] * Sv_o
    inward_proj = projection.proj_vec(o_uv, inward_3d)

    # dir' = perpendicular to OTHER edge's projected tangent, oriented toward
    # surface via inward_proj.
    dir_perp = np.array([-T_other_proj[1], T_other_proj[0]], dtype=float)
    if float(np.dot(dir_perp, inward_proj)) < 0.0:
        dir_perp = -dir_perp

    axis = projection._axis
    a = float(axis @ SN_other) * float(Tb @ SN_other)
    b = float(np.dot(Tb_proj, dir_perp))

    if a * b > 0.0:
        return 0
    return 1 if a > 0.0 else -1


def _s_of_uv_on_edge(uv_on_edge: np.ndarray, edge: np.void) -> float:
    """Parametric position s ∈ [0, 1] of `uv_on_edge` on `edge`, via projection
    onto `edge["pq"]`. Used to convert a DP's uv1/uv2 into a bary for the BE.
    """
    p = np.asarray(edge["p"], dtype=float).reshape(2)
    pq = np.asarray(edge["pq"], dtype=float).reshape(2)
    d = np.asarray(uv_on_edge, dtype=float).reshape(2) - p
    denom = float(np.dot(pq, pq))
    if denom == 0.0:
        return 0.5
    return float(np.dot(d, pq) / denom)


def split_bcs_at_bdps(
    mesh: "Mesh",
    bcs: "list[BoundaryCurve]",
    sics: list,
    sis_pairs: np.ndarray,
    dps: np.ndarray,
    splits: SplitArrays,
    surface: "SurfaceParams",
    projection: "Projection",
) -> None:
    """Split BCs and SISs at every DP with on_boundary=True.

    For each BDP: create one SP of type 'bdp', attach one SPT on the unique
    SIS containing the DP (vis_chge=0), and one SPT per incident BE
    (`E1` always, plus `E2` if EE-type DP with both edges boundary). The
    BE-side vis_chge follows spec case 1 (two boundary edges → vis ∈ {-1, 0, +1})
    or case 2 (one boundary edge → vis ∈ {-1, +1}).

    `bcs` and `sics` are unused by the algorithm (incident BEs come from
    `dp.E1/E2` + `mesh.edges["g"]`; the SIS is found by scanning sis_pairs);
    kept in the signature for API parity with the roadmap.
    """
    del bcs, sics

    if len(dps) == 0:
        return

    bdp_idx = np.nonzero(dps["on_boundary"])[0]
    if len(bdp_idx) == 0:
        return

    sis_p = sis_pairs["p_dp"] if len(sis_pairs) else np.empty(0, dtype=np.int32)
    sis_q = sis_pairs["q_dp"] if len(sis_pairs) else np.empty(0, dtype=np.int32)

    edges_g = mesh.edges["g"]

    for i in bdp_idx:
        dp = dps[i]
        xyz = np.asarray(dp["xyz"], dtype=float)
        uv = np.asarray(dp["uv1"], dtype=float)  # sheet-1 uv as the SP's anchor
        xy = projection.XY(xyz)
        sp_idx = splits.add_sp(uv=uv, xyz=xyz, xy=xy, sp_type="bdp")

        # SIS containing this DP (spec line 305: exactly one on the boundary).
        if len(sis_pairs) > 0:
            on_p = sis_p == int(i)
            on_q = sis_q == int(i)
            sis_hits = np.flatnonzero(on_p | on_q)
        else:
            sis_hits = np.empty(0, dtype=np.int32)
        if len(sis_hits) > 1:
            raise RuntimeError(
                f"O8: boundary DP {int(i)} appears in {len(sis_hits)} SISs; "
                f"expected at most 1 (spec line 305)"
            )
        if len(sis_hits) == 1:
            sis_idx = int(sis_hits[0])
            bary_sis = 0.0 if bool(on_p[sis_idx]) else 1.0
            spt_sis = splits.add_spt(sp_idx=sp_idx, bary=bary_sis, vis_chge=0)
            splits.attach_to_segment(sis_pairs, sis_idx, spt_sis,
                                     segment_label="SIS")

        # Identify incident boundary edges (E1 always; E2 if EE-type).
        e1_idx = int(dp["E1"])
        e1_is_bnd = e1_idx >= 0 and int(edges_g[e1_idx]) == -1
        e2_idx = int(dp["E2"])
        e2_is_bnd = (
            str(dp["type"]) == "EE"
            and e2_idx >= 0
            and int(edges_g[e2_idx]) == -1
        )

        incident: list[tuple[int, np.ndarray, np.ndarray]] = []
        if e1_is_bnd:
            incident.append((
                e1_idx,
                np.asarray(dp["uv1"], dtype=float),
                np.asarray(dp["uv2"], dtype=float),
            ))
        if e2_is_bnd:
            incident.append((
                e2_idx,
                np.asarray(dp["uv2"], dtype=float),
                np.asarray(dp["uv1"], dtype=float),
            ))
        if not incident:
            continue

        case1 = len(incident) == 2

        for k, (e_idx, this_uv, other_uv) in enumerate(incident):
            edge = mesh.edges[e_idx]
            s_here = _s_of_uv_on_edge(this_uv, edge)
            if case1:
                # The "other" incident is the other entry.
                o_e_idx, o_this_uv, _ = incident[1 - k]
                other_edge = mesh.edges[o_e_idx]
                s_other = _s_of_uv_on_edge(o_this_uv, other_edge)
                vis = _bdp_vis_chge_case1(
                    surface, projection, edge, s_here,
                    other_edge, s_other, other_uv,
                )
            else:
                vis = _bdp_vis_chge_case2(
                    surface, projection, edge, s_here, other_uv,
                )

            spt_be = splits.add_spt(sp_idx=sp_idx, bary=s_here, vis_chge=vis)
            splits.attach_to_segment(mesh.edges, e_idx, spt_be,
                                     segment_label="BE")
