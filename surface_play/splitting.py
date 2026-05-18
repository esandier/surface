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

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from surface_play.contour import ContourCurve
    from surface_play.curves import BoundaryCurve
    from surface_play.mesh import Mesh
    from surface_play.projection import Projection
    from surface_play.surface import SurfaceParams


@dataclass
class SubCurve:
    """One piece of a curve between two SPs (Split Points).

    Defined here so O12 (HC construction) and O13 (BC/CC/SIC splitting) share
    the same dataclass without a later move. HCs emit `SubCurve(kind="HC",
    internal=[], is_closed=False)` straight from O12 — the two-SP definition
    (spec line 408). O13 will produce kind ∈ {"BC", "CC", "SIC"} by cutting
    parent chains at SPTs, populating `internal` with the segment_idx/bary
    points between the two endpoints.

    G5: `start`/`end` are SP indices into `SplitArrays.sps`. Identity sharing
    is integer equality of indices — adjacent SubCurves at the same SP have
    `==` ints, which lookup the same Python tuple in `splits.sps[idx]`.
    """

    kind: Literal["BC", "CC", "SIC", "HC"]
    is_closed: bool
    start: int                                       # SP index
    end: int                                         # SP index
    internal: list[tuple[int, float]] = field(default_factory=list)
    vc_in: int = 0                                   # ∈ {-1, 0}
    vc_out: int = 0                                  # ∈ {-1, 0}
    parent_idx: int = -1                             # index into bcs/ccs/sics; -1 for HC


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


# ── O9 ────────────────────────────────────────────────────────────────────────

def _sis_preimage_faces(sis_pair_row: np.void, dps: np.ndarray) -> tuple[int, int]:
    """Two faces on which the SIS's preimages live.

    Mirrors `intersections.py:_first_shared_face` usage at lines 800-801 / 805-806
    of `find_triple_points`: the SIS has two preimage segments, one on
    `_first_shared_face(A_p.A1, A_q.A1|A2)` and one on `_first_shared_face(
    A_p.A2, A_q.A2|A1)`, where the |A2/|A1 alternative is selected by `flip`.
    """
    from surface_play.intersections import _first_shared_face

    p = int(sis_pair_row["p_dp"])
    q = int(sis_pair_row["q_dp"])
    f = int(sis_pair_row["flip"])
    f_A = _first_shared_face(
        dps["A1"][p], dps["A1"][q] if f == 1 else dps["A2"][q],
    )
    f_B = _first_shared_face(
        dps["A2"][p], dps["A2"][q] if f == 1 else dps["A1"][q],
    )
    return f_A, f_B


def split_sics_at_tps(
    tps: np.ndarray,
    sis_pairs: np.ndarray,
    dps: np.ndarray,
    splits: SplitArrays,
    surface: "SurfaceParams",
    projection: "Projection",
) -> None:
    """Split SISs at every triple point. One 'tp' SP and 3 SPTs (one per SIS).

    For each TP and each `i ∈ {0, 1, 2}` (local index over `tp["faces"]`):
      - `Fl = tp["faces"][i]`, `Pl = tp["uv"][i]` (preimage on Fl).
      - Identify the SIS "Si" in `tp["sis_indices"]` whose preimages do NOT
        live on Fl (the OTHER two SISs share Fl).
      - `N = SN(Pl)`, flipped so `axis · N > 0` (spec line 345 says `Z(N) > 0`
        but Projection.Z subtracts an anchor and is wrong for normals — see
        the O7 `_bvis_chge` precedent and `resume_O7_DONE.md` point 2).
      - `T_Si = dps[q_dp].xyz - dps[p_dp].xyz` (chord — adequate for sign).
      - `vis_chge = +1 if T_Si · N > 0 else -1` (G10 — TP splits never 0).
      - `bary = ‖TP - p_xyz‖ / ‖q_xyz - p_xyz‖`, clipped to [0, 1].

    Mutates `sis_pairs` (split1/split2 on the 3 involved SISs) and `splits`.
    """
    if projection.mode != "ortho":
        raise NotImplementedError("split_sics_at_tps: perspective not yet supported")

    if len(tps) == 0:
        return

    axis = projection._axis

    for tp in tps:
        tp_xyz = np.asarray(tp["xyz"], dtype=float).reshape(3)
        tp_xy = projection.XY(tp_xyz)
        # SP uv anchor: pick the first preimage uv (any of the three would do;
        # it is only metadata since the TP is a single 3D point).
        sp_uv = np.asarray(tp["uv"][0], dtype=float).reshape(2)
        sp_idx = splits.add_sp(uv=sp_uv, xyz=tp_xyz, xy=tp_xy, sp_type="tp")

        sis_idx_arr = np.asarray(tp["sis_indices"], dtype=np.int64).reshape(3)
        face_arr = np.asarray(tp["faces"], dtype=np.int64).reshape(3)
        uv_arr = np.asarray(tp["uv"], dtype=float).reshape(3, 2)

        # Precompute (f_A, f_B) for each of the 3 SISs.
        sis_pre_faces = [
            _sis_preimage_faces(sis_pairs[int(sis_idx_arr[k])], dps)
            for k in range(3)
        ]

        for i in range(3):
            Fl = int(face_arr[i])
            Pl = uv_arr[i]

            # Find Si: the SIS not preimaging on Fl.
            si_local = -1
            for k in range(3):
                f_A, f_B = sis_pre_faces[k]
                if Fl != f_A and Fl != f_B:
                    si_local = k
                    break
            if si_local < 0:
                raise RuntimeError(
                    f"O9: no SIS in tp.sis_indices={sis_idx_arr.tolist()} avoids "
                    f"face Fl={Fl}; preimage-faces={sis_pre_faces}"
                )
            s_idx = int(sis_idx_arr[si_local])

            # Normal at Pl, oriented so axis · N > 0.
            N = np.asarray(surface.SN(float(Pl[0]), float(Pl[1])),
                           dtype=float).reshape(3)
            if float(axis @ N) < 0.0:
                N = -N

            # Chord tangent of Si.
            row = sis_pairs[s_idx]
            p_dp = int(row["p_dp"])
            q_dp = int(row["q_dp"])
            p_xyz = np.asarray(dps["xyz"][p_dp], dtype=float).reshape(3)
            q_xyz = np.asarray(dps["xyz"][q_dp], dtype=float).reshape(3)
            chord = q_xyz - p_xyz
            chord_norm = float(np.linalg.norm(chord))

            vis_chge = 1 if float(np.dot(chord, N)) > 0.0 else -1

            if chord_norm == 0.0:
                bary = 0.5
            else:
                bary = float(np.linalg.norm(tp_xyz - p_xyz) / chord_norm)
                if bary < 0.0:
                    bary = 0.0
                elif bary > 1.0:
                    bary = 1.0

            spt_idx = splits.add_spt(sp_idx=sp_idx, bary=bary, vis_chge=vis_chge)
            splits.attach_to_segment(sis_pairs, s_idx, spt_idx,
                                     segment_label="SIS")


# ── O10 ───────────────────────────────────────────────────────────────────────

def _vp_uv_for_sp(vp: np.void) -> np.ndarray:
    """The VP's uv in the domain — used as the SP's uv anchor."""
    return np.asarray(vp["uv"], dtype=float).reshape(2)


def split_ccs_at_vps(
    vps: np.ndarray,
    css: np.ndarray,
    ccs: list,
    splits: SplitArrays,
    surface: "SurfaceParams",
    projection: "Projection",
) -> None:
    """Split CSs at every cusp (VP). One 'vp' SP and one SPT per VP.

    For each VP:
      - Create SP with `type="vp"`, `xyz=vp.xyz`, `xy=projection.XY(xyz)`,
        `uv=vp.uv`.
      - Attach a single SPT on the host CS (index `vp["cs"]`) with
        `bary = vp.s` and `vis_chge = vp.vis_change` (G10 — already ±1 from O4).

    `ccs` and `surface` are unused (no geometry — VPs already carry vis_change);
    kept in the signature for API parity with the roadmap.
    """
    del ccs, surface  # API parity only; algorithm reads vp + writes css.

    if len(vps) == 0:
        return

    for vp in vps:
        uv = _vp_uv_for_sp(vp)
        xyz = np.asarray(vp["xyz"], dtype=float).reshape(3)
        xy = projection.XY(xyz)
        sp_idx = splits.add_sp(uv=uv, xyz=xyz, xy=xy, sp_type="vp")

        cs_idx = int(vp["cs"])
        bary = float(vp["s"])
        vis = int(vp["vis_change"])
        spt_idx = splits.add_spt(sp_idx=sp_idx, bary=bary, vis_chge=vis)
        splits.attach_to_segment(css, cs_idx, spt_idx, segment_label="CS")


# ── O11 ───────────────────────────────────────────────────────────────────────

def _build_sis_preimage_segments(
    sis_pairs: np.ndarray,
    dps: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """For each SIS, build its two preimage segments (mirrors C12).

    Returns six (2K,) / (2K, 2) arrays (K = len(sis_pairs)):
      - seg_uv0, seg_uv1 : segment endpoints in uv (this preimage).
      - other_uv0, other_uv1 : endpoints of the OTHER preimage (same SIS).
      - sis_idx_of_seg : (2K,) — which sis_pair each segment belongs to.
      - valid_mask : (2K,) bool — whether the segment is well-formed (skip
        segments where flip / A-set determination fails).

    Both preimages of a SIS walk from p_dp to q_dp with consistent
    parameterization: segment k=2*pair is the "A" preimage, k=2*pair+1 is "B".
    The OTHER preimage of segment k can therefore be linearly interpolated at
    the same parameter t to recover an approximate uv on the other sheet.
    """
    K = len(sis_pairs)
    if K == 0:
        empty2 = np.empty((0, 2), dtype=float)
        empty1 = np.empty(0, dtype=np.int32)
        empty_b = np.empty(0, dtype=bool)
        return empty2, empty2, empty2, empty2, empty1, empty_b

    seg_uv0 = np.empty((2 * K, 2), dtype=float)
    seg_uv1 = np.empty((2 * K, 2), dtype=float)
    other_uv0 = np.empty((2 * K, 2), dtype=float)
    other_uv1 = np.empty((2 * K, 2), dtype=float)
    sis_idx_of_seg = np.empty(2 * K, dtype=np.int32)
    valid = np.zeros(2 * K, dtype=bool)

    for k in range(K):
        p = int(sis_pairs["p_dp"][k])
        q = int(sis_pairs["q_dp"][k])
        f = int(sis_pairs["flip"][k])

        A_p_uv = dps["uv1"][p]
        A_q_uv = dps["uv1"][q] if f == 1 else dps["uv2"][q]
        B_p_uv = dps["uv2"][p]
        B_q_uv = dps["uv2"][q] if f == 1 else dps["uv1"][q]

        # Segment 2k (A preimage), other is B.
        seg_uv0[2 * k] = A_p_uv
        seg_uv1[2 * k] = A_q_uv
        other_uv0[2 * k] = B_p_uv
        other_uv1[2 * k] = B_q_uv
        sis_idx_of_seg[2 * k] = k
        valid[2 * k] = True

        # Segment 2k+1 (B preimage), other is A.
        seg_uv0[2 * k + 1] = B_p_uv
        seg_uv1[2 * k + 1] = B_q_uv
        other_uv0[2 * k + 1] = A_p_uv
        other_uv1[2 * k + 1] = A_q_uv
        sis_idx_of_seg[2 * k + 1] = k
        valid[2 * k + 1] = True

    return seg_uv0, seg_uv1, other_uv0, other_uv1, sis_idx_of_seg, valid


def _cs_vis_chge_at_cdp(
    surface: "SurfaceParams",
    projection: "Projection",
    uv_hit: np.ndarray,
    cs_dir_uv: np.ndarray,
    other_uv: np.ndarray,
) -> int:
    """CS-side vis_chge at a CDP. Spec line 404 (= case-2 BDP formula reused).

    `vis = +1 if (axis·SN') * (T·SN') > 0 else -1`,
    where SN' is the OTHER sheet's 3D normal (at the SIS's other preimage uv)
    and T is the 3D CS tangent at the CDP, computed as `Su*cs_dir[0] +
    Sv*cs_dir[1]`. `cs_dir_uv` is the uv chord direction of the CS (close-
    adjusted by the caller).
    """
    if projection.mode != "ortho":
        raise NotImplementedError(
            "_cs_vis_chge_at_cdp: perspective not yet supported"
        )

    u, v = float(uv_hit[0]), float(uv_hit[1])
    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    T = float(cs_dir_uv[0]) * Su + float(cs_dir_uv[1]) * Sv

    ou, ov = float(other_uv[0]), float(other_uv[1])
    SN_other = np.asarray(surface.SN(ou, ov), dtype=float).reshape(3)

    axis = projection._axis
    factor = float(axis @ SN_other) * float(T @ SN_other)
    return 1 if factor > 0.0 else -1


def _sis_vis_chge_at_cdp(
    surface: "SurfaceParams",
    projection: "Projection",
    uv_hit: np.ndarray,
    cs_dir_uv: np.ndarray,
    sis_dir_uv: np.ndarray,
) -> int:
    """SIS-side vis_chge at a CDP. Spec line 408.

    `vis = +1 if inner(T, N) > 0 else -1`, where:
      - T is the 2D SIS preimage tangent in uv (close-adjusted by caller).
      - N is the CS normal in uv, oriented toward the front sheet. "Front sheet"
        per spec is interpretation (1) — intrinsic to THIS sheet: N points
        into the half of the CS where THIS sheet curves toward the viewer.
        Under "axis points toward viewer", that is the half where moving in
        +N direction increases axis·S in 3D.

    Implementation: take N0 = (-cs_dir_v, cs_dir_u) (perpendicular to CS
    chord in uv), then orient N0 so that
       `axis · (Su*N0[0] + Sv*N0[1]) > 0`
    (a small step in +N0 moves the surface point toward the viewer).
    """
    if projection.mode != "ortho":
        raise NotImplementedError(
            "_sis_vis_chge_at_cdp: perspective not yet supported"
        )

    u, v = float(uv_hit[0]), float(uv_hit[1])
    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    axis = projection._axis

    N = np.array([-float(cs_dir_uv[1]), float(cs_dir_uv[0])], dtype=float)
    # Orient N toward the front sheet (axis·dS_N > 0 under axis-toward-viewer).
    dS_N = N[0] * Su + N[1] * Sv
    if float(axis @ dS_N) < 0.0:
        N = -N

    T = np.asarray(sis_dir_uv, dtype=float).reshape(2)
    return 1 if float(T @ N) > 0.0 else -1


def split_at_cdps(
    mesh,
    css: np.ndarray,
    cps: np.ndarray,
    sis_pairs: np.ndarray,
    dps: np.ndarray,
    splits: SplitArrays,
    surface: "SurfaceParams",
    projection: "Projection",
) -> None:
    """Split CSs and SISs at every CS×SIS domain intersection (CDP).

    Pipeline:
      1. Build CS preimage segments (1 per CS) from `cps[p_cp].uv → cps[q_cp].uv`.
      2. Build SIS preimage segments (2 per SIS) — mirrors C12's
         `find_triple_points` preimage-segment build (intersections.py:782-807).
      3. Cross-mode sweep via `sweep_segments` (G3 close-aware via `mesh.domain`).
      4. For each hit: create SP type='cdp'; attach SPT on CS with vis_chge
         per spec line 404 (`_cs_vis_chge_at_cdp`); attach SPT on SIS with
         vis_chge per spec line 408 (`_sis_vis_chge_at_cdp`).

    The `cps` parameter is added to the signature (roadmap signature omits
    it — but the CS preimage segments need `cps[*].uv`).

    Mutates `css` (split1/split2), `sis_pairs` (split1/split2), and `splits`.
    """
    from surface_play.intersections import sweep_segments

    if projection.mode != "ortho":
        raise NotImplementedError("split_at_cdps: perspective not yet supported")

    if len(css) == 0 or len(sis_pairs) == 0:
        return

    domain = getattr(mesh, "domain", None)

    # Step 1: CS preimage segments.
    cs_uv0 = cps["uv"][css["p_cp"]]      # (M, 2)
    cs_uv1 = cps["uv"][css["q_cp"]]      # (M, 2)

    # Step 2: SIS preimage segments (2 per SIS).
    (sis_uv0, sis_uv1, other_uv0, other_uv1,
     sis_idx_of_seg, sis_seg_valid) = _build_sis_preimage_segments(sis_pairs, dps)
    if not sis_seg_valid.any():
        return
    # Drop invalid segments (none expected with current C10 output, but be safe).
    if not sis_seg_valid.all():
        sis_uv0 = sis_uv0[sis_seg_valid]
        sis_uv1 = sis_uv1[sis_seg_valid]
        other_uv0 = other_uv0[sis_seg_valid]
        other_uv1 = other_uv1[sis_seg_valid]
        sis_idx_of_seg = sis_idx_of_seg[sis_seg_valid]

    # Step 3: cross-mode close-aware sweep.
    hits = sweep_segments(cs_uv0, cs_uv1, sis_uv0, sis_uv1, domain)
    if len(hits) == 0:
        return

    # Step 4: process hits.
    for h in hits:
        cs_idx = int(h["a"])
        sis_seg_idx = int(h["b"])
        sis_idx = int(sis_idx_of_seg[sis_seg_idx])
        uv_hit = np.asarray(h["uv"], dtype=float).reshape(2)
        t_a = float(h["t_a"])
        t_b = float(h["t_b"])

        # 3D position + screen XY of the CDP.
        xyz_hit = np.asarray(
            surface.S(float(uv_hit[0]), float(uv_hit[1])), dtype=float,
        ).reshape(3)
        xy_hit = projection.XY(xyz_hit)
        sp_idx = splits.add_sp(uv=uv_hit, xyz=xyz_hit, xy=xy_hit, sp_type="cdp")

        # OTHER preimage uv: lerp on the other preimage's segment at t_b.
        # Both preimages walk from p_dp to q_dp (consistent flip-handled
        # parametrization built in `_build_sis_preimage_segments`).
        o_uv0 = np.asarray(other_uv0[sis_seg_idx], dtype=float).reshape(2)
        o_uv1 = np.asarray(other_uv1[sis_seg_idx], dtype=float).reshape(2)
        # Close-aware lerp: shift o_uv1 toward o_uv0 if domain identifies axes.
        if domain is not None and getattr(domain, "type", None) == "rect":
            o_uv1 = domain.close(o_uv0, o_uv1)
        other_uv = o_uv0 + t_b * (o_uv1 - o_uv0)

        # Direction vectors (close-adjusted) for the vis_chge formulas.
        cs_p = np.asarray(cs_uv0[cs_idx], dtype=float).reshape(2)
        cs_q = np.asarray(cs_uv1[cs_idx], dtype=float).reshape(2)
        if domain is not None and getattr(domain, "type", None) == "rect":
            cs_q = domain.close(cs_p, cs_q)
        cs_dir_uv = cs_q - cs_p

        sis_p = np.asarray(sis_uv0[sis_seg_idx], dtype=float).reshape(2)
        sis_q = np.asarray(sis_uv1[sis_seg_idx], dtype=float).reshape(2)
        if domain is not None and getattr(domain, "type", None) == "rect":
            sis_q = domain.close(sis_p, sis_q)
        sis_dir_uv = sis_q - sis_p

        # CS-side SPT.
        vis_cs = _cs_vis_chge_at_cdp(
            surface, projection, uv_hit, cs_dir_uv, other_uv,
        )
        spt_cs = splits.add_spt(sp_idx=sp_idx, bary=t_a, vis_chge=vis_cs)
        splits.attach_to_segment(css, cs_idx, spt_cs, segment_label="CS")

        # SIS-side SPT.
        vis_sis = _sis_vis_chge_at_cdp(
            surface, projection, uv_hit, cs_dir_uv, sis_dir_uv,
        )
        spt_sis = splits.add_spt(sp_idx=sp_idx, bary=t_b, vis_chge=vis_sis)
        splits.attach_to_segment(sis_pairs, sis_idx, spt_sis,
                                 segment_label="SIS")
