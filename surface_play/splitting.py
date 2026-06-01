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

    `internal[i] = (abs_seg_idx, segment_local_bary)`. `bary` is in the segment's
    own parameterization (0 or 1 for join-vertex points); O14 reads the uv as
    `seg.p + bary * seg.pq` without needing the chain direction.
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
    ("uv",     "2f8"),
    ("xyz",    "3f8"),
    ("xy",     "2f8"),
    ("type",   "U3"),
    # Alternate preimage uv for SPs that sit on a SIC chain (bdp, cdp, and ha
    # anchored on a SIC). Two sheets meet at a SIC's 3D self-intersection;
    # `uv` records the canonical sheet, `uv_alt` records the other sheet's
    # uv. NaN for SPs that don't sit on a SIC (cn, bcp, vp, ha-on-CC).
    # TPs have 3 preimages — not handled by uv_alt; treated as NaN for now.
    ("uv_alt", "2f8"),
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

    def add_sp(self, uv, xyz, xy, sp_type, uv_alt=None) -> int:
        if uv_alt is None:
            uv_alt = (float("nan"), float("nan"))
        self.sps.append(
            (tuple(uv), tuple(xyz), tuple(xy), sp_type, tuple(uv_alt))
        )
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
    front_uv: np.ndarray | None = None,
) -> int:
    """Visibility change on a boundary edge at parametric position s ∈ (0, 1).

    Adopted 2026-05-24 (Etienne): the "alternative formula" first proposed
    on 2026-05-19. Two changes vs spec:

      1. `Np = ∇_uv(axis·SN)`, not the spec's `∇_uv(axis·S)` (depth gradient).
         The silhouette curve is `axis·SN = 0`, so its gradient is the
         geometrically correct quantity for sign-of-silhouette-crossing.
         The depth-gradient version gave the same sign at both BCPs of the
         same CC on paraboloid trial 2 (cycle did not close).

      2. The "no visibility change" gate uses the projected boundary and
         silhouette tangents in the image (`f3 = Tb_proj · Tp_proj`)
         instead of the brittle `ker · Nb` test.

    The alternative formula computes:

        ker   = kernel of dS-projection in uv, oriented so axis·dS(ker) > 0.
        Tb    = dp / |dp|  (boundary tangent in uv).
        Np    = ∇_uv(axis · SN)
               = ( (Suu × Sv + Su × Suv) · axis,
                   (Suv × Sv + Su × Svv) · axis ),
                oriented so Np · ker > 0.
        Tp    = (-Np[1], Np[0])    (uv-perpendicular to Np = silhouette tangent in uv)
                oriented inward (so Tp · Nb > 0).
        Tb_proj, Tp_proj  = projections of dS(Tb), dS(Tp) onto the screen plane.

        f1 = Tb · Np         (sign of silhouette-function change along Tb)
        f2 = Np · ker        (always +1 after orientation; kept for symmetry)
        f3 = Tb_proj · Tp_proj   (whether projected boundary tangent and
                                  silhouette tangent agree in image)

        if f1·f2·f3 > 0:  return 0    (no observable change; "exit" BCP)
        return +1 if f1·f2 > 0 else -1

    Known degeneracy: disk_paraboloid_ca trial 3 (nearly edge-on view) is
    INFEASIBLE at certain mesh resolutions because `|Tb_proj|` becomes
    tiny and f3 is FP-noise. A `axis · dS(Nb) ≤ 0` gate was tried as a
    replacement but introduced new persistent failures on the paraboloid
    (likely mathematically wrong); reverted 2026-05-24.

    `Nb = edge["dir"]` is the INWARD boundary normal in uv (mesh.py flips
    `d` to point toward the boundary triangle's apex).

    Persp (added 2026-05-27): `axis` is the per-point viewer direction
    `eye − S(p0)`. The Np formula gains terms `−Su·SN` and `−Sv·SN` which
    vanish identically (`Su · (Su × Sv) = 0`), so the symbolic Np
    expression is unchanged. Screen-space tangents Tb_proj / Tp_proj
    switch from the 3D-residual form `T − (axis·T)·axis` (valid only when
    axis is the constant image-plane normal) to `projection.proj_vec(uv,
    T_3d)` (2D) — ortho-equivalent in sign, persp-correct.
    """
    p = np.asarray(edge["p"], dtype=float).reshape(2)
    dp = np.asarray(edge["pq"], dtype=float).reshape(2)
    Nb = np.asarray(edge["dir"], dtype=float).reshape(2)

    norm_dp = float(np.linalg.norm(dp))
    if norm_dp == 0.0:
        return 0

    p0 = p + float(s) * dp
    u, v = float(p0[0]), float(p0[1])

    # Tb: analytic boundary tangent at the BCP, sign-matched to the edge's
    # chord direction. The discretized chord can deviate from the smooth
    # tangent by O(1/res), and in near-edge-on views that deviation can flip
    # f3's sign (whose magnitude is set by the screen-projection of the
    # boundary's image-space tangent, which is itself small in such views).
    Tb = surface.domain.boundary_tangent(p0, dp)

    ker = projection.ker_param(p0)

    S_p = np.asarray(surface.S(u, v), dtype=float).reshape(3)
    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    axis = projection.viewer_direction(S_p).reshape(3)

    dS_ker = ker[0] * Su + ker[1] * Sv
    if float(axis @ dS_ker) < 0.0:
        ker = -ker

    # Np = gradient of axis·SN in uv, oriented so Np·ker > 0. If the caller
    # has the precomputed `front_uv` from the CP that lives at this BCP,
    # use it directly (Spec change 2026-05-27 — one source of truth for the
    # front-sheet vector). Otherwise re-derive locally.
    if front_uv is not None:
        Np = np.asarray(front_uv, dtype=float).reshape(2)
    else:
        Suu = np.asarray(surface.Suu(u, v), dtype=float).reshape(3)
        Suv = np.asarray(surface.Suv(u, v), dtype=float).reshape(3)
        Svv = np.asarray(surface.Svv(u, v), dtype=float).reshape(3)
        Np = np.array([
            float(np.cross(Suu, Sv) @ axis + np.cross(Su, Suv) @ axis),
            float(np.cross(Suv, Sv) @ axis + np.cross(Su, Svv) @ axis),
        ])
        if float(np.dot(Np, ker)) < 0.0:
            Np = -Np

    Tp = np.array([-Np[1], Np[0]])
    if float(np.dot(Tp, Nb)) < 0.0:
        Tp = -Tp

    dS_Tb = Tb[0] * Su + Tb[1] * Sv
    dS_Tp = Tp[0] * Su + Tp[1] * Sv
    Tb_proj = projection.proj_vec(p0, dS_Tb)
    Tp_proj = projection.proj_vec(p0, dS_Tp)

    f1 = float(np.dot(Tb, Np))
    f2 = float(np.dot(Np, ker))
    f3 = float(np.dot(Tb_proj, Tp_proj))

    if f1 * f2 * f3 > 0.0:
        return 0
    return 1 if f1 * f2 > 0.0 else -1


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
        vis = _bvis_chge(surface, projection, mesh.edges[e_idx], bary_be,
                         front_uv=cp["front_uv"])
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

    Persp (2026-05-27): `axis = viewer_direction(S(uv_DP))`. Both sheets
    share the same 3D point at a DP, so any uv mapping to it gives the
    same xyz.
    """
    p = np.asarray(edge["p"], dtype=float).reshape(2)
    pq = np.asarray(edge["pq"], dtype=float).reshape(2)
    uv = p + float(s) * pq
    u, v = float(uv[0]), float(uv[1])
    S_p = np.asarray(surface.S(u, v), dtype=float).reshape(3)
    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    Tb = pq[0] * Su + pq[1] * Sv

    ou, ov = float(other_uv[0]), float(other_uv[1])
    SN_other = np.asarray(surface.SN(ou, ov), dtype=float).reshape(3)

    axis = projection.viewer_direction(S_p).reshape(3)
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

    Persp (2026-05-27): `axis = viewer_direction(S(uv_DP))`. proj_vec
    already handles both modes.
    """
    # OUR edge geometry.
    p = np.asarray(edge_self["p"], dtype=float).reshape(2)
    pq = np.asarray(edge_self["pq"], dtype=float).reshape(2)
    uv = p + float(s_self) * pq
    u, v = float(uv[0]), float(uv[1])
    S_p = np.asarray(surface.S(u, v), dtype=float).reshape(3)
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

    axis = projection.viewer_direction(S_p).reshape(3)
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
        uv_alt = np.asarray(dp["uv2"], dtype=float)  # sheet-2 uv for SIC walks
        xy = projection.XY(xyz)
        sp_idx = splits.add_sp(uv=uv, xyz=xyz, xy=xy, sp_type="bdp",
                               uv_alt=uv_alt)

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

    Persp (2026-05-27): `axis = viewer_direction(tp_xyz)` — same 3D point
    for all three preimages, so axis is per-TP not per-preimage.
    """
    if len(tps) == 0:
        return

    for tp in tps:
        tp_xyz = np.asarray(tp["xyz"], dtype=float).reshape(3)
        tp_xy = projection.XY(tp_xyz)
        axis = projection.viewer_direction(tp_xyz).reshape(3)
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

    Persp (2026-05-27): `axis = viewer_direction(S(uv_hit))`.
    """
    u, v = float(uv_hit[0]), float(uv_hit[1])
    S_p = np.asarray(surface.S(u, v), dtype=float).reshape(3)
    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
    T = float(cs_dir_uv[0]) * Su + float(cs_dir_uv[1]) * Sv

    ou, ov = float(other_uv[0]), float(other_uv[1])
    SN_other = np.asarray(surface.SN(ou, ov), dtype=float).reshape(3)

    axis = projection.viewer_direction(S_p).reshape(3)
    factor = float(axis @ SN_other) * float(T @ SN_other)
    return 1 if factor > 0.0 else -1


def _sis_vis_chge_at_cdp(
    surface: "SurfaceParams",
    projection: "Projection",
    uv_hit: np.ndarray,
    cs_dir_uv: np.ndarray,
    sis_dir_uv: np.ndarray,
    front_uv: np.ndarray | None = None,
) -> int:
    """SIS-side vis_chge at a CDP. Spec line 408.

    `vis = +1 if inner(T, N) > 0 else -1`, where:
      - T is the 2D SIS preimage tangent in uv (close-adjusted by caller).
      - N is the CS normal in uv, oriented toward the front sheet. "Front sheet"
        per spec is interpretation (1) — intrinsic to THIS sheet: N points
        into the half of the CS where THIS sheet curves toward the viewer.

    If the caller provides `front_uv` (the per-CP front-sheet vector
    interpolated at the CDP from the CS's two endpoints), use it directly as
    N. This is the spec change of 2026-05-27 (one source of truth for the
    front-sheet vector). Otherwise fall back to the local construction:
    N = (-cs_dir_v, cs_dir_u), oriented so that `axis · dS(N) > 0`.

    For a CS lying on the silhouette `axis·SN = 0`, `front_uv = grad_uv
    (axis·SN)` is perpendicular to the contour tangent in uv, exactly aligned
    with the local "perpendicular to cs_dir_uv toward front sheet" — only the
    sign matters for the final dot, so the substitution is sign-equivalent.

    Persp (2026-05-27): `axis = viewer_direction(S(uv_hit))`.
    """
    if front_uv is not None:
        N = np.asarray(front_uv, dtype=float).reshape(2)
    else:
        u, v = float(uv_hit[0]), float(uv_hit[1])
        S_p = np.asarray(surface.S(u, v), dtype=float).reshape(3)
        Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
        Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)
        axis = projection.viewer_direction(S_p).reshape(3)
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

    # Step 3: cross-mode seam-aware sweep. `sweep_segments` localizes each
    # segment via `domain.interpolate` (the seam-aware primitive), so a seam
    # face's CS / seam-spanning SIS preimage is brought into its own local frame
    # rather than appearing as a spurious diameter-long chord.
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

        # OTHER preimage uv: lerp on the other preimage's segment at t_b.
        # Both preimages walk from p_dp to q_dp (consistent flip-handled
        # parametrization built in `_build_sis_preimage_segments`).
        # Stored as the SP's uv_alt so SIC polyline walks can land on the
        # matching sheet at the chain endpoint.
        o_uv0 = np.asarray(other_uv0[sis_seg_idx], dtype=float).reshape(2)
        o_uv1 = np.asarray(other_uv1[sis_seg_idx], dtype=float).reshape(2)
        # Seam-aware lerp at t_b along the other preimage segment.
        if domain is not None:
            other_uv = domain.interpolate(o_uv0, o_uv1, t_b)
        else:
            other_uv = o_uv0 + t_b * (o_uv1 - o_uv0)

        sp_idx = splits.add_sp(uv=uv_hit, xyz=xyz_hit, xy=xy_hit, sp_type="cdp",
                               uv_alt=other_uv)

        # Direction vectors (seam-aware) for the vis_chge formulas.
        cs_p = np.asarray(cs_uv0[cs_idx], dtype=float).reshape(2)
        cs_q = np.asarray(cs_uv1[cs_idx], dtype=float).reshape(2)
        cs_q_loc = domain.localize(cs_p, cs_q) if domain is not None else cs_q
        cs_dir_uv = cs_q_loc - cs_p

        sis_p = np.asarray(sis_uv0[sis_seg_idx], dtype=float).reshape(2)
        sis_q = np.asarray(sis_uv1[sis_seg_idx], dtype=float).reshape(2)
        sis_q_loc = domain.localize(sis_p, sis_q) if domain is not None else sis_q
        sis_dir_uv = sis_q_loc - sis_p

        # CS-side SPT.
        vis_cs = _cs_vis_chge_at_cdp(
            surface, projection, uv_hit, cs_dir_uv, other_uv,
        )
        spt_cs = splits.add_spt(sp_idx=sp_idx, bary=t_a, vis_chge=vis_cs)
        splits.attach_to_segment(css, cs_idx, spt_cs, segment_label="CS")

        # SIS-side SPT.
        # Interpolate front_uv at the CDP from the CS endpoints — used as the
        # CS-normal-toward-front-sheet (spec change 2026-05-27).
        cs_row = css[cs_idx]
        front_p = np.asarray(cps[int(cs_row["p_cp"])]["front_uv"], dtype=float)
        front_q = np.asarray(cps[int(cs_row["q_cp"])]["front_uv"], dtype=float)
        front_uv_cdp = (1.0 - t_a) * front_p + t_a * front_q
        vis_sis = _sis_vis_chge_at_cdp(
            surface, projection, uv_hit, cs_dir_uv, sis_dir_uv,
            front_uv=front_uv_cdp,
        )
        spt_sis = splits.add_spt(sp_idx=sp_idx, bary=t_b, vis_chge=vis_sis)
        splits.attach_to_segment(sis_pairs, sis_idx, spt_sis,
                                 segment_label="SIS")


# ── O13 ───────────────────────────────────────────────────────────────────────

def _split_one_chain(
    kind: str,
    parent_idx: int,
    chain: np.ndarray,
    seg_array: np.ndarray,
    splits: SplitArrays,
    is_closed: bool,
    chain_to_abs,
) -> list[SubCurve]:
    """Cut a single BC/CC/SIC chain at its SPTs into SubCurves.

    `chain` is the parent's signed token array (BC.edge_indices /
    CC.cs_indices / SIC.sis_indices). `chain_to_abs(k)` maps `abs(chain[k])-1`
    to the absolute index into `seg_array` (identity for CC/SIC; for BC it
    is `mesh.boundary_edge_idx[abs(chain[k])-1]`).

    Returns SubCurves; mutates nothing.
    """
    n_chain = len(chain)
    if n_chain == 0:
        return []

    # Unique segments along the chain (closed chains repeat the first token at the end).
    if is_closed and n_chain >= 2:
        N_seg = n_chain - 1
    else:
        N_seg = n_chain
    if N_seg == 0:
        return []

    abs_segs = np.empty(N_seg, dtype=np.int64)
    signs = np.empty(N_seg, dtype=np.int8)
    for k in range(N_seg):
        signed = int(chain[k])
        signs[k] = 1 if signed > 0 else -1
        abs_segs[k] = int(chain_to_abs(abs(signed) - 1))

    # Walk chain, collecting SPTs in traversal order.
    # all_spts[i] = (chain_pos, slot, sp_idx, bary_chain, vis_chge).
    all_spts: list[tuple[int, int, int, float, int]] = []
    for cp in range(N_seg):
        abs_idx = int(abs_segs[cp])
        sign = int(signs[cp])
        seg = seg_array[abs_idx]
        local: list[tuple[float, int, int, int]] = []
        for slot_name in ("split1", "split2"):
            slot = int(seg[slot_name])
            if slot < 0:
                continue
            spt = splits.spts[slot]
            sp_idx = int(spt[0])
            bary_native = float(spt[1])
            vis_chge = int(spt[2])
            bary_chain = bary_native if sign > 0 else 1.0 - bary_native
            local.append((bary_chain, slot, sp_idx, vis_chge))
        local.sort(key=lambda x: x[0])
        for bary_chain, slot, sp_idx, vis_chge in local:
            all_spts.append((cp, slot, sp_idx, bary_chain, vis_chge))

    # Merge consecutive same-SP SPTs (corner case: 2 SPTs at one corner SP, one
    # on each incident BE, with bary_chain 1.0 then 0.0 across the chain join).
    # Marker fields:
    #   sp_idx
    #   out_cp, out_bc, out_vc   — arc ARRIVING at this SP reads these
    #   in_cp,  in_bc,  in_vc    — arc LEAVING from this SP reads these
    Marker = tuple  # (sp_idx, out_cp, out_bc, out_vc, in_cp, in_bc, in_vc)
    markers: list[tuple] = []
    for (cp, slot, sp_idx, bary_chain, vis_chge) in all_spts:
        if markers and markers[-1][0] == sp_idx:
            prev = markers[-1]
            # Update IN side to the newer SPT.
            markers[-1] = (prev[0], prev[1], prev[2], prev[3], cp, bary_chain, vis_chge)
        else:
            markers.append((sp_idx, cp, bary_chain, vis_chge, cp, bary_chain, vis_chge))

    # For closed chains, the last marker can wrap into the first if they share
    # the same SP across the chain seam.
    if is_closed and len(markers) >= 2 and markers[-1][0] == markers[0][0]:
        last = markers.pop()
        first = markers[0]
        # The last marker's OUT-side becomes the first marker's OUT-side.
        markers[0] = (first[0], last[1], last[2], last[3], first[4], first[5], first[6])

    def _internal_forward(cp_a: int, cp_b: int, wrap: bool,
                          in_bc: float | None = None,
                          out_bc: float | None = None) -> list[tuple[int, float]]:
        """Segment-end intermediate points for an arc starting in segment cp_a
        and ending in segment cp_b. Each entry is (abs_seg_idx, segment_local_bary)
        — the end of segment j in chain-traversal direction, expressed in the
        segment's NATIVE parameterization (1.0 if the segment is forward in the
        chain, 0.0 if reversed). O14 reads uv as `seg.p + bary * seg.pq` directly.

        When the start marker sits at the END of cp_a in chain direction
        (in_bc == 1.0), the j=cp_a entry would duplicate the start SP — drop it.
        Symmetrically, when the end marker sits at the START of cp_b in chain
        direction (out_bc == 0.0), the j=cp_b-1 entry would duplicate the end
        SP — drop it. Without this, `_build_polyline` produces zero-length
        tail/head segments at vertex-located SPTs (e.g., HAs).
        """
        pts: list[tuple[int, float]] = []
        def _end_bary(j: int) -> float:
            return 1.0 if int(signs[j]) > 0 else 0.0
        skip_first = (in_bc is not None and in_bc == 1.0)
        skip_last = (out_bc is not None and out_bc == 0.0)
        if not wrap:
            if cp_a == cp_b:
                return pts
            j_lo = cp_a + 1 if skip_first else cp_a
            j_hi = cp_b - 1 if skip_last else cp_b
            for j in range(j_lo, j_hi):
                pts.append((int(abs_segs[j]), _end_bary(j)))
            return pts
        # Wrap: walk forward (mod N_seg) from cp_a, stop when we'd re-enter cp_b.
        j = cp_a
        is_first = True
        while True:
            will_be_last = (((j + 1) % N_seg) == cp_b)
            skip_this = (is_first and skip_first) or (will_be_last and skip_last)
            if not skip_this:
                pts.append((int(abs_segs[j]), _end_bary(j)))
            j = (j + 1) % N_seg
            is_first = False
            if j == cp_b:
                break
        return pts

    out: list[SubCurve] = []

    if not markers:
        # No SPTs anywhere on the chain.
        internal = [
            (int(abs_segs[j]), 1.0 if int(signs[j]) > 0 else 0.0)
            for j in range(N_seg)
        ]
        out.append(SubCurve(
            kind=kind, is_closed=bool(is_closed),
            start=-1, end=-1, internal=internal,
            vc_in=0, vc_out=0, parent_idx=parent_idx,
        ))
        return out

    if is_closed:
        n = len(markers)
        for k in range(n):
            cur = markers[k]
            nxt = markers[(k + 1) % n]
            sp_a = int(cur[0])
            in_cp_a, in_bc_a, in_vc_a = int(cur[4]), float(cur[5]), int(cur[6])
            sp_b = int(nxt[0])
            out_cp_b, out_bc_b, out_vc_b = int(nxt[1]), float(nxt[2]), int(nxt[3])

            # Wrap detection: does this arc cross the chain seam?
            if in_cp_a < out_cp_b:
                wrap = False
            elif in_cp_a > out_cp_b:
                wrap = True
            else:  # in_cp_a == out_cp_b: same chain segment
                if n == 1:
                    wrap = True   # single SP closed loop — full loop
                else:
                    wrap = in_bc_a > out_bc_b

            internal = _internal_forward(in_cp_a, out_cp_b, wrap,
                                         in_bc=in_bc_a, out_bc=out_bc_b)

            sign_in = int(signs[in_cp_a])
            sign_out = int(signs[out_cp_b])
            vc_in = min(0, sign_in * in_vc_a)
            vc_out = min(0, -sign_out * out_vc_b)

            sub_closed = (n == 1)
            out.append(SubCurve(
                kind=kind, is_closed=sub_closed,
                start=sp_a, end=sp_b, internal=internal,
                vc_in=vc_in, vc_out=vc_out, parent_idx=parent_idx,
            ))
        return out

    # Open chain: n-1 SubCurves between consecutive markers (markers act as endpoints).
    n = len(markers)
    if n < 2:
        return []
    for k in range(n - 1):
        cur = markers[k]
        nxt = markers[k + 1]
        sp_a = int(cur[0])
        in_cp_a, in_bc_a, in_vc_a = int(cur[4]), float(cur[5]), int(cur[6])
        sp_b = int(nxt[0])
        out_cp_b, out_bc_b, out_vc_b = int(nxt[1]), float(nxt[2]), int(nxt[3])

        internal = _internal_forward(in_cp_a, out_cp_b, wrap=False,
                                     in_bc=in_bc_a, out_bc=out_bc_b)
        sign_in = int(signs[in_cp_a])
        sign_out = int(signs[out_cp_b])
        vc_in = min(0, sign_in * in_vc_a)
        vc_out = min(0, -sign_out * out_vc_b)

        out.append(SubCurve(
            kind=kind, is_closed=False,
            start=sp_a, end=sp_b, internal=internal,
            vc_in=vc_in, vc_out=vc_out, parent_idx=parent_idx,
        ))
    return out


def assemble_subcurves(
    bcs: "list[BoundaryCurve]",
    ccs: "list[ContourCurve]",
    sics: list,
    hcs: list[SubCurve],
    mesh: "Mesh",
    css: np.ndarray,
    sis_pairs: np.ndarray,
    splits: SplitArrays,
) -> list[SubCurve]:
    """Cut every BC/CC/SIC at its SPTs into SubCurves; append the HCs verbatim.

    Spec: §"Construction of Split Curves (SCs)" (lines 442-450).

    HCs from O12 are already SubCurve(kind="HC", internal=[]) — they are NOT
    re-split here (spec: "Apart from HCs, ..."), simply concatenated.
    """
    out: list[SubCurve] = []

    bnd_idx = mesh.boundary_edge_idx
    for bc_idx, bc in enumerate(bcs):
        out.extend(_split_one_chain(
            "BC", bc_idx, bc.edge_indices, mesh.edges, splits,
            bool(bc.is_closed),
            chain_to_abs=lambda k: int(bnd_idx[k]),
        ))

    for cc_idx, cc in enumerate(ccs):
        out.extend(_split_one_chain(
            "CC", cc_idx, cc.cs_indices, css, splits,
            bool(cc.is_closed),
            chain_to_abs=lambda k: k,
        ))

    for sic_idx, sic in enumerate(sics):
        out.extend(_split_one_chain(
            "SIC", sic_idx, sic.sis_indices, sis_pairs, splits,
            bool(sic.is_closed),
            chain_to_abs=lambda k: k,
        ))

    out.extend(hcs)
    return out
