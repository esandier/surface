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
