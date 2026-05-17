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
    from surface_play.curves import BoundaryCurve
    from surface_play.mesh import Mesh
    from surface_play.projection import Projection


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
