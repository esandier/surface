"""helpers.py — O12: `build_helper_curves`.

Bridges disconnected curve components (BCs / CCs / SICs) into one connected
graph by emitting "helper curves" (HCs) — straight lines in the domain
between the closest sample-cloud points of two components.

Pipeline:
  1. Build a Sample for every segment endpoint of every BC/CC/SIC.
  2. Union-find over (BC, CC, SIC) nodes by SP-identity (G5): two nodes
     are united iff they share an SP via any SPT.
  3. Loop until 1 component remains:
       a. Compute min Euclidean uv-distance between every active pair.
       b. Pick the closest pair (i, j) and its argmin samples (qi, qj).
       c. Emit two HA SPs + SPTs (G6 dedup: reuse SP at same target sample).
       d. Compute vc_in / vc_out per spec lines 432-440 (BC → 0;
          CC → 0/-1 by T·N; SIC → 0/-1 by (axis·SN')·(T·SN')).
       e. Append SubCurve(kind="HC", ...).
       f. Merge components.

G19: distances are plain Euclidean in (uv) space; identification's
`close()` is NOT applied — HC straight lines run through the original
rectangle, not via seam wrap-around.

Spec source: `Reorganization_of_the_surface_app.md` §"Construction of
helper curves" (lines 410-440) and roadmap O12 (lines 1268-1304).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from surface_play.splitting import SplitArrays, SubCurve

if TYPE_CHECKING:
    from surface_play.contour import ContourCurve
    from surface_play.curves import BoundaryCurve
    from surface_play.intersections import SelfIntersectingCurve
    from surface_play.mesh import Mesh
    from surface_play.projection import Projection
    from surface_play.surface import SurfaceParams


@dataclass
class _Sample:
    """One candidate attachment point for an HA.

    `uv`           : 2D position in domain (G19: original-rect coords).
    `parent_kind`  : "BC" / "CC" / "SIC".
    `parent_idx`   : index into bcs / ccs / sics.
    `seg_global`   : index into mesh.edges (BC) / css (CC) / sis_pairs (SIC).
    `end`          : 0 (p-end of segment) or 1 (q-end).
    `dedup_key`    : hashable per-point identity for G6 dedup —
                       BC : mesh vertex index.
                       CC : cp index.
                       SIC: (dp_idx, sheet) — sheet 1 = uv1, 2 = uv2.
                     Two samples representing the same physical point on
                     the same parent share this key.
    `other_uv`     : SIC-only. uv of the OTHER preimage at this DP; None
                     for BC/CC. Used by `_vc_for_sic` to evaluate SN' at
                     the other sheet.
    `cs_dir_uv`    : CC-only. The CS chord direction in uv at this sample,
                     used to derive the front-sheet normal N for vc.
    """

    uv: np.ndarray
    parent_kind: str
    parent_idx: int
    seg_global: int
    end: int
    dedup_key: tuple
    other_uv: np.ndarray | None = None
    cs_dir_uv: np.ndarray | None = None


# ── sample collection ────────────────────────────────────────────────────────

def _bc_samples(
    bc_idx: int, bc: "BoundaryCurve", mesh: "Mesh"
) -> list[_Sample]:
    """One sample per endpoint of every BE in the BC chain.

    The BC chain is signed tokens over `mesh.boundary_edge_idx` (P4
    convention); we map back to absolute `mesh.edges` indices so the
    segment_global field matches what `splits.attach_to_segment(
    mesh.edges, ...)` writes into.
    """
    out: list[_Sample] = []
    if len(bc.edge_indices) == 0:
        return out
    bnd_idx = mesh.boundary_edge_idx
    for signed in bc.edge_indices:
        i = abs(int(signed)) - 1
        e_idx = int(bnd_idx[i])
        edge = mesh.edges[e_idx]
        p_v = int(edge["p_idx"])
        q_v = int(edge["q_idx"])
        out.append(_Sample(
            uv=np.asarray(mesh.uv[p_v], dtype=float).copy(),
            parent_kind="BC", parent_idx=bc_idx,
            seg_global=e_idx, end=0,
            dedup_key=("BC", bc_idx, p_v),
        ))
        out.append(_Sample(
            uv=np.asarray(mesh.uv[q_v], dtype=float).copy(),
            parent_kind="BC", parent_idx=bc_idx,
            seg_global=e_idx, end=1,
            dedup_key=("BC", bc_idx, q_v),
        ))
    return out


def _cc_samples(
    cc_idx: int, cc: "ContourCurve",
    css: np.ndarray, cps: np.ndarray,
) -> list[_Sample]:
    """One sample per endpoint of every CS in the CC chain.

    `cs_dir_uv` is the chord direction `cps[q].uv - cps[p].uv` (NOT
    close-adjusted; on identified domains the CS lives inside a single
    face whose uv anchors are in the same rectangle copy).
    """
    out: list[_Sample] = []
    for signed in cc.cs_indices:
        cs_idx = abs(int(signed)) - 1
        cs = css[cs_idx]
        p_cp = int(cs["p_cp"])
        q_cp = int(cs["q_cp"])
        p_uv = np.asarray(cps[p_cp]["uv"], dtype=float).reshape(2)
        q_uv = np.asarray(cps[q_cp]["uv"], dtype=float).reshape(2)
        chord = q_uv - p_uv
        out.append(_Sample(
            uv=p_uv.copy(),
            parent_kind="CC", parent_idx=cc_idx,
            seg_global=cs_idx, end=0,
            dedup_key=("CC", cc_idx, p_cp),
            cs_dir_uv=chord,
        ))
        out.append(_Sample(
            uv=q_uv.copy(),
            parent_kind="CC", parent_idx=cc_idx,
            seg_global=cs_idx, end=1,
            dedup_key=("CC", cc_idx, q_cp),
            cs_dir_uv=chord,
        ))
    return out


def _sic_samples(
    sic_idx: int, sic: "SelfIntersectingCurve",
    sis_pairs: np.ndarray, dps: np.ndarray,
) -> list[_Sample]:
    """Both preimages of every SIS contribute samples (one per endpoint).

    Per spec line 415, an SIC's two preimages are ONE curve (one node).
    For each SIS in the chain:
      - A preimage : (dps[p].uv1, dps[q].uv1 if flip=+1 else dps[q].uv2)
      - B preimage : (dps[p].uv2, dps[q].uv2 if flip=+1 else dps[q].uv1)
    `other_uv` at each endpoint is the uv on the OPPOSITE preimage at
    the same DP (i.e., the other sheet's representation of this 3D
    intersection point). `flip` re-orders only at the q-end.
    """
    out: list[_Sample] = []
    for signed in sic.sis_indices:
        sis_idx = abs(int(signed)) - 1
        row = sis_pairs[sis_idx]
        p_dp = int(row["p_dp"])
        q_dp = int(row["q_dp"])
        f = int(row["flip"])

        A_p = np.asarray(dps[p_dp]["uv1"], dtype=float).reshape(2)
        B_p = np.asarray(dps[p_dp]["uv2"], dtype=float).reshape(2)
        if f == 1:
            A_q = np.asarray(dps[q_dp]["uv1"], dtype=float).reshape(2)
            B_q = np.asarray(dps[q_dp]["uv2"], dtype=float).reshape(2)
        else:
            A_q = np.asarray(dps[q_dp]["uv2"], dtype=float).reshape(2)
            B_q = np.asarray(dps[q_dp]["uv1"], dtype=float).reshape(2)

        # A preimage endpoints (sheet 1 at p_dp; sheet depends on flip at q).
        out.append(_Sample(
            uv=A_p.copy(), parent_kind="SIC", parent_idx=sic_idx,
            seg_global=sis_idx, end=0,
            dedup_key=("SIC", sic_idx, p_dp, 1),
            other_uv=B_p.copy(),
        ))
        out.append(_Sample(
            uv=A_q.copy(), parent_kind="SIC", parent_idx=sic_idx,
            seg_global=sis_idx, end=1,
            dedup_key=("SIC", sic_idx, q_dp, 1 if f == 1 else 2),
            other_uv=B_q.copy(),
        ))
        # B preimage endpoints.
        out.append(_Sample(
            uv=B_p.copy(), parent_kind="SIC", parent_idx=sic_idx,
            seg_global=sis_idx, end=0,
            dedup_key=("SIC", sic_idx, p_dp, 2),
            other_uv=A_p.copy(),
        ))
        out.append(_Sample(
            uv=B_q.copy(), parent_kind="SIC", parent_idx=sic_idx,
            seg_global=sis_idx, end=1,
            dedup_key=("SIC", sic_idx, q_dp, 2 if f == 1 else 1),
            other_uv=A_q.copy(),
        ))
    return out


# ── union-find over SP-identity (G5) ─────────────────────────────────────────

def _build_components(
    bcs: "list[BoundaryCurve]",
    ccs: "list[ContourCurve]",
    sics: "list[SelfIntersectingCurve]",
    mesh: "Mesh",
    css: np.ndarray,
    sis_pairs: np.ndarray,
    splits: SplitArrays,
) -> list[int]:
    """Return a parent[i] array of length n_curves = |bcs|+|ccs|+|sics|.

    Two curves are united iff they share an SP via any SPT (G5).
    """
    n_bc = len(bcs)
    n_cc = len(ccs)
    n_sic = len(sics)
    n = n_bc + n_cc + n_sic

    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    # SP -> set of curve_idx (across all kinds, with kind-prefix offsets).
    sp_to_curves: dict[int, list[int]] = {}

    def _note(curve_idx: int, spt_idx: int) -> None:
        if spt_idx < 0:
            return
        sp_idx = int(splits.spts[spt_idx][0])
        sp_to_curves.setdefault(sp_idx, []).append(curve_idx)

    # BC chains → mesh.edges[mesh.boundary_edge_idx[i]].split1/split2
    bnd_idx = mesh.boundary_edge_idx
    for bc_idx, bc in enumerate(bcs):
        cidx = bc_idx
        for signed in bc.edge_indices:
            i = abs(int(signed)) - 1
            e_idx = int(bnd_idx[i])
            edge = mesh.edges[e_idx]
            _note(cidx, int(edge["split1"]))
            _note(cidx, int(edge["split2"]))

    # CC chains → css[cs_idx].split1/split2
    for cc_idx, cc in enumerate(ccs):
        cidx = n_bc + cc_idx
        for signed in cc.cs_indices:
            cs_idx = abs(int(signed)) - 1
            cs = css[cs_idx]
            _note(cidx, int(cs["split1"]))
            _note(cidx, int(cs["split2"]))

    # SIC chains → sis_pairs[sis_idx].split1/split2
    for sic_idx, sic in enumerate(sics):
        cidx = n_bc + n_cc + sic_idx
        for signed in sic.sis_indices:
            sis_idx = abs(int(signed)) - 1
            row = sis_pairs[sis_idx]
            _note(cidx, int(row["split1"]))
            _note(cidx, int(row["split2"]))

    for curves in sp_to_curves.values():
        if len(curves) < 2:
            continue
        for c in curves[1:]:
            union(curves[0], c)

    return [find(i) for i in range(n)]


# ── vc_in / vc_out per spec lines 436-440 ────────────────────────────────────

def _vc_for_endpoint(
    sample: _Sample, T_proj_uv: np.ndarray,
    surface: "SurfaceParams", projection: "Projection",
) -> int:
    """Visibility change ∈ {-1, 0} for an HC endpoint.

    Spec lines 436-440: dispatch by `sample.parent_kind`.

    - BC : 0 (spec line 436).
    - CC : N = CS uv-normal oriented toward front sheet (axis·dS_N > 0
           under axis-toward-viewer). 0 if T·N > 0 else -1.
    - SIC: 0 if (axis·SN') * (T_3d·SN') > 0 else -1, with T_3d the
           differential of S applied to T_proj at sample.uv.
    """
    if sample.parent_kind == "BC":
        return 0

    if projection.mode != "ortho":
        raise NotImplementedError(
            "_vc_for_endpoint: perspective mode not yet supported"
        )

    axis = projection._axis
    u, v = float(sample.uv[0]), float(sample.uv[1])
    Su = np.asarray(surface.Su(u, v), dtype=float).reshape(3)
    Sv = np.asarray(surface.Sv(u, v), dtype=float).reshape(3)

    if sample.parent_kind == "CC":
        cs_dir = sample.cs_dir_uv
        if cs_dir is None:
            raise RuntimeError("_vc_for_endpoint: CC sample missing cs_dir_uv")
        N = np.array([-float(cs_dir[1]), float(cs_dir[0])], dtype=float)
        # Orient N toward the front sheet under axis-toward-viewer:
        # axis · (Su*N[0] + Sv*N[1]) > 0  ⇒  +N moves S closer to viewer.
        dS_N = N[0] * Su + N[1] * Sv
        if float(axis @ dS_N) < 0.0:
            N = -N
        return 0 if float(T_proj_uv @ N) > 0.0 else -1

    if sample.parent_kind == "SIC":
        other = sample.other_uv
        if other is None:
            raise RuntimeError("_vc_for_endpoint: SIC sample missing other_uv")
        T_3d = float(T_proj_uv[0]) * Su + float(T_proj_uv[1]) * Sv
        ou, ov = float(other[0]), float(other[1])
        SN_other = np.asarray(surface.SN(ou, ov), dtype=float).reshape(3)
        factor = float(axis @ SN_other) * float(T_3d @ SN_other)
        return 0 if factor > 0.0 else -1

    raise RuntimeError(f"_vc_for_endpoint: unknown kind {sample.parent_kind!r}")


# ── component pairwise distance ──────────────────────────────────────────────

def _argmin_pair(
    samples_i: list[_Sample], samples_j: list[_Sample],
) -> tuple[float, int, int]:
    """Return (min_dist, idx_i, idx_j) over the cartesian product.

    Plain Euclidean (G19) — no `close()`.
    """
    uv_i = np.stack([s.uv for s in samples_i], axis=0)   # (Ni, 2)
    uv_j = np.stack([s.uv for s in samples_j], axis=0)   # (Nj, 2)
    diff = uv_i[:, None, :] - uv_j[None, :, :]
    d2 = (diff * diff).sum(axis=-1)                       # (Ni, Nj)
    k = int(np.argmin(d2))
    Ni, Nj = d2.shape
    ki, kj = divmod(k, Nj)
    return float(np.sqrt(d2[ki, kj])), ki, kj


# ── O12 main entry ───────────────────────────────────────────────────────────

def build_helper_curves(
    bcs: "list[BoundaryCurve]",
    ccs: "list[ContourCurve]",
    sics: "list[SelfIntersectingCurve]",
    css: np.ndarray,
    sis_pairs: np.ndarray,
    cps: np.ndarray,
    dps: np.ndarray,
    mesh: "Mesh",
    splits: SplitArrays,
    projection: "Projection",
    surface: "SurfaceParams",
    domain,
) -> list[SubCurve]:
    """Bridge disconnected curve components with helper curves.

    Returns one `SubCurve(kind="HC", internal=[], is_closed=False)` per
    emitted HC; mutates `splits` (appends HA SPs and SPTs) and writes
    split1/split2 on the parent segments.

    `domain` is accepted for API symmetry with the roadmap signature; the
    mesh already carries it.
    """
    del domain  # mesh.domain is the source of truth; kept for API parity.

    n_bc = len(bcs)
    n_cc = len(ccs)
    n_sic = len(sics)
    n_curves = n_bc + n_cc + n_sic
    if n_curves == 0:
        return []

    # 1. Samples per curve.
    curve_samples: list[list[_Sample]] = []
    for bc_idx, bc in enumerate(bcs):
        curve_samples.append(_bc_samples(bc_idx, bc, mesh))
    for cc_idx, cc in enumerate(ccs):
        curve_samples.append(_cc_samples(cc_idx, cc, css, cps))
    for sic_idx, sic in enumerate(sics):
        curve_samples.append(_sic_samples(sic_idx, sic, sis_pairs, dps))

    # 2. Components via SP-identity (G5).
    roots = _build_components(bcs, ccs, sics, mesh, css, sis_pairs, splits)

    # Aggregate samples per component root, skipping curves with no samples.
    comp_samples: dict[int, list[_Sample]] = {}
    for c_idx, samples in enumerate(curve_samples):
        if not samples:
            continue
        r = roots[c_idx]
        comp_samples.setdefault(r, []).extend(samples)

    if len(comp_samples) <= 1:
        return []

    # 3. Bridging loop. Active components keyed by representative root.
    hcs: list[SubCurve] = []
    ha_dedup: dict[tuple, int] = {}   # sample.dedup_key -> SP index

    def _ensure_ha(sample: _Sample) -> int:
        key = sample.dedup_key
        if key in ha_dedup:
            return ha_dedup[key]
        uv = sample.uv
        xyz = np.asarray(
            surface.S(float(uv[0]), float(uv[1])), dtype=float,
        ).reshape(3)
        xy = projection.XY(xyz)
        sp_idx = splits.add_sp(uv=uv, xyz=xyz, xy=xy, sp_type="ha")
        ha_dedup[key] = sp_idx

        bary = 0.0 if sample.end == 0 else 1.0
        spt_idx = splits.add_spt(sp_idx=sp_idx, bary=bary, vis_chge=0)
        if sample.parent_kind == "BC":
            splits.attach_to_segment(
                mesh.edges, sample.seg_global, spt_idx, segment_label="BE",
            )
        elif sample.parent_kind == "CC":
            splits.attach_to_segment(
                css, sample.seg_global, spt_idx, segment_label="CS",
            )
        else:  # SIC
            splits.attach_to_segment(
                sis_pairs, sample.seg_global, spt_idx, segment_label="SIS",
            )
        return sp_idx

    while len(comp_samples) > 1:
        roots_list = list(comp_samples.keys())
        best: tuple[float, int, int, _Sample, _Sample] | None = None
        for i in range(len(roots_list)):
            for j in range(i + 1, len(roots_list)):
                ri, rj = roots_list[i], roots_list[j]
                d, ki, kj = _argmin_pair(comp_samples[ri], comp_samples[rj])
                if best is None or d < best[0]:
                    best = (d, ri, rj, comp_samples[ri][ki], comp_samples[rj][kj])
        assert best is not None  # len(comp_samples) > 1 ⇒ at least one pair
        _, ri, rj, qi, qj = best

        sp_qi = _ensure_ha(qi)
        sp_qj = _ensure_ha(qj)

        T_in = qj.uv - qi.uv
        T_out = qi.uv - qj.uv
        vc_in = _vc_for_endpoint(qi, T_in, surface, projection)
        vc_out = _vc_for_endpoint(qj, T_out, surface, projection)

        hcs.append(SubCurve(
            kind="HC", is_closed=False,
            start=sp_qi, end=sp_qj,
            internal=[], vc_in=vc_in, vc_out=vc_out, parent_idx=-1,
        ))

        # Merge component ri into rj.
        comp_samples[rj].extend(comp_samples[ri])
        del comp_samples[ri]

    return hcs
