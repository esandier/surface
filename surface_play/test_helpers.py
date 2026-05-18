"""Tests for surface_play.helpers (Layer O, step O12)."""

import numpy as np

from surface_play.contour import (
    build_contour_curves,
    build_contour_segments,
    find_contour_points,
)
from surface_play.curves import BoundaryCurve, build_bcs
from surface_play.domain import Domain
from surface_play.helpers import build_helper_curves
from surface_play.intersections import dp_dtype, sis_dtype
from surface_play.mesh import Mesh, build_mesh, edge_dtype, face_dtype
from surface_play.projection import Projection
from surface_play.splitting import SplitArrays, SubCurve, split_bcs_at_bcps
from surface_play.test_fixtures import paraboloid, torus


def _ortho_proj(surface, I=(1.0, 0.0, 0.0), J=(0.0, 1.0, 0.0)):
    return Projection(surface, I=list(I), J=list(J))


def _build_layer_c_o(surf, J, resolution=15):
    """Mesh + BCs + CPs + CSs + CCs at the given view."""
    mesh = build_mesh(surf.domain, surf, resolution=resolution, jitter=True, seed=42)
    bcs = build_bcs(mesh)
    proj = _ortho_proj(surf, J=J)
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    return mesh, bcs, proj, cps, css, ccs


def _empty_dps_sis():
    return (
        np.zeros(0, dtype=dp_dtype),
        np.zeros(0, dtype=sis_dtype),
    )


# ── Sub-assert 1: single-component starting state → empty HC list ────────────

def test_build_helper_curves_single_component_paraboloid_side():
    """Paraboloid side view (axis = -Y): silhouette is an open arc with both
    endpoints on boundary edges. After `split_bcs_at_bcps`, the BC and the
    CC share SPs at the 2 BCPs → single component → empty HC list.
    """
    surf = paraboloid(perturb=False)
    # I=(1,0,0), J=(0,0,1) → axis = (0,-1,0). Silhouette at v=0.
    mesh, bcs, proj, cps, css, ccs = _build_layer_c_o(surf, J=(0.0, 0.0, 1.0))

    assert len(bcs) >= 1
    assert len(ccs) >= 1
    assert int((cps["ptype"] == 4).sum()) >= 2  # ≥2 BCPs at ±u boundaries

    splits = SplitArrays()
    split_bcs_at_bcps(mesh, bcs, ccs, css, cps, splits, surf, proj)

    dps, sis_pairs = _empty_dps_sis()
    hcs = build_helper_curves(
        bcs, ccs, [], css, sis_pairs, cps, dps,
        mesh, splits, proj, surf, mesh.domain,
    )
    assert hcs == [], f"expected 0 HCs (single component), got {len(hcs)}"


# ── Sub-asserts 2, 3, 4, 6, 7: torus Z-view → 2 disconnected CCs → 1 HC ──────

def test_build_helper_curves_torus_multi_component():
    """Torus Z view: 2 closed CCs (at v=0 and v=π), no BC, no SISs.

    Components share no SPTs → 2 components → exactly 1 HC.
    """
    surf = torus(perturb=False)
    # I=(1,0,0), J=(0,1,0) → axis = (0,0,1). Silhouette at v=0 and v=π.
    mesh, bcs, proj, cps, css, ccs = _build_layer_c_o(surf, J=(0.0, 1.0, 0.0))

    assert len(bcs) == 0, "torus has no BC"
    assert len(ccs) >= 2, f"expected ≥2 CCs, got {len(ccs)}"

    n_ccs_before = len(ccs)

    splits = SplitArrays()
    # No splits to run pre-O12: no corners/BCPs/BDPs/CDPs/TPs/VPs are
    # forced for this view; CCs are bare closed loops.

    dps, sis_pairs = _empty_dps_sis()
    hcs = build_helper_curves(
        bcs, ccs, [], css, sis_pairs, cps, dps,
        mesh, splits, proj, surf, mesh.domain,
    )

    # (2) n_HCs == n_components_before - 1; here n_components_before == n_ccs.
    assert len(hcs) == n_ccs_before - 1, (
        f"expected {n_ccs_before - 1} HCs to bridge {n_ccs_before} CCs, "
        f"got {len(hcs)}"
    )

    sps = splits.sps_array()
    spts = splits.spts_array()

    # (3) Every freshly-added SP has type=='ha' (no other splits were run).
    assert (sps["type"] == "ha").all(), (
        f"expected all SPs to be 'ha', got types {set(sps['type'].tolist())}"
    )

    # (4) Every HC has internal == [] and is_closed == False.
    for k, hc in enumerate(hcs):
        assert hc.kind == "HC"
        assert hc.is_closed is False, f"HC {k} is_closed should be False"
        assert hc.internal == [], f"HC {k} internal should be empty"

    # (6) vc_in, vc_out ∈ {-1, 0} for every HC.
    for k, hc in enumerate(hcs):
        assert hc.vc_in in (-1, 0), f"HC {k}.vc_in = {hc.vc_in}"
        assert hc.vc_out in (-1, 0), f"HC {k}.vc_out = {hc.vc_out}"

    # (7) G5 identity: each HC's start SP equals the sp_idx of an SPT that
    #     is attached to a segment of the parent curve of qi (similarly end).
    #     For HCs whose qi is on a CC, the SPT must reside in css[*].split*.
    for k, hc in enumerate(hcs):
        # The SP at `hc.start` must be referenced by at least one SPT
        # attached to a CS (since on torus all parents are CCs).
        spt_for_start = [
            i for i, t in enumerate(splits.spts) if int(t[0]) == hc.start
        ]
        assert spt_for_start, f"HC {k} start SP has no SPT (G5 broken)"
        attached_on_css = False
        for spt_idx in spt_for_start:
            in_css = (
                (css["split1"] == spt_idx) | (css["split2"] == spt_idx)
            ).any()
            if in_css:
                attached_on_css = True
                break
        assert attached_on_css, (
            f"HC {k} start SPT not attached to any CS (G5: SP index integrity)"
        )


# ── Sub-assert 5: G6 dedup — synthetic 3-BC fixture ──────────────────────────

def _make_synthetic_mesh_3_bcs():
    """Build a Mesh with 3 disjoint boundary edges (no shared vertices).

    Vertex layout (uv):
        0 = (0, 0)      = hub (shared argmin from both other BCs to BC0)
        1 = (0, 10)
        2 = (1, 0)
        3 = (1.5, 0)
        4 = (0, -1)
        5 = (0, -2)

    BEs:
        e0 : vertices 0 → 1   (BC0)
        e1 : vertices 2 → 3   (BC1)
        e2 : vertices 4 → 5   (BC2)

    Distances:
        BC0↔BC1: argmin (0,0)↔(1,0)   = 1.0      → uses hub on BC0.
        BC0↔BC2: argmin (0,0)↔(0,-1)  = 1.0      → uses hub on BC0.
        BC1↔BC2: argmin (1,0)↔(0,-1)  = sqrt(2)  ≈ 1.414.

    Iteration 1 picks (BC0↔BC1, hub-on-BC0, (1,0)-on-BC1). Merge.
    Iteration 2 picks merged ↔ BC2 — argmin still includes hub-on-BC0
    ↔ (0,-1) at distance 1 (smaller than the (1,0)↔(0,-1) = sqrt(2)).
    So hub is reused → G6 dedup fires.
    """
    uv = np.array(
        [
            [0.0,  0.0],
            [0.0, 10.0],
            [1.0,  0.0],
            [1.5,  0.0],
            [0.0, -1.0],
            [0.0, -2.0],
        ],
        dtype=float,
    )
    edges = np.zeros(3, dtype=edge_dtype)
    for i, (p, q) in enumerate([(0, 1), (2, 3), (4, 5)]):
        edges[i]["p_idx"] = p
        edges[i]["q_idx"] = q
        edges[i]["p"] = uv[p]
        edges[i]["pq"] = uv[q] - uv[p]
        edges[i]["f"] = 0
        edges[i]["g"] = -1
        edges[i]["dir"] = [0.0, 0.0]
        edges[i]["flip"] = 1
        edges[i]["split1"] = -1
        edges[i]["split2"] = -1

    return Mesh(
        domain=Domain(type="rect", bounds=(-0.5, 2.0, -2.5, 11.0)),
        surface=paraboloid(perturb=False),
        uv=uv,
        tris=np.zeros((0, 3), dtype=np.int32),
        edges=edges,
        faces=np.zeros(0, dtype=face_dtype),
        SN=np.zeros((6, 3), dtype=float),
        xyz=np.zeros((6, 3), dtype=float),
        boundary_edge_idx=np.array([0, 1, 2], dtype=np.int32),
        corner_idx=np.array([], dtype=np.int32),
    )


def test_build_helper_curves_g6_dedup_three_bcs_share_hub():
    """G6: hub vertex on BC0 is the argmin target for both HCs → one HA SP
    on BC0 (referenced by both HC.start indices), not two."""
    mesh = _make_synthetic_mesh_3_bcs()

    # Token convention: signed (i+1) where i is index into boundary_edge_idx.
    bcs = [
        BoundaryCurve(edge_indices=np.array([1], dtype=np.intp), is_closed=False),
        BoundaryCurve(edge_indices=np.array([2], dtype=np.intp), is_closed=False),
        BoundaryCurve(edge_indices=np.array([3], dtype=np.intp), is_closed=False),
    ]

    surf = mesh.surface
    proj = _ortho_proj(surf)
    splits = SplitArrays()
    dps, sis_pairs = _empty_dps_sis()
    # Empty cps/css/sics — no CCs / SICs in this fixture.
    cps = np.zeros(0, dtype=np.dtype([
        ("e", "i4"), ("s", "f8"), ("uv", "2f8"), ("xyz", "3f8"),
        ("d", "2f8"), ("ptype", "u1"),
    ]))
    css = np.zeros(0, dtype=np.dtype([
        ("p_cp", "i4"), ("q_cp", "i4"), ("face", "i4"),
        ("split1", "i4"), ("split2", "i4"),
    ]))

    hcs = build_helper_curves(
        bcs, [], [], css, sis_pairs, cps, dps,
        mesh, splits, proj, surf, mesh.domain,
    )

    # 3 components → exactly 2 HCs.
    assert len(hcs) == 2, f"expected 2 HCs, got {len(hcs)}"

    # G6: both HC.start indices on the BC0 side must point to the SAME SP
    # (the hub at vertex 0). Equivalently: only ONE 'ha' SP sits on BC0's
    # hub vertex.
    sps = splits.sps_array()
    hub_uv = mesh.uv[0]
    hub_hits = [
        i for i in range(len(sps))
        if sps[i]["type"] == "ha" and np.allclose(sps[i]["uv"], hub_uv)
    ]
    assert len(hub_hits) == 1, (
        f"G6 dedup: hub should produce exactly 1 HA SP, got {len(hub_hits)} "
        f"(sps at hub: {hub_hits})"
    )
    hub_sp = hub_hits[0]

    # Both HCs must have their hub-side endpoint pointing to that one SP.
    hub_starts = [hc for hc in hcs if hc.start == hub_sp or hc.end == hub_sp]
    assert len(hub_starts) == 2, (
        f"G6 dedup: expected 2 HCs to reference the shared hub SP={hub_sp}, "
        f"got {len(hub_starts)} (HC starts/ends: "
        f"{[(hc.start, hc.end) for hc in hcs]})"
    )

    # The hub SP should be attached to e0 (BC0's edge) via exactly ONE SPT
    # (not duplicated despite being referenced by 2 HCs).
    e0 = mesh.edges[0]
    slots = [int(e0["split1"]), int(e0["split2"])]
    spt_for_hub = [
        s for s in slots if s >= 0 and int(splits.spts[s][0]) == hub_sp
    ]
    assert len(spt_for_hub) == 1, (
        f"G6 dedup: BC0's hub edge should have exactly 1 SPT pointing to "
        f"the hub SP, got {len(spt_for_hub)} (slots: {slots})"
    )
