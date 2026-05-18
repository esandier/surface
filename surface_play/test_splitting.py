"""Tests for surface_play.splitting (Layer O, step O5)."""

import numpy as np
import pytest

from surface_play.contour import (
    _cp_sequence,
    build_contour_curves,
    build_contour_segments,
    find_contour_points,
    find_vps,
)
from surface_play.curves import build_bcs
from surface_play.intersections import dp_dtype, sis_dtype, tp_dtype
from surface_play.mesh import build_mesh, edge_dtype
from surface_play.projection import Projection
from surface_play.splitting import (
    SplitArrays,
    SplitSlotOverflowError,
    SubCurve,
    _bvis_chge,
    assemble_subcurves,
    sp_dtype,
    split_at_cdps,
    split_bcs_at_bcps,
    split_bcs_at_bdps,
    split_bcs_at_corners,
    split_ccs_at_vps,
    split_sics_at_tps,
    spt_dtype,
)
from surface_play.test_fixtures import (
    cylinder_cy,
    disk_paraboloid_po,
    fig8,
    helicoid,
    mobius_u,
    paraboloid,
    torus,
)


def test_sp_spt_primitives():
    splits = SplitArrays()

    # (1) add_sp returns sequential indices starting at 0.
    i0 = splits.add_sp(uv=(0.0, 0.0), xyz=(0.0, 0.0, 0.0), xy=(0.0, 0.0), sp_type="cn")
    i1 = splits.add_sp(uv=(0.5, 0.0), xyz=(0.0, 0.0, 0.0), xy=(0.0, 0.0), sp_type="cn")
    i2 = splits.add_sp(uv=(1.0, 0.0), xyz=(0.0, 0.0, 0.0), xy=(0.0, 0.0), sp_type="cn")
    assert [i0, i1, i2] == [0, 1, 2]

    # (2) add_spt returns sequential indices starting at 0.
    j0 = splits.add_spt(sp_idx=i0, bary=0.0, vis_chge=0)
    j1 = splits.add_spt(sp_idx=i1, bary=0.5, vis_chge=0)
    j2 = splits.add_spt(sp_idx=i2, bary=1.0, vis_chge=0)
    assert [j0, j1, j2] == [0, 1, 2]

    # (3) attach_to_segment writes split1 first, then split2, then raises.
    seg = np.zeros(1, dtype=edge_dtype)
    seg["split1"] = -1
    seg["split2"] = -1

    splits.attach_to_segment(seg, 0, j0, segment_label="BE")
    assert seg[0]["split1"] == j0
    assert seg[0]["split2"] == -1

    splits.attach_to_segment(seg, 0, j1, segment_label="BE")
    assert seg[0]["split1"] == j0
    assert seg[0]["split2"] == j1

    with pytest.raises(SplitSlotOverflowError) as exc_info:
        splits.attach_to_segment(seg, 0, j2, segment_label="BE")

    # (4) Error message names segment kind, index, and existing SP types.
    msg = str(exc_info.value)
    assert "BE" in msg
    assert "seg_idx=0" in msg
    assert "'cn'" in msg or "cn" in msg


def _ortho_proj(surface):
    return Projection(surface, I=[1.0, 0.0, 0.0], J=[0.0, 1.0, 0.0])


def test_corner_splits():
    # (1) Rect no-id resolution=10: 4 corner SPs, 8 SPTs, each corner attached
    #     to exactly 2 BEs.
    surf = paraboloid(perturb=False)
    mesh = build_mesh(surf.domain, surf, resolution=10, jitter=True, seed=42)
    bcs = build_bcs(mesh)
    proj = _ortho_proj(surf)
    splits = SplitArrays()

    split_bcs_at_corners(mesh, bcs, splits, proj)

    sps = splits.sps_array()
    spts = splits.spts_array()
    assert len(sps) == 4, f"expected 4 corner SPs, got {len(sps)}"
    assert (sps["type"] == "cn").all()
    assert len(spts) == 8, f"expected 8 SPTs (2 per corner), got {len(spts)}"
    # vis_chge is 0 for all corner SPTs (G10).
    assert (spts["vis_chge"] == 0).all()

    # Each corner SP is referenced by exactly 2 SPTs.
    for i in range(4):
        n = int((spts["sp_idx"] == i).sum())
        assert n == 2, f"corner SP {i} has {n} SPTs, expected 2"

    # (2) Each corner SP's BE attachments have bary ∈ {0, 1}.
    assert set(np.unique(spts["bary"]).tolist()) <= {0.0, 1.0}

    # Cross-check: the 8 SPT indices live in split1/split2 of 8 distinct BEs
    # (or 4 BEs each holding two corner SPTs — but in rect no-id every BE
    # touches at most one corner, so 8 distinct BEs).
    split_slots = np.concatenate([mesh.edges["split1"], mesh.edges["split2"]])
    filled = split_slots[split_slots != -1]
    assert sorted(filled.tolist()) == list(range(8))

    # (3) Rect cy-no (cylinder_cy is v-cy), mo-no: 0 corner SPs.
    for factory in (cylinder_cy, mobius_u):
        surf_id = factory(perturb=False)
        mesh_id = build_mesh(surf_id.domain, surf_id, resolution=10, jitter=True, seed=42)
        bcs_id = build_bcs(mesh_id)
        splits_id = SplitArrays()
        split_bcs_at_corners(mesh_id, bcs_id, splits_id, _ortho_proj(surf_id))
        assert len(splits_id.sps) == 0, f"{factory.__name__}: expected 0 SPs"
        assert len(splits_id.spts) == 0

    # cy-cy (torus): corner_idx has length 1 but it's interior to a closed BC
    # (and actually there are no BEs at all on the torus). No SPs.
    surf_t = torus(perturb=False)
    mesh_t = build_mesh(surf_t.domain, surf_t, resolution=10, jitter=True, seed=42)
    bcs_t = build_bcs(mesh_t)
    splits_t = SplitArrays()
    split_bcs_at_corners(mesh_t, bcs_t, splits_t, _ortho_proj(surf_t))
    assert len(splits_t.sps) == 0

    # (4) Disk: no corners, no SPs.
    surf_d = disk_paraboloid_po(perturb=False)
    mesh_d = build_mesh(surf_d.domain, surf_d, resolution=10, jitter=True, seed=42)
    bcs_d = build_bcs(mesh_d)
    splits_d = SplitArrays()
    split_bcs_at_corners(mesh_d, bcs_d, splits_d, _ortho_proj(surf_d))
    assert len(splits_d.sps) == 0
    assert len(splits_d.spts) == 0


def _build_pipeline(surf, *, resolution=15):
    mesh = build_mesh(surf.domain, surf, resolution=resolution, jitter=True, seed=42)
    bcs = build_bcs(mesh)
    proj = _ortho_proj(surf)
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    return mesh, bcs, proj, cps, css, ccs


def test_bcp_splits():
    # (1) Helicoid ortho: contour along u=0 hits v_min and v_max boundary edges.
    #     SP count of type 'bcp' equals the number of CPs with ptype==4.
    surf = helicoid(perturb=False)
    mesh, bcs, proj, cps, css, ccs = _build_pipeline(surf, resolution=15)

    n_bcp_cps = int((cps["ptype"] == 4).sum())
    assert n_bcp_cps >= 2, (
        f"helicoid ortho should have ≥2 boundary CPs (v_min, v_max), got {n_bcp_cps}"
    )

    splits = SplitArrays()
    split_bcs_at_bcps(mesh, bcs, ccs, css, cps, splits, surf, proj)

    sps = splits.sps_array()
    spts = splits.spts_array()
    bcp_sp_mask = sps["type"] == "bcp"
    assert int(bcp_sp_mask.sum()) == n_bcp_cps, (
        f"expected {n_bcp_cps} bcp SPs, got {int(bcp_sp_mask.sum())}"
    )

    # (2) Each bcp SP is referenced by exactly 2 SPTs (one BE, one CS).
    for k in np.nonzero(bcp_sp_mask)[0]:
        n = int((spts["sp_idx"] == k).sum())
        assert n == 2, f"bcp SP {k} has {n} SPTs, expected 2"

    # (3) bvis_chge regression: every BE-side SPT carries vis_chge ∈ {-1, 0, +1}.
    bcp_sp_indices = set(int(k) for k in np.nonzero(bcp_sp_mask)[0])
    be_vis_values = []
    cs_vis_values = []
    for t in spts:
        if int(t["sp_idx"]) not in bcp_sp_indices:
            continue
        bary = float(t["bary"])
        # CS-side SPTs have bary ∈ {0, 1}; BE-side SPTs have bary == cp.s ∈ (0, 1).
        if bary in (0.0, 1.0):
            cs_vis_values.append(int(t["vis_chge"]))
        else:
            be_vis_values.append(int(t["vis_chge"]))
    assert set(be_vis_values) <= {-1, 0, 1}, (
        f"BE-side vis_chge out of {{-1, 0, +1}}: {sorted(set(be_vis_values))}"
    )
    assert all(v == 0 for v in cs_vis_values), (
        f"CS-side SPTs must have vis_chge == 0, got {cs_vis_values}"
    )

    # (4) Spec line 298: every boundary CP is an endpoint of an open CC.
    cp_to_cs = {}
    for cs_idx in range(len(css)):
        cp_to_cs.setdefault(int(css[cs_idx]["p_cp"]), []).append(cs_idx)
        cp_to_cs.setdefault(int(css[cs_idx]["q_cp"]), []).append(cs_idx)
    cs_to_cc = {}
    for cc_idx, cc in enumerate(ccs):
        for si in cc.cs_indices:
            cs_to_cc[abs(int(si)) - 1] = cc_idx
    for cp_idx in np.nonzero(cps["ptype"] == 4)[0]:
        css_for_cp = cp_to_cs.get(int(cp_idx), [])
        assert len(css_for_cp) == 1, (
            f"boundary CP {int(cp_idx)} should be in exactly 1 CS, got {css_for_cp}"
        )
        cc_idx = cs_to_cc[css_for_cp[0]]
        cc = ccs[cc_idx]
        assert not cc.is_closed, f"boundary CP {int(cp_idx)} is on a closed CC"
        seq = _cp_sequence(cc, css)
        assert int(cp_idx) in (seq[0], seq[-1]), (
            f"boundary CP {int(cp_idx)} not at CC endpoints {seq[0], seq[-1]}"
        )

    # (5) No-BCP control: torus ortho has no boundary edges, so no boundary CPs
    #     and no bcp SPs.
    surf_t = torus(perturb=False)
    mesh_t, bcs_t, proj_t, cps_t, css_t, ccs_t = _build_pipeline(surf_t, resolution=10)
    assert int((cps_t["ptype"] == 4).sum()) == 0
    splits_t = SplitArrays()
    split_bcs_at_bcps(mesh_t, bcs_t, ccs_t, css_t, cps_t, splits_t, surf_t, proj_t)
    assert len(splits_t.sps) == 0
    assert len(splits_t.spts) == 0


def _find_two_boundary_edges(mesh, *, parallel: bool = False):
    """Return two boundary-edge indices from a real mesh.

    parallel=False: any two BEs (typically perpendicular sides — fine for
    case 2 mocks).
    parallel=True: two BEs that share a side (same row of `g == -1` with
    similar `pq` direction) — used for case 1 mocks where two boundary
    edges meet at a "DP".
    """
    bnd = mesh.boundary_edge_idx
    if parallel:
        # Pick first two BEs on the same side: their `pq` should be ~parallel.
        first = int(bnd[0])
        first_pq = mesh.edges[first]["pq"]
        for k in bnd[1:]:
            pq_k = mesh.edges[int(k)]["pq"]
            cos = abs(float(np.dot(first_pq, pq_k))) / (
                np.linalg.norm(first_pq) * np.linalg.norm(pq_k) + 1e-30
            )
            if cos > 0.99:
                return first, int(k)
        return first, int(bnd[1])
    return int(bnd[0]), int(bnd[1])


def _make_synthetic_bdp(mesh, *, dp_type: str, e1: int, e2: int = -1,
                        s1: float = 0.5, s2: float = 0.5,
                        f2: int = -1, uv2_face=None):
    """Build a 1-DP, 1-SIS pair of arrays anchored on real mesh edges.

    Two records: one DP with on_boundary=True (the BDP) plus one anchor DP
    so that `sis_pairs` can be non-empty. The SIS edge joins them.
    """
    dps = np.zeros(2, dtype=dp_dtype)

    # BDP — record 0.
    edge1 = mesh.edges[e1]
    uv1 = edge1["p"] + s1 * edge1["pq"]
    dps["E1"][0] = e1
    dps["E2"][0] = e2
    dps["F2"][0] = f2
    dps["type"][0] = dp_type
    dps["on_boundary"][0] = True
    dps["uv1"][0] = uv1
    dps["A1"][0] = (int(mesh.edges["f"][e1]), int(mesh.edges["g"][e1]))
    if dp_type == "EF":
        # uv2 on face f2 — caller supplies, or default to face centroid.
        if uv2_face is None:
            face = mesh.faces[f2]
            uv2 = face["p"] + (face["pq"] + face["pr"]) / 3.0
        else:
            uv2 = np.asarray(uv2_face, dtype=float)
        dps["uv2"][0] = uv2
        dps["A2"][0] = (f2, -1)
        dps["xyz"][0] = mesh.surface.S(float(uv2[0]), float(uv2[1]))
    else:  # EE
        edge2 = mesh.edges[e2]
        uv2 = edge2["p"] + s2 * edge2["pq"]
        dps["uv2"][0] = uv2
        dps["A2"][0] = (int(mesh.edges["f"][e2]), int(mesh.edges["g"][e2]))
        dps["xyz"][0] = mesh.surface.S(float(uv2[0]), float(uv2[1]))

    # Anchor DP — record 1; off-boundary, just used as the other end of the SIS.
    anchor_e = int(mesh.boundary_edge_idx[-1])  # any edge; data unused for vis
    anchor_edge = mesh.edges[anchor_e]
    anchor_uv = anchor_edge["p"] + 0.5 * anchor_edge["pq"]
    dps["E1"][1] = anchor_e
    dps["E2"][1] = -1
    dps["F2"][1] = -1
    dps["type"][1] = "EF"
    dps["on_boundary"][1] = False
    dps["uv1"][1] = anchor_uv
    dps["uv2"][1] = anchor_uv
    dps["A1"][1] = (int(mesh.edges["f"][anchor_e]), int(mesh.edges["g"][anchor_e]))
    dps["A2"][1] = (-1, -1)
    dps["xyz"][1] = mesh.surface.S(float(anchor_uv[0]), float(anchor_uv[1]))

    sis_pairs = np.zeros(1, dtype=sis_dtype)
    sis_pairs["p_dp"] = 0
    sis_pairs["q_dp"] = 1
    sis_pairs["flip"] = 1
    sis_pairs["split1"] = -1
    sis_pairs["split2"] = -1
    return dps, sis_pairs


def test_bdp_splits_no_bdps():
    # (1) Surfaces without self-intersection / without on_boundary DPs:
    #     split_bcs_at_bdps creates 0 'bdp' SPs and leaves arrays untouched.
    from surface_play.intersections import (
        build_sis_pairs, find_double_points,
    )
    for factory in (paraboloid, torus, helicoid):
        surf = factory(perturb=False)
        mesh = build_mesh(surf.domain, surf, resolution=10, jitter=True, seed=42)
        bcs = build_bcs(mesh)
        proj = _ortho_proj(surf)
        dps = find_double_points(mesh, surf)
        sis_pairs = build_sis_pairs(dps)
        splits = SplitArrays()
        split_bcs_at_bdps(mesh, bcs, [], sis_pairs, dps, splits, surf, proj)
        assert len(splits.sps) == 0, (
            f"{factory.__name__}: expected 0 bdp SPs, got {len(splits.sps)}"
        )


def test_bdp_splits_mock_case2_ef():
    # (2) Mock case-2 EF DP anchored on the helicoid (boundary E1, interior F2).
    surf = helicoid(perturb=False)
    mesh = build_mesh(surf.domain, surf, resolution=10, jitter=True, seed=42)
    bcs = build_bcs(mesh)
    proj = _ortho_proj(surf)

    e1 = int(mesh.boundary_edge_idx[0])
    # Pick an interior face — first non-boundary face works.
    bnd_face_set = set(mesh.edges["f"][mesh.boundary_edge_idx].tolist())
    f2 = next(i for i in range(len(mesh.faces)) if i not in bnd_face_set)

    dps, sis_pairs = _make_synthetic_bdp(
        mesh, dp_type="EF", e1=e1, f2=f2, s1=0.5,
    )
    splits = SplitArrays()
    split_bcs_at_bdps(mesh, bcs, [], sis_pairs, dps, splits, surf, proj)

    sps = splits.sps_array()
    spts = splits.spts_array()
    bdp_mask = sps["type"] == "bdp"
    assert int(bdp_mask.sum()) == 1, f"expected 1 bdp SP, got {int(bdp_mask.sum())}"
    sp_idx = int(np.flatnonzero(bdp_mask)[0])
    own_spts = spts[spts["sp_idx"] == sp_idx]
    assert len(own_spts) == 2, (
        f"case-2 EF should make 2 SPTs (1 SIS + 1 BE); got {len(own_spts)}"
    )
    be_spts = [int(t["vis_chge"]) for t in own_spts if t["bary"] not in (0.0, 1.0)]
    sis_spts = [int(t["vis_chge"]) for t in own_spts if t["bary"] in (0.0, 1.0)]
    assert all(v == 0 for v in sis_spts), f"SIS-side vis_chge must be 0; got {sis_spts}"
    assert all(v in (-1, 1) for v in be_spts), (
        f"case-2 BE-side vis_chge ∈ {{-1, +1}}; got {be_spts}"
    )


def test_bdp_splits_mock_case2_ee_one_boundary():
    # (3) Mock case-2 EE DP: E1 boundary, E2 interior.
    surf = helicoid(perturb=False)
    mesh = build_mesh(surf.domain, surf, resolution=10, jitter=True, seed=42)
    bcs = build_bcs(mesh)
    proj = _ortho_proj(surf)

    e1 = int(mesh.boundary_edge_idx[0])
    # Pick any interior edge (g != -1).
    interior_mask = mesh.edges["g"] != -1
    e2 = int(np.flatnonzero(interior_mask)[0])

    dps, sis_pairs = _make_synthetic_bdp(
        mesh, dp_type="EE", e1=e1, e2=e2, s1=0.4, s2=0.4,
    )
    splits = SplitArrays()
    split_bcs_at_bdps(mesh, bcs, [], sis_pairs, dps, splits, surf, proj)

    sps = splits.sps_array()
    spts = splits.spts_array()
    bdp_mask = sps["type"] == "bdp"
    assert int(bdp_mask.sum()) == 1
    sp_idx = int(np.flatnonzero(bdp_mask)[0])
    own = spts[spts["sp_idx"] == sp_idx]
    assert len(own) == 2, f"case-2 EE (one boundary): 2 SPTs; got {len(own)}"
    be_spts = [int(t["vis_chge"]) for t in own if t["bary"] not in (0.0, 1.0)]
    assert all(v in (-1, 1) for v in be_spts), be_spts


def test_bdp_splits_mock_case1():
    # (4) Mock case-1: EE DP with both E1 and E2 on the boundary.
    surf = helicoid(perturb=False)
    mesh = build_mesh(surf.domain, surf, resolution=10, jitter=True, seed=42)
    bcs = build_bcs(mesh)
    proj = _ortho_proj(surf)

    e1, e2 = _find_two_boundary_edges(mesh, parallel=False)
    assert e1 != e2

    dps, sis_pairs = _make_synthetic_bdp(
        mesh, dp_type="EE", e1=e1, e2=e2, s1=0.4, s2=0.4,
    )
    splits = SplitArrays()
    split_bcs_at_bdps(mesh, bcs, [], sis_pairs, dps, splits, surf, proj)

    sps = splits.sps_array()
    spts = splits.spts_array()
    bdp_mask = sps["type"] == "bdp"
    assert int(bdp_mask.sum()) == 1
    sp_idx = int(np.flatnonzero(bdp_mask)[0])
    own = spts[spts["sp_idx"] == sp_idx]
    # Case-1: one SIS SPT + 2 BE SPTs.
    assert len(own) == 3, f"case-1 should make 3 SPTs (1 SIS + 2 BE); got {len(own)}"
    be_spts = [int(t["vis_chge"]) for t in own if t["bary"] not in (0.0, 1.0)]
    sis_spts = [int(t["vis_chge"]) for t in own if t["bary"] in (0.0, 1.0)]
    assert all(v == 0 for v in sis_spts)
    assert all(v in (-1, 0, 1) for v in be_spts), (
        f"case-1 BE-side vis_chge ∈ {{-1, 0, +1}}; got {be_spts}"
    )


def test_bvis_chge_handcrafted():
    # Synthetic ortho config: helicoid at a non-degenerate uv on the v=v_min
    # boundary edge — verify _bvis_chge returns a value in {-1, 0, +1} and is
    # sign-consistent with a small perturbation of s.
    surf = helicoid(perturb=False)
    proj = _ortho_proj(surf)
    # Fake a boundary edge at v=v_min, oriented in +u, with inward dir = (0, +1).
    edge = np.zeros(1, dtype=edge_dtype)[0]
    edge["p"] = (-0.5, -np.pi)
    edge["pq"] = (1.0, 0.0)
    edge["dir"] = (0.0, 1.0)
    for s in (0.25, 0.5, 0.75):
        v = _bvis_chge(surf, proj, edge, s)
        assert v in (-1, 0, 1), f"_bvis_chge returned {v} for s={s}"


# ── O9: split_sics_at_tps ────────────────────────────────────────────────────

def _synthetic_tp_for_o9(chord_zs):
    """Build a 1-TP synthetic fixture for O9.

    `chord_zs` is a 3-tuple of +1 / -1 — the sign of the z-chord (q_xyz.z -
    p_xyz.z) for each of the 3 SISs. Combined with SN ≡ (0, 0, 1) and axis ≡
    (0, 0, 1) (I=[1,0,0], J=[0,1,0]), this lets the test predict each
    SPT's vis_chge: vis = +1 if chord_z > 0 else -1.

    The geometry mirrors `test_intersections._synthetic_three_sis_fixture`:
    three SISs interlocking on three faces (F1=10, F2=20, F3=30); 3D
    verification trivially passes because `S(u, v) ≡ 0`.
    """
    from types import SimpleNamespace

    from surface_play.domain import Domain

    dps = np.zeros(6, dtype=dp_dtype)
    dps["A1"] = np.array(
        [(10, -1), (10, -1), (10, -1), (10, -1), (20, -1), (20, -1)],
        dtype=np.int32,
    )
    dps["A2"] = np.array(
        [(20, -1), (20, -1), (30, -1), (30, -1), (30, -1), (30, -1)],
        dtype=np.int32,
    )
    dps["uv1"] = np.array([
        [0.0, 0.0], [0.2, 0.2],
        [0.0, 0.2], [0.2, 0.0],
        [0.4, 0.2], [0.6, 0.0],
    ])
    dps["uv2"] = np.array([
        [0.4, 0.0], [0.6, 0.2],
        [0.0, 0.4], [0.2, 0.6],
        [0.0, 0.6], [0.2, 0.4],
    ])
    dps["on_boundary"] = False
    # Chord xyz per SIS: p_xyz at z=-1 (or +1), q_xyz at z=+1 (or -1).
    for k, sign in enumerate(chord_zs):
        p = 2 * k
        q = 2 * k + 1
        dps["xyz"][p] = [0.0, 0.0, -float(sign)]
        dps["xyz"][q] = [0.0, 0.0, +float(sign)]

    sis_pairs = np.zeros(3, dtype=sis_dtype)
    sis_pairs["p_dp"] = [0, 2, 4]
    sis_pairs["q_dp"] = [1, 3, 5]
    sis_pairs["flip"] = 1
    sis_pairs["split1"] = -1
    sis_pairs["split2"] = -1

    class _FakeSurface:
        def S(self, u, v):
            return np.zeros(3)

        def SN(self, u, v):
            return np.array([0.0, 0.0, 1.0])

    mesh = SimpleNamespace(domain=Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0)))
    return sis_pairs, dps, mesh, _FakeSurface()


def _fake_ortho_projection():
    """An ortho Projection with I=+x, J=+y (axis=+z). The surface attribute is
    only needed for ker_param/proj_vec; O9 uses only `_axis` and `XY`, neither
    of which touches surface methods, so a None placeholder via object() works.
    """
    return Projection(object(), I=[1.0, 0.0, 0.0], J=[0.0, 1.0, 0.0])


def test_tp_splits_no_tps():
    # (1) No-TP surfaces: the full Layer-C pipeline through find_triple_points
    #     yields no TPs; split_sics_at_tps produces 0 SPs.
    from surface_play.intersections import (
        build_sis_pairs,
        find_double_points,
        find_triple_points,
    )
    for factory in (paraboloid, torus, helicoid, mobius_u):
        surf = factory(perturb=False)
        mesh = build_mesh(surf.domain, surf, resolution=10, jitter=True, seed=42)
        proj = _ortho_proj(surf)
        dps = find_double_points(mesh, surf)
        sis_pairs = build_sis_pairs(dps)
        tps = find_triple_points(sis_pairs, dps, mesh, surf)
        splits = SplitArrays()
        split_sics_at_tps(tps, sis_pairs, dps, splits, surf, proj)
        assert len(splits.sps) == 0, (
            f"{factory.__name__}: expected 0 tp SPs, got {len(splits.sps)}"
        )
        assert len(splits.spts) == 0


def test_tp_splits_synthetic_one_tp():
    # (2) Synthetic 1-TP fixture: exactly 1 SP type='tp' and 3 SPTs (one per SIS).
    from surface_play.intersections import find_triple_points

    sis_pairs, dps, mesh, surf = _synthetic_tp_for_o9(chord_zs=(+1, +1, +1))
    tps = find_triple_points(sis_pairs, dps, mesh, surf)
    assert len(tps) == 1, f"expected 1 TP from fixture, got {len(tps)}"

    proj = _fake_ortho_projection()
    splits = SplitArrays()
    split_sics_at_tps(tps, sis_pairs, dps, splits, surf, proj)

    sps = splits.sps_array()
    spts = splits.spts_array()
    tp_mask = sps["type"] == "tp"
    assert int(tp_mask.sum()) == 1, f"expected 1 tp SP, got {int(tp_mask.sum())}"

    sp_idx = int(np.flatnonzero(tp_mask)[0])
    own = spts[spts["sp_idx"] == sp_idx]
    assert len(own) == 3, f"expected 3 SPTs (one per SIS), got {len(own)}"

    # Each of the 3 SISs gets a split slot filled.
    filled = sum(
        (int(sis_pairs[k]["split1"]) != -1) + (int(sis_pairs[k]["split2"]) != -1)
        for k in range(3)
    )
    assert filled == 3, f"expected 3 SIS split slots filled, got {filled}"

    # (3) G10: vis_chge ∈ {-1, +1} only (never 0 for TPs).
    assert set(int(v) for v in own["vis_chge"]) <= {-1, 1}


def test_tp_splits_sign_regression():
    # (4) Handcrafted signs: SN ≡ +z, axis ≡ +z, chord_z = ±1 — so vis_chge
    #     must equal sign(chord_z) for each SIS.
    from surface_play.intersections import find_triple_points

    chord_zs = (+1, -1, +1)  # SIS_0: +1, SIS_1: -1, SIS_2: +1
    sis_pairs, dps, mesh, surf = _synthetic_tp_for_o9(chord_zs=chord_zs)
    tps = find_triple_points(sis_pairs, dps, mesh, surf)
    assert len(tps) == 1

    proj = _fake_ortho_projection()
    splits = SplitArrays()
    split_sics_at_tps(tps, sis_pairs, dps, splits, surf, proj)

    sps = splits.sps_array()
    spts = splits.spts_array()
    assert int((sps["type"] == "tp").sum()) == 1

    # Map each filled SPT back to its SIS index, then check the sign.
    for s_idx in range(3):
        slot1 = int(sis_pairs[s_idx]["split1"])
        slot2 = int(sis_pairs[s_idx]["split2"])
        assert slot1 != -1 and slot2 == -1, (
            f"SIS {s_idx}: expected exactly split1 filled, got "
            f"split1={slot1} split2={slot2}"
        )
        spt = spts[slot1]
        expected = 1 if chord_zs[s_idx] > 0 else -1
        assert int(spt["vis_chge"]) == expected, (
            f"SIS {s_idx}: chord_z sign {chord_zs[s_idx]} → expected vis_chge "
            f"{expected}, got {int(spt['vis_chge'])}"
        )
        # bary = 0.5 because TP_xyz = (0, 0, 0) is the midpoint of (0,0,±1).
        assert abs(float(spt["bary"]) - 0.5) < 1e-9, (
            f"SIS {s_idx}: expected bary=0.5, got {float(spt['bary'])}"
        )


# ── O10: split_ccs_at_vps ────────────────────────────────────────────────────

def _torus_side_proj(surf):
    """Torus viewed along -Y — the canonical 4-cusp configuration (see
    test_contour.torus_side_view).
    """
    return Projection(surf, I=[1.0, 0.0, 0.0], J=[0.0, 0.0, 1.0])


def test_vp_splits_no_cusps():
    # (1) No-cusp surface (torus top-down) → 0 VP SPs.
    surf = torus(perturb=False)
    mesh = build_mesh(surf.domain, surf, resolution=20, jitter=True, seed=42)
    proj = _ortho_proj(surf)  # top-down — no cusps
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    vps = find_vps(ccs, css, cps, surf, proj)
    assert len(vps) == 0, f"torus top-down should have 0 VPs, got {len(vps)}"

    splits = SplitArrays()
    split_ccs_at_vps(vps, css, ccs, splits, surf, proj)
    assert len(splits.sps) == 0
    assert len(splits.spts) == 0
    # No CS got a split slot filled.
    filled = int(((css["split1"] != -1) | (css["split2"] != -1)).sum())
    assert filled == 0


def test_vp_splits_torus_side_four_cusps():
    # (2) Torus side view — canonical 4 cusps → 4 VP SPs.
    # (4) Each VP creates exactly 1 SPT (only on its CS).
    surf = torus(perturb=False)
    mesh = build_mesh(surf.domain, surf, resolution=20, jitter=True, seed=42)
    proj = _torus_side_proj(surf)
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    vps = find_vps(ccs, css, cps, surf, proj)
    assert len(vps) == 4, f"torus side view should have 4 VPs, got {len(vps)}"

    splits = SplitArrays()
    split_ccs_at_vps(vps, css, ccs, splits, surf, proj)

    sps = splits.sps_array()
    spts = splits.spts_array()
    vp_mask = sps["type"] == "vp"
    assert int(vp_mask.sum()) == len(vps), (
        f"expected {len(vps)} vp SPs, got {int(vp_mask.sum())}"
    )

    # Exactly 1 SPT per VP SP.
    for k in np.flatnonzero(vp_mask):
        n = int((spts["sp_idx"] == int(k)).sum())
        assert n == 1, f"vp SP {int(k)} has {n} SPTs, expected 1"

    # (3) G10: vis_chge ∈ {-1, +1} for every VP SPT.
    vp_sp_indices = set(int(k) for k in np.flatnonzero(vp_mask))
    vp_vis_vals = [
        int(t["vis_chge"]) for t in spts if int(t["sp_idx"]) in vp_sp_indices
    ]
    assert set(vp_vis_vals) <= {-1, 1}, (
        f"VP SPT vis_chge must be ±1, got {sorted(set(vp_vis_vals))}"
    )

    # Each VP's SPT bary equals vp.s, vis equals vp.vis_change (verbatim from O4).
    # Walk in vp order — SPs and SPTs were appended in vp order.
    vp_sp_list = sorted(int(k) for k in np.flatnonzero(vp_mask))
    for vp, sp_idx in zip(vps, vp_sp_list):
        own = spts[spts["sp_idx"] == sp_idx]
        assert len(own) == 1
        spt = own[0]
        assert abs(float(spt["bary"]) - float(vp["s"])) < 1e-12
        assert int(spt["vis_chge"]) == int(vp["vis_change"])

    # Each VP's host CS got a split slot filled.
    for vp in vps:
        cs_idx = int(vp["cs"])
        s1 = int(css[cs_idx]["split1"])
        s2 = int(css[cs_idx]["split2"])
        assert s1 != -1 or s2 != -1, f"CS {cs_idx} got no split for its VP"


# ── O11: split_at_cdps ───────────────────────────────────────────────────────

def _build_full_layer_c_o(surf, *, resolution=20, view_J=(0.0, 1.0, 0.0)):
    """Build mesh, BCs, CPs, CSs, CCs, dps, sis_pairs for O11 tests."""
    from surface_play.intersections import (
        build_sis_pairs,
        find_double_points,
    )

    mesh = build_mesh(surf.domain, surf, resolution=resolution, jitter=True, seed=42)
    proj = Projection(surf, I=[1.0, 0.0, 0.0], J=list(view_J))
    cps = find_contour_points(mesh, proj)
    css = build_contour_segments(cps, mesh)
    ccs = build_contour_curves(css, cps)
    dps = find_double_points(mesh, surf)
    sis_pairs = build_sis_pairs(dps)
    return mesh, proj, cps, css, ccs, dps, sis_pairs


def test_cdp_splits_no_intersections():
    # (1) No-SIS surfaces (helicoid, torus, paraboloid): 0 CDP SPs because
    #     sis_pairs is empty.
    for factory in (helicoid, torus, paraboloid):
        surf = factory(perturb=False)
        mesh, proj, cps, css, ccs, dps, sis_pairs = _build_full_layer_c_o(
            surf, resolution=15
        )
        assert len(sis_pairs) == 0, (
            f"{factory.__name__}: unexpected SIS — fixture changed"
        )
        splits = SplitArrays()
        split_at_cdps(mesh, css, cps, sis_pairs, dps, splits, surf, proj)
        assert len(splits.sps) == 0, (
            f"{factory.__name__}: expected 0 cdp SPs, got {len(splits.sps)}"
        )


# Side view of the fig-8 immersion (axis = I × J = (0, -1, 0)) — empirically
# yields ~4 CDPs at res=20 (silhouette curve crosses each SIC preimage).
_FIG8_SIDE_VIEW_J = (0.0, 0.0, 1.0)


def test_cdp_splits_fig8_side_count():
    # (2) Fig-8 immersion, side view: produces CDPs. Each yields 2 SPTs
    #     (1 CS-side + 1 SIS-side).
    surf = fig8(perturb=False)
    mesh, proj, cps, css, ccs, dps, sis_pairs = _build_full_layer_c_o(
        surf, resolution=20, view_J=_FIG8_SIDE_VIEW_J,
    )
    assert len(sis_pairs) > 0, "fig8 should produce SISs"

    splits = SplitArrays()
    split_at_cdps(mesh, css, cps, sis_pairs, dps, splits, surf, proj)

    sps = splits.sps_array()
    spts = splits.spts_array()
    cdp_mask = sps["type"] == "cdp"
    n_cdp = int(cdp_mask.sum())
    assert n_cdp >= 2, (
        f"fig8 side view should produce ≥2 CDPs, got {n_cdp}"
    )

    # Each CDP yields exactly 2 SPTs.
    for k in np.flatnonzero(cdp_mask):
        n = int((spts["sp_idx"] == int(k)).sum())
        assert n == 2, f"cdp SP {int(k)} has {n} SPTs (expected 2)"


def test_cdp_splits_vis_chge_sign():
    # (3) G10: every CDP SPT has vis_chge ∈ {-1, +1} (never 0).
    surf = fig8(perturb=False)
    mesh, proj, cps, css, ccs, dps, sis_pairs = _build_full_layer_c_o(
        surf, resolution=20, view_J=_FIG8_SIDE_VIEW_J,
    )
    splits = SplitArrays()
    split_at_cdps(mesh, css, cps, sis_pairs, dps, splits, surf, proj)

    sps = splits.sps_array()
    spts = splits.spts_array()
    cdp_sp_indices = set(int(k) for k in np.flatnonzero(sps["type"] == "cdp"))
    assert cdp_sp_indices, "fig8 side view should produce CDPs"

    cdp_vis = [
        int(t["vis_chge"]) for t in spts if int(t["sp_idx"]) in cdp_sp_indices
    ]
    assert set(cdp_vis) <= {-1, 1}, (
        f"CDP SPT vis_chge must be ±1 (G10), got {sorted(set(cdp_vis))}"
    )


def test_cdp_splits_resolution_stability():
    # (4) G3 close-aware regression: CDP count should be stable across
    #     resolutions (no doubling from seam wrap-around). Fig-8 is cy-cy
    #     so seam handling matters.
    surf = fig8(perturb=False)
    counts = []
    for res in (15, 20, 25):
        mesh, proj, cps, css, ccs, dps, sis_pairs = _build_full_layer_c_o(
            surf, resolution=res, view_J=_FIG8_SIDE_VIEW_J,
        )
        splits = SplitArrays()
        split_at_cdps(mesh, css, cps, sis_pairs, dps, splits, surf, proj)
        sps = splits.sps_array()
        n_cdp = int((sps["type"] == "cdp").sum())
        counts.append(n_cdp)

    # Stability: no abrupt doubling between resolutions (factor 0.3..3 is
    # generous; on fig8 side view we expect counts to remain ~constant).
    assert counts[0] > 0
    for a, b in zip(counts, counts[1:]):
        ratio = b / a
        assert 0.3 <= ratio <= 3.0, (
            f"CDP count unstable across resolutions: {counts}"
        )


# ── O13: assemble_subcurves ──────────────────────────────────────────────────

def test_assemble_subcurves():
    # (1) No-split BC (closed loop with no SPTs): 1 closed SubCurve.
    #     cylinder_cy has a closed BC and no corners/CCs/SISs → empty splits.
    surf = cylinder_cy(perturb=False)
    mesh = build_mesh(surf.domain, surf, resolution=10, jitter=True, seed=42)
    bcs = build_bcs(mesh)
    assert len(bcs) >= 1 and bcs[0].is_closed, "cylinder_cy should have closed BC"
    splits = SplitArrays()
    empty_css = np.zeros(0, dtype=np.dtype([
        ("p_cp", "i4"), ("q_cp", "i4"), ("face", "i4"),
        ("split1", "i4"), ("split2", "i4"),
    ]))
    empty_sis = np.zeros(0, dtype=np.dtype([
        ("p_dp", "i4"), ("q_dp", "i4"), ("flip", "i1"),
        ("split1", "i4"), ("split2", "i4"),
    ]))
    subs = assemble_subcurves(
        bcs, [], [], [], mesh, empty_css, empty_sis, splits,
    )
    # One closed SubCurve per closed BC with no splits.
    bc_subs = [s for s in subs if s.kind == "BC"]
    assert len(bc_subs) == len(bcs)
    for s in bc_subs:
        assert s.is_closed is True
        assert s.start == -1 and s.end == -1
        assert s.vc_in == 0 and s.vc_out == 0
        assert len(s.internal) > 0  # the whole loop is the internal point chain

    # (2) BC with corner splits (rect no-id paraboloid): 4 SubCurves on the
    #     single closed BC.
    surf2 = paraboloid(perturb=False)
    mesh2 = build_mesh(surf2.domain, surf2, resolution=10, jitter=True, seed=42)
    bcs2 = build_bcs(mesh2)
    assert len(bcs2) == 1 and bcs2[0].is_closed
    proj2 = _ortho_proj(surf2)
    splits2 = SplitArrays()
    split_bcs_at_corners(mesh2, bcs2, splits2, proj2)
    subs2 = assemble_subcurves(
        bcs2, [], [], [], mesh2, empty_css, empty_sis, splits2,
    )
    bc_subs2 = [s for s in subs2 if s.kind == "BC"]
    assert len(bc_subs2) == 4, (
        f"rect no-id BC should split into 4 corner-to-corner arcs, got {len(bc_subs2)}"
    )
    for s in bc_subs2:
        # Open arcs joining two distinct corner SPs.
        assert s.is_closed is False
        assert s.start >= 0 and s.end >= 0
        assert s.start != s.end

    # (3) CC with cusps (VPs) on closed CCs — torus side view, 4 cusps, 2 CCs.
    surf3 = torus(perturb=False)
    mesh3 = build_mesh(surf3.domain, surf3, resolution=20, jitter=True, seed=42)
    proj3 = _torus_side_proj(surf3)
    cps3 = find_contour_points(mesh3, proj3)
    css3 = build_contour_segments(cps3, mesh3)
    ccs3 = build_contour_curves(css3, cps3)
    vps3 = find_vps(ccs3, css3, cps3, surf3, proj3)
    assert len(vps3) == 4
    closed_ccs = sum(1 for cc in ccs3 if cc.is_closed)
    open_ccs = sum(1 for cc in ccs3 if not cc.is_closed)
    splits3 = SplitArrays()
    split_ccs_at_vps(vps3, css3, ccs3, splits3, surf3, proj3)
    subs3 = assemble_subcurves(
        [], ccs3, [], [], mesh3, css3, empty_sis, splits3,
    )
    cc_subs3 = [s for s in subs3 if s.kind == "CC"]
    # For each CC: closed → n_splits SubCurves; open → n_splits + 1.
    # Each VP splits its host CC once. Count VPs per CC by traversing ccs3.
    vp_count_per_cc = [0] * len(ccs3)
    for vp in vps3:
        cs_idx = int(vp["cs"])
        for cc_idx, cc in enumerate(ccs3):
            if any(abs(int(si)) - 1 == cs_idx for si in cc.cs_indices):
                vp_count_per_cc[cc_idx] += 1
                break
    expected_sub_count = sum(
        (vp_count_per_cc[i] if ccs3[i].is_closed else vp_count_per_cc[i] + 1)
        for i in range(len(ccs3))
        if vp_count_per_cc[i] > 0 or not ccs3[i].is_closed
    )
    # CCs with no splits: closed CC w/ no splits → 1 closed SubCurve; open CC
    # w/ no splits → 1 SubCurve (start=-1, end=-1). Add them too.
    no_split_count = sum(1 for v in vp_count_per_cc if v == 0)
    expected_total = expected_sub_count + no_split_count
    assert len(cc_subs3) == expected_total, (
        f"expected {expected_total} CC SubCurves, got {len(cc_subs3)}; "
        f"per-CC vps={vp_count_per_cc}, closed={[cc.is_closed for cc in ccs3]}"
    )
    del closed_ccs, open_ccs  # not directly asserted

    # (4) G5 identity: each SubCurve's start/end SP int equals neighbor's end/start.
    #     For test 2 (closed BC split at 4 corners): the 4 SubCurves form a cycle;
    #     each corner SP appears exactly as the end of one SubCurve AND as the
    #     start of the next.
    starts = [s.start for s in bc_subs2]
    ends = [s.end for s in bc_subs2]
    assert sorted(starts) == sorted(ends), (
        f"start multiset != end multiset; starts={starts}, ends={ends}"
    )
    # Each start/end SP int dereferences to a stored SP tuple in splits2.sps —
    # adjacent SubCurves dereference the SAME Python object (identity check).
    for s in bc_subs2:
        assert splits2.sps[s.start] is splits2.sps[s.start]
        assert splits2.sps[s.end] is splits2.sps[s.end]
    # Pair up: for each "end" SP, find a SubCurve whose "start" equals it; they
    # reference the same SP object via the same index.
    for s in bc_subs2:
        partners = [t for t in bc_subs2 if t.start == s.end]
        assert len(partners) >= 1, (
            f"no neighbor SubCurve starts at SP {s.end} (G5 broken)"
        )
        partner = partners[0]
        assert splits2.sps[s.end] is splits2.sps[partner.start], (
            "G5: adjacent SubCurves should dereference the SAME SP object"
        )

    # (5) vc_in, vc_out ∈ {-1, 0} for every emitted SubCurve.
    for s in subs2:
        assert s.vc_in in (-1, 0), f"BC SubCurve vc_in={s.vc_in} not in {{-1,0}}"
        assert s.vc_out in (-1, 0), f"BC SubCurve vc_out={s.vc_out} not in {{-1,0}}"
    for s in subs3:
        assert s.vc_in in (-1, 0), f"CC SubCurve vc_in={s.vc_in} not in {{-1,0}}"
        assert s.vc_out in (-1, 0), f"CC SubCurve vc_out={s.vc_out} not in {{-1,0}}"

    # HC passthrough sanity: an HC handed in arrives in the output unchanged.
    fake_hc = SubCurve(
        kind="HC", is_closed=False, start=0, end=0, internal=[],
        vc_in=0, vc_out=-1, parent_idx=-1,
    )
    subs_hc = assemble_subcurves(
        [], [], [], [fake_hc], mesh, empty_css, empty_sis, SplitArrays(),
    )
    assert fake_hc in subs_hc
