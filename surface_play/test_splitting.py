"""Tests for surface_play.splitting (Layer O, step O5)."""

import numpy as np
import pytest

from surface_play.curves import build_bcs
from surface_play.mesh import build_mesh, edge_dtype
from surface_play.projection import Projection
from surface_play.splitting import (
    SplitArrays,
    SplitSlotOverflowError,
    sp_dtype,
    split_bcs_at_corners,
    spt_dtype,
)
from surface_play.test_fixtures import (
    cylinder_cy,
    disk_paraboloid_po,
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
