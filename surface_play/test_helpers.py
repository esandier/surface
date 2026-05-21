"""Tests for surface_play.helpers (Layer O, step O12)."""

import numpy as np

from surface_play.contour import (
    build_contour_curves,
    build_contour_segments,
    find_contour_points,
)
from surface_play.curves import BoundaryCurve, build_bcs
from surface_play.helpers import build_helper_curves
from surface_play.intersections import dp_dtype, sis_dtype
from surface_play.mesh import build_mesh
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

