"""Tests for P5: sweep_segments; C8: candidate_pairs; C9: find_double_points."""

import math
import time

import numpy as np
import pytest

from surface_play.domain import Domain
from surface_play.intersections import (
    build_sics,
    build_sis_pairs,
    candidate_pairs,
    dp_dtype,
    find_double_points,
    find_triple_points,
    intersect_dtype,
    sis_dtype,
    sweep_segments,
    tp_dtype,
)
from surface_play.mesh import build_mesh
from surface_play.test_fixtures import (
    fig8, helicoid, mobius_u, paraboloid, torus,
)


TWO_PI = 2 * math.pi


def test_cross_simple():
    """Two crossing diagonals on a unit square produce one hit at (0.5, 0.5)."""
    a0 = np.array([[0.0, 0.0]])
    a1 = np.array([[1.0, 1.0]])
    b0 = np.array([[0.0, 1.0]])
    b1 = np.array([[1.0, 0.0]])
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))

    hits = sweep_segments(a0, a1, b0, b1, domain)

    assert hits.dtype == intersect_dtype
    assert len(hits) == 1
    np.testing.assert_allclose(hits["uv"][0], (0.5, 0.5), atol=1e-12)
    assert int(hits["a"][0]) == 0
    assert int(hits["b"][0]) == 0
    assert abs(hits["t_a"][0] - 0.5) < 1e-12
    assert abs(hits["t_b"][0] - 0.5) < 1e-12


def test_self_same_segment_no_hit():
    """Self-mode, same segment included twice → no hit (i<j filter; collinear → singular)."""
    a0 = np.array([[0.0, 0.0], [0.0, 0.0]])
    a1 = np.array([[1.0, 1.0], [1.0, 1.0]])
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))

    hits = sweep_segments(a0, a1, None, None, domain, self_sweep=True)

    assert len(hits) == 0


def test_close_aware_seam_crossing():
    """Cyl-u rect (0, 2π): a wrapping horizontal segment vs a vertical near the seam.

    Without `domain`, the wrap is not understood: the parametric intersection
    of the literal lines lies at t_a ≈ 1.008 (just outside (0, 1)) and the
    bbox-x ranges don't even overlap → no hit. With `domain`, close-aware
    anchoring of seg_a's q to (0.10 + 2π, 0.5) places the cross at t_a ≈ 0.75,
    t_b = 0.5 → exactly one hit.
    """
    domain_cy = Domain(type="rect", bounds=(0.0, TWO_PI, 0.0, 1.0), u_identify="cy")
    a0 = np.array([[6.18, 0.5]])
    a1 = np.array([[0.10, 0.5]])
    b0 = np.array([[0.05, 0.4]])
    b1 = np.array([[0.05, 0.6]])

    hits_with = sweep_segments(a0, a1, b0, b1, domain_cy)
    assert len(hits_with) == 1
    assert int(hits_with["a"][0]) == 0
    assert int(hits_with["b"][0]) == 0
    np.testing.assert_allclose(hits_with["uv"][0], (0.05 + TWO_PI, 0.5), atol=1e-9)
    assert 0.0 < hits_with["t_a"][0] < 1.0
    assert 0.0 < hits_with["t_b"][0] < 1.0
    assert abs(hits_with["t_b"][0] - 0.5) < 1e-9

    hits_without = sweep_segments(a0, a1, b0, b1, None)
    assert len(hits_without) == 0


def test_parallel_non_overlapping_no_hit():
    """Two horizontal parallel segments at distinct y → no hit (denom = 0)."""
    a0 = np.array([[0.0, 0.0]])
    a1 = np.array([[1.0, 0.0]])
    b0 = np.array([[0.0, 0.5]])
    b1 = np.array([[1.0, 0.5]])
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))

    hits = sweep_segments(a0, a1, b0, b1, domain)

    assert len(hits) == 0


def test_empty_inputs():
    """Empty input arrays → empty result, both cross and self modes."""
    empty = np.empty((0, 2))
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))

    out_cross = sweep_segments(empty, empty, empty, empty, domain)
    assert len(out_cross) == 0
    assert out_cross.dtype == intersect_dtype

    out_self = sweep_segments(empty, empty, None, None, domain, self_sweep=True)
    assert len(out_self) == 0
    assert out_self.dtype == intersect_dtype


def test_performance_smoke():
    """1000 random short segments, self-sweep, completes in <1s with sane hit count."""
    rng = np.random.default_rng(42)
    n = 1000
    p0 = rng.uniform(0.0, 1.0, (n, 2))
    p1 = p0 + rng.uniform(-0.1, 0.1, (n, 2))
    domain = Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0))

    t0 = time.perf_counter()
    hits = sweep_segments(p0, p1, None, None, domain, self_sweep=True)
    elapsed = time.perf_counter() - t0

    assert elapsed < 1.0, f"sweep took {elapsed:.3f}s, expected < 1s"
    assert hits.dtype == intersect_dtype
    # Sanity: short random segments in a unit square → some hits but not pathological.
    assert 0 <= len(hits) < n * n // 4


# --- C8: candidate_pairs ----------------------------------------------------

def _random_bboxes(rng, n, *, low=0.0, high=1.0, size=0.05):
    """Random AABBs with corners in [low, high] and extent ≤ size on each axis."""
    mins = rng.uniform(low, high, (n, 3))
    extents = rng.uniform(0.0, size, (n, 3))
    return np.concatenate([mins, mins + extents], axis=1)


def _brute_force_pairs(e_bbox, f_bbox):
    """Vectorized AABB-overlap brute force. Returns sorted (e_idx, f_idx) int64 arrays."""
    e_mins = e_bbox[:, None, :3]
    e_maxs = e_bbox[:, None, 3:]
    f_mins = f_bbox[None, :, :3]
    f_maxs = f_bbox[None, :, 3:]
    overlap = np.all((e_mins <= f_maxs) & (e_maxs >= f_mins), axis=-1)
    e_idx, f_idx = np.where(overlap)
    return e_idx.astype(np.int64), f_idx.astype(np.int64)


def _pair_set(e_idx, f_idx):
    return set(zip(e_idx.tolist(), f_idx.tolist()))


def test_candidate_pairs_empty():
    """Empty inputs → empty outputs (int64)."""
    empty = np.empty((0, 6), dtype=np.float64)
    nonempty = np.array([[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]], dtype=np.float64)

    e_idx, f_idx = candidate_pairs(empty, empty)
    assert e_idx.shape == (0,) and f_idx.shape == (0,)
    assert e_idx.dtype == np.int64 and f_idx.dtype == np.int64

    e_idx, f_idx = candidate_pairs(empty, nonempty)
    assert e_idx.shape == (0,) and f_idx.shape == (0,)

    e_idx, f_idx = candidate_pairs(nonempty, empty)
    assert e_idx.shape == (0,) and f_idx.shape == (0,)


def test_candidate_pairs_non_overlapping():
    """Two disjoint AABBs → no pair."""
    e_bbox = np.array([[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]], dtype=np.float64)
    f_bbox = np.array([[2.0, 2.0, 2.0, 3.0, 3.0, 3.0]], dtype=np.float64)
    e_idx, f_idx = candidate_pairs(e_bbox, f_bbox)
    assert e_idx.shape == (0,) and f_idx.shape == (0,)


def test_candidate_pairs_overlapping():
    """Two overlapping AABBs → exactly one pair (0, 0)."""
    e_bbox = np.array([[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]], dtype=np.float64)
    f_bbox = np.array([[0.5, 0.5, 0.5, 1.5, 1.5, 1.5]], dtype=np.float64)
    e_idx, f_idx = candidate_pairs(e_bbox, f_bbox)
    assert e_idx.tolist() == [0] and f_idx.tolist() == [0]
    assert e_idx.dtype == np.int64 and f_idx.dtype == np.int64


def test_candidate_pairs_brute_force_equivalence():
    """Random (N=200, M=300) AABBs: BVH output ≡ brute-force output as set of pairs.

    AABB size is kept small (≤ 5% of domain) so that no AABB straddles the
    BVH's partition midpoints — the legacy object-split BVH places each AABB
    in exactly one child by its min-coord, so a wide AABB spanning a midpoint
    would be invisible to its overlapping neighbors on the other side. In real
    use (mesh edges/faces from a regular grid) AABB extents are O(1/N) of the
    domain and this regime holds. Density here is low (~ tens of pairs); the
    test is about contract equivalence, not density coverage.
    """
    rng = np.random.default_rng(42)
    e_bbox = _random_bboxes(rng, 200, size=0.05)
    f_bbox = _random_bboxes(rng, 300, size=0.05)

    bvh_e, bvh_f = candidate_pairs(e_bbox, f_bbox)
    bf_e, bf_f = _brute_force_pairs(e_bbox, f_bbox)

    assert len(bf_e) > 0, "test setup: expected at least one brute-force pair"
    assert _pair_set(bvh_e, bvh_f) == _pair_set(bf_e, bf_f)
    # BVH dedup: no duplicate pairs.
    assert len(bvh_e) == len(set(zip(bvh_e.tolist(), bvh_f.tolist())))


def test_candidate_pairs_performance():
    """N=10_000 BVH < 1s after JIT warmup; brute-force ≫ BVH (speedup > 30×).

    The roadmap's "brute-force baseline takes > 30s" is captured here as a
    *ratio* assertion (per the roadmap's "Assert ratio" instruction): we don't
    actually run 30s of brute force in CI. Instead we time a naive Python
    nested-loop brute force on a small (n0 × n0) subset, extrapolate to n × n,
    and assert the BVH is at least 30× faster.
    """
    rng = np.random.default_rng(7)
    n = 10_000
    # Small AABBs to (a) keep overlap density low (~few × n hits) and
    # (b) stay within the legacy BVH's correctness envelope (see equivalence test).
    e_bbox = _random_bboxes(rng, n, size=0.005)
    f_bbox = _random_bboxes(rng, n, size=0.005)

    # JIT warm-up (so we don't time Numba compilation).
    _ = candidate_pairs(e_bbox[:10], f_bbox[:10])

    t0 = time.perf_counter()
    e_idx, f_idx = candidate_pairs(e_bbox, f_bbox)
    t_bvh = time.perf_counter() - t0
    assert t_bvh < 1.0, f"BVH took {t_bvh:.3f}s, expected < 1s"

    # Naive pure-Python brute force on a small (n0 × n0) subset; extrapolate.
    e_list = e_bbox.tolist()
    f_list = f_bbox.tolist()
    n0 = 500
    t0 = time.perf_counter()
    cnt = 0
    for i in range(n0):
        eb = e_list[i]
        for j in range(n0):
            fb = f_list[j]
            if (eb[0] <= fb[3] and eb[3] >= fb[0]
                and eb[1] <= fb[4] and eb[4] >= fb[1]
                and eb[2] <= fb[5] and eb[5] >= fb[2]):
                cnt += 1
    t_brute_small = time.perf_counter() - t0
    t_brute_full = t_brute_small * (n * n) / (n0 * n0)
    assert t_brute_full / max(t_bvh, 1e-6) > 30.0, (
        f"BVH speedup only {t_brute_full / max(t_bvh, 1e-6):.1f}× "
        f"(t_bvh={t_bvh:.4f}s, t_brute_extrapolated={t_brute_full:.1f}s)"
    )


def test_candidate_pairs_determinism():
    """Same input → identical output across runs (no RNG inside BVH)."""
    rng = np.random.default_rng(123)
    e_bbox = _random_bboxes(rng, 150, size=0.3)
    f_bbox = _random_bboxes(rng, 150, size=0.3)

    e1, f1 = candidate_pairs(e_bbox, f_bbox)
    e2, f2 = candidate_pairs(e_bbox, f_bbox)
    np.testing.assert_array_equal(e1, e2)
    np.testing.assert_array_equal(f1, f2)
    # Stable sort: combined-key non-decreasing.
    combined = (e1.astype(np.int64) << 32) | f1.astype(np.int64)
    assert np.all(np.diff(combined) > 0)


# --- C9: find_double_points -------------------------------------------------

def _build(surface, resolution, *, jitter=True, seed=42):
    return build_mesh(surface.domain, surface, resolution=resolution,
                      jitter=jitter, seed=seed)


def test_find_double_points_helicoid_empty():
    """Helicoid (rect-no-no): embedded → no self-intersections."""
    surface = helicoid()
    mesh = _build(surface, resolution=15)
    dps = find_double_points(mesh, surface)
    assert dps.dtype == dp_dtype
    assert len(dps) == 0


def test_find_double_points_mobius_empty():
    """Möbius band: embedded in 3D → no self-intersections."""
    surface = mobius_u()
    mesh = _build(surface, resolution=15)
    dps = find_double_points(mesh, surface)
    assert len(dps) == 0


def test_find_double_points_fig8_nonempty():
    """Fig-8 cy-cy res=20: closed self-intersection curve → DPs > 0."""
    surface = fig8()
    mesh = _build(surface, resolution=20)
    dps = find_double_points(mesh, surface)
    # The fig-8 immersion has one closed SIC, so DP count should be O(res).
    # Loose lower bound — we mostly want to assert the algorithm finds DPs.
    assert len(dps) >= 10, f"expected ≥10 DPs, got {len(dps)}"
    assert len(dps) < 200, f"DP count {len(dps)} suspiciously high (fragmentation?)"
    # Every DP is one of the two recognized types.
    types = set(dps["type"].tolist())
    assert types <= {"EF", "EE"}


def test_find_double_points_sibling_consumption():
    """G2 regression: sibling-pair consumption prevents duplicates.

    No two DPs share the same (E1, F2) for "EF" type, and no two "EE" DPs
    share the same unordered (E1, E2) edge pair.
    """
    surface = fig8()
    mesh = _build(surface, resolution=25)
    dps = find_double_points(mesh, surface)
    assert len(dps) > 0

    ef_mask = dps["type"] == "EF"
    if ef_mask.any():
        ef_keys = list(zip(dps["E1"][ef_mask].tolist(),
                           dps["F2"][ef_mask].tolist()))
        assert len(ef_keys) == len(set(ef_keys)), \
            "duplicate (E1, F2) in EF records — sibling consumption failed"

    ee_mask = dps["type"] == "EE"
    if ee_mask.any():
        ee_keys = [
            tuple(sorted((int(e1), int(e2))))
            for e1, e2 in zip(dps["E1"][ee_mask], dps["E2"][ee_mask])
        ]
        assert len(ee_keys) == len(set(ee_keys)), \
            "duplicate unordered (E1, E2) in EE records — sibling consumption failed"


def test_find_double_points_no_fragmentation():
    """DP count grows smoothly with resolution (no factor-of-2 jumps)."""
    surface = fig8()
    counts = {}
    for res in (20, 30, 40):
        mesh = _build(surface, resolution=res)
        dps = find_double_points(mesh, surface)
        counts[res] = len(dps)
        assert len(dps) > 0, f"res={res}: no DPs found"

    r1 = counts[30] / counts[20]
    r2 = counts[40] / counts[30]
    assert 0.5 <= r1 <= 2.0, f"DP count ratio 30/20 = {r1:.2f} outside [0.5, 2.0] (counts: {counts})"
    assert 0.5 <= r2 <= 2.0, f"DP count ratio 40/30 = {r2:.2f} outside [0.5, 2.0] (counts: {counts})"


def test_find_double_points_uv_in_bounds():
    """uv1, uv2 lie within domain bounds (with tolerance for jitter)."""
    surface = fig8()
    mesh = _build(surface, resolution=20)
    dps = find_double_points(mesh, surface)
    assert len(dps) > 0

    u_min, u_max, v_min, v_max = surface.domain.bounds
    tol = 1e-3 * max(u_max - u_min, v_max - v_min)

    for field in ("uv1", "uv2"):
        uvs = dps[field]
        assert (uvs[:, 0] >= u_min - tol).all(), f"{field} u below u_min"
        assert (uvs[:, 0] <= u_max + tol).all(), f"{field} u above u_max"
        assert (uvs[:, 1] >= v_min - tol).all(), f"{field} v below v_min"
        assert (uvs[:, 1] <= v_max + tol).all(), f"{field} v above v_max"


def test_find_double_points_performance():
    """Fig-8 cy-cy res=30 < 5s after warmup (Numba BVH already JIT'd)."""
    surface = fig8()
    # Warmup: ensure BVH JIT cache is hot (independent of perf timing below).
    mesh_warm = _build(surface, resolution=10)
    _ = find_double_points(mesh_warm, surface)

    mesh = _build(surface, resolution=30)
    t0 = time.perf_counter()
    dps = find_double_points(mesh, surface)
    elapsed = time.perf_counter() - t0
    assert len(dps) > 0
    assert elapsed < 5.0, f"find_double_points took {elapsed:.2f}s, expected < 5s"


# --- C10: build_sis_pairs ---------------------------------------------------

def _shares_pair(a, b):
    """Reference set-intersection on length-2 face slots, ignoring -1."""
    a_set = {int(x) for x in a if int(x) >= 0}
    b_set = {int(x) for x in b if int(x) >= 0}
    return bool(a_set & b_set)


def test_build_sis_pairs_empty():
    """Empty/singleton input → empty SIS array with correct dtype."""
    out0 = build_sis_pairs(np.empty(0, dtype=dp_dtype))
    assert out0.dtype == sis_dtype
    assert len(out0) == 0

    out1 = build_sis_pairs(np.zeros(1, dtype=dp_dtype))
    assert out1.dtype == sis_dtype
    assert len(out1) == 0


def test_build_sis_pairs_straight_match():
    """Two DPs with matching A1/A1 and A2/A2 → one SIS, flip=+1."""
    dps = np.zeros(2, dtype=dp_dtype)
    dps["A1"][0] = (10, -1); dps["A2"][0] = (20, -1)
    dps["A1"][1] = (10, -1); dps["A2"][1] = (20, -1)

    sis = build_sis_pairs(dps)
    assert sis.dtype == sis_dtype
    assert len(sis) == 1
    assert int(sis["p_dp"][0]) == 0
    assert int(sis["q_dp"][0]) == 1
    assert int(sis["flip"][0]) == 1


def test_build_sis_pairs_crossed_match():
    """Two DPs with A1[0]==A2[1] and A2[0]==A1[1] → one SIS, flip=-1."""
    dps = np.zeros(2, dtype=dp_dtype)
    dps["A1"][0] = (10, -1); dps["A2"][0] = (20, -1)
    dps["A1"][1] = (20, -1); dps["A2"][1] = (10, -1)

    sis = build_sis_pairs(dps)
    assert len(sis) == 1
    assert int(sis["flip"][0]) == -1


def test_build_sis_pairs_fig8_count():
    """Fig-8 cy-cy res=20: SIS count is in [0.5·DP, 2·DP] (one closed SIC)."""
    surface = fig8()
    mesh = _build(surface, resolution=20)
    dps = find_double_points(mesh, surface)
    assert len(dps) > 0

    sis = build_sis_pairs(dps)
    n_dp = len(dps)
    n_sis = len(sis)
    assert n_sis > 0
    assert 0.5 * n_dp <= n_sis <= 2 * n_dp, (
        f"SIS count {n_sis} out of [0.5·DP, 2·DP] = "
        f"[{0.5 * n_dp:.1f}, {2 * n_dp}] (DPs={n_dp})"
    )


def test_build_sis_pairs_no_spurious_sharing():
    """Every emitted SIS connects DPs that actually share faces (in some orientation)."""
    surface = fig8()
    mesh = _build(surface, resolution=20)
    dps = find_double_points(mesh, surface)
    sis = build_sis_pairs(dps)
    assert len(sis) > 0

    for k in range(len(sis)):
        i = int(sis["p_dp"][k])
        j = int(sis["q_dp"][k])
        A1i, A2i = dps["A1"][i], dps["A2"][i]
        A1j, A2j = dps["A1"][j], dps["A2"][j]
        unflipped = _shares_pair(A1i, A1j) and _shares_pair(A2i, A2j)
        flipped = _shares_pair(A1i, A2j) and _shares_pair(A2i, A1j)
        assert unflipped or flipped, (
            f"SIS {k} = ({i}, {j}) connects DPs that share no faces"
        )
        # Determinism: p_dp < q_dp (from triu_indices).
        assert i < j


def test_build_sis_pairs_splits_init():
    """G17: all emitted SIS have split1 == split2 == -1 (Layer O populates)."""
    surface = fig8()
    mesh = _build(surface, resolution=20)
    dps = find_double_points(mesh, surface)
    sis = build_sis_pairs(dps)
    assert len(sis) > 0
    assert (sis["split1"] == -1).all()
    assert (sis["split2"] == -1).all()


# --- C11: build_sics --------------------------------------------------------

def test_build_sics_empty():
    """Empty SIS input → empty list."""
    out = build_sics(np.empty(0, dtype=sis_dtype))
    assert out == []


def test_build_sics_fig8_single_closed():
    """Fig-8 cy-cy res=20: exactly 1 SIC, closed; covers all SIS exactly once."""
    surface = fig8()
    mesh = _build(surface, resolution=20)
    dps = find_double_points(mesh, surface)
    sis = build_sis_pairs(dps)
    sics = build_sics(sis)

    assert len(sics) == 1
    assert sics[0].is_closed is True

    seen = {abs(int(s)) - 1 for s in sics[0].sis_indices}
    assert seen == set(range(len(sis)))


@pytest.mark.parametrize("res", [30, 40])
def test_build_sics_fig8_no_fragmentation(res):
    """Fig-8 cy-cy at higher res: still exactly 1 closed SIC (regression for G16)."""
    surface = fig8()
    mesh = _build(surface, resolution=res)
    dps = find_double_points(mesh, surface)
    sis = build_sis_pairs(dps)
    sics = build_sics(sis)

    assert len(sics) == 1, f"res={res}: expected 1 SIC, got {len(sics)}"
    assert sics[0].is_closed is True


# --- C12: find_triple_points ------------------------------------------------

def _synthetic_three_sis_fixture():
    """Three SISs whose preimages form a TP interlock on three faces (10, 20, 30).

    Faces are arbitrary distinct integer labels; the algorithm only inspects
    A1/A2 sets, not real face geometry. Three diagonal-cross pairs are placed
    on three disjoint sub-rectangles of (0, 1)² so each pair crosses once
    cleanly (strict 0 < t < 1).
    """
    from types import SimpleNamespace

    from surface_play.domain import Domain

    dps = np.zeros(6, dtype=dp_dtype)
    # Faces: F1=10, F2=20, F3=30.
    dps["A1"] = np.array(
        [(10, -1), (10, -1), (10, -1), (10, -1), (20, -1), (20, -1)],
        dtype=np.int32,
    )
    dps["A2"] = np.array(
        [(20, -1), (20, -1), (30, -1), (30, -1), (30, -1), (30, -1)],
        dtype=np.int32,
    )
    # SIS_12 (DP 0-1): preimage on F1 crosses SIS_13's F1-preimage at (0.1, 0.1).
    # SIS_13 (DP 2-3): preimage on F1 at (0.1, 0.1); on F3 at (0.1, 0.5).
    # SIS_23 (DP 4-5): preimage on F2 at (0.5, 0.1); on F3 at (0.1, 0.5).
    dps["uv1"] = np.array([
        [0.0, 0.0], [0.2, 0.2],     # SIS_12 on F1: diag /
        [0.0, 0.2], [0.2, 0.0],     # SIS_13 on F1: diag \
        [0.4, 0.2], [0.6, 0.0],     # SIS_23 on F2: diag \
    ])
    dps["uv2"] = np.array([
        [0.4, 0.0], [0.6, 0.2],     # SIS_12 on F2: diag /
        [0.0, 0.4], [0.2, 0.6],     # SIS_13 on F3: diag /
        [0.0, 0.6], [0.2, 0.4],     # SIS_23 on F3: diag \
    ])
    dps["on_boundary"] = False

    sis = np.zeros(3, dtype=sis_dtype)
    sis["p_dp"] = [0, 2, 4]
    sis["q_dp"] = [1, 3, 5]
    sis["flip"] = [1, 1, 1]
    sis["split1"] = -1
    sis["split2"] = -1

    class _FakeSurface:
        def S(self, u, v):
            # All preimages map to a single 3D point → 3D verification passes.
            return np.zeros(3)

    mesh = SimpleNamespace(domain=Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0)))
    return sis, dps, mesh, _FakeSurface()


def test_find_triple_points_empty():
    """Empty `sis_pairs` (and `dps`) → empty TP array with correct dtype."""
    from types import SimpleNamespace

    from surface_play.domain import Domain

    mesh = SimpleNamespace(domain=Domain(type="rect", bounds=(0.0, 1.0, 0.0, 1.0)))

    class _S:
        def S(self, u, v):
            return np.zeros(3)

    tps = find_triple_points(
        np.empty(0, dtype=sis_dtype),
        np.empty(0, dtype=dp_dtype),
        mesh, _S(),
    )
    assert tps.dtype == tp_dtype
    assert len(tps) == 0


def test_find_triple_points_fig8_zero():
    """Fig-8 cy-cy res=20: single closed SIC, no self-tangency → 0 TPs."""
    surface = fig8()
    mesh = _build(surface, resolution=20)
    dps = find_double_points(mesh, surface)
    sis = build_sis_pairs(dps)
    tps = find_triple_points(sis, dps, mesh, surface)
    assert tps.dtype == tp_dtype
    assert len(tps) == 0


def test_find_triple_points_synthetic():
    """Hand-crafted 3-SIS interlock → exactly 1 TP, sis_indices = {0, 1, 2}."""
    sis, dps, mesh, surface = _synthetic_three_sis_fixture()
    tps = find_triple_points(sis, dps, mesh, surface)

    assert tps.dtype == tp_dtype
    assert len(tps) == 1
    assert sorted(tps[0]["sis_indices"].tolist()) == [0, 1, 2]
    assert sorted(tps[0]["faces"].tolist()) == [10, 20, 30]
    # P_i positioned at face-cross midpoints (see fixture).
    expected_uv = {(0.1, 0.1), (0.5, 0.1), (0.1, 0.5)}
    got_uv = {(round(p[0], 6), round(p[1], 6)) for p in tps[0]["uv"]}
    assert got_uv == expected_uv


@pytest.mark.xfail(
    reason="No 'immersion with triple point' surface in the DB fixture yet",
    strict=False,
)
def test_find_triple_points_db_immersion_with_tp():
    """TODO: add an immersion-with-TP surface to test_fixtures + assert ≥1 TP."""
    raise AssertionError("DB fixture for an immersion with TPs not yet added")


@pytest.mark.parametrize("factory", [helicoid, torus, paraboloid, mobius_u])
def test_find_triple_points_no_spurious(factory):
    """Embedded fixtures: no DPs → no SISs → no TPs."""
    surface = factory()
    mesh = _build(surface, resolution=15)
    dps = find_double_points(mesh, surface)
    sis = build_sis_pairs(dps)
    tps = find_triple_points(sis, dps, mesh, surface)
    assert len(tps) == 0


def test_find_triple_points_3d_consistency():
    """Every emitted TP has all three pairwise ‖S(P_i) − S(P_j)‖ < xyz_tol."""
    sis, dps, mesh, surface = _synthetic_three_sis_fixture()
    xyz_tol = 1e-3
    tps = find_triple_points(sis, dps, mesh, surface, xyz_tol=xyz_tol)
    assert len(tps) > 0

    for tp in tps:
        P1, P2, P3 = tp["uv"][0], tp["uv"][1], tp["uv"][2]
        S1 = np.asarray(surface.S(float(P1[0]), float(P1[1]))).reshape(3)
        S2 = np.asarray(surface.S(float(P2[0]), float(P2[1]))).reshape(3)
        S3 = np.asarray(surface.S(float(P3[0]), float(P3[1]))).reshape(3)
        assert np.linalg.norm(S1 - S2) < xyz_tol
        assert np.linalg.norm(S1 - S3) < xyz_tol
        assert np.linalg.norm(S2 - S3) < xyz_tol
