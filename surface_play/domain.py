import math
from dataclasses import dataclass, field
from typing import Literal

import numpy as np


@dataclass
class Domain:
    type: Literal["rect", "disk", "annulus"]
    bounds: tuple  # (u_min, u_max, v_min, v_max); for disk/annulus: (r_min, r_max, 0, 2π)
    u_identify: Literal["no", "cy", "mo"] = "no"
    v_identify: Literal["no", "cy", "mo"] = "no"
    coord_type: Literal["ca", "po"] = "ca"
    # Boundary identification for disk domains: "antipodal" glues the boundary
    # circle (u, v) ~ (-u, -v), turning the disk into ℝP² (orientation-
    # reversing, like Möbius). "no" leaves the disk boundary a real edge.
    boundary_identify: Literal["no", "antipodal"] = "no"

    def __post_init__(self) -> None:
        if self.boundary_identify == "antipodal":
            # Antipodal gluing identifies the OUTER boundary circle |z| = r_max
            # by (u,v) ~ (-u,-v). It works at ANY outer radius R = r_max: the
            # involution is σ_R(z) = -R²/z̄ (= -z on |z| = R), so the surface
            # just has to satisfy S(z) = S(σ_R(z)) there (the user's
            # responsibility, like periodicity for cy/mo). The inner boundary is
            # NOT identified — r_min is free: a disk (r_min=0) glues to ℝP², an
            # annulus (r_min>0) to ℝP² minus a disk (a Möbius band with the
            # inner circle left as a real boundary).
            if self.type not in ("disk", "annulus"):
                raise ValueError(
                    f"antipodal identification is only valid on a disk/annulus, "
                    f"not type={self.type!r}."
                )

    @property
    def is_antipodal(self) -> bool:
        return self.type in ("disk", "annulus") and self.boundary_identify == "antipodal"

    @property
    def needs_close(self) -> bool:
        """True iff `close()` can actually move a point — i.e. the domain has an
        identification (rect cy/mo, or antipodal disk). For unidentified
        rect/disk/annulus `close()` is a no-op, so callers skip it."""
        if self.type == "rect":
            return self.u_identify != "no" or self.v_identify != "no"
        return self.is_antipodal

    @property
    def period_u(self) -> float:
        if self.type == "rect" and self.u_identify in ("cy", "mo"):
            return self.bounds[1] - self.bounds[0]
        return float("nan")

    @property
    def period_v(self) -> float:
        if self.type == "rect" and self.v_identify in ("cy", "mo"):
            return self.bounds[3] - self.bounds[2]
        return float("nan")

    def close(self, p: np.ndarray, q: np.ndarray) -> np.ndarray:
        """Shift q toward p by integer multiples of the period on identified axes.

        Accepts (2,) or (N, 2) arrays. Disk/annulus domains return q unchanged.
        G3: all 2D sweeps must use this before computing segment directions.

        Möbius (mo) identification: (u_min, v) ~ (u_max, -v) — when an odd
        number of seam crossings is applied to one axis, the OTHER axis must
        flip sign accordingly. mo on u flips v on each odd u-wrap; mo on v
        flips u on each odd v-wrap.
        """
        p = np.asarray(p, dtype=float)
        q = np.asarray(q, dtype=float).copy()

        # Antipodal disks are handled by `localize` / `interpolate` (via the σ_R
        # involution), NOT here — `close` is the rect cy/mo period-shift only.
        if self.type != "rect":
            return q

        scalar = q.ndim == 1
        u_min, u_max, v_min, v_max = self.bounds
        v_center = 0.5 * (v_min + v_max)
        u_center = 0.5 * (u_min + u_max)

        if self.u_identify in ("cy", "mo"):
            pu = self.period_u
            if scalar:
                diff = q[0] - p[0]
                k = int(round(diff / pu))
                q[0] -= k * pu
                if self.u_identify == "mo" and (k % 2) != 0:
                    q[1] = 2.0 * v_center - q[1]
            else:
                diff = q[:, 0] - p[:, 0]
                k = np.round(diff / pu).astype(int)
                q[:, 0] -= k * pu
                if self.u_identify == "mo":
                    flip = (k % 2) != 0
                    q[flip, 1] = 2.0 * v_center - q[flip, 1]

        if self.v_identify in ("cy", "mo"):
            pv = self.period_v
            if scalar:
                diff = q[1] - p[1]
                k = int(round(diff / pv))
                q[1] -= k * pv
                if self.v_identify == "mo" and (k % 2) != 0:
                    q[0] = 2.0 * u_center - q[0]
            else:
                diff = q[:, 1] - p[:, 1]
                k = np.round(diff / pv).astype(int)
                q[:, 1] -= k * pv
                if self.v_identify == "mo":
                    flip = (k % 2) != 0
                    q[flip, 0] = 2.0 * u_center - q[flip, 0]

        return q

    def _sigma(self, x: np.ndarray) -> np.ndarray:
        """ℝP² involution at outer radius R = r_max: σ_R(z) = -R²/z̄ = -R² x/|x|².

        `S(z) = S(σ_R(z))`, so `x` and `σ_R(x)` are the SAME surface point; σ_R
        equals `-z` on |z| = R (the gluing), swaps the disk interior (|z|<R) with
        its exterior, and σ_R∘σ_R = id. Near the origin (|x|² ≈ 0) σ_R blows up —
        return `x` there (degenerate guard; never the closer rep anyway).
        """
        r2 = float(self.bounds[1]) ** 2
        n2 = (x * x).sum(axis=1, keepdims=True)
        safe = np.where(n2 > 0.0, n2, 1.0)
        return np.where(n2 > 1e-15, -r2 * x / safe, x)

    def localize(self, p: np.ndarray, q: np.ndarray) -> np.ndarray:
        """The representative of `q` in `p`'s frame — NO map-back.

        This is what the SWEEP and the `vis_chge` direction vectors need: the
        endpoint written so the segment is local to `p`. For the antipodal disk
        that is the closer of `{q, σ(q)}`, which may legitimately lie OUTSIDE the
        unit disk (an honest inside→outside segment). For rect cy/mo it is the
        period-shift; otherwise `q` unchanged. `p`, `q` are `(2,)` or `(N, 2)`.
        """
        p = np.asarray(p, dtype=float)
        q = np.asarray(q, dtype=float)
        if not self.is_antipodal:
            return self.close(p, q) if self.needs_close else np.array(q, copy=True)
        scalar = (q.ndim == 1)
        p2 = np.atleast_2d(p)
        q2 = np.atleast_2d(q)
        if p2.shape[0] != q2.shape[0]:
            if p2.shape[0] == 1:
                p2 = np.broadcast_to(p2, q2.shape)
            elif q2.shape[0] == 1:
                q2 = np.broadcast_to(q2, p2.shape)
        qp = self._sigma(q2)
        d_q = ((q2 - p2) ** 2).sum(axis=1)
        d_qp = ((qp - p2) ** 2).sum(axis=1)
        closer = np.where((d_qp < d_q)[:, None], qp, q2)
        return closer[0] if scalar else closer

    def interpolate(self, p: np.ndarray, q: np.ndarray, s) -> np.ndarray:
        """uv at parameter `s` ∈ [0, 1] along `p → q`, ALWAYS inside the disk.

        Localizes `q` to `p`'s frame, lerps, and — for the antipodal disk — maps
        an outside result back into the disk via `σ` (so the returned point is a
        valid fundamental-domain coordinate to feed `S`). Rendering/resampling
        call this. The sweep and `vis_chge` directions use `localize` instead,
        because they need the endpoint left outside.

        `p`, `q` are `(2,)` or `(N, 2)`; `s` is a scalar or `(N,)`.
        """
        p = np.asarray(p, dtype=float)
        q = np.asarray(q, dtype=float)
        q_loc = self.localize(p, q)
        s_arr = np.asarray(s, dtype=float)
        coef = float(s_arr) if s_arr.ndim == 0 else s_arr.reshape(-1, 1)
        res = p + coef * (q_loc - p)
        if self.is_antipodal:
            R2 = float(self.bounds[1]) ** 2
            scalar = (res.ndim == 1)
            r2 = np.atleast_2d(res)
            out = (r2 * r2).sum(axis=1) > R2
            if np.any(out):
                r2 = r2.copy()
                r2[out] = self._sigma(r2[out])
            res = r2[0] if scalar else r2
        return res

    def bary_edge(self, p: np.ndarray, pq: np.ndarray, s: float) -> np.ndarray:
        """Return p + s·pq."""
        return np.asarray(p, dtype=float) + s * np.asarray(pq, dtype=float)

    def bary_face(
        self,
        p: np.ndarray,
        pq: np.ndarray,
        pr: np.ndarray,
        s: float,
        t: float,
    ) -> np.ndarray:
        """Return p + s·pq + t·pr."""
        return (
            np.asarray(p, dtype=float)
            + s * np.asarray(pq, dtype=float)
            + t * np.asarray(pr, dtype=float)
        )

    def boundary_tangent(self, uv: np.ndarray, edge_dp: np.ndarray) -> np.ndarray:
        """Unit tangent of the smooth boundary curve at `uv`, sign-matched to `edge_dp`.

        At a BCP the discretized edge chord can deviate from the analytic boundary
        tangent by O(1/res); in near-edge-on views this can flip the sign of
        `_bvis_chge`'s f3 test (whose magnitude scales with the angular gap
        between chord and tangent). Using the analytic tangent makes the result
        resolution-independent.

        - rect: each side is axis-aligned, so the chord already IS the analytic
          tangent — return the normalized `edge_dp` directly.
        - disk/annulus: the boundary is a circle, tangent at (u, v) is
          (-v, u)/sqrt(u²+v²). Sign chosen to match `edge_dp`'s direction.
        """
        uv = np.asarray(uv, dtype=float).reshape(2)
        edge_dp = np.asarray(edge_dp, dtype=float).reshape(2)

        if self.type == "rect":
            n = float(np.linalg.norm(edge_dp))
            return edge_dp / n if n > 0.0 else np.zeros(2)

        # disk / annulus
        u, v = float(uv[0]), float(uv[1])
        r = math.sqrt(u * u + v * v)
        if r == 0.0:
            n = float(np.linalg.norm(edge_dp))
            return edge_dp / n if n > 0.0 else np.zeros(2)
        Tan = np.array([-v / r, u / r])
        if float(Tan @ edge_dp) < 0.0:
            Tan = -Tan
        return Tan

    def bbox_diag(self) -> float:
        """Diagonal of the parameter-domain bounding box."""
        u_min, u_max, v_min, v_max = self.bounds
        if self.type == "rect":
            return math.sqrt((u_max - u_min) ** 2 + (v_max - v_min) ** 2)
        # disk/annulus: spatial diameter
        r_max = u_max
        return 2.0 * r_max

    def width(self) -> float:
        """rect diagonal or 2·r_max for disk/annulus."""
        return self.bbox_diag()
