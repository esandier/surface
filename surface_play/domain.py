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
            # Antipodal gluing (u,v)~(-u,-v) closes the FULL unit disk into ℝP².
            # It requires r_min = 0 (an annulus would leave an inner boundary —
            # a Möbius band, not ℝP²) and r_max = 1 (the simple point-reflection
            # coincides with a surface's antipodal involution only on the unit
            # circle |z| = 1; e.g. Bryant–Kustner Boy uses z ↦ -1/z̄ = -z there).
            r_min, r_max = self.bounds[0], self.bounds[1]
            if self.type != "disk":
                raise ValueError(
                    f"antipodal identification is only valid on a disk, "
                    f"not type={self.type!r}."
                )
            if abs(r_min) > 1e-9 or abs(r_max - 1.0) > 1e-9:
                raise ValueError(
                    f"antipodal identification requires the full unit disk "
                    f"(r_min=0, r_max=1); got r_min={r_min}, r_max={r_max}. "
                    f"An annulus glued on one boundary is a Möbius band, not "
                    f"ℝP², and the gluing is only geometrically valid at |z|=1."
                )

    @property
    def is_antipodal(self) -> bool:
        return self.type in ("disk", "annulus") and self.boundary_identify == "antipodal"

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

        if self.is_antipodal:
            # Boundary glued (u, v) ~ (-u, -v): the only alias of q is its
            # point-reflection -q. Choose whichever is closer to p — for
            # interior points q wins, only seam-spanning pairs flip to -q.
            if q.ndim == 1:
                if np.sum((-q - p) ** 2) < np.sum((q - p) ** 2):
                    return -q
                return q
            qr = -q
            flip = np.sum((qr - p) ** 2, axis=-1) < np.sum((q - p) ** 2, axis=-1)
            q[flip] = qr[flip]
            return q

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
