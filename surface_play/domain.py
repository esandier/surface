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
