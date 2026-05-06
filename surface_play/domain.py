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
        """
        p = np.asarray(p, dtype=float)
        q = np.asarray(q, dtype=float).copy()

        if self.type != "rect":
            return q

        scalar = q.ndim == 1

        if self.u_identify in ("cy", "mo"):
            pu = self.period_u
            if scalar:
                diff = q[0] - p[0]
                q[0] -= round(diff / pu) * pu
            else:
                diff = q[:, 0] - p[:, 0]
                q[:, 0] -= np.round(diff / pu) * pu

        if self.v_identify in ("cy", "mo"):
            pv = self.period_v
            if scalar:
                diff = q[1] - p[1]
                q[1] -= round(diff / pv) * pv
            else:
                diff = q[:, 1] - p[:, 1]
                q[:, 1] -= np.round(diff / pv) * pv

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
