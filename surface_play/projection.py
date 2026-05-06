"""Projection: ortho/perspective view of a parametric surface.

Spec: Reorganization_of_the_surface_app.md §"Projection of the surface" (lines 155-159);
kerdS formula in spec line 273.
Gotchas: G8 (ker_param via SVD).

P3 sub-step status:
- Ortho mode: full implementation (XY, Z, viewer_direction, per_vertex_viewer_dot,
  ker_param, kerdS, proj_vec).
- Perspective mode: XY, Z, viewer_direction, per_vertex_viewer_dot implemented.
  ker_param / kerdS / proj_vec raise NotImplementedError — the chain-rule
  contribution from the perspective divide is deferred to a follow-up sub-step
  inside P3 (G8 TODO).

API note: there is no public `axis` attribute. Use `viewer_direction(xyz=None)`:
in ortho the no-arg form returns the constant (3,) view direction; in persp
xyz must be supplied (raises otherwise) and the result is `xyz - eye`
(unnormalized — every consumer is sign-only on dot products with SN, which is
itself unnormalized). The internal `_axis` (= I × J normalized) is the
image-plane normal used by XY's perspective divide and by Z's depth formula.
"""

import numpy as np

from surface_play.surface import SurfaceParams


class Projection:
    """Ortho or perspective projection wired to a parametric surface."""

    def __init__(
        self,
        surface: SurfaceParams,
        I: list[float] | np.ndarray,
        J: list[float] | np.ndarray,
        O: list[float] | np.ndarray | None = None,
        eye: list[float] | np.ndarray | None = None,
    ):
        self.surface = surface
        self.I = np.asarray(I, dtype=float).reshape(3)
        self.J = np.asarray(J, dtype=float).reshape(3)

        axis = np.cross(self.I, self.J)
        n = float(np.linalg.norm(axis))
        if n == 0.0:
            raise ValueError("I × J must be non-zero")
        self._axis = axis / n  # image-plane normal, unit length

        self.mode = "ortho" if eye is None else "persp"

        if self.mode == "ortho":
            self.eye = None
            self.O = (
                np.zeros(3) if O is None
                else np.asarray(O, dtype=float).reshape(3)
            )
        else:
            self.eye = np.asarray(eye, dtype=float).reshape(3)
            if O is None:
                self.O = self.eye.copy()
            else:
                O_arr = np.asarray(O, dtype=float).reshape(3)
                if not np.allclose(O_arr, self.eye):
                    raise ValueError(
                        "In perspective mode, O must equal eye "
                        "(frontend contract from templates/play.html:330)"
                    )
                self.O = O_arr

    # ── XY / Z ────────────────────────────────────────────────────────────────
    def XY(self, xyz: np.ndarray) -> np.ndarray:
        """3D → 2D view plane. Accepts (3,) or (N, 3); returns (2,) or (N, 2)."""
        xyz = np.asarray(xyz, dtype=float)
        anchor = self.O if self.mode == "ortho" else self.eye
        d = xyz - anchor
        if self.mode == "ortho":
            if d.ndim == 1:
                return np.array([self.I @ d, self.J @ d])
            return np.stack([d @ self.I, d @ self.J], axis=-1)
        # perspective: standard perspective divide along the image-plane normal;
        # eye projects to (0, 0) by definition (zero-vector special case).
        if d.ndim == 1:
            z = float(self._axis @ d)
            if z == 0.0:
                return np.zeros(2)
            return np.array([self.I @ d, self.J @ d]) / z
        z = d @ self._axis
        out = np.zeros((d.shape[0], 2))
        nz = z != 0.0
        if np.any(nz):
            out[nz, 0] = (d[nz] @ self.I) / z[nz]
            out[nz, 1] = (d[nz] @ self.J) / z[nz]
        return out

    def Z(self, xyz: np.ndarray) -> float:
        """Depth scalar along image-plane normal. Accepts (3,) or (N, 3)."""
        xyz = np.asarray(xyz, dtype=float)
        anchor = self.O if self.mode == "ortho" else self.eye
        d = xyz - anchor
        if d.ndim == 1:
            return float(self._axis @ d)
        return d @ self._axis

    # ── Viewer direction ──────────────────────────────────────────────────────
    def viewer_direction(self, xyz: np.ndarray | None = None) -> np.ndarray:
        """Viewer direction.

        Ortho: xyz may be omitted; returns the constant (3,) unit axis.
               If xyz is given (N, 3), broadcast-returns (N, 3) for convenience.
        Persp: xyz is required; returns `xyz - eye` (unnormalized — sign-only
               consumers don't need unit length, matching SN's convention).
               Shape mirrors xyz.
        """
        if self.mode == "ortho":
            if xyz is None:
                return self._axis.copy()
            xyz = np.asarray(xyz, dtype=float)
            return np.broadcast_to(self._axis, xyz.shape).copy()
        if xyz is None:
            raise ValueError(
                "viewer_direction(xyz) requires xyz in perspective mode"
            )
        xyz = np.asarray(xyz, dtype=float)
        return xyz - self.eye

    def per_vertex_viewer_dot(self, mesh) -> np.ndarray:
        """SN · viewer_direction at every vertex, vectorized via einsum."""
        SN = np.asarray(mesh.SN, dtype=float)
        if SN.ndim != 2 or SN.shape[1] != 3:
            raise ValueError(f"mesh.SN must be (N, 3), got {SN.shape}")
        if self.mode == "ortho":
            return np.einsum("ij,j->i", SN, self._axis)
        xyz = np.asarray(mesh.xyz, dtype=float)
        d = xyz - self.eye  # unnormalized; sign-only downstream
        return np.einsum("ij,ij->i", SN, d)

    # ── Kernel direction (G8) ─────────────────────────────────────────────────
    def ker_param(self, p: np.ndarray) -> np.ndarray:
        """2D kernel direction in (u, v): right-singular vector of d(XY∘S) at p
        for the smallest singular value (G8). Ortho only — persp TODO."""
        if self.mode == "persp":
            raise NotImplementedError(
                "ker_param for perspective mode is a P3 follow-up sub-step (G8 TODO)"
            )
        p = np.asarray(p, dtype=float).reshape(2)
        u, v = float(p[0]), float(p[1])
        Su = np.asarray(self.surface.Su(u, v), dtype=float).reshape(3)
        Sv = np.asarray(self.surface.Sv(u, v), dtype=float).reshape(3)
        # O drops out of the Jacobian (constant translation).
        J_mat = np.array([
            [self.I @ Su, self.I @ Sv],
            [self.J @ Su, self.J @ Sv],
        ])
        _U, _S, Vh = np.linalg.svd(J_mat)
        return Vh[-1]

    def kerdS(self, uv: np.ndarray) -> np.ndarray:
        """3D image of ker_param via dS at uv. Ortho only — persp TODO."""
        uv = np.asarray(uv, dtype=float).reshape(2)
        kp = self.ker_param(uv)  # raises in persp mode
        u, v = float(uv[0]), float(uv[1])
        Su = np.asarray(self.surface.Su(u, v), dtype=float).reshape(3)
        Sv = np.asarray(self.surface.Sv(u, v), dtype=float).reshape(3)
        return kp[0] * Su + kp[1] * Sv

    def proj_vec(self, uv: np.ndarray, dir3d: np.ndarray) -> np.ndarray:
        """Map a 3D vector to its 2D image in the view plane.
        Ortho: linear (I·v, J·v). Perspective: TODO (chain rule through divide)."""
        if self.mode == "persp":
            raise NotImplementedError(
                "proj_vec for perspective mode is a P3 follow-up sub-step"
            )
        d = np.asarray(dir3d, dtype=float).reshape(3)
        return np.array([self.I @ d, self.J @ d])
