"""SurfaceParams: parametric surface with perturbation and lambdified callables.

Spec: Reorganization_of_the_surface_app.md §"Parametrization of the surface".
Gotchas: G1 (perturbation phases), G11 (cse), G12 (vectorized eval).
"""

import math
from typing import Literal

import numpy as np
import sympy as sp

from surface_play.domain import Domain


# G1 perturbation constants
_PHASE_U_DEFAULT = 0.543367
_PHASE_V_DEFAULT = 0.172935
_PERT_AMPLITUDE = 0.0005
_PERT_WEIGHTS = (0.346, 0.632, 0.693)

_BBOX_SAMPLE_N = 25  # 25×25 grid for sampling bbox of unperturbed surface


class SurfaceParams:
    """Parametric surface with cse'd lambdified evaluators for S and its derivatives."""

    def __init__(
        self,
        X: str,
        Y: str,
        Z: str,
        parameter_names: str,
        domain: Domain,
        perturb: bool = True,
        output_type: Literal["ca", "cy"] = "ca",
    ):
        self.domain = domain
        self.output_type = output_type

        names = parameter_names.split()
        if len(names) != 2:
            raise ValueError(
                f"parameter_names must contain 2 names, got {parameter_names!r}"
            )
        u_name, v_name = names
        u_sym, v_sym = sp.symbols(f"{u_name} {v_name}", real=True)
        self._u, self._v = u_sym, v_sym

        local = {u_name: u_sym, v_name: v_sym}
        X_sym = sp.sympify(X, locals=local)
        Y_sym = sp.sympify(Y, locals=local)
        Z_sym = sp.sympify(Z, locals=local)

        # coord_type == "po": user wrote in polar (r, θ); substitute polar→cartesian
        # so that internal _u, _v are cartesian disk coords (x, y) and downstream
        # callers always see a unified (x, y) evaluation interface for disk/annulus.
        if domain.coord_type == "po":
            if domain.type not in ("disk", "annulus"):
                raise ValueError(
                    f"coord_type='po' requires disk/annulus domain, got {domain.type!r}"
                )
            x_sym, y_sym = sp.symbols("_x _y", real=True)
            sub = {
                u_sym: sp.sqrt(x_sym**2 + y_sym**2),
                v_sym: sp.atan2(y_sym, x_sym),
            }
            X_sym = sp.simplify(X_sym.subs(sub))
            Y_sym = sp.simplify(Y_sym.subs(sub))
            Z_sym = sp.simplify(Z_sym.subs(sub))
            self._u, self._v = x_sym, y_sym
            u_sym, v_sym = x_sym, y_sym

        if output_type == "cy":
            r_e, theta_e, z_e = X_sym, Y_sym, Z_sym
            X_sym = r_e * sp.cos(theta_e)
            Y_sym = r_e * sp.sin(theta_e)
            Z_sym = z_e

        S_pure = sp.Matrix([X_sym, Y_sym, Z_sym])

        dX, dY, dZ, bbox = self._sample_bbox(S_pure)
        self.bbox_diag = bbox

        # G1 perturbation only defined on rect domains (where identification
        # seams need invariance). Disk/annulus perturbation is a future concern.
        if perturb and domain.type == "rect":
            S_sym = self._add_perturbation(S_pure, dX, dY, dZ)
            self.perturb_flag = True
        else:
            S_sym = S_pure
            self.perturb_flag = False

        Su_sym = sp.diff(S_sym, u_sym)
        Sv_sym = sp.diff(S_sym, v_sym)
        Suu_sym = sp.diff(Su_sym, u_sym)
        Suv_sym = sp.diff(Su_sym, v_sym)
        Svv_sym = sp.diff(Sv_sym, v_sym)
        SN_sym = Su_sym.cross(Sv_sym)

        joined = (
            list(S_sym) + list(Su_sym) + list(Sv_sym)
            + list(Suu_sym) + list(Suv_sym) + list(Svv_sym)
            + list(SN_sym)
        )

        # G11 + G12: single cse'd lambdified callable returning all 21 components.
        self._eval_fn = sp.lambdify(
            (u_sym, v_sym), joined, modules="numpy", cse=True
        )

        # Naive per-function lambdified callables (no cse, separate calls) are
        # used ONLY by the cse-speedup benchmark in test_surface. Building them
        # eagerly here lambdifies 7 large expressions with no CSE — very
        # expensive for heavy formulas (e.g. the Boy surface). Defer them: store
        # the symbolic lists and lambdify on first access via __getattr__ so the
        # production path (which only uses the cse'd `_eval_fn`) never pays.
        self._naive_syms = {
            "_naive_S": list(S_sym),    "_naive_Su": list(Su_sym),
            "_naive_Sv": list(Sv_sym),  "_naive_Suu": list(Suu_sym),
            "_naive_Suv": list(Suv_sym), "_naive_Svv": list(Svv_sym),
            "_naive_SN": list(SN_sym),
        }
        self._naive_args = (u_sym, v_sym)

    def __getattr__(self, name):
        # Lazily build the test-only naive (no-cse) callables on first access.
        # __getattr__ runs only for attributes missing from __dict__, so this
        # never interferes with the eagerly-set production attributes.
        if name.startswith("_naive_"):
            syms = self.__dict__.get("_naive_syms")
            if syms is not None and name in syms:
                fn = sp.lambdify(self._naive_args, syms[name], modules="numpy")
                setattr(self, name, fn)  # cache for subsequent accesses
                return fn
        raise AttributeError(name)

    # ── Public callables ──────────────────────────────────────────────────────
    def _eval_all(self, u, v):
        u_arr = np.asarray(u, dtype=float)
        v_arr = np.asarray(v, dtype=float)
        u_arr, v_arr = np.broadcast_arrays(u_arr, v_arr)
        target_shape = u_arr.shape

        out = self._eval_fn(u_arr, v_arr)
        out = [np.broadcast_to(np.asarray(x, dtype=float), target_shape) for x in out]
        out_arr = np.stack(out)  # (21, *target_shape)

        return (
            out_arr[0:3], out_arr[3:6], out_arr[6:9],
            out_arr[9:12], out_arr[12:15], out_arr[15:18],
            out_arr[18:21],
        )

    def S(self, u, v):
        return self._eval_all(u, v)[0]

    def Su(self, u, v):
        return self._eval_all(u, v)[1]

    def Sv(self, u, v):
        return self._eval_all(u, v)[2]

    def Suu(self, u, v):
        return self._eval_all(u, v)[3]

    def Suv(self, u, v):
        return self._eval_all(u, v)[4]

    def Svv(self, u, v):
        return self._eval_all(u, v)[5]

    def SN(self, u, v):
        return self._eval_all(u, v)[6]

    # ── Internals ─────────────────────────────────────────────────────────────
    def _sample_bbox(self, S_pure):
        """Sample S over a 25×25 grid to compute (dX, dY, dZ) extents and bbox diagonal.

        For disk/annulus domains the internal eval params (_u, _v) are cartesian
        (x, y), so we sample (r, θ) over the polar bounds and convert before eval.
        """
        u_min, u_max, v_min, v_max = self.domain.bounds
        N = _BBOX_SAMPLE_N

        if self.domain.type == "rect":
            u_grid = np.linspace(u_min, u_max, N)
            v_grid = np.linspace(v_min, v_max, N)
            UU, VV = np.meshgrid(u_grid, v_grid, indexing="ij")
        else:
            rs = np.linspace(u_min, u_max, N)
            thetas = np.linspace(v_min, v_max, N)
            RR, TT = np.meshgrid(rs, thetas, indexing="ij")
            UU = RR * np.cos(TT)
            VV = RR * np.sin(TT)

        fn = sp.lambdify((self._u, self._v), list(S_pure), modules="numpy")
        out = fn(UU.ravel(), VV.ravel())
        out = [np.broadcast_to(np.asarray(x, dtype=float), (UU.size,)) for x in out]
        Xs, Ys, Zs = out
        dX = float(np.ptp(Xs))
        dY = float(np.ptp(Ys))
        dZ = float(np.ptp(Zs))
        return dX, dY, dZ, math.sqrt(dX**2 + dY**2 + dZ**2)

    def _add_perturbation(self, S_pure, dX, dY, dZ):
        """G1: phase-offset perturbation; vanishes on mo seams, repeats on cy seams."""
        u, v = self._u, self._v
        u_min, u_max, v_min, v_max = self.domain.bounds
        d_u = u_max - u_min
        d_v = v_max - v_min

        u_id = self.domain.u_identify
        v_id = self.domain.v_identify
        freq_u = 2 if u_id != "mo" else 1
        freq_v = 2 if v_id != "mo" else 1
        phi_u = sp.Float(0.0) if u_id == "mo" else sp.Float(_PHASE_U_DEFAULT)
        phi_v = sp.Float(0.0) if v_id == "mo" else sp.Float(_PHASE_V_DEFAULT)

        wave = (
            sp.sin((u - u_min) * freq_u * sp.pi / d_u + phi_u)
            * sp.sin((v - v_min) * freq_v * sp.pi / d_v + phi_v)
        )
        amp = _PERT_AMPLITUDE
        wx, wy, wz = _PERT_WEIGHTS
        delta = sp.Matrix([
            wave * amp * dX * wx,
            wave * amp * dY * wy,
            wave * amp * dZ * wz,
        ])
        return S_pure + delta
