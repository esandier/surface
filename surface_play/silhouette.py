# Todo:
# Corriger la connection qui traverse une ligne de contour, probablement dû au raccord cylindre qui fout la merde. visible sur le 
# vase bulle sous certains angles. 

# Corriger l'exclusion d'intersections (en prenant comme règle que les segments ne doivent pas appartenir à une même face ?)
# Corriger le bug dans simplify: cf message ci-dessous, plot de (sin(u)+sin(v)^2 
# Ajouter la perspective (axis dépend du point) ?
# Corriger pour avoir un résultat correct en cas d'échelle non-isotrope ?
# automatiser la résolution. Soit en fonction de la courbure max, soit en fonction des incohérences de visibilité.

# ajouter les auto intersections
# Ajouter des hachures décoratives
# ajouter un bord donné par une équation
# optimiser breaklines: d'abord, dans la triangulation, couper les lignes de bord aux coins, puis vectoriser.. Pour dans très longtemps, prendre les lignes domaine, diviser pour garantir que chacune ne s'autointersecte pas, regarder
# les couples qui s'intersectent. Diviser en deux, boucler. attention autour des cusps et bords-contours, il faut prendre les courbes dans un
# voisinage autour et étudier les auto intersections. Autre possibilité: utiliser shapely. Optimiser connect.

import sympy as sp
import numpy as np
import math
import triangle as tr

from numpy.linalg import norm
from numpy.core.records import fromrecords
from numpy.core.records import fromarrays
from numpy.linalg import svd
from scipy.optimize import newton, linprog
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree

from sympy.utilities import lambdify
from collections import namedtuple
from collections import deque
from collections import defaultdict
import bisect

import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import time

colormap = [
    "navy",
    "indianred",
    "slateblue",
    "firebrick",
    "mediumpurple",
    "blue",
    "red",
    "green",
    "black",
    "gold",
    "orange",
    "navajowhite",
    "aquamarine",
    "greenyellow",
    "lime",
]
linestyles = [(0, (3, 1)), (0, (3, 2)), (0, (2, 2)), (0, (1, 2)), (0, (1, 3))]

# faces, les indices des 3 points dans la grid
f_type = np.dtype([("p", "uint", 2), ("q", "uint", 2), ("r", "uint", 2)])
# arêtes, les indices des 2 points, des faces (g=-1 s'il n'y a qu'une face), coord dans le domaine du point, du vecteur, et de la direction intérieure
e_type = np.dtype(
    [
        ("p", "uint"),
        ("q", "uint"),
        ("f", "i4"),
        ("g", "i4"),
        ("fp", "f8", 2),
        ("fdp", "f8", 2),
        ("dir", "f8", 2),
        ("flip", "i2"),
    ]
)
# point de contour, indice de son arête, coord bary sur l'arête, normale au contour dans le plan de vision dirigée vers S, coord dans le domaine
c_type = np.dtype(
    [("e", "uint"), ("s", "f8"), ("d", "f8", 2), ("fp", "f8", 2)]
)  # point de contour
# cusp
cusp_type = np.dtype(
    [("p", "int"), ("q", "int"), ("s", "f8"), ("fp", "f8", 2)]
)  # point de cusp
# point de ligne type a 2 car: 1er pour le point ('c' pour contour 'v' pour cusp), le 2ème pour la ligne ('b' pour bord, 'c' pour contour),
# puis arête qui contient le point, point domaine, et direction surface, puis changement de visibilité.
# v est le changement de visibilité en entrant  dans la ligne depuis le point. toujours négatif ou nul
# fp change avec le relèvement d'une ligne pour éviter les saut dans le cas de quotient
# ixfp est dans le carré de coordonnées et ne change pas, permet de tester le contact de 2 lignes.
line_pt_type = np.dtype(
    [
        ("type", "<U2"),
        ("e", "int"),
        ("v", "int"),
        ("d", "f8", 2),
        ("fp", "f8", 2),
        ("ixfp", "f8", 2),
    ]
)
# stockage des arêtes pour calculer leurs intersections
bag_type = np.dtype(
    [
        ("x", "f8"),
        ("type", "<U1"),
        ("key", "int"),
        ("n", int),
        ("vec", "f8", 2),
        ("p", "f8", 2),
        ("z", "f8"),
        ("dp", "f8", 2),
        ("dz", "f8"),
        ("s0", "f8"),
        ("s1", "f8"),
    ]
)
bagtype_bis = np.dtype(
    [
        ("x", "f8"),
        ("l_idx", "int"),
        ("e_idx", "int"),
        ("e", "int"),
        ("n", "int"),
        ("d", "f8", 2),
        ("p", "f8", 2),
        ("z", "f8"),
        ("dp", "f8", 2),
        ("dz", "f8"),
        ("fp", "f8", 2),
        ("fdp", "f8", 2),
    ]
)

def make_lines(
    ps, qs, blocks=set()
):  # ps et qs sont les listes des extrémités des arêtes, ou des faces qui touchent,
    # selon que c'est le bord ou les contours. renvoie une liste de np.array.
    # pour les courbes fermées, le dernier élément de la liste est égal au premier.
    t0 = time.perf_counter()

    if isinstance(
        ps[0], np.ndarray
    ):  # si ce sont des points de contours, on transforme en tuples pour en faire des clés de dictionnaire
        pts, qts = [tuple(p) for p in ps], [tuple(q) for q in qs]
    else:
        pts, qts = ps, qs

    def pt(i):
        res = (
            pts[i] if i >= 0 else qts[i]
        )  # indice négatif = arête renversée dans la ligne
        return res

    def qt(i):
        res = qts[i] if i >= 0 else pts[i]
        return res

    def rev(i):
        j = i - len(pts) if i >= 0 else i + len(pts)
        return j

    dic = defaultdict(list)
    # clé =  point, valeur =  arêtes auxquelles il appartient, index >0 ou <0, selon que c'est le 'p' de l'arête ou le 'q'
    for i, (p, q) in enumerate(zip(pts, qts)):
        dic[p].append(i)
        if (
            q != -1
        ):  # q = -1 quand c'est une face absente (pour une point de contour au bord)
            dic[q].append(i - len(pts))

    lines = []

    while dic:  # while dic is not empty
        p, l = dic.popitem()  # un point, et le s arêtes auxquelles il appartient
        line = deque([l[0]])
        if len(l) > 1:
            line.appendleft(rev(l[1]))  # p est le 'q' de line[0]

        while pt(line[0]) in dic:
            l = dic.pop(pt(line[0]))
            if qt(l[0]) != qt(line[0]):
                line.appendleft(rev(l[0]))
            elif len(l) > 1:
                line.appendleft(rev(l[1]))
        while qt(line[-1]) in dic:
            l = dic.pop(qt(line[-1]))
            if qt(l[0]) != pt(line[-1]):
                line.append(l[0])
            elif len(l) > 1:
                line.append(l[1])
        lines.append(np.array(line))

    # print('[%0.3fs] %s' % (time.perf_counter()-t0, 'make lines'))
    return lines

def vect_prod(p, q):
    return np.array(
        [p[i % 3] * q[(i + 1) % 3] - p[(i + 1) % 3] * q[i % 3] for i in range(1, 4)]
    )

def _eval_vec(fn, u, v):
    """Call a lambdified vector function and return a (3, N) float64 array.
    Handles the case where one or more components are Python scalars (constant
    expressions) by broadcasting them to the size of u."""
    raw = fn(u, v)
    n = u.size
    return np.row_stack([
        np.full(n, c, dtype=np.float64) if np.ndim(c) == 0 else np.asarray(c, dtype=np.float64)
        for c in raw
    ])

def _center_radius(S_vals):
    """Bounding-sphere center and radius from a (3, N) array of surface points."""
    xmin, xmax = float(S_vals[0].min()), float(S_vals[0].max())
    ymin, ymax = float(S_vals[1].min()), float(S_vals[1].max())
    zmin, zmax = float(S_vals[2].min()), float(S_vals[2].max())
    center = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2]
    radius = np.sqrt((xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2) / 2
    return center, float(radius)


class Surface:
    def print(self, s):
        # pass
        print(s)

    def __init__(self, X="x", Y="y", Z="x*y", param_names="x y", bounds=(0,1,0,1),
                 quotient=("no", "no"), domain_type="rect", coord_type="ca", r_min=0.0, r_max=1.0,
                 output_type="ca"):
        t0 = time.perf_counter()
        self.low_res = 80
        self.domain_type = domain_type
        self.r_min = r_min
        self.r_max = r_max

        self.u_min, self.u_max, self.v_min, self.v_max = (
            bounds[0], bounds[1], bounds[2], bounds[3],
        )
        self.u_identify = quotient[0]
        self.v_identify = quotient[1]

        if domain_type == 'disk':
            # For disk domains, override bounds to bounding square and clear identifications
            self.u_min, self.u_max = -r_max, r_max
            self.v_min, self.v_max = -r_max, r_max
            self.u_identify = 'no'
            self.v_identify = 'no'

        # Polar → cartesian substitution for disk domain
        if domain_type == 'disk' and coord_type == 'po':
            var_syms = sp.symbols(param_names)   # (r_sym, t_sym)
            x_s, y_s = sp.symbols('_x_ _y_')
            subs = [(var_syms[0], sp.sqrt(x_s**2 + y_s**2)),
                    (var_syms[1], sp.atan2(y_s, x_s))]
            X = str(sp.sympify(X).subs(subs))
            Y = str(sp.sympify(Y).subs(subs))
            Z = str(sp.sympify(Z).subs(subs))
            param_names = '_x_ _y_'

        # Cylindrical output → cartesian: X=R, Y=Θ, Z stays
        if output_type == 'cy':
            R_expr = sp.sympify(X)
            T_expr = sp.sympify(Y)
            Z_expr = sp.sympify(Z)
            X = str(R_expr * sp.cos(T_expr))
            Y = str(R_expr * sp.sin(T_expr))
            Z = str(Z_expr)

        u, v = sp.symbols(param_names)

        # S_pur: unperturbed surface expression (truc_nul trick forces SymPy to keep array structure)
        with sp.evaluate(False) :
            truc_nul = u*sp.sin(sp.pi)
        S_pur = sp.Array(sp.sympify([X, Y, Z])) + sp.Array([truc_nul,truc_nul,truc_nul])

        # Coarse sampling to estimate bounding-box extents for the perturbation amplitude.
        # (for_3js and triangulate do their own full sampling later.)
        _Slambda = lambdify([u, v], S_pur, "numpy")
        _bb_res = 20
        if domain_type == 'disk':
            _r = np.linspace(max(r_min, r_max / 100), r_max, _bb_res)
            _t = np.linspace(0, 2 * np.pi, _bb_res)
            _rg, _tg = np.meshgrid(_r, _t, indexing='ij')
            _ug, _vg = _rg * np.cos(_tg), _rg * np.sin(_tg)
        else:
            _ug, _vg = np.meshgrid(
                np.linspace(self.u_min, self.u_max, _bb_res),
                np.linspace(self.v_min, self.v_max, _bb_res),
                indexing='ij'
            )
        _Sg = _eval_vec(_Slambda, _ug.flatten(), _vg.flatten())
        dX = float(_Sg[0].max() - _Sg[0].min()) or 1.0
        dY = float(_Sg[1].max() - _Sg[1].min()) or 1.0
        dZ = float(_Sg[2].max() - _Sg[2].min()) or 1.0

        # On perturbe S, et on fabrique les fonctions lambdifiées
        u, v = sp.symbols(param_names)
        d_u, d_v = self.u_max - self.u_min, self.v_max - self.v_min
        ax_x, ax_y, ax_z = sp.symbols("ax_x ax_y ax_z")
        ax = sp.Array([ax_x, ax_y, ax_z])
        du, dv = sp.symbols("du dv")

        # Le choix de la perturbation est délicat. En gros il faut une dérivée suffisament grande
        # (de l'ordre de dérivée de X divisé par resolution), amplitude petite pour que ça se voie pas, oscillation pas trop
        # grande pour pas perturber les cusps.
        freq_u = 2 if self.u_identify != "mo" else 1
        freq_v = 2 if self.v_identify != "mo" else 1

        S_sp = S_pur + (sp.sin((u - self.u_min) * freq_u * sp.pi / d_u) * sp.sin((v - self.v_min) * freq_v * sp.pi / d_v)) * sp.Array(
            [
                0.0005 * dX * .346,
                0.0005 * dY * .632,
                0.0005 * dZ * .693,
            ]
        )

        # dérivée de S.N dans la direction (du, dv)
        # DS_axis_sp = sp.diff(S_sp,u).dot([ax_x,ax_y,ax_z]) * du + sp.diff(S_sp,v).dot([ax_x,ax_y,ax_z]) * dv
        Su_sp = sp.diff(S_sp, u)  # dérivée par rapport à u, comme expression sympy
        Sv_sp = sp.diff(S_sp, v)
        Suu_sp = sp.diff(Su_sp, u)  # dérivée par rapport à u, comme expression sympy
        Suv_sp = sp.diff(Su_sp, v)
        Svv_sp = sp.diff(Sv_sp, v)  # dérivée par rapport à u, comme expression sympy

        def cross(X, Y):
            return sp.Array(
                [
                    X[1] * Y[2] - X[2] * Y[1],
                    X[2] * Y[0] - X[0] * Y[2],
                    X[0] * Y[1] - X[1] * Y[0],
                ]
            )

        def dot(X, Y):
            return X[0] * Y[0] + X[1] * Y[1] + X[2] * Y[2]

        SN_sp = cross(
            Su_sp, Sv_sp
        )  # vecteur normal, non normé, comme expression sympy,
        # Orthographic silhouette condition: N · axis = 0
        SN_axis_sp = dot(SN_sp, ax)
        SD_jac_sp = sp.Array(
            [
                dot(
                    cross(Suu_sp, Sv_sp) + cross(Su_sp, Suv_sp), ax
                ),  # grad du jacobien, orth aux contours
                dot(cross(Suv_sp, Sv_sp) + cross(Su_sp, Svv_sp), ax),
            ]
        )

        # Perspective silhouette condition: N · (S - E) = 0  (ax = eye position E)
        # Note: N · Su = N · Sv = 0 always, so the gradient simplifies to the same structure
        SN_axis_persp_sp = dot(SN_sp, S_sp - ax)
        SD_jac_persp_sp = sp.Array(
            [
                dot(cross(Suu_sp, Sv_sp) + cross(Su_sp, Suv_sp), S_sp - ax),
                dot(cross(Suv_sp, Sv_sp) + cross(Su_sp, Svv_sp), S_sp - ax),
            ]
        )

        self.S = lambdify([u, v], S_sp, "numpy")# applied to a numpy array, creates a list of three np.arrays of the same dimension
        self.Su = lambdify([u, v], Su_sp, "numpy")
        self.Sv = lambdify([u, v], Sv_sp, "numpy")
        self.Suu = lambdify([u, v], Suu_sp, "numpy")
        self.Suv = lambdify([u, v], Suv_sp, "numpy")
        self.Svv = lambdify([u, v], Svv_sp, "numpy")
        self.SN = lambdify([u, v], SN_sp, "numpy")
        self.SN_axis_ortho = lambdify([u, v, ax_x, ax_y, ax_z], SN_axis_sp, "numpy")
        self.SD_jac_ortho  = lambdify([u, v, ax_x, ax_y, ax_z], SD_jac_sp, "numpy")
        self.SN_axis_persp = lambdify([u, v, ax_x, ax_y, ax_z], SN_axis_persp_sp, "numpy")
        self.SD_jac_persp  = lambdify([u, v, ax_x, ax_y, ax_z], SD_jac_persp_sp, "numpy")
        self.SN_axis = self.SN_axis_ortho  # default; overridden by set_axis
        self.SD_jac  = self.SD_jac_ortho

        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "init surface"))

    def set_axis(
        self, I=[1,0,0], J=[0,0,1], eye=None
    ):  # definit l'axe et les fonctions self.XY et self.Z
        vecI, vecJ = np.array(I), np.array(J)
        vecI, vecJ = vecI/norm(vecI), vecJ/norm(vecJ)
        self.axis = np.cross(vecI, vecJ)
        self.axis = self.axis/norm(self.axis) # coordonnées x,y,z de la normale au plan de vision

        if eye is None:
            # Orthographic
            self.SN_axis  = self.SN_axis_ortho
            self.SD_jac   = self.SD_jac_ortho
            self.ax_param = self.axis  # passed as (ax_x, ax_y, ax_z) to lambdify calls
            def XY(vec):  # renvoie coordonnées d'un point de l'espace dans le repère I,J
                return np.array([np.inner(vecI, vec.T), np.inner(vecJ, vec.T)])
            def proj_vec(u, v, T):  # project a tangent vector T at (u,v) — linear for ortho
                return np.array([np.inner(vecI, T.T), np.inner(vecJ, T.T)])
            def ker_param(u, v):  # kernel of orthographic projection Jacobian in parameter space
                du = vect_prod(self.Su(u, v), self.axis)
                dv = vect_prod(self.Sv(u, v), self.axis)
                A = np.stack((du.T, dv.T), axis=-1)
                _, _, vh = svd(A)
                return vh[..., -1, :].T
        else:
            # Perspective: silhouette is N·(S−E)=0, XY uses perspective division
            eye = np.array(eye)
            self.SN_axis  = self.SN_axis_persp
            self.SD_jac   = self.SD_jac_persp
            self.ax_param = eye  # passed as (ax_x, ax_y, ax_z) = eye position
            def XY(vec):
                vec = np.asarray(vec)
                d = vec - (eye[:, np.newaxis] if vec.ndim == 2 else eye)
                depth = -np.inner(self.axis, d.T)  # axis points from surface toward eye, so d·axis < 0; negate for positive depth
                return np.array([np.inner(vecI, d.T) / depth, np.inner(vecJ, d.T) / depth])
            def proj_vec(u, v, T):  # project a tangent vector T at (u,v) using the perspective Jacobian
                # d(XY)/dt = (inner(vecI/J, T) + xy * inner(axis, T)) / depth
                T = np.asarray(T)
                S_pt = np.asarray(self.S(u, v))
                d = S_pt - (eye[:, np.newaxis] if S_pt.ndim == 2 else eye)
                depth = -np.inner(self.axis, d.T)
                xy = np.array([np.inner(vecI, d.T) / depth, np.inner(vecJ, d.T) / depth])
                axis_T = np.inner(self.axis, T.T)
                return np.array([(np.inner(vecI, T.T) + xy[0] * axis_T) / depth,
                                 (np.inner(vecJ, T.T) + xy[1] * axis_T) / depth])
            def ker_param(u, v):  # kernel of perspective projection Jacobian in parameter space (2×2 SVD)
                Su_a = np.array(self.Su(u, v)); Sv_a = np.array(self.Sv(u, v))
                S_pt = np.asarray(self.S(u, v))
                scalar = S_pt.ndim == 1
                if scalar:
                    Su_a = Su_a[:, np.newaxis]; Sv_a = Sv_a[:, np.newaxis]; S_pt = S_pt[:, np.newaxis]
                d = S_pt - eye[:, np.newaxis]
                depth = -np.inner(self.axis, d.T)
                xy = np.array([np.inner(vecI, d.T) / depth, np.inner(vecJ, d.T) / depth])
                aSu = np.inner(self.axis, Su_a.T); aSv = np.inner(self.axis, Sv_a.T)
                # 2×2 perspective Jacobian: rows are vecI, vecJ; cols are u, v directions
                J = np.array([[np.inner(vecI, Su_a.T) + xy[0] * aSu, np.inner(vecI, Sv_a.T) + xy[0] * aSv],
                              [np.inner(vecJ, Su_a.T) + xy[1] * aSu, np.inner(vecJ, Sv_a.T) + xy[1] * aSv]])
                # J shape (2, 2, N) → (N, 2, 2) for batched SVD
                _, _, vh = svd(J.transpose(2, 0, 1))  # transpose axes to (N, 2, 2) without transposing the matrix itself
                ker = vh[:, -1, :].T  # (2, N)
                return ker[:, 0] if scalar else ker

        def Z(vec):  # depth along view axis — relative ordering valid for both ortho and persp
            return np.inner(self.axis, vec.T)

        def XYZ(vec):  # renvoie le vecteur dont les coord. sur (I,J) sont 'vec'
            return vec[0] * vecI + vec[1] * vecJ

        self.XY = XY
        self.proj_vec = proj_vec
        self.ker_param = ker_param
        self.Z = Z
        self.XYZ = XYZ

    def for_3js(self) :
        if self.domain_type == 'disk':
            return self._for_3js_disk()
        res = self.low_res
        u_grid, v_grid = np.meshgrid(
            np.linspace(self.u_min, self.u_max, res + 1),
            np.linspace(self.v_min, self.v_max, res + 1),
            indexing="ij"
        )
        u_flat, v_flat = u_grid.flatten(), v_grid.flatten()
        S_vals  = _eval_vec(self.S,  u_flat, v_flat)   # (3, N)
        SN_vals = _eval_vec(self.SN, u_flat, v_flat)   # (3, N)
        positions = S_vals.T.flatten()
        normals   = SN_vals.T.flatten()
        center, radius = _center_radius(S_vals)
        index = np.reshape(np.arange((res + 1) ** 2), (res + 1, res + 1))
        a = index[:-1, :-1]; b = index[:-1, 1:]
        c = index[1:, :-1];  d = index[1:, 1:]
        faces = np.concatenate((np.transpose(np.array([a,c,b])).flatten(),
                                np.transpose(np.array([c,d,b])).flatten()))
        return positions.tolist(), normals.tolist(), faces.tolist(), center, radius

    def _for_3js_disk(self):
        res = self.low_res
        n_bnd = max(64, res)
        theta_bnd = np.linspace(0, 2 * np.pi, n_bnd, endpoint=False)
        outer = np.column_stack([self.r_max * np.cos(theta_bnd),
                                  self.r_max * np.sin(theta_bnd)])
        outer_segs = np.column_stack([np.arange(n_bnd), (np.arange(n_bnd) + 1) % n_bnd])
        if self.r_min > 0:
            inner = np.column_stack([self.r_min * np.cos(theta_bnd),
                                      self.r_min * np.sin(theta_bnd)])
            inner_segs = np.column_stack([n_bnd + np.arange(n_bnd),
                                           n_bnd + (np.arange(n_bnd) + 1) % n_bnd])
            verts = np.vstack([outer, inner])
            segs  = np.vstack([outer_segs, inner_segs])
            mesh_in = {'vertices': verts, 'segments': segs, 'holes': np.array([[0.0, 0.0]])}
        else:
            mesh_in = {'vertices': outer, 'segments': outer_segs}
        disk_area = np.pi * (self.r_max**2 - self.r_min**2)
        target_area = disk_area / (res * res)
        mesh = tr.triangulate(mesh_in, f'pYq30a{target_area:.15f}')
        u_flat = mesh['vertices'][:, 0]
        v_flat = mesh['vertices'][:, 1]
        tris   = mesh['triangles'].astype(np.int32)
        S_vals  = _eval_vec(self.S,  u_flat, v_flat)   # (3, N)
        SN_vals = _eval_vec(self.SN, u_flat, v_flat)   # (3, N)
        positions = S_vals.T.flatten()
        normals   = SN_vals.T.flatten()
        faces     = tris.flatten()
        center, radius = _center_radius(S_vals)
        return positions.tolist(), normals.tolist(), faces.tolist(), center, radius

    def triangulate(self, res):
        t0 = time.perf_counter()
        self.res = res

        if self.domain_type == 'disk':
            self._triangulate_disk(res, t0)
        else:
            self._triangulate_rect(res, t0)

    def _triangulate_disk(self, res, t0):
        """Triangulate a disk/annulus domain using the triangle library."""
        n_bnd = max(64, res)
        theta_bnd = np.linspace(0, 2 * np.pi, n_bnd, endpoint=False)

        outer = np.column_stack([self.r_max * np.cos(theta_bnd),
                                  self.r_max * np.sin(theta_bnd)])
        outer_segs = np.column_stack([np.arange(n_bnd), (np.arange(n_bnd) + 1) % n_bnd])

        if self.r_min > 0:
            inner = np.column_stack([self.r_min * np.cos(theta_bnd),
                                      self.r_min * np.sin(theta_bnd)])
            inner_segs = np.column_stack([n_bnd + np.arange(n_bnd),
                                           n_bnd + (np.arange(n_bnd) + 1) % n_bnd])
            verts = np.vstack([outer, inner])
            segs  = np.vstack([outer_segs, inner_segs])
            mesh_in = {'vertices': verts, 'segments': segs, 'holes': np.array([[0.0, 0.0]])}
        else:
            mesh_in = {'vertices': outer, 'segments': outer_segs}

        disk_area = np.pi * (self.r_max**2 - self.r_min**2)
        target_area = disk_area / (res * res)
        mesh = tr.triangulate(mesh_in, f'pYq30a{target_area:.15f}')

        u_flat = mesh['vertices'][:, 0]
        v_flat = mesh['vertices'][:, 1]
        tris   = mesh['triangles'].astype(np.int32)   # (M, 3)
        N_tris = len(tris)

        self.S_grid  = _eval_vec(self.S,  u_flat, v_flat)
        self.SN_grid = _eval_vec(self.SN, u_flat, v_flat)
        self.center, self.radius = _center_radius(self.S_grid)
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "calcul grid (disk)"))
        t0 = time.perf_counter()

        mean_edge = np.sqrt(disk_area / N_tris)
        self.uv_vector = np.array([5.0 / mean_edge, 5.0 / mean_edge])

        # Build edge → face map
        edge_to_faces = defaultdict(list)
        for f_idx, (va, vb, vc) in enumerate(tris):
            for p_v, q_v in [(va, vb), (vb, vc), (vc, va)]:
                edge_to_faces[(min(p_v, q_v), max(p_v, q_v))].append((f_idx, p_v, q_v))

        int_ps, int_qs, int_fs, int_gs = [], [], [], []
        int_fps, int_fdps, int_dirs, int_flips = [], [], [], []
        bnd_ps, bnd_qs, bnd_fs, bnd_gs = [], [], [], []
        bnd_fps, bnd_fdps, bnd_dirs, bnd_flips = [], [], [], []

        for (mn, mx), face_list in edge_to_faces.items():
            fp_c  = np.array([u_flat[mn], v_flat[mn]])
            fdp_c = np.array([u_flat[mx] - u_flat[mn], v_flat[mx] - v_flat[mn]])
            if len(face_list) == 2:
                int_ps.append(mn); int_qs.append(mx)
                int_fs.append(face_list[0][0]); int_gs.append(face_list[1][0])
                int_fps.append(fp_c); int_fdps.append(fdp_c)
                int_dirs.append(np.array([1.0, 0.0])); int_flips.append(1)
            else:
                f_idx = face_list[0][0]
                edge_len = norm(fdp_c)
                edge_dir = fdp_c / (edge_len + 1e-30)
                inward = np.array([-edge_dir[1], edge_dir[0]])
                centroid = np.array([u_flat[tris[f_idx]].mean(), v_flat[tris[f_idx]].mean()])
                if np.dot(inward, centroid - fp_c) < 0:
                    inward = -inward
                bnd_ps.append(mn); bnd_qs.append(mx)
                bnd_fs.append(f_idx); bnd_gs.append(-1)
                bnd_fps.append(fp_c); bnd_fdps.append(fdp_c)
                bnd_dirs.append(inward); bnd_flips.append(1)

        def _make_eds(ps, qs, fs, gs, fps, fdps, dirs, flips):
            if not ps:
                return np.array([], dtype=e_type)
            return fromarrays(
                (np.array(ps, dtype=np.uint32),
                 np.array(qs, dtype=np.uint32),
                 np.array(fs, dtype=np.int32),
                 np.array(gs, dtype=np.int32),
                 np.array(fps, dtype=np.float64),
                 np.array(fdps, dtype=np.float64),
                 np.array(dirs, dtype=np.float64),
                 np.array(flips, dtype=np.int16)),
                dtype=e_type,
            )

        i_eds = _make_eds(int_ps, int_qs, int_fs, int_gs, int_fps, int_fdps, int_dirs, int_flips)
        b_eds = _make_eds(bnd_ps, bnd_qs, bnd_fs, bnd_gs, bnd_fps, bnd_fdps, bnd_dirs, bnd_flips)
        self.b_index = len(i_eds)
        self.eds  = np.concatenate([i_eds, b_eds]).view(e_type)
        self.ieds = self.eds[:self.b_index]
        self.beds = self.eds[self.b_index:]
        self.corners = set()
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "triangulation (disk)"))

        t0 = time.perf_counter()
        self.b_lines = make_lines(self.beds["p"], self.beds["q"], self.corners)
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "make boundary lines"))

        def simple_close(pt, qt):
            return qt
        self.close = simple_close
        self.generator_u = lambda pt, i: pt
        self.generator_v = lambda pt, i: pt

    def _triangulate_rect(self, res, t0):
        """Triangulate a rectangular domain using a structured grid (flat vertex indices)."""
        u_list = np.linspace(self.u_min, self.u_max, res + 1)
        v_list = np.linspace(self.v_min, self.v_max, res + 1)

        self.uv_vector = 5 * np.array(
            [res / (self.u_max - self.u_min), res / (self.v_max - self.v_min)]
        )

        u_grid, v_grid = np.meshgrid(u_list, v_list, indexing="ij")
        u_flat, v_flat = u_grid.flatten(), v_grid.flatten()
        self.S_grid  = _eval_vec(self.S,  u_flat, v_flat)
        self.SN_grid = _eval_vec(self.SN, u_flat, v_flat)
        self.center, self.radius = _center_radius(self.S_grid)
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "calcul grid"))
        t0 = time.perf_counter()

        # Flat vertex index: vertex (i,j) → i*(res+1)+j
        points = np.arange((res + 1) * (res + 1)).reshape(res + 1, res + 1)
        fpoints = np.array(np.meshgrid(u_list, v_list)).T  # (res+1,res+1,2), u=row, v=col

        if self.u_identify == "cy":
            points[res, :] = points[0, :]
        if self.v_identify == "cy":
            points[:, res] = points[:, 0]
        if self.u_identify == "mo":
            points[res, :] = points[0, ::-1]
        if self.v_identify == "mo":
            points[:, res] = points[::-1, 0]

        a = points[:res, :res]
        b = points[:res, 1:res + 1]
        c = points[1:res + 1, :res]
        d = points[1:res + 1, 1:res + 1]

        fa = fpoints[:res, :res]
        fb = fpoints[:res, 1:res + 1]
        fc = fpoints[1:res + 1, :res]
        fd = fpoints[1:res + 1, 1:res + 1]

        f = np.reshape(np.arange(res * res), (res, res))          # up faces
        g = np.reshape(np.arange(res * res, 2 * res * res), (res, res))  # dn faces
        o = np.full((res,), -1)

        drte = np.array([1, 0])
        gche = np.array([-1, 0])
        haut = np.array([0, 1])
        bas  = np.array([0, -1])

        flip_square_no = np.full((res, res), 1)
        flip_line_no   = np.full((res,), 1)

        dn_eds = fromarrays((a[:, 0], c[:, 0], g[:, 0], o,
                              fa[:, 0], fc[:, 0] - fa[:, 0],
                              np.tile(haut, (res, 1)), flip_line_no), dtype=e_type)  # left
        up_eds = fromarrays((b[:, -1], d[:, -1], f[:, -1], o,
                              fb[:, -1], fd[:, -1] - fb[:, -1],
                              np.tile(bas, (res, 1)), flip_line_no), dtype=e_type)   # right
        lf_eds = fromarrays((a[0, :], b[0, :], f[0, :], o,
                              fa[0, :], fb[0, :] - fa[0, :],
                              np.tile(drte, (res, 1)), flip_line_no), dtype=e_type)  # bottom
        rt_eds = fromarrays((c[-1, :], d[-1, :], g[-1, :], o,
                              fc[-1, :], fd[-1, :] - fc[-1, :],
                              np.tile(gche, (res, 1)), flip_line_no), dtype=e_type)  # top
        dg_eds = fromarrays((a, d, f, g, fa, fd - fa,
                              np.tile(drte, (res, res, 1)), flip_square_no), dtype=e_type)
        hr_eds = fromarrays((a[:, 1:], c[:, 1:], g[:, 1:], f[:, :-1],
                              fa[:, 1:], fc[:, 1:] - fa[:, 1:],
                              np.tile(drte, (res, res - 1, 1)), flip_square_no[:, 1:]), dtype=e_type)
        vr_eds = fromarrays((a[1:, :], b[1:, :], f[1:, :], g[:-1, :],
                              fa[1:, :], fb[1:, :] - fa[1:, :],
                              np.tile(drte, (res - 1, res, 1)), flip_square_no[1:, :]), dtype=e_type)

        if self.u_identify == "cy":
            lf_eds["g"] = g[-1, :]
        if self.v_identify == "cy":
            dn_eds["g"] = f[:, -1]
        if self.u_identify == "mo":
            lf_eds["g"] = np.flip(g[-1, :], 0)
            dg_eds["flip"][-1, :] = -1
            hr_eds["flip"][-1, :] = -1
            dn_eds["flip"][-1] = up_eds["flip"][-1] = -1
        if self.v_identify == "mo":
            dg_eds["flip"][:, -1] = -1
            dn_eds["g"] = np.flip(f[:, -1], 0)
            vr_eds["flip"][:, -1] = -1
            lf_eds["flip"][-1] = rt_eds["flip"][-1] = -1

        if (self.u_identify in ("cy", "mo")) and self.v_identify == "no":
            self.eds = np.array(np.concatenate((dg_eds.flatten(), hr_eds.flatten(),
                                                 vr_eds.flatten(), lf_eds, dn_eds, up_eds)), dtype=e_type)
            self.b_index = res * res + 2 * (res - 1) * res + res
        elif self.u_identify == "no" and self.v_identify in ("cy", "mo"):
            self.eds = np.array(np.concatenate((dg_eds.flatten(), hr_eds.flatten(),
                                                 vr_eds.flatten(), dn_eds, lf_eds, rt_eds)), dtype=e_type)
            self.b_index = res * res + 2 * (res - 1) * res + res
        elif self.u_identify != "no" and self.v_identify != "no":
            self.eds = np.array(np.concatenate((dg_eds.flatten(), hr_eds.flatten(),
                                                 vr_eds.flatten(), lf_eds, dn_eds)), dtype=e_type)
            self.b_index = res * res + 2 * (res - 1) * res + 2 * res
        else:
            self.eds = np.array(np.concatenate((dg_eds.flatten(), hr_eds.flatten(),
                                                 vr_eds.flatten(), lf_eds, rt_eds, up_eds, dn_eds)), dtype=e_type)
            self.b_index = res * res + 2 * (res - 1) * res

        self.ieds = self.eds[:self.b_index]
        self.beds = self.eds[self.b_index:]
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "triangulation"))

        # Corners: flat integer indices (hashable directly as ints)
        if self.u_identify == self.v_identify == "no":
            self.corners = {int(points[0, 0]), int(points[0, res]),
                            int(points[res, 0]), int(points[res, res])}
        else:
            self.corners = set()

        t0 = time.perf_counter()
        if len(self.beds) > 0:
            self.b_lines = make_lines(self.beds["p"], self.beds["q"], self.corners)
        else:
            self.b_lines = []
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "make boundary lines"))

        # les generateurs du groupe par lequel on quotiente R^2.

        def f_u(pt, i):  # itéré i fois, où i est un entier relatif
            if self.u_identify == "cy":
                return pt + i * np.array([self.u_max - self.u_min, 0])
            elif self.u_identify == "mo":
                y = pt[1] if i % 2 == 0 else self.v_min + self.v_max - pt[1]
                return np.array([pt[0] + i * (self.u_max - self.u_min), y])
            else:
                return pt

        self.generator_u = f_u

        def f_v(pt, i):
            if self.v_identify == "cy":
                return pt + i * np.array([0, self.v_max - self.v_min])
            elif self.v_identify == "mo":
                x = pt[0] if i % 2 == 0 else self.u_max + self.u_min - pt[0]
                return np.array([x, pt[1] + i * (self.v_max - self.v_min)])
            else:
                return pt

        self.generator_v = f_v

        def close(pt, qt):  # renvoie le représentant du point qt le plus proche de pt
            m = math.floor((pt[0] - qt[0]) / (self.u_max - self.u_min))
            n = math.floor((pt[1] - qt[1]) / (self.v_max - self.v_min))
            four_pts = [
                self.generator_u(self.generator_v(qt, j), i)
                for i in range(m, m + 2)
                for j in range(n, n + 2)
            ]
            return min(four_pts, key=lambda q: norm(q - pt))

        def simple_close(pt, qt):
            return qt

        if self.u_identify != "no" or self.v_identify != "no":
            self.close = close
        else:
            self.close = simple_close  # s'il n'y a pas de raccord, il faut éviter de ralentir le programme avec la fonction close.

    def traitement(self, use_lp=True, use_newton=True, simplify_pts=None):
        self._use_newton = use_newton
        self.breaks = {
            "b": {},
            "c": {},
            "f": {},
            "g": {},
        }  # b=bord, c = contour, f,g = fictifs bord, contour (arêtes qui s'échappent vers un point de visi = 0)
        self.find_silhouette()  # trouve les lignes de plis et les cusps
        self.lines = []
        self.break_lines()  # découpe les bords et plis à chaque bord-pli ou cusp.(long)
        self.connect()  # ajoute des segments pour que l'ensemble des bords et plis soit connexe
        self.simplify_lines(n_points=simplify_pts)  # espace les points régulièrement (à l'écran) pour éviter les artefacts de discrétisation.
        self.line_bks = [[[] for j in range(len(l) - 1)] for l in self.lines]
        line, idx = self.intersections() # calcule les intersections/chgement de visi. Renvoie: point visible (long)
        # line, idx = 0, 0
        self.visibilities = {}
        if use_lp:
            self.visibilite(line, idx)   # LP visibility
        else:
            self._bfs_visibility(line, idx)  # BFS propagation (fast, no LP)
        #for p, v in self.visibilities.items(): # debug
        #    print("point:[%0.2f,%0.2f] visibilité:%i"%(*self.XY(np.array(self.S(*p))), v))

        self.origine = self.lines[line][idx]['fp'] #debug

    def find_silhouette(self):
        t0 = time.perf_counter()

        ######## Calcul des points de contour #########

        #### sélection des arêtes
        if self.SN_axis is self.SN_axis_persp:
            # Perspective: silhouette condition is N·(S−eye)=0
            eye = self.ax_param
            NZ_grid = (self.SN_grid * (self.S_grid - eye[:, np.newaxis])).sum(axis=0)
        else:
            # Orthographic: silhouette condition is N·axis=0
            NZ_grid = (
                self.SN_grid[0] * self.axis[0]
                + self.SN_grid[1] * self.axis[1]
                + self.SN_grid[2] * self.axis[2]
            )

        vec_i = np.where(
            NZ_grid[self.eds["p"]]
            * NZ_grid[self.eds["q"]]
            * self.eds["flip"]
            < 0
        )
        vec_i = vec_i[0]  # np.where renvoie un tuple, on prend le 1er élément.
        self.print(
            "[%0.3fs] %s"
            % (time.perf_counter() - t0, "sélection des arêtes avec signe changeant")
        )

        #### newton
        t0 = time.perf_counter()

        def NZbis(s, u, v, du, dv):
            return self.SN_axis(u + s * du, v + s * dv, *self.ax_param)

        init = np.full(len(vec_i), 0.5)
        if len(vec_i) > 0 and getattr(self, '_use_newton', True):
            vec_s = newton(
                NZbis,
                init,
                args=(
                    self.eds["fp"][vec_i, 0],
                    self.eds["fp"][vec_i, 1],
                    self.eds["fdp"][vec_i, 0],
                    self.eds["fdp"][vec_i, 1],
                ),
            )
            # Filter out diverged Newton iterates (NaN or outside [0,1]).
            # This can happen on disk boundary edges when the secant step
            # overshoots into a region where the surface is undefined (e.g.
            # sqrt of a negative number), propagating NaN through the
            # subsequent SVD call in kerdS.
            valid = np.isfinite(vec_s) & (vec_s >= 0) & (vec_s <= 1)
            vec_i = vec_i[valid]
            vec_s = vec_s[valid]
        elif len(vec_i) > 0:
            vec_s = np.full(len(vec_i), 0.5)  # midpoint when Newton disabled
        else:
            vec_s = np.array([])
        self.print(
            "[%0.3fs] %s"
            % (time.perf_counter() - t0, "calcul des points de contour avec newton")
        )

        ######## Ajout de la dir. orth. au contour dirigée vers la surf. aux pts du contour #########
        t0 = time.perf_counter()
        vec = fromarrays(
            (self.eds["fp"][vec_i, :], self.eds["fdp"][vec_i, :]),
            dtype=[("fp", "f", 2), ("fdp", "f", 2)],
        )

        pts_num = np.array(vec["fp"] + vec_s.reshape(len(vec_s), 1) * vec["fdp"])

        u_vec, v_vec = pts_num[:, 0], pts_num[:, 1]
        _, dir_vec = self.kerdS(u_vec, v_vec)

        self.print(
            "[%0.3fs] %s" % (time.perf_counter() - t0, "direction orth au contour")
        )

        ######## Calcul des lignes de contour
        t0 = time.perf_counter()

        self.c_pts = fromarrays((vec_i, vec_s, dir_vec.T, pts_num), dtype=c_type)
        if len(self.c_pts) != 0 :
            self.c_lines = make_lines(
                self.eds[self.c_pts["e"]]["f"], self.eds[self.c_pts["e"]]["g"]
            )
            # self.c_lines = self.make_c_lines(self.cpoints)
        else:
            self.c_lines = []
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "make contour lines"))

        ######## breaks au  bord #########
        t0 = time.perf_counter()

        for i, idx in enumerate(vec_i):
            if idx < self.b_index:  # idx ne correspond pas à une arête de bord
                continue
            v = self.bvis_chge(self.eds[idx], vec_s[i])
            self.breaks["b"][idx - self.b_index] = (vec_s[i], v, i)
            # print('visi change ',v)
            # self.add_break('b', i, idx - self.b_index,vec_s[i],v) # index dans c_pts, index dans self.beds
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "breaks au bord"))

        ######## Breaks de cusps
        t0 = time.perf_counter()
        cusps = []

        for l in self.c_lines:
            # ancienne version
            directions = dir_vec[:, l]
            signes = np.sum(directions[:, 0:-1] * directions[:, 1:], axis=0)

            indices = np.where(signes < 0)[0]
            for i in indices:
                p = np.array([u_vec[l[i]], v_vec[l[i]]])
                q = self.close(p, [u_vec[l[i + 1]], v_vec[l[i + 1]]])
                r = (p + q) / 2
                vis = 1 if np.inner(self.dS(*p, *(q - p)), self.axis) > 0 else -1
                cusps.append((l[i], l[i + 1], 1 / 2, r))
                self.breaks["c"][l[i]] = (0.5, vis, len(cusps) - 1)

        self.cusps = fromrecords(cusps, dtype=cusp_type)
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "breaks de cusps"))
        print('c_lines=%i b_lines=%i cusps=%i' % (len(self.c_lines), len(self.b_lines), len(self.cusps)))
        for ci, cl in enumerate(self.c_lines):
            print('  c_line %i: len=%i closed=%s' % (ci, len(cl), cl[-1]==cl[0] if len(cl)>0 else 'empty'))
        # print(len(self.breaks['b']),len(self.breaks['c']))

    def break_lines(self):
        t0 = time.perf_counter()

        # Precompute dirint for all beds at s=0 (p endpoint) and s=1 (q endpoint), vectorized.
        # dirint_pre[s, k] = 2D direction for beds[k] evaluated at parameter s.
        if len(self.beds) > 0:
            dirint_pre = np.empty((2, len(self.beds), 2))
            nu_u  = self.beds["dir"][:, 0];  nu_v  = self.beds["dir"][:, 1]
            fdp_u = self.beds["fdp"][:, 0];  fdp_v = self.beds["fdp"][:, 1]
            for si, fp_arr in enumerate([self.beds["fp"],
                                         self.beds["fp"] + self.beds["fdp"]]):
                u, v   = fp_arr[:, 0], fp_arr[:, 1]
                Su_v   = np.array(self.Su(u, v))             # (3, N)
                Sv_v   = np.array(self.Sv(u, v))             # (3, N)
                dS_nu  = Su_v * nu_u  + Sv_v * nu_v          # (3, N)
                dS_fdp = Su_v * fdp_u + Sv_v * fdp_v         # (3, N)
                vec_xy = self.proj_vec(u, v, dS_nu)          # (2, N)
                t_xy   = self.proj_vec(u, v, dS_fdp)         # (2, N)
                t_n    = t_xy / np.sqrt((t_xy**2).sum(axis=0))
                inner  = (vec_xy * t_n).sum(axis=0)
                vp     = vec_xy - inner * t_n
                dirint_pre[si] = (vp / np.sqrt((vp**2).sum(axis=0))).T  # (N, 2)

        for l in self.b_lines:
            lines = [[]]  # chaque b_line est découpée
            for index, i in enumerate(l):
                e = self.beds[i]
                e_idx = i if i >= 0 else i + len(self.beds)
                p = e["q"] if i < 0 else e["p"]
                q = e["q"] if i >= 0 else e["p"]
                sens = (
                    0 if i >= 0 else 1
                )  # indique si il y a inversion du sens de parcours.
                fp = e["fp"] if i >= 0 else e["fp"] + e["fdp"]
                fq = e["fp"] if i < 0 else e["fp"] + e["fdp"]

                lines[-1].append(
                    ("bb", e_idx + self.b_index, 0, dirint_pre[sens, e_idx], fp, fp)
                )  # la direction intérieure est calculée au point.
                if index == len(l) - 1:
                    continue  # le dernier point est à part, il sert juste à fermer

                if e_idx in self.breaks["b"]:
                    s, v, c_pt = self.breaks["b"][e_idx]
                    v = v if i >= 0 else -v
                    # la visib est le changement de  vis quand on entre dans la ligne, ça ne peut pas être > 0
                    # quand il y a un break au bord, la direction intérieure change, donc on la met à 0 au break
                    lines[-1].append(
                        (
                            "cb",
                            e_idx + self.b_index,
                            min(-v, 0),
                            0,
                            self.c_pts[c_pt]["fp"],
                            self.c_pts[c_pt]["fp"],
                        )
                    )
                    lines.append(
                        [
                            (
                                "cb",
                                e_idx + self.b_index,
                                min(v, 0),
                                0,
                                self.c_pts[c_pt]["fp"],
                                self.c_pts[c_pt]["fp"],
                            )
                        ]
                    )
                if q in self.corners:
                    # en cas de coin, on casse, car la direction intérieure est discontinue aux coins.
                    lines[-1].append(
                        (
                            "bb",
                            e_idx + self.b_index,
                            0,
                            dirint_pre[1 - sens, e_idx],
                            fq,
                            fq,
                        )
                    )
                    # lines[-1].append(('bb', e_idx + self.b_index, 0, 0, fq, fq))
                    lines.append([])
                    # print('corner')

            if lines[-1] == []:  # peut arriver si le dernier point était un coin
                print("pop")
                lines.pop()
            if (
                len(lines) > 1
            ):  # si il y a eu un break, le premier bout et le dernier sont ressoudés
                lines[0] = lines.pop() + lines[0][1:]
            for i, l in enumerate(lines):
                np_l = fromrecords(l, dtype=line_pt_type)
                self.lines.append(np_l)

        for l in self.c_lines:
            # Fast path: no cusps on this line → direct vectorized construction
            if not any(i in self.breaks["c"] for i in l[:-1]):
                idxs = np.asarray(l, dtype=np.intp)
                pts  = self.c_pts[idxs]
                np_l = np.empty(len(idxs), dtype=line_pt_type)
                np_l["type"]  = "cc"
                np_l["e"]     = pts["e"]
                np_l["v"]     = 0
                np_l["d"]     = pts["d"]
                np_l["fp"]    = pts["fp"]
                np_l["ixfp"]  = pts["fp"]
                self.lines.append(np_l)
                continue

            lines = [[]]  # chaque b_line est découpée (cusps present)

            for index, i in enumerate(l):
                point = self.c_pts[i]
                lines[-1].append(
                    ("cc", point["e"], 0, point["d"], point["fp"], point["fp"])
                )
                if index == len(l) - 1:
                    continue  # le dernier point est à part, il sert juste à fermer

                if i in self.breaks["c"]:
                    s, v, p = self.breaks["c"][i]  # nouveau self.breaks
                    # s, v, p = self.breaks['c'][i][0] # à ce stade, il y a au plus un break par arête.
                    pt = self.c_pts[self.cusps[p].p]
                    # la direction de la surface n'est pas définie au cusp. On la met à 0.
                    lines[-1].append(
                        (
                            "vc",
                            pt.e,
                            min(-v, 0),
                            0,
                            self.cusps[p]["fp"],
                            self.cusps[p]["fp"],
                        )
                    )

                    lines.append(
                        [
                            (
                                "vc",
                                pt.e,
                                min(v, 0),
                                0,
                                self.cusps[p]["fp"],
                                self.cusps[p]["fp"],
                            )
                        ]
                    )

            if (
                len(lines) > 1 and lines[0][0][1] == lines[-1][-1][1]
            ):  # si il y a eu un break, et que l était un loop, le premier bout et le dernier sont ressoudés
                lines[0] = lines.pop() + lines[0][1:]
            for i, l in enumerate(lines):
                np_l = fromrecords(l, dtype=line_pt_type)
                self.lines.append(np_l)
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "breaks lines"))

    def connect(self):
        t0 = time.perf_counter()

        def trim(
            l,
        ):  # on enlève les cusps de l pour ne pas connecter à un cusp. on renvoie les indices dans la ligne initiale, et la ligne tronquée.
            i, j = 0, len(l)
            if l["type"][0][0] == "v":
                i = min(5,len(l)//3) # on évite les extrêmités des lignes, dans la mesure du possible. Ici: marge de 10, à trois ça déconne parfois
            if l["type"][-1][0] == "v":
                j = len(l)-min(5, len(l)//3)
            return np.arange(i, j), l[i:j]

        self.cc_s = (
            self.connected_components()
        )  # liste de listes d'indices de lignes, qui sont connectées.
        print('nombre de composantes connexes %i'%(len(self.cc_s),))
        for ci, cc in enumerate(self.cc_s):
            types = [self.lines[l][0]['type'] for l in cc]
            print('  composante %i: %i lignes, types=%s' % (ci, len(cc), types))

        if len(self.cc_s) < 2:
            return
        print('connect !!')
        # distances entre les composantes connexes, links dit quels points sont reliés : indice de ligne, indice dans la ligne
        trees = {}
        l_index = {}
        e_index = {}
        links = np.empty(
            (len(self.cc_s), len(self.cc_s)),
            dtype=[("il", "int"), ("ie", "int"), ("jl", "int"), ("je", "int")],
        )
        distances = np.empty((len(self.cc_s), len(self.cc_s)), "f")
        for i, cc in enumerate(self.cc_s):
            distances[i, i] = 0
            l_index[i] = np.concatenate(
                [np.full((len(trim(self.lines[l])[0]),), l) for l in cc]
            )
            e_index[i] = np.concatenate([trim(self.lines[l])[0] for l in cc])
            pts = np.concatenate([trim(self.lines[l])[1]["ixfp"] for l in cc])
            trees[i] = cKDTree(pts)
            for j in range(i):
                i_point, j_point = len(e_index[i]) // 2, len(e_index[j]) // 2
                for k in range(2):
                    j_point = trees[j].query(
                        self.lines[l_index[i][i_point]][e_index[i][i_point]]["ixfp"]
                    )[1]
                    i_point = trees[i].query(
                        self.lines[l_index[j][j_point]][e_index[j][j_point]]["ixfp"]
                    )[1]
                # print('ipoint %i jpoint %i' % (i_point, j_point))
                il, ie, jl, je = (
                    l_index[i][i_point],
                    e_index[i][i_point],
                    l_index[j][j_point],
                    e_index[j][j_point],
                )
                links[i, j], links[j, i] = (il, ie, jl, je), (jl, je, il, ie)
                distances[i, j] = distances[j, i] = norm(
                    self.lines[il][ie]["ixfp"] - self.lines[jl][je]["ixfp"]
                )
                distances[i, i] = distances[j, j] = max(
                    distances[i, i], distances[j, j], distances[i, j] + 1
                )  # il faut le max sur la diagonale

        # recherche de la façon la plus courte de relier les composantes connexes
        i0, j0 = np.argwhere(distances == np.min(distances))[0]
        big_cc = set([i0, j0])
        E = set(range(len(self.cc_s)))
        cuts = [(i0, j0)]
        while E - big_cc:
            A = list(big_cc)
            B = list(E - big_cc)
            i, j = np.argwhere(
                distances[np.ix_(A, B)] == np.min(distances[np.ix_(A, B)])
            )[0]
            big_cc.add(B[j])
            cuts.append((A[i], B[j]))

        # Création des lignes de connection dic contient les indices des points de connection sur chaque ligne
        dic = defaultdict(list)
        for i, j in cuts:
            il, ie, jl, je = links[i, j]
            print('  cut: comp%i ligne%i[%i](type=%s) <-> comp%i ligne%i[%i](type=%s)' % (
                i, il, ie, self.lines[il][ie]['type'],
                j, jl, je, self.lines[jl][je]['type']))
            bisect.insort(dic[il], ie)
            bisect.insort(dic[jl], je)
            p, q = self.lines[il][ie], self.lines[jl][je]
            self.lines.append(self.connection(p, q))
        # découpage des lignes pour que les points de connection soient des extrêmités.
        for l in dic:
            # print(' line ', l, '   cut at edges  ', [self.lines[l][i].e for i in dic[l]])
            if dic[l][0] != 0:
                dic[l] = [0] + dic[l]
            if dic[l][-1] != len(self.lines[l]) - 1:
                dic[l].append(len(self.lines[l]) - 1)

            lines = [
                self.lines[l][dic[l][i] : dic[l][i + 1] + 1]
                for i in range(len(dic[l]) - 1)
            ]
            self.lines[l] = lines[0]
            self.lines = self.lines + lines[1:]
            # print('results in lines ', [[m[0]['e'], m[-1]['e']] for m in lines])

        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "connect"))

    def simplify_lines(self, n_points=None):
        t0 = time.perf_counter()

        # define sample
        def sample_fn(x):
            return x * x / (1 + x)

        p = []  # list of lists of image points of the lines
        # n_points is "output points per unit of screen width"; multiply by 50 to get sample count.
        # Old default: 50 * max(res // 5, 100).  With n_points=80 → 50*80 = 4000 ≈ same order.
        _n = (n_points * 50) if n_points is not None else 50 * max(self.res // 5, 100)
        sample = sample_fn(np.linspace(0, 50, _n))

        im_min = np.inf
        im_max = -np.inf

        for l in self.lines:
            p.append(self.XY(np.array(self.S(l["fp"][:, 0], l["fp"][:, 1]))).T)
            lmin, lmax = np.nanmin(p[-1]), np.nanmax(p[-1])
            if np.isfinite(lmin): im_min = min(lmin, im_min)
            if np.isfinite(lmax): im_max = max(lmax, im_max)

        # vector which scales XY coordinates to a square of sidelength 1
        scale = np.array([1, 1]) / (im_max - im_min)

        for i, l in enumerate(self.lines):
            origin, end = l[0], l[-1]
            intervals = norm(scale * (p[i][1:] - p[i][:-1]), axis=1) # interval[i] = dist(p[i], p[i+1])
            lengths = np.cumsum(intervals)  # lengths[i+1] - lengths[i] = interval[i+1], lengths[i] = dist(p[0], p[i+1])
            index = np.searchsorted(sample, lengths[-1] / 2) # sample[index - 1] < l/2 <= sample[index]
            points = np.concatenate(
                (sample[:index], lengths[-1] - sample[index - 1 :: -1]) 
            )[
                1:-1
            ]  # points[i] = sample[i+1]

            j = np.searchsorted(lengths, points) # lengths[j[k]-1] < points[k] <= lengths[j[k]]  
            line = np.empty((len(j) + 2,), dtype=line_pt_type)
            l_fp = self.relevement(l["fp"]) # comme on interpole, il faut un relèvement dans le cas 'cy' ou 'mo'

            # p[k] = ((length[j[k]] - points[k]) * p[j[k]] + (point[k] - lenght[j[k] - 1]) * p[j[k]])/interval[j[k]]
            fp_s = (
                (points - lengths[j] + intervals[j]) * l_fp[j + 1].T
                + (lengths[j] - points) * l_fp[j].T
            ) / intervals[j]
            dirint_s = (
                (points - lengths[j] + intervals[j]) * l[j + 1]["d"].T
                + (lengths[j] - points) * l[j]["d"].T
            ) / intervals[j]
            fp_s = fp_s.T
            dirint_s = (
                dirint_s / (0.1 + norm(dirint_s, axis=0))
            ).T  # .1 pour éviter div par 0

            line[0] = origin
            line[-1] = end
            line["type"][1:-1] = l[1]["type"]
            line["v"][1:-1] = 0
            line["d"][1:-1] = dirint_s
            line["fp"][1:-1] = fp_s
            self.lines[i] = line
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "simplify lines"))

    def intersections(self):
        t0 = time.perf_counter()
        # print(len(self.lines))

        bag_bis = np.empty((sum([len(l) - 1 for l in self.lines]),), dtype=bagtype_bis)
        # index = 0
        line_indices = [0]
        for i, l in enumerate(self.lines):
            fp, fq = l["fp"][:-1], l["fp"][1:]  # les points dans le domaine
            p = self.XY(np.array(self.S(fp[:, 0], fp[:, 1]))).T  # les points à l'écran
            q = self.XY(np.array(self.S(fq[:, 0], fq[:, 1]))).T
            zp = self.Z(np.array(self.S(fp[:, 0], fp[:, 1])))
            zq = self.Z(np.array(self.S(fq[:, 0], fq[:, 1])))
            d = (l["d"][:-1] + l["d"][1:]) / 2
            # d = l['d'][:-1] # c'est le premier point de l'arête qui donne la direction intérieure
            e = l["e"][:-1]
            length = len(l) - 1
            n = (
                np.full((length,), 1)
                if l[0]["type"][1] == "b"
                else np.full((length,), 2)
                if l[0]["type"][1] == "c"
                else np.full((length,), 0)
            )
            l_idx = np.full((length,), i)
            e_idx = np.arange(length)
            line_indices.append(line_indices[i] + length)
            bag_bis[line_indices[i] : line_indices[i + 1]] = fromarrays(
                (
                    np.minimum(p[:, 0], q[:, 0]), # le min de la coord 'x' sur l'arête
                    l_idx,
                    e_idx,
                    e,
                    n,
                    d,
                    p,
                    zp,
                    q - p,
                    zq - zp,
                    fp,
                    fq - fp,
                ),
                dtype=bagtype_bis,
            )

        self.print(
            "[%0.3fs] %s"
            % (time.perf_counter() - t0, "Intersections bis - préparation")
        )

        t0 = time.perf_counter()
        u_tol, v_tol = (self.u_max - self.u_min) / self.res, (
            self.v_max - self.v_min
        ) / self.res
        tol = 0.5 * np.array( # sert pour exclure certaines intersections de segments, s'ils sont trop proches dans le domaine
            [u_tol, v_tol]
        )  # pas clair par quoi il faut multiplier, avec .O1 pose des problèmes en raccord cylindre

        bb_im = np.concatenate(
            (
                np.minimum(bag_bis["p"], bag_bis["p"] + bag_bis["dp"]),
                np.maximum(bag_bis["p"], bag_bis["p"] + bag_bis["dp"]),
            ),
            axis=1,
        )

        use_batch = (self.u_identify == "no" and self.v_identify == "no")
        if use_batch:
            bag_p_c   = np.ascontiguousarray(bag_bis["p"])
            bag_dp_c  = np.ascontiguousarray(bag_bis["dp"])
            bag_fp_c  = np.ascontiguousarray(bag_bis["fp"])
            bag_fdp_c = np.ascontiguousarray(bag_bis["fdp"])
            bag_d_c   = np.ascontiguousarray(bag_bis["d"])
            bag_z     = bag_bis["z"]
            bag_dz    = bag_bis["dz"]
            bag_n     = bag_bis["n"]
            bag_li    = bag_bis["l_idx"]
            bag_ei    = bag_bis["e_idx"]

            # Fully vectorized sweep-line pair collection — no Python loop over edges.
            # Sort edges by x_min; for each sorted edge i find all j>i with x_min[j]<=x_max[i]
            # (x-overlap guaranteed by construction), then filter y-overlap.
            n = len(bag_bis)
            order   = np.argsort(bb_im[:, 0])
            xmin_s  = bb_im[order, 0];  xmax_s = bb_im[order, 2]
            ymin_s  = bb_im[order, 1];  ymax_s = bb_im[order, 3]
            end_idx = np.searchsorted(xmin_s, xmax_s, side='right')
            counts  = np.maximum(0, end_idx - np.arange(1, n + 1)).astype(np.intp)
            total_pairs = int(counts.sum())

            if total_pairs > 0:
                offsets = np.empty(n + 1, dtype=np.intp)
                offsets[0] = 0
                np.cumsum(counts, out=offsets[1:])
                i_s = np.repeat(np.arange(n, dtype=np.intp), counts)
                j_s = i_s + 1 + (np.arange(total_pairs, dtype=np.intp) - offsets[i_s])

                # y-overlap filter
                y_ok = (ymin_s[j_s] <= ymax_s[i_s]) & (ymin_s[i_s] <= ymax_s[j_s])
                i_s = i_s[y_ok];  j_s = j_s[y_ok]

                if len(i_s) > 0:
                    # Map sorted → original indices, enforce lo < hi
                    i_orig = order[i_s];  j_orig = order[j_s]
                    lo = np.minimum(i_orig, j_orig);  hi = np.maximum(i_orig, j_orig)

                    # Domain filter: keep pairs whose param-space intervals are disjoint
                    e_fp_d = bag_fp_c[lo];  e_fq_d = e_fp_d + bag_fdp_c[lo]
                    f_fp_d = bag_fp_c[hi];  f_fq_d = f_fp_d + bag_fdp_c[hi]
                    emax_d = np.maximum(e_fp_d, e_fq_d)
                    emin_d = np.minimum(e_fp_d, e_fq_d)
                    fmax_d = np.maximum(f_fp_d, f_fq_d)
                    fmin_d = np.minimum(f_fp_d, f_fq_d)
                    dom_ok = (
                        (fmin_d[:, 0] > tol[0] + emax_d[:, 0]) |
                        (fmin_d[:, 1] > tol[1] + emax_d[:, 1]) |
                        (emin_d[:, 0] > tol[0] + fmax_d[:, 0]) |
                        (emin_d[:, 1] > tol[1] + fmax_d[:, 1])
                    )
                    ip = lo[dom_ok];  jp = hi[dom_ok]

                    if len(ip) > 0:
                        # Batch intersection arithmetic
                        e_pxy  = bag_p_c[ip];  e_dpxy = bag_dp_c[ip]
                        f_pxy  = bag_p_c[jp];  f_dpxy = bag_dp_c[jp]
                        D   = e_dpxy[:, 0] * f_dpxy[:, 1] - e_dpxy[:, 1] * f_dpxy[:, 0]
                        nz  = D != 0
                        D_s = np.where(nz, D, 1.0)
                        uv  = e_pxy - f_pxy
                        s   = (f_dpxy[:, 0] * uv[:, 1] - f_dpxy[:, 1] * uv[:, 0]) / D_s
                        t   = (e_dpxy[:, 0] * uv[:, 1] - e_dpxy[:, 1] * uv[:, 0]) / D_s
                        hit = nz & (s > 0) & (s < 1) & (t > 0) & (t < 1)

                        if np.any(hit):
                            ip_h = ip[hit];  jp_h = jp[hit]
                            s_h  = s[hit];   t_h  = t[hit]
                            ze   = bag_z[ip_h] + s_h * bag_dz[ip_h]
                            zf   = bag_z[jp_h] + t_h * bag_dz[jp_h]
                            maxz = np.maximum(np.abs(ze),
                                   np.maximum(np.abs(zf),
                                   np.maximum(np.abs(bag_dz[ip_h]),
                                   np.maximum(np.abs(bag_dz[jp_h]), 1e-10))))
                            near_eq = np.abs(ze - zf) <= 1e-4 * maxz
                            e_gt_f  = (~near_eq) & (ze > zf)
                            f_gt_e  = (~near_eq) & (zf > ze)

                            if np.any(e_gt_f):  # ip above jp: breakpoint on jp-line at t
                                ip_egt = ip_h[e_gt_f];  jp_egt = jp_h[e_gt_f]
                                tt_egt = t_h[e_gt_f]
                                inner = (bag_dp_c[jp_egt, 0] * bag_d_c[ip_egt, 0] +
                                         bag_dp_c[jp_egt, 1] * bag_d_c[ip_egt, 1])
                                vs = np.where(inner > 0, -bag_n[ip_egt], bag_n[ip_egt])
                                for k_idx in range(len(jp_egt)):
                                    if vs[k_idx] != 0:
                                        self.add_line_bk(int(bag_li[jp_egt[k_idx]]),
                                                         int(bag_ei[jp_egt[k_idx]]),
                                                         tt_egt[k_idx], int(vs[k_idx]))

                            if np.any(f_gt_e):  # jp above ip: breakpoint on ip-line at s
                                ip_fgt = ip_h[f_gt_e];  jp_fgt = jp_h[f_gt_e]
                                ss_fgt = s_h[f_gt_e]
                                inner2 = (bag_dp_c[ip_fgt, 0] * bag_d_c[jp_fgt, 0] +
                                          bag_dp_c[ip_fgt, 1] * bag_d_c[jp_fgt, 1])
                                vs2 = np.where(inner2 > 0, -bag_n[jp_fgt], bag_n[jp_fgt])
                                for k_idx in range(len(ip_fgt)):
                                    if vs2[k_idx] != 0:
                                        self.add_line_bk(int(bag_li[ip_fgt[k_idx]]),
                                                         int(bag_ei[ip_fgt[k_idx]]),
                                                         ss_fgt[k_idx], int(vs2[k_idx]))

                            if np.any(near_eq):  # rare: near-equal Z, re-evaluate with S
                                for k_idx in np.where(near_eq)[0]:
                                    self.edge_intersect(bag_bis[ip_h[k_idx]],
                                                        bag_bis[jp_h[k_idx]])
        else:
            # Periodic surfaces: numpy bbox overlap check (close() needed for wrap-around)
            xmin = bb_im[:, 0]; xmax = bb_im[:, 2]
            ymin = bb_im[:, 1]; ymax = bb_im[:, 3]
            for i in range(len(bag_bis) - 1):
                js = np.where(
                    (xmin[i+1:] <= xmax[i]) & (xmin[i] <= xmax[i+1:]) &
                    (ymin[i+1:] <= ymax[i]) & (ymin[i] <= ymax[i+1:])
                )[0] + (i + 1)
                if not len(js):
                    continue
                e = bag_bis[i]
                e_p_fp = e["fp"]
                e_q_fp = e_p_fp + e["fdp"]
                efmax = np.maximum(e_p_fp, e_q_fp)
                efmin = np.minimum(e_p_fp, e_q_fp)
                for j in js:
                    f = bag_bis[j]
                    f_p = self.close(e_p_fp, f["fp"])
                    f_q = f_p + f["fdp"]
                    if np.any(np.minimum(f_p, f_q) > tol + efmax) or np.any(
                        efmin > tol + np.maximum(f_p, f_q)
                    ):
                        self.edge_intersect(e, f)

        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "Intersections bis"))

        bag_bis.sort(order=("x"))
        e = bag_bis[0]  # arête où le min de x est atteint
        if e["x"] == e["p"][0]:  # c'est p qui a le min de x
            return e["l_idx"], e["e_idx"]
        else:  # sinon c'est q qui a le min de x.
            return e["l_idx"], e["e_idx"] + 1

    def _bfs_visibility(self, line, idx):
        """BFS propagation of visibility (pre-LP implementation)."""
        t0 = time.perf_counter()
        # Reset any stale LP results so plot_for_browser uses current line_bks and l[0]["v"]
        self.opt_line_bks = self.line_bks
        self.opt_cp = None
        visited = set()
        unseen = set()
        lines_dic = defaultdict(set)
        for i, l in enumerate(self.lines):
            p, q = tuple(l[0]["ixfp"]), tuple(l[-1]["ixfp"])
            lines_dic[p].add((i, 0))
            lines_dic[q].add((i, -1))
            unseen.update([p, q])

        v = (
            self.lines[line][idx]["v"]
            if (idx == 0 or idx == len(self.lines[line]) - 1)
            else 0
        )

        vis = self.propagation(line, idx, v)
        p, q = self.lines[line][0], self.lines[line][-1]
        point = q if idx == 0 else p
        self.visibilities[tuple(point["ixfp"])] = vis
        visited.add(tuple(point["ixfp"]))
        unseen.discard(tuple(point["ixfp"]))

        while visited:
            pt = visited.pop()
            for i, j in lines_dic[pt]:
                l = self.lines[i]
                qt = tuple(l[-1 - j]["ixfp"])
                if qt in unseen:
                    self.visibilities[qt] = self.propagation(
                        i, j, self.visibilities[pt] + self.lines[i][j]["v"]
                    )
                    unseen.discard(qt)
                    visited.add(qt)
            visited.discard(pt)

        # _bfs_visibility doesn't set opt_cp / opt_line_bks — plot_for_browser uses fallbacks
        self.print("[%0.3fs] Visibilité BFS" % (time.perf_counter() - t0))

    def visibilite(self, line, idx):
        # LP complet: V[j], cp[i], cq[i], b[k] variables libres; tcp, tcq, tb slacks L1.
        # Toujours faisable (solution nulle satisfait toutes les contraintes).
        # Minimize Σtcp + Σtcq + Σtb
        t0 = time.perf_counter()
        self.visibilities = {}
        m = len(self.lines)
        self.opt_line_bks = self.line_bks
        self.opt_cp = [float(self.lines[i][0]["v"]) for i in range(m)]
        if m == 0:
            return

        # Node index
        nodes = {}
        for l in self.lines:
            for pt in [tuple(l[0]["ixfp"]), tuple(l[-1]["ixfp"])]:
                if pt not in nodes:
                    nodes[pt] = len(nodes)
        n = len(nodes)

        p_idx = [nodes[tuple(self.lines[i][0]["ixfp"])] for i in range(m)]
        q_idx = [nodes[tuple(self.lines[i][-1]["ixfp"])] for i in range(m)]

        orig_cp = [float(self.lines[i][0]["v"]) for i in range(m)]
        orig_cq = [float(self.lines[i][-1]["v"]) for i in range(m)]

        # Énumération globale des breakpoints
        line_bps = [[] for _ in range(m)]
        bp_vals = []
        for i in range(m):
            for j in range(len(self.lines[i]) - 1):
                for s, v in self.line_bks[i][j]:
                    line_bps[i].append(len(bp_vals))
                    bp_vals.append(float(v))
        B = len(bp_vals)

        # Segments: (ligne i, nb de breakpoints dans le préfixe)
        seg_info = []
        for i in range(m):
            for j in range(len(line_bps[i]) + 1):
                seg_info.append((i, j))
        nsegs = len(seg_info)

        # Composantes connexes
        adj = defaultdict(list)
        for i in range(m):
            adj[p_idx[i]].append(i)
            adj[q_idx[i]].append(i)
        seen_lines = set()
        components = []
        for start in range(m):
            if start in seen_lines:
                continue
            comp, stack = [], [start]
            while stack:
                li = stack.pop()
                if li in seen_lines:
                    continue
                seen_lines.add(li)
                comp.append(li)
                for nb in adj[p_idx[li]] + adj[q_idx[li]]:
                    if nb not in seen_lines:
                        stack.append(nb)
            components.append(comp)

        # Layout: [V(n), cp(m), cq(m), b(B), tcp(m), tcq(m), tb(B)]
        i_V   = 0
        i_cp  = n
        i_cq  = n + m
        i_b   = n + 2*m
        i_tcp = n + 2*m + B
        i_tcq = n + 3*m + B
        i_tb  = n + 4*m + B
        nvars = n + 4*m + 2*B

        c_obj = np.zeros(nvars)
        c_obj[i_tcp:i_tcp+m] = 1.0
        c_obj[i_tcq:i_tcq+m] = 1.0
        c_obj[i_tb:i_tb+B]   = 1.0

        # Égalités: delta + ancre
        eq_rows, eq_cols, eq_vals, b_eq_list = [], [], [], []
        eq_row = 0

        # (1) V[q] - V[p] - cp[i] + cq[i] - Σb[k on i] = 0
        for i in range(m):
            bk_cols = [i_b+k for k in line_bps[i]]
            cr = [i_V+q_idx[i], i_V+p_idx[i], i_cp+i, i_cq+i] + bk_cols
            vr = [1., -1., -1., 1.] + [-1.]*len(bk_cols)
            eq_rows += [eq_row]*len(cr); eq_cols += cr; eq_vals += vr
            b_eq_list.append(0.); eq_row += 1

        # (2) Ancre: cp[line] + V[p_line] + Σ_prefix b[k] = 0
        main_comp_lines = next(comp for comp in components if line in comp)
        prefix_count = sum(len(self.line_bks[line][j]) for j in range(idx))
        anchor_prefix_bks = [i_b+k for k in line_bps[line][:prefix_count]]
        for comp in components:
            if comp is main_comp_lines:
                cr = [i_cp+line, i_V+p_idx[line]] + anchor_prefix_bks
                vr = [1., 1.] + [1.]*len(anchor_prefix_bks)
                eq_rows += [eq_row]*len(cr); eq_cols += cr; eq_vals += vr
                b_eq_list.append(0.)
            else:
                eq_rows.append(eq_row); eq_cols.append(i_V+p_idx[comp[0]]); eq_vals.append(1.)
                b_eq_list.append(0.)
            eq_row += 1

        A_eq = csc_matrix((eq_vals, (eq_rows, eq_cols)), shape=(eq_row, nvars))
        b_eq = np.array(b_eq_list)

        # Inégalités: vis ≤ 0 (dures) + slacks L1
        rows, cols, vals, b_ub_list = [], [], [], []
        row = 0

        # (3) vis ≤ 0 par segment: cp[i] + V[p_i] + Σ_prefix b[k] ≤ 0
        for seg_i, seg_j in seg_info:
            prefix_bk = [i_b+k for k in line_bps[seg_i][:seg_j]]
            cr = [i_cp+seg_i, i_V+p_idx[seg_i]] + prefix_bk
            rows += [row]*len(cr); cols += cr; vals += [1., 1.] + [1.]*len(prefix_bk)
            b_ub_list.append(0.); row += 1

        # (4) Slacks L1 pour cp, cq, b
        for i in range(m):
            rows += [row,row]; cols += [i_cp+i, i_tcp+i]; vals += [1., -1.]
            b_ub_list.append(orig_cp[i]); row += 1
            rows += [row,row]; cols += [i_cp+i, i_tcp+i]; vals += [-1., -1.]
            b_ub_list.append(-orig_cp[i]); row += 1
            rows += [row,row]; cols += [i_cq+i, i_tcq+i]; vals += [1., -1.]
            b_ub_list.append(orig_cq[i]); row += 1
            rows += [row,row]; cols += [i_cq+i, i_tcq+i]; vals += [-1., -1.]
            b_ub_list.append(-orig_cq[i]); row += 1
        for k in range(B):
            rows += [row,row]; cols += [i_b+k, i_tb+k]; vals += [1., -1.]
            b_ub_list.append(bp_vals[k]); row += 1
            rows += [row,row]; cols += [i_b+k, i_tb+k]; vals += [-1., -1.]
            b_ub_list.append(-bp_vals[k]); row += 1

        A_ub = csc_matrix((vals, (rows, cols)), shape=(row, nvars)) if row > 0 else None
        b_ub = np.array(b_ub_list) if row > 0 else None

        # V, cp, cq, b libres; tcp, tcq, tb ≥ 0
        bounds = [(None, None)]*(n + 2*m + B) + [(0, None)]*(2*m + B)

        result = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                         bounds=bounds, method='highs')

        if not result.success:
            self.print("L1 visibility LP failed: " + result.message)
            return

        V      = result.x[i_V:i_V+n]
        opt_cp = result.x[i_cp:i_cp+m]
        opt_b  = result.x[i_b:i_b+B]

        for pt, j in nodes.items():
            self.visibilities[pt] = int(round(V[j]))
        self.opt_cp = [round(float(opt_cp[i])) for i in range(m)]

        # Reconstruire line_bks avec valeurs optimisées
        self.opt_line_bks = [[[] for _ in range(len(self.lines[i]) - 1)] for i in range(m)]
        bp_counter = 0
        for i in range(m):
            for j in range(len(self.lines[i]) - 1):
                for s, v in self.line_bks[i][j]:
                    self.opt_line_bks[i][j].append((s, round(opt_b[bp_counter])))
                    bp_counter += 1

        self.print("[%0.3fs] Visibilité LP (n=%d, m=%d, B=%d, segs=%d, %d comp)" %
                   (time.perf_counter() - t0, n, m, B, nsegs, len(components)))
        # DIAG: per-line visibility
        for i in range(m):
            eff_vis = int(round(V[p_idx[i]])) + int(round(float(opt_cp[i])))
            # collect breakpoints for this line
            bk_list = []
            ki = 0
            for j in range(len(self.lines[i]) - 1):
                for s, _ in self.line_bks[i][j]:
                    bk_list.append("(j=%i,s=%.2f,b=%+.0f)" % (j, s, round(float(opt_b[line_bps[i][ki]]))))
                    ki += 1
            bk_str = " bks=[%s]" % ",".join(bk_list) if bk_list else ""
            print("  line%i type=%s len=%i  V[p]=%+.0f cp=%+.0f cq=%+.0f  eff_vis=%i orig_cp=%+.0f orig_cq=%+.0f%s" % (
                i, self.lines[i][0]["type"], len(self.lines[i]),
                round(V[p_idx[i]]), round(float(opt_cp[i])), round(float(result.x[i_cq+i])),
                eff_vis, orig_cp[i], orig_cq[i], bk_str))

    def plot_for_browser(self):
        t0 = time.perf_counter()

        bks    = getattr(self, 'opt_line_bks', self.line_bks)
        opt_cp = getattr(self, 'opt_cp', None)
        lines  = defaultdict(list)

        # --- Phase 1: batch-evaluate S and XY for all stored line points ---
        # Concatenate all l["fp"] arrays; line i starts at line_offsets[i].
        all_fp = np.concatenate([l["fp"] for l in self.lines], axis=0)  # (total, 2)
        line_offsets = np.zeros(len(self.lines) + 1, dtype=np.intp)
        for i, l in enumerate(self.lines):
            line_offsets[i + 1] = line_offsets[i] + len(l)
        xy_pts = self.XY(np.array(self.S(all_fp[:, 0], all_fp[:, 1])))  # (2, total)

        # Collect all interpolated breakpoint positions (non-?? lines only)
        bk_pts_list = []
        for i, l in enumerate(self.lines):
            if l[0]["type"] == "??":
                continue
            for j in range(len(l) - 1):
                p, q = l[j]["fp"], l[j + 1]["fp"]
                for s, v in bks[i][j]:
                    bk_pts_list.append((1.0 - s) * p + s * q)

        if bk_pts_list:
            bk_pts = np.array(bk_pts_list)  # (K, 2)
            xy_bk  = self.XY(np.array(self.S(bk_pts[:, 0], bk_pts[:, 1])))  # (2, K)
        bk_idx = 0

        # --- Phase 2: reconstruct lines dict using precomputed projections ---
        for i, l in enumerate(self.lines):
            if l[0]["type"] == "??":
                continue
            base = int(line_offsets[i])
            cp   = float(opt_cp[i]) if opt_cp is not None else float(l[0]["v"])
            vis  = cp
            if tuple(l[0]["ixfp"]) in self.visibilities:
                vis += self.visibilities[tuple(l[0]["ixfp"])]
            line = []
            for j in range(len(l) - 1):
                line.append(xy_pts[:, base + j])
                for s, v in bks[i][j]:
                    line.append(xy_bk[:, bk_idx]);  bk_idx += 1
                    lines[min(vis, 0)].append(line)
                    vis  = vis + v
                    line = [line[-1]]
            line.append(xy_pts[:, base + len(l) - 1])
            lines[min(vis, 0)].append(line)

        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "préparation du dessin"))
        return lines

    def plot_for_terminal(self, T0):

        t0 = time.perf_counter()
        fig = plt.figure(dpi=100)

        ax = fig.add_subplot(projection="3d")
        ax.set(proj_type="ortho")
        # axis = (  cos azim cos elev, sin azim cos elev, sin elev)
        self.elev = np.arcsin(self.axis[2])
        self.azim = np.arctan2(self.axis[1]/np.cos(self.elev), self.axis[0]/np.cos(self.elev))  * 180 / np.pi # atan2(y,x)
        self.elev = self.elev  * 180 / np.pi
        ax.view_init(elev=self.elev, azim=self.azim )
        ax.set_axis_off()
        fig.tight_layout()

        ptref = [0, 0]
        global point_mark
        point_mark = ax.scatter(
            *self.S(*self.lines[ptref[0]][ptref[1]]["fp"]), marker=""
        )
        ax.scatter(*self.S(*self.origine))# debug

        ## print visibilities
        # for i, l in enumerate(self.lines):
        #    print(' line number : ', i, 'visibilities : ', self.visibilities[tuple(l[0]['ixfp'])], ', ', self.visibilities[tuple(l[-1]['ixfp'])])

        def zoom(factor):
            limits = np.array([ax.get_xlim(), ax.get_ylim(), ax.get_zlim()])
            average = np.column_stack(
                ((limits[:, 0] + limits[:, 1]) / 2, (limits[:, 0] + limits[:, 1]) / 2)
            )
            new = (limits - average) / factor + average
            ax.set(xlim=new[0], ylim=new[1], zlim=new[2])

        def zoom_point():
            res = self.res
            p = self.lines[ptref[0]][ptref[1]]["fp"]
            pt = self.S(*p)
            limits = np.array([ax.get_xlim(), ax.get_ylim(), ax.get_zlim()])
            average = np.column_stack(
                ((limits[:, 0] + limits[:, 1]) / 2, (limits[:, 0] + limits[:, 1]) / 2)
            )
            new = (limits - average) / res + np.column_stack((pt, pt))
            ax.set(xlim=new[0], ylim=new[1], zlim=new[2])
            plt.show()

        def pt_change():
            global point_mark
            point_mark.remove()
            point_mark = ax.scatter(
                *self.S(*self.lines[ptref[0]][ptref[1]]["fp"]), marker="x"
            )

        def onkey(event):
            # print(event.key)
            if event.key == "u":
                zoom(0.9)
            if event.key == "y":
                zoom(1.1)
            if event.key == "ù":
                ptref[0] = (ptref[0] + 1) % len(self.lines)
                print(ptref[0], " ", ptref[1])
                pt_change()
            if event.key == "m":
                ptref[0] = (ptref[0] - 1) % len(self.lines)
                print(ptref[0], " ", ptref[1])
                pt_change()
            if event.key == " ":
                ptref[1] = -1 - ptref[1]
                print(ptref[0], " ", ptref[1])
                pt_change()
            if event.key == "&":
                ptref[1] = max(ptref[1] - 1, -len(self.lines[ptref[0]]))
                print(ptref[0], " ", ptref[1])
                pt_change()
            if event.key == "é":
                ptref[1] = min(ptref[1] + 1, len(self.lines[ptref[0]])-1)
                print(ptref[0], " ", ptref[1])
                pt_change()                
            if event.key == "x":
                zoom_point()
            if event.key == "!":
                theta = ax.azim * np.pi / 180
                phi = ax.elev * np.pi / 180
                self.set_axis([-np.sin(theta), np.cos(theta), 0], [-np.cos(theta)*np.sin(phi), -np.sin(theta)*np.sin(phi), np.cos(phi)])
                self.traitement()
                ax.clear()
                ax.set(proj_type="ortho")
                ax.view_init(elev=ax.elev, azim=ax.azim)
                ax.set_axis_off()
                self.plot_lines(ax)
            if event.key == "d":
                try:
                    self.plot_domain()
                except Exception as e:
                    print("plot_domain error:", e)
            plt.show()

        # plt.rcParams['keymap.left'].remove('left')
        cid = fig.canvas.mpl_connect("key_press_event", onkey)

        self.plot_lines(ax)


        limits = np.array([ax.get_xlim(), ax.get_ylim(), ax.get_zlim()])
        ax.set_box_aspect(aspect=np.ptp(limits, axis=1))
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "préparation du dessin"))
        self.print("[%0.3fs] %s" % (time.perf_counter() - T0, "Total"))
        plt.show()

    def plot_domain(self):
        from matplotlib.collections import LineCollection
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_aspect('equal')
        ax.set_title('Domaine (u,v)  —  rouge=contour  bleu=bord  gris=connexion  ★=cusp  ●=breakpoint')

        # Maillage en gris clair (LineCollection pour la performance)
        segs = np.stack([self.eds["fp"], self.eds["fp"] + self.eds["fdp"]], axis=1)
        ax.add_collection(LineCollection(segs, colors='lightgray', linewidths=0.4, zorder=1))

        # Lignes silhouette
        line_colors = {'b': 'steelblue', 'c': 'crimson', '?': 'dimgray'}
        for i, l in enumerate(self.lines):
            fp = l["fp"]
            t = l[0]["type"]
            col = line_colors.get(t[1] if len(t) > 1 else '?', 'dimgray')
            ax.plot(fp[:, 0], fp[:, 1], color=col, lw=1.5, zorder=3)

            # Indice de ligne + changements de vis aux extrémités
            mid = fp[len(fp) // 2]
            ax.annotate(f'{i} ({l[0]["v"]},{l[-1]["v"]})', mid,
                        fontsize=6, color=col, zorder=5,
                        bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.5, ec='none'))

            # Breakpoints
            for j in range(len(l) - 1):
                p, q = l[j]["fp"], l[j + 1]["fp"]
                for s, v in self.line_bks[i][j]:
                    bk = (1 - s) * p + s * q
                    ax.scatter(bk[0], bk[1], s=25, color='darkorange', zorder=4)
                    ax.annotate(f'{v:+d}', bk, fontsize=6, color='darkorange', zorder=5)

        # Cusps
        if hasattr(self, 'cusps') and len(self.cusps) > 0:
            ax.scatter(self.cusps["fp"][:, 0], self.cusps["fp"][:, 1],
                       s=80, color='magenta', marker='*', zorder=6, label='cusps')
            ax.legend(fontsize=8)

        # Extrémités des lignes (noeuds du graphe)
        for l in self.lines:
            for pt in [l[0]["fp"], l[-1]["fp"]]:
                ax.scatter(pt[0], pt[1], s=15, color='black', zorder=5)

        ax.set_xlabel('u'); ax.set_ylabel('v')
        plt.tight_layout()
        plt.show(block=False)

    def domain_data(self):
        # Données domaine (u,v) pour la vue de débogage dans le navigateur
        border_edges = []
        for e in self.beds:
            border_edges.append([e["fp"].tolist(), (e["fp"] + e["fdp"]).tolist()])

        lines = []
        for i, l in enumerate(self.lines):
            if l[0]["type"] == "??":
                continue
            pts = l["fp"].tolist()
            bks = []
            for j in range(len(l) - 1):
                p, q = l[j]["fp"], l[j + 1]["fp"]
                for s, v in self.line_bks[i][j]:
                    bks.append({"pt": ((1 - s) * p + s * q).tolist(), "v": int(v)})
            lines.append({
                "idx": i,
                "type": l[0]["type"],
                "v0": int(l[0]["v"]),
                "v1": int(l[-1]["v"]),
                "pts": pts,
                "bks": bks,
            })

        cusps = self.cusps["fp"].tolist() if hasattr(self, "cusps") and len(self.cusps) > 0 else []

        return {
            "bounds": [self.u_min, self.u_max, self.v_min, self.v_max],
            "res": self.res,
            "border_edges": border_edges,
            "lines": lines,
            "cusps": cusps,
        }

    def json_out(self):
        pass

    def relevement(self, fp):
        # relève la courbe fp
        line = np.array(fp)
        for i in range(1, len(line)):
            line[i] = self.close(line[i - 1], line[i])
        return line

    def dS(self, u, v, du, dv):
        return np.array(self.Su(u, v)) * du + np.array(self.Sv(u, v)) * dv

    def d2S(self, u, v, du, dv):
        return (
            np.array(self.Suu(u, v)) * du * du
            + np.array(self.Suv(u, v)) * du * dv * 2
            + np.array(self.Svv(u, v)) * dv * dv
        )

    def dirint(
        self, e, s=0.5
    ):  # Pour une arête de bord, normale à S dans le plan de vision, dirigée vers S, calculée au point de coord bary s sur l'arête
        nu = e["dir"]
        fp = e["fp"] + s * e["fdp"]
        vec = self.XY(self.dS(*fp, *nu))
        t = self.XY(self.dS(*fp, *e["fdp"]))
        t = t / norm(t)
        vec = vec - np.inner(vec, t) * t
        return vec / norm(vec)

    def dirvec(self, p):  # direction de la surface dans l'image au point de contour p
        return self.kerdS(*p)[1]

    def bvis_chge(self, e, s):  # e est une arête de bord
        p, dp, Nb = e["fp"], e["fdp"], e["dir"]  # segment dans le domaine

        p0 = p + s * dp  # le point où la visibilité peut changer

        ker = self.kerdS(*p0)[0]
        # on le dirige vers le pli supérieur, càd les z croissants
        if self.Z(self.dS(*p0, *ker)) < 0:
            ker = -ker

        if np.inner(ker, Nb) <= 0:
            return 0  # silhouette sortant du domaine, pas de changement de visibilité

        Tb = dp / norm(dp)  # tangent au bord

        Np = np.array(self.SD_jac(*p0, *self.ax_param))  # normale au pli
        if np.inner(Np, ker) < 0:
            Np = -Np  # dirigée vers le pli supérieur

        inner = np.inner(Tb, Np)
        if abs(inner) > 0.05:  # cas non dégénéré : formule standard
            return 1 if inner > 0 else -1

        # Cas dégénéré (cusp proche du bord) : ker ≈ Tb donc Np ≈ ±Nb et Tb·Np ≈ 0.
        # On évalue SN_axis en un point légèrement décalé le long de Tb pour capturer
        # les termes d'ordre supérieur qui déterminent le signe.
        eps = (self.u_max - self.u_min + self.v_max - self.v_min) / self.res
        sn = self.SN_axis(*(p0 + eps * Tb), *self.ax_param)
        self.print("Warning: near-degenerate border-contour crossing at [%0.4f,%0.4f] (Tb·Np=%.4f), using SN_axis offset (sn=%.4f)" % (p0[0], p0[1], inner, sn))
        return 1 if sn > 0 else -1

    def vis_chge(
        self, p, dir
    ):  # changement de visi quand part d'un pli en p dans la direction dir
        ker = self.kerdS(*p)[0]
        # on le dirige vers le pli supérieur, càd les z croissants
        if np.inner(self.dS(*p, *ker), self.axis) < 0:
            ker = -ker
        Np = np.array(self.SD_jac(*p, *self.ax_param))  # normale au pli
        if np.inner(Np, ker) < 0:
            Np = -Np  # dirigée vers le pli supérieur


        return 0 if np.inner(Np, dir) >= 0 else -1

    def kerdS(self, u, v):  # u et v peuvent être des float ou des vecteurs de floats
        ker = self.ker_param(u, v)  # kernel of projection Jacobian (ortho or perspective)

        diff2 = self.proj_vec(u, v, self.d2S(u, v, ker[0], ker[1]))
        im = self.proj_vec(u, v, self.dS(u, v, -ker[1], ker[0]))
        diff2 = diff2 - np.sum(diff2 * im, axis=0) * im / np.sum(im * im, axis=0)
        dir_vec = diff2 / norm(diff2, axis=0)

        # ker = noyau de (dS projeté sur le plan de vision), dir_vec = dérivée seconde dans la direction du noyau, projetée sur la coimage
        return ker, dir_vec

    def edge_intersect(
        self, e, f
    ):  # calcule les coord barycentriques de l'intersection sur chacune des arêtes, renseigne self.breaks

        p, dp, zp, dzp, q, dq, zq, dzq = (
            e["p"],
            e["dp"],
            e["z"],
            e["dz"],
            f["p"],
            f["dp"],
            f["z"],
            f["dz"],
        )
        D = np.cross(dp, dq)
        u = p - q
        if D != 0:
            s = np.cross(dq, u) / D
            t = (
                np.cross(dp, u) / D
            )  # coord barycentriques de l'intersection sur chacune des arêtes
            if s > 0 and s < 1 and t > 0 and t < 1:
                ze_lin = zp + s * dzp
                zf_lin = zq + t * dzq
                z_tol = 1e-4 * max(abs(ze_lin), abs(zf_lin), abs(dzp), abs(dzq), 1e-10)
                if abs(ze_lin - zf_lin) <= z_tol:
                    # Z-difference too small for reliable comparison — evaluate actual surface Z
                    ze = float(self.Z(np.array(self.S(*(e["fp"] + s * e["fdp"])))))
                    zf = float(self.Z(np.array(self.S(*(f["fp"] + t * f["fdp"])))))
                else:
                    ze, zf = ze_lin, zf_lin
                if ze > zf:
                    v = -e["n"] if np.inner(dq, e["d"]) > 0 else e["n"]
                    if v != 0:  # v= 0 pour les arêtes fictives
                        self.add_line_bk(
                            f["l_idx"], f["e_idx"], t, v
                        )  # e est audessus de f
                elif ze < zf:
                    v = -f["n"] if np.inner(dp, f["d"]) > 0 else f["n"]
                    if v != 0:
                        self.add_line_bk(
                            e["l_idx"], e["e_idx"], s, v
                        )  # f est audessus de e
                return True
            else:
                return False

    def add_line_bk(
        self, l_idx, e_idx, s, v
    ):  # l_idx: indice de la ligne dans self.lines, e_idx: indice du point p de l'arête dans la ligne, s: coord bary, v:chgt de visi
        bisect.insort(self.line_bks[l_idx][e_idx], (s, v))

    def propagation(self, l_idx, idx, visib):
        # propage la visibilité sur une ligne l_idx à partir du point idx de visibilité visib
        # dans le sens + si idx est <> -1, dans le sens - sinon. Renvoie la visibilité à l'extrêmité
        # attention, visib est la visibilité quand on est déjà dans la ligne
        l = self.lines[l_idx]
        vis = 0
        if idx == 0 :
            for bks in self.line_bks[l_idx][
                idx:
            ]:  # attention, idx peut-être négatif, y faut pas utiliser range
                for s, v in bks:
                    vis = vis + v
            # on soustrait le changement de visi à l'extrêmité
            return visib + vis - l[-1]["v"]
        else:
            for bks in self.line_bks[l_idx][0:idx]:
                for s, v in bks:
                    vis = vis + v
            return visib - vis - l[0]["v"]

    def connected_components(self):
        dic = defaultdict(list)
        for i, l in enumerate(self.lines):
            dic[tuple(l[0]["fp"])].append(i)
            dic[tuple(l[-1]["fp"])].append(i - len(self.lines))

        components = []

        lines = set(range(len(self.lines)))
        while lines:
            l = lines.pop()
            component = set([l])
            nodes = set(
                [tuple(self.lines[l][0]["ixfp"]), tuple(self.lines[l][-1]["ixfp"])]
            )
            already_seen = set()
            while nodes:
                node = nodes.pop()
                already_seen.add(node)
                for i in dic[node]:
                    iplus = i if i >= 0 else i + len(self.lines)
                    lines -= set([iplus])
                    component.add(iplus)
                    j = 0 if i >= 0 else -1  # node est au début ou à la fin
                    pt = self.lines[i][-1 - j]  # l'autre extremité de la ligne
                    if tuple(pt["ixfp"]) not in already_seen:
                        nodes.add(tuple(pt["ixfp"]))
            components.append(list(component))
        return components

    def connection(self, p, q):
        delta = min(
            (self.u_max - self.u_min) / self.res, (self.v_max - self.v_min) / self.res
        )
        #n = max(int(norm(q["ixfp"] - p["ixfp"]) / delta),2)
        n = 10 # a priori ne change pas grand chose puisque simplify_lines le redécoupe
        l = np.empty((n,), dtype=line_pt_type)
        l["fp"] = l["ixfp"] = np.column_stack(
            (
                np.linspace(p["ixfp"][0], q["ixfp"][0], n),
                np.linspace(p["ixfp"][1], q["ixfp"][1], n),
            )
        )
        l[0], l[-1] = p, q # p.copy(), q.copy()  # le copy est il nécessaire??
        if p["type"][0] == "v" or q["type"][0] == "v":
            print("pappa mia")
        if p["type"][0] == "c":
            l[0]["v"] = self.vis_chge(p["ixfp"], q["ixfp"] - p["ixfp"])
        if q["type"][0] == "c":
            l[-1]["v"] = self.vis_chge(q["ixfp"], p["ixfp"] - q["ixfp"])
        l["type"] = "??"
        return l

    def plot_lines(self, ax): # utilisé uniquement dans 'plot_for_terminal'
        bks = getattr(self, 'opt_line_bks', self.line_bks)
        opt_cp = getattr(self, 'opt_cp', None)
        visible_lines = []  # tracées en dernier
        for i, l in enumerate(self.lines):
            type = "?"
            # if l[0]["type"] == "??":
            #     continue
            cp = float(opt_cp[i]) if opt_cp is not None else float(l[0]["v"])
            vis = cp
            if tuple(l[0]["ixfp"]) in self.visibilities:
                type = l[0]["type"][1]
                vis += self.visibilities[tuple(l[0]["ixfp"])]
            line = []
            for j in range(len(l) - 1):
                p, q = l[j]["fp"], l[j + 1]["fp"]
                line.append(p)
                for s, v in bks[i][j]:
                    line.append((1 - s) * p + s * q)
                    if vis == 0:
                        visible_lines.append(line)
                    else:
                        self.plot_line(
                            ax, line, vis
                        )
                    vis = vis + v
                    line = [line[-1]]
            line.append(l[-1]["fp"])
            if vis == 0:
                visible_lines.append(line)
            else:
                self.plot_line(ax, line, vis)
        for line in visible_lines:
            self.plot_line(ax, line, 0)

    def plot_line(self, ax, line, vis):
        l = np.array(line)
        l = self.S(l[:, 0], l[:, 1])
        if vis == 0:
            ax.plot(l[0], l[1], l[2], color="black")
        elif vis != 0:

            # col = "black"
            col = colormap[min(14, max(vis+8, 0))]
            # col = 'powderblue' if type == 'b' else 'thistle' if type =='c' else 'navajowhite'
            # col = (1+.5/min(vis, -1),1+.5/min(vis, -1),1+.5/min(vis, -1)) # le min evite une erreur si vis > 0

            linestyle = 'solid'
            #linestyle = linestyles[min(4, max(-vis, 0))]
            # linestyle = (0,(2,1))

            ax.plot(l[0], l[1], l[2], linestyle=linestyle, color=col)

# T0 = time.perf_counter()

