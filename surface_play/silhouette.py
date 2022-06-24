# Todo:
# Corriger le problème su simplify lines aux points de cusps qui a l'air de déconnecter la ligne.

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
import rtree as rt

from numpy.linalg import norm
from numpy.core.records import fromrecords
from numpy.core.records import fromarrays
from numpy.linalg import svd
from scipy.optimize import newton

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
        ("p", "uint", 2),
        ("q", "uint", 2),
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

def generator_function(data):
    for i, obj in enumerate(data):
        yield (i, (obj[0], obj[1], obj[2], obj[3]), 42)

class Surface:
    def print(self, s):
        # pass
        print(s)

    def __init__(self, X="x", Y="y", Z="x*y", param_names="x y", bounds = (0,1,0,1) , quotient=("no", "no")):
        t0 = time.perf_counter()
        self.low_res = 80

        self.u_min, self.u_max, self.v_min, self.v_max = (
            bounds[0],
            bounds[1],
            bounds[2],
            bounds[3],
        )
        self.u_identify = quotient[0]
        self.v_identify = quotient[1]

        u, v = sp.symbols(param_names)

        # On détermine la bounding box (avec résolution low_res) pour perturber S, ça sert aussi dans for_3js
        with sp.evaluate(False) :
            truc_nul = u*sp.sin(sp.pi)
        S_pur = sp.Array(sp.sympify([X, Y, Z])) + sp.Array([truc_nul,truc_nul,truc_nul])
        Slambda = lambdify([u, v], S_pur, "numpy")
        u_grid, v_grid = np.meshgrid(
            np.linspace(self.u_min, self.u_max, self.low_res+1), 
            np.linspace(self.v_min, self.v_max, self.low_res+1),
            indexing="ij"
        )
        self.X_grid, self.Y_grid, self.Z_grid = Slambda(u_grid, v_grid)
        self.bbox = (
            np.amin(self.X_grid), 
            np.amax(self.X_grid), 
            np.amin(self.Y_grid), 
            np.amax(self.Y_grid), 
            np.amin(self.Z_grid), 
            np.amax(self.Z_grid)
        )
        self.center = [(self.bbox[0]+self.bbox[1])/2, (self.bbox[2]+self.bbox[3])/2, (self.bbox[4]+self.bbox[5])/2]
        dX, dY, dZ = self.bbox[1] - self.bbox[0], self.bbox[3] - self.bbox[2], self.bbox[5] - self.bbox[4]
        self.radius = np.sqrt(dX**2 + dY**2 + dZ**2)/2

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
                0.005 * dX * .346,  
                0.005 * dY * .632,  
                0.005 * dZ * .693,  
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
        SN_axis_sp = dot(SN_sp, ax)  # produit scalaire du précédent avec un axe
        SD_jac_sp = sp.Array(
            [
                dot(
                    cross(Suu_sp, Sv_sp) + cross(Su_sp, Suv_sp), ax
                ),  # grad du jacobien, orth aux contours
                dot(cross(Suv_sp, Sv_sp) + cross(Su_sp, Svv_sp), ax),
            ]
        )

        self.S = lambdify([u, v], S_sp, "numpy")# applied to a numpy array, creates a list of three np.arrays of the same dimension
        self.Su = lambdify([u, v], Su_sp, "numpy")
        self.Sv = lambdify([u, v], Sv_sp, "numpy")
        self.Suu = lambdify([u, v], Suu_sp, "numpy")
        self.Suv = lambdify([u, v], Suv_sp, "numpy")
        self.Svv = lambdify([u, v], Svv_sp, "numpy")
        self.SN = lambdify([u, v], SN_sp, "numpy")
        self.SN_axis = lambdify([u, v, ax_x, ax_y, ax_z], SN_axis_sp, "numpy")
        self.SD_jac = lambdify([u, v, ax_x, ax_y, ax_z], SD_jac_sp, "numpy")

        # finalement, on fabrique une grille de normales pour for_3js
        self.NX_grid, self.NY_grid, self.NZ_grid = self.SN(u_grid, v_grid)

        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "init surface"))

    def set_axis(
        self, I = [1,0,0], J = [0,0,1]
    ):  # definit l'axe et les fonctions self.XY et self.Z
        vecI, vecJ = np.array(I), np.array(J)
        vecI, vecJ = vecI/norm(vecI), vecJ/norm(vecJ)
        self.axis = np.cross(vecI, vecJ)  
        self.axis = self.axis/norm(self.axis) # coordonnées x,y,z de la normale au plan de vision
        def XY(vec):  # renvoie coordonnées d'un point de l'espace dans le repère I,J
            return np.array([np.inner(vecI, vec.T), np.inner(vecJ, vec.T)]) # le .T est pour le cas d'appel vectoriel. Pour 1 vecteur, ça change rien

        def Z(vec):  # renvoie coordonnées d'un point de l'espace sur l'axe 'axis'
            return np.inner(self.axis, vec.T)

        def XYZ(vec):  # renvoie le vecteur dont les coord. sur (I,J) sont 'vec'
            return vec[0] * vecI + vec[1] * vecJ

        self.XY = XY
        self.Z = Z
        self.XYZ = XYZ

    def for_3js(self) :
        res = self.low_res
        # liste 1d  x_1,y_1,z_1,x_2,y_2,z_2,... des positions/normales. l'indice du point (i,j) de la grille est i*(res+1) + j
        positions = np.transpose(np.array([self.X_grid.flatten(), self.Y_grid.flatten(), self.Z_grid.flatten()])).flatten() 
        normals = np.transpose(np.array([self.NX_grid.flatten(), self.NY_grid.flatten(), self.NZ_grid.flatten()])).flatten()
        # index(i,j) = indice du point p(i,j)
        index = np.reshape(np.arange((res+1)*(res+1)), (res+1, res+1))
        a = index[:-1,:-1]
        b = index[:-1,1:]
        c = index[1:,:-1]
        d = index[1:,1:]
        faces = np.concatenate((np.transpose(np.array([a,c,b])).flatten(), np.transpose(np.array([c,d,b])).flatten()))
        
        return positions.tolist(), normals.tolist(), faces.tolist(), self.center, self.radius    

    def triangulate(
        self, res
    ):  
        t0 = time.perf_counter()
        self.res = res
        
        d_u, d_v = self.u_max - self.u_min, self.v_max - self.v_min

        u_list = np.linspace(self.u_min, self.u_max, res + 1)  # res = number of squares
        v_list = np.linspace(self.v_min, self.v_max, res + 1)

        self.uv_vector = 5 * np.array(
            [res / (self.u_max - self.u_min), res / (self.v_max - self.v_min)]
        )  # multiplié par la longueur, donne le nombre de subdivision. Grand = plus de subdivisions.

        u_grid, v_grid = np.meshgrid(
            u_list, v_list, indexing="ij"
        )  # matrix indexing: u est l'indice de ligne, v l'indice de colonne
        self.S_grid = self.S(u_grid, v_grid)
        self.SN_grid = self.SN(u_grid, v_grid)
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "calcul grid"))

        t0 = time.perf_counter()

        ############ listes d'indices
        points = np.array(np.meshgrid(np.arange(res + 1), np.arange(res + 1))).T
        fpoints = np.array(
            np.meshgrid(u_list, v_list)
        ).T  # transposé fait que u = ligne, v = colonne, et que le dernier axe = point

        if self.u_identify == "cy":
            points[res, :] = points[0, :]
        if self.v_identify == "cy":
            points[:, res] = points[:, 0]
        if self.u_identify == "mo":
            points[res, :] = points[0, ::-1]  # reverse order
        if self.v_identify == "mo":
            points[:, res] = points[::-1, 0]

        # les sommets du carré: b au dessus de a, d au dessus de c, c à droite de a
        a = points[:res, :res]
        b = points[:res, 1 : res + 1]
        c = points[1 : res + 1, :res]
        d = points[1 : res + 1, 1 : res + 1]

        # en coordonnées u,v
        fa = fpoints[:res, :res]
        fb = fpoints[:res, 1 : res + 1]
        fc = fpoints[1 : res + 1, :res]
        fd = fpoints[1 : res + 1, 1 : res + 1]
        # on ne fait pas d'identifications pour les raccord cyl et moebius, sinon les arêtes ont des extremités trop éloignées.

        # Les faces : devenu inutile depuis qu'on utilise plus "other", la face est juste un numéro
        # up_faces = fromarrays((a,b,d), dtype = f_type).flatten()
        # dn_faces = fromarrays((a,d,c), dtype = f_type).flatten()
        # self.faces = np.array(np.concatenate((up_faces, dn_faces)), dtype = f_type)

        # les indexs des faces (down et up), dans un tableau indicé par les coord. du point bas gauche du carré
        f = np.reshape(np.arange(res * res), (res, res))  # f = up faces
        g = np.reshape(np.arange(res * res, 2 * res * res), (res, res))  # g = dn faces
        o = np.full(
            (res,), -1
        )  # l'indice -1 signifie qu'il n'y a pas de face (pour le bord)

        # les directions (n'ont un sens que pour les arêtes de bord)
        drte = np.array([1, 0])
        gche = np.array([-1, 0])
        haut = np.array([0, 1])
        bas = np.array([0, -1])
        # haut = direction vers le haut. C'est la normale aux arêtes en bas. Idem pour les autres

        # le champ 'flip'
        flip_square_no = np.full((res, res), 1)
        flip_line_no = np.full((res,), 1)

        # arêtes de bord
        dn_eds = fromarrays(
            (
                a[:, 0],
                c[:, 0],
                g[:, 0],
                o,
                fa[:, 0],
                fc[:, 0] - fa[:, 0],
                np.tile(haut, (res, 1)),
                flip_line_no,
            ),
            dtype=e_type,
        )  # left
        up_eds = fromarrays(
            (
                b[:, -1],
                d[:, -1],
                f[:, -1],
                o,
                fb[:, -1],
                fd[:, -1] - fb[:, -1],
                np.tile(bas, (res, 1)),
                flip_line_no,
            ),
            dtype=e_type,
        )  # right
        lf_eds = fromarrays(
            (
                a[0, :],
                b[0, :],
                f[0, :],
                o,
                fa[0, :],
                fb[0, :] - fa[0, :],
                np.tile(drte, (res, 1)),
                flip_line_no,
            ),
            dtype=e_type,
        )  # up
        rt_eds = fromarrays(
            (
                c[-1, :],
                d[-1, :],
                g[-1, :],
                o,
                fc[-1, :],
                fd[-1, :] - fc[-1, :],
                np.tile(gche, (res, 1)),
                flip_line_no,
            ),
            dtype=e_type,
        )  # down
        # arêtes intérieures
        dg_eds = fromarrays(
            (a, d, f, g, fa, fd - fa, np.tile(drte, (res, res, 1)), flip_square_no),
            dtype=e_type,
        )  # diagonales
        hr_eds = fromarrays(
            (
                a[:, 1:],
                c[:, 1:],
                g[:, 1:],
                f[:, :-1],
                fa[:, 1:],
                fc[:, 1:] - fa[:, 1:],
                np.tile(drte, (res, res - 1, 1)),
                flip_square_no[:, 1:],
            ),
            dtype=e_type,
        )  # horiz
        vr_eds = fromarrays(
            (
                a[1:, :],
                b[1:, :],
                f[1:, :],
                g[:-1, :],
                fa[1:, :],
                fb[1:, :] - fa[1:, :],
                np.tile(drte, (res - 1, res, 1)),
                flip_square_no[1:, :],
            ),
            dtype=e_type,
        )  # vertic
        if self.u_identify == "cy":
            lf_eds["g"] = g[-1, :]
        if self.v_identify == "cy":
            dn_eds["g"] = f[:, -1]
        if self.u_identify == "mo":
            lf_eds["g"] = np.flip(g[-1, :], 0)
            dg_eds["flip"][-1, :] = -1
            hr_eds["flip"][-1, :] = -1
            dn_eds["flip"][-1] = up_eds["flip"][
                -1
            ] = (
                -1
            )  # indique que la normale calculée avec NZ grid donne une orientation incohérente aux extremites de l'arete
        if self.v_identify == "mo":
            dg_eds["flip"][:, -1] = -1
            dn_eds["g"] = np.flip(f[:, -1], 0)
            vr_eds["flip"][:, -1] = -1
            lf_eds["flip"][-1] = rt_eds["flip"][-1] = -1

        # tous ensemble
        if (self.u_identify == "cy" or self.u_identify == "mo") and self.v_identify == "no":
            self.eds = np.array(
                np.concatenate(
                    (
                        dg_eds.flatten(),
                        hr_eds.flatten(),
                        vr_eds.flatten(),
                        lf_eds,
                        dn_eds,
                        up_eds,
                    )
                ),
                dtype=e_type,
            )
            self.b_index = (
                res * res + 2 * (res - 1) * res + res
            )  # diagonales + 2* verticales intérieures + verticale gauche
            self.ieds = self.eds[
                : self.b_index
            ]  # a view of the edges containing only the interior edges
            self.beds = self.eds[
                self.b_index :
            ]  # a view of the edges containing only the boundary edges
        elif self.u_identify == "no" and (self.v_identify == "cy" or self.v_identify == "mo"):
            self.eds = np.array(
                np.concatenate(
                    (
                        dg_eds.flatten(),
                        hr_eds.flatten(),
                        vr_eds.flatten(),
                        dn_eds,
                        lf_eds,
                        rt_eds,
                    )
                ),
                dtype=e_type,
            )
            self.b_index = res * res + 2 * (res - 1) * res + res
            self.ieds = self.eds[
                : self.b_index
            ]  # a view of the edges containing only the interior edges
            self.beds = self.eds[
                self.b_index :
            ]  # a view of the edges containing only the boundary edges
        elif self.u_identify != "no" and self.v_identify != "no":
            self.eds = np.array(
                np.concatenate(
                    (
                        dg_eds.flatten(),
                        hr_eds.flatten(),
                        vr_eds.flatten(),
                        lf_eds,
                        dn_eds,
                    )
                ),
                dtype=e_type,
            )
            self.b_index = res * res + 2 * (res - 1) * res + 2 * res
            self.ieds = self.eds[
                : self.b_index
            ]  # a view of the edges containing only the interior edges
            self.beds = self.eds[
                self.b_index :
            ]  # a view of the edges containing only the boundary edges
        else:
            self.eds = np.array(
                np.concatenate(
                    (
                        dg_eds.flatten(),
                        hr_eds.flatten(),
                        vr_eds.flatten(),
                        lf_eds,
                        rt_eds,
                        up_eds,
                        dn_eds,
                    )
                ),
                dtype=e_type,
            )
            self.b_index = res * res + 2 * (res - 1) * res
            self.ieds = self.eds[
                : self.b_index
            ]  # a view of the edges containing only the interior edges
            self.beds = self.eds[
                self.b_index :
            ]  # a view of the edges containing only the boundary edges
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "triangulation"))


        # les coins. Points transformés en tuples pour être hashables
        if self.u_identify == self.v_identify == "no":
            self.corners = set(
                [
                    tuple(points[0, 0]),
                    tuple(points[0, res]),
                    tuple(points[res, 0]),
                    tuple(points[res, res]),
                ]
            )
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

    def traitement(self):
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
        self.simplify_lines()  # espace les points régulièrement (à l'écran) pour éviter les artefacts de discrétisation.
        self.line_bks = [[[] for j in range(len(l) - 1)] for l in self.lines]
        line, idx = self.intersections() # calcule les intersections/chgement de visi. Renvoie: point visible (long)
        # line, idx = 0, 0
        for i, l in enumerate(self.lines): #debug
           self.print("ligne:%i type:%s len:%i debut:[%0.2f,%0.2f] fin:[%0.2f,%0.2f]  chgt visi début:%i chgt visi fin:%i"%(i,l[0]['type'], 
                                                                                                            len(l), *self.XY(np.array(self.S(*l['fp'][0]))), 
                                                                                                            *self.XY(np.array(self.S(*l['fp'][-1]))), l['v'][0],  l['v'][-1]))
           string = 'breaks : ['
           for j in range(len(l) - 1):
               for (s,v) in self.line_bks[i][j]:
                   string+= "(%i, %i) "%(j,v)
           self.print(string+']')
        self.visibilities = {}
        # self.visibility(line, idx)  # propage la visibilité à toutes les lignes
        self.visibilite(line, idx)  # propage la visibilité à toutes les lignes
        #for p, v in self.visibilities.items(): # debug
        #    print("point:[%0.2f,%0.2f] visibilité:%i"%(*self.XY(np.array(self.S(*p))), v))

        self.origine = self.lines[line][idx]['fp'] #debug

    def find_silhouette(self):
        t0 = time.perf_counter()

        ######## Calcul des points de contour #########

        #### sélection des arêtes
        NZ_grid = (
            self.SN_grid[0] * self.axis[0]
            + self.SN_grid[1] * self.axis[1]
            + self.SN_grid[2] * self.axis[2]
        )

        vec_i = np.where(
            NZ_grid[self.eds["p"][:, 0], self.eds["p"][:, 1]]
            * NZ_grid[self.eds["q"][:, 0], self.eds["q"][:, 1]]
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
            return self.SN_axis(u + s * du, v + s * dv, *self.axis)

        init = np.full(len(vec_i), 0.5)
        if len(vec_i) > 0:
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
        dir_vec = self.kerdS(u_vec, v_vec)[1]
        # print(dir_vec)

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
            directions = dir_vec[:, l]
            signes = np.sum(directions[:, 0:-1] * directions[:, 1:], axis=0)
            indices = np.where(signes < 0)[0]
            for i in indices:
                p, q = np.array([u_vec[l[i]], v_vec[l[i]]]), (
                    [u_vec[l[i + 1]], v_vec[l[i + 1]]]
                )
                q = self.close(p, q)
                vis = (
                    1 if np.inner(self.dS(*p, *(q - p)), self.axis) > 0 else -1
                )  # ça monte dans le sens p, q. donc la visibilité augmente
                r = (p + q) / 2
                # nu = np.array([-(q-p)[1], (q-p)[0]])
                # for j in range(20): # bisection pour trouver le point précis, avec newton pour rester sur le contour
                #    a, b = (p+q)/2 - nu/2, (p+q)/2 + nu/2
                #    s = newton(NZbis, 0.5, args=(a[0], a[1], nu[0], nu[1]))
                #    r = a + s * (b-a)
                #    if np.inner(self.dirvec(p), self.dirvec(r)) < 0:
                #        q = r
                #    elif np.inner(self.dirvec(q), self.dirvec(r)) < 0:
                #        p = r
                #    else:
                #        print('problème de cusp!!!')
                cusps.append((l[i], l[i + 1], 1 / 2, r))
                # print(vis)
                self.breaks["c"][l[i]] = (0.5, vis, len(cusps) - 1)

        self.cusps = fromrecords(cusps, dtype=cusp_type)
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "breaks de cusps"))
        # print(len(self.c_lines), len(self.b_lines), len(self.cusps))
        # print(len(self.breaks['b']),len(self.breaks['c']))

    def break_lines(self):
        t0 = time.perf_counter()

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
                fp = e["fp"] if i >= 0 else e["fp"] + e.fdp
                fq = e["fp"] if i < 0 else e["fp"] + e.fdp

                lines[-1].append(
                    ("bb", e_idx + self.b_index, 0, self.dirint(e, sens), fp, fp)
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
                if tuple(q) in self.corners:
                    # en cas de coin, on casse, car la direction intérieure est discontinue aux coins.
                    lines[-1].append(
                        (
                            "bb",
                            e_idx + self.b_index,
                            0,
                            self.dirint(e, 1 - sens),
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
            lines = [[]]  # chaque b_line est découpée

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
            #i, j = 0, len(l)
            #if l["type"][0][0] == "v":
            #    i = 1
            #if l["type"][-1][0] == "v":
            #    j = len(l) - 1
            i, j = 0, len(l)
            if l["type"][0][0] == "v":
                i = min(5,len(l)//3) # on évite les extrêmités des lignes, dans la mesure du possible. Ici: marge de 10, à trois ça déconne parfois
            if l["type"][-1][0] == "v":
                j = len(l)-min(5, len(l)//3)
            return np.arange(i, j), l[i:j]

        self.cc_s = (
            self.connected_components()
        )  # liste de listes d'indices de lignes, qui sont connectées.
        # print('nombre de composantes connexes %i'%(len(self.cc_s),))

        if len(self.cc_s) < 2:
            return
        # print('connect !!')
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
            # e_index[i] = np.concatenate([np.arange(len(self.lines[l])) for l in cc])
            e_index[i] = np.concatenate([trim(self.lines[l])[0] for l in cc])
            bb = np.concatenate(
                (
                    np.concatenate([trim(self.lines[l])[1]["ixfp"] for l in cc]),
                    np.concatenate([trim(self.lines[l])[1]["ixfp"] for l in cc]),
                ),
                axis=1,
            )
            trees[i] = rt.index.Index(
                generator_function(bb), properties=rt.index.Property()
            )
            for j in range(i):
                i_point, j_point = len(e_index[i]) // 2, len(e_index[j]) // 2
                for k in range(2):
                    j_point = list(
                        trees[j].nearest(
                            tuple(
                                self.lines[l_index[i][i_point]][e_index[i][i_point]][
                                    "ixfp"
                                ]
                            )
                        )
                    )[0]
                    i_point = list(
                        trees[i].nearest(
                            tuple(
                                self.lines[l_index[j][j_point]][e_index[j][j_point]][
                                    "ixfp"
                                ]
                            )
                        )
                    )[0]
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

    def simplify_lines(self):
        t0 = time.perf_counter()

        # define sample
        def sample_fn(x):
            #return x
            #return x + np.sqrt(x)
            #return np.sqrt(self.res) * x * x / (1 + np.sqrt(self.res) * x) # moins cher que celui du dessous mais moins bien
            return x * x / (1 + x)
            #return x*np.sqrt(x)/(1+np.sqrt(x)) # semble être le bonne fonction, sais pas pourquoi

        p = []  # list of lists of image points of the lines
        # pas clair ce qu'il faut. le 50 représente la longueur max d'une ligne (largeur écran = 1), 100 = nombre de points par unité de longueur
        # sample = sample_fn(np.linspace(0,50, 50*self.res))
        # sample = sample_fn(np.linspace(0,50, 50*min(self.res, 200))) # pas clair ce qu'il faut.
        sample = sample_fn(
            np.linspace(0, 50, 50 * max(self.res // 5, 100))
        )  # pas clair ce qu'il faut.

        # compute bounds of image self.im_bounds = ([xmin, ymin], [xmax, ymax])
        im_min = np.min(self.XY(np.array(self.S(self.u_min, self.v_min))))
        im_max = np.max(self.XY(np.array(self.S(self.u_min, self.v_min))))
        for l in self.lines:
            # p.append((self.S(l['fp'][:,0], l['fp'][:,1])[:,0]).T)
            p.append(self.XY(np.array(self.S(l["fp"][:, 0], l["fp"][:, 1]))).T)

            lmin, lmax = np.min(p[-1]), np.max(p[-1])
            im_min, im_max = min(lmin, im_min), max(lmax, im_max)

        # vector which scales XY coordinates to a square of sidelength 1
        scale = np.array([1, 1]) / (im_max - im_min)

        for i, l in enumerate(self.lines):
            origin, end = l[0], l[-1]
            l["fp"] = self.relevement(l["fp"])
            intervals = norm(scale * (p[i][1:] - p[i][:-1]), axis=1) # interval[i] = dist(p[i], p[i+1])
            # indices = np.where(intervals==0)[0]
            # if indices.size >0:
            #    print(len(l))
            #    print(l['fp'])
            #    print("l(i) %0.3f, %0.3f  l(i+1) %0.3f, %0.3f"%tuple(list(l['fp'][i])+list(l['fp'][i+1]))
            lengths = np.cumsum(intervals)  # lengths[i+1] - lengths[i] = interval[i+1], lengths[i] = dist(p[0], p[i+1])
            index = np.searchsorted(sample, lengths[-1] / 2) # sample[index - 1] < l/2 <= sample[index]
            points = np.concatenate(
                (sample[:index], lengths[-1] - sample[index - 1 :: -1]) 
            )[
                1:-1
            ]  # points[i] = sample[i+1]
            j = np.searchsorted(lengths, points) # lengths[j[k]-1] < points[k] <= lengths[j[k]]  
            # p[k] = ((length[j[k]] - points[k]) * p[j[k]] + (point[k] - lenght[j[k] - 1]) * p[j[k]])/interval[j[k]]
            fp_s = (
                (points - lengths[j] + intervals[j]) * l[j + 1]["fp"].T
                + (lengths[j] - points) * l[j]["fp"].T
            ) / intervals[j]
            dirint_s = (
                (points - lengths[j] + intervals[j]) * l[j + 1]["d"].T
                + (lengths[j] - points) * l[j]["d"].T
            ) / intervals[j]
            fp_s = fp_s.T
            dirint_s = (
                dirint_s / (0.1 + norm(dirint_s, axis=0))
            ).T  # .1 pour éviter div par 0
            line = np.empty((len(j) + 2,), dtype=line_pt_type)
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

        # bag_bis.sort(order = ('x')) # incompatible avec la suite, où les indices de ligne sont utilisés
        # j = np.searchsorted(bag_bis.x, np.maximum(bag_bis.p[:,0], bag_bis.p[:,0] + bag_bis.dp[:,0]))
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

        p = rt.index.Property()
        # p.leaf_capacity = 50
        # p.fill_factor = 0.7
        bb_im = np.concatenate(
            (
                np.minimum(bag_bis["p"], bag_bis["p"] + bag_bis["dp"]),
                np.maximum(bag_bis["p"], bag_bis["p"] + bag_bis["dp"]),
            ),
            axis=1,
        )
        index_im = rt.index.Index(generator_function(bb_im), properties=p)
        # bb_dom = np.concatenate((np.minimum(bag_bis['fp'], bag_bis['fp'] + bag_bis.fdp) - tol,
        #                         np.maximum(bag_bis['fp'], bag_bis['fp'] + bag_bis.fdp) + tol), axis = 1)
        # index_dom = rt.index.Index(generator_function(bb_dom), properties = p)

        for i, e in enumerate(bag_bis):
            box_im = [j for j in index_im.intersection(bb_im[i]) if j > i]
            # box_im = set([j for j in index_im.intersection(bb_im[i]) if j > i])
            # box_dom = set([j for j in index_dom.intersection(bb_dom[i]) if j > i])
            # emax, emin = np.maximum(e.p, e.p + e.dp), np.minimum(e.p, e.p + e.dp)
            e_p, e_q = e["fp"], e["fp"] + e["fdp"]
            efmax, efmin = np.maximum(e_p, e_q), np.minimum(e_p, e_q)
            # for j in box_im - box_dom:
            #    self.edge_inter_bis(e, bag_bis[j])
            for j in box_im:
                f = bag_bis[j]
                f_p = self.close(e_p, f["fp"])
                f_q = f_p + f["fdp"]
                if np.any(np.minimum(f_p, f_q) > tol + efmax) or np.any(
                    efmin > tol + np.maximum(f_p, f_q)
                ):
                    self.edge_intersect(e, f)

        # for i in range(len(bag_bis)):
        #    e = bag_bis[i]
        #    emax, emin = np.maximum(e.p, e.p + e.dp), np.minimum(e.p, e.p + e.dp)
        #    efmax, efmin = np.maximum(e['fp'], e['fp'] + e.fdp), np.minimum(e['fp'], e['fp'] + e.fdp)
        #    sub = bag_bis[i+1:j[i]]
        #    subq = sub.p + sub.dp
        #    subj = sub[
        #        np.all(np.minimum(sub.p, subq) < emax, axis = 1) &
        #        np.all(emin < np.maximum(sub.p, subq), axis = 1)
        #        ]
        #    for f in subj:
        #        if (np.any(np.minimum(f['fp'], f['fp'] + f.fdp) > tol + efmax) or
        #        np.any(efmin > tol + np.maximum(f['fp'], f['fp'] + f.fdp))):
        #            self.edge_inter_bis(e, f)
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "Intersections bis"))

        bag_bis.sort(order=("x"))
        e = bag_bis[0]  # arête où le min de x est atteint
        if e["x"] == e["p"][0]:  # c'est p qui a le min de x
            return e["l_idx"], e["e_idx"]
        else:  # sinon c'est q qui a le min de x.
            return e["l_idx"], e["e_idx"] + 1

    def visibility(self, line, idx): # remplacé par "visibilite"
        t0 = time.perf_counter()

        # calcul de la visibilité des extrêmités des lignes, sachant que
        # le point idx de la ligne line a la visibilité 0
        t0 = time.perf_counter()
        v = (
            self.lines[line][idx]["v"]
            if (idx == 0 or idx == len(self.lines[line]) - 1)
            else 0
        )
        # print(' visibilities  ', v)
        self.propagate(
            line, idx, v
        )  # inscrit dans visibilities la visibilité des extrêmités de line.

        dic = defaultdict(set)
        # dic[point] = tuple(indices des lignes dont le point est l'extrêmité, 0 si point est le début, -1 si point est la fin)
        for i, l in enumerate(self.lines):
            if i != line:  # on exclut la ligne initiale qui est déjà faite.
                dic[tuple(l[0]["ixfp"])].add((i, 0))
                dic[tuple(l[-1]["ixfp"])].add((i, -1))

        # print(self.visibilities)
        nodes = set([next(iter(self.visibilities.keys()))])
        already_seen = set()
        while nodes:
            node = nodes.pop()
            # print(dic[node])
            already_seen.add(node)
            for i, j in dic[node]:
                pt = self.lines[i][-1 - j]["ixfp"]  # l'autre extremité de la ligne
                if tuple(pt) not in already_seen:
                    # j = 0 if i>= 0 else -1 # node est au début ou à la fin
                    self.propagate(
                        i, j, self.visibilities[node] + self.lines[i][j]["v"]
                    )
                    nodes.add(tuple(pt))
                    dic[tuple(pt)].discard((i, j))
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "Visibilité"))

    def visibilite(self, line, idx):
        # calcul de la visibilité des extrêmités des lignes, sachant que
        # le point idx de la ligne line a la visibilité 0
        t0 = time.perf_counter()

        ####### Initialisation
        visited = set()
        unseen = set()
        lines_dic = defaultdict(set) # lignes par sommet, indice de ligne, et 0 si le sommet est l[0], -1 si sommet = l[-1]
        for i, l in enumerate(self.lines):
            p, q = tuple(l[0]["ixfp"]), tuple(l[-1]["ixfp"])
            lines_dic[p].add((i, 0))
            lines_dic[q].add((i, -1))
            unseen.update([p,q])

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

        ####### Propagation
        dic = defaultdict(set)
        while visited :
            pt = visited.pop()
            for i, j in lines_dic[pt] :
                l = self.lines[i]
                qt = tuple(l[-1 - j]["ixfp"])
                if qt in unseen :
                    self.visibilities[qt] = self.propagation(i, j, self.visibilities[pt] + self.lines[i][j]["v"])
                    if self.visibilities[qt]>0:
                        self.print("visibilité > 0 !!!") # debug
                    unseen.discard(qt)
                    visited.add(qt)
                    # lines_dic[qt].discard((i,-1-j))
                elif self.visibilities[qt] != self.propagation(i, j, self.visibilities[pt] + self.lines[i][j]["v"]):
                    self.print("inconsistency!!!!") # debug
            visited.discard(pt)
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "Visibilité"))

    def plot_for_browser(self):
        t0 = time.perf_counter()

        lines = defaultdict(list)  # clé = visibilité, valeur = liste de lignes. ligne = [x_1,y_1,x_2,y_2...]
        for i, l in enumerate(self.lines):
            type = "?"
            if l[0]["type"] == "??":
                continue
            vis = l[0]["v"]
            if tuple(l[0]["ixfp"]) in self.visibilities:
                type = l[0]["type"][1]
                vis += self.visibilities[tuple(l[0]["ixfp"])]
            line = []
            for j in range(len(l) - 1):
                p, q = l[j]["fp"], l[j + 1]["fp"]
                line.append(self.XY(np.array(self.S(*p))))
                for s, v in self.line_bks[i][j]:
                    line.append(self.XY(np.array(self.S(*((1 - s) * p + s * q)))))
                    lines[vis].append(line)
                    vis = vis + v
                    line = [line[-1]] # coordonnées x,y du dernier point
            last_pt = l[-1]["fp"]
            line.append(self.XY(np.array(self.S(*last_pt))))
            lines[vis].append(line)
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
                #ax.scatter(*self.S(*self.origine))# debug
                self.plot_lines(ax)
                #for i, l in enumerate(self.lines): # debug
                #    for j in range(len(l) - 1):
                #        for (s,v) in self.line_bks[i][j]:
                #            p, q = l[j]["fp"], l[j + 1]["fp"]
                #            ax.scatter(*self.S(*((1 - s) * p + s * q)))
            plt.show()

        # plt.rcParams['keymap.left'].remove('left')
        cid = fig.canvas.mpl_connect("key_press_event", onkey)

        self.plot_lines(ax)

        #visible_lines = []  # tracées en dernier

        #for i, l in enumerate(self.lines):
        #    type = "?"
        #    if l[0]["type"] == "??":
        #        continue
        #    # if l[0].type[0] == 'v':
        #    #    ax.scatter(*self.S(*l[0]['fp']), marker = 'x')
        #    # if l[-1].type[0] == 'v':
        #    #    ax.scatter(*self.S(*l[-1]['fp']), marker = 'x')
        #    vis = l[0]["v"]
        #    if tuple(l[0]["ixfp"]) in self.visibilities:
        #        type = l[0]["type"][1]
        #        vis += self.visibilities[tuple(l[0]["ixfp"])]

        #    # des flèches a chaque début et fin de ligne qui indiquent un éventuel changement de visibilité
        #    # if l[0]['v'] != 0:
        #    #    color = 'green' if l[0]['v'] > 0 else 'red'
        #    #    p, q = self.S(*l[0]['fp']), self.S(*l[1]['fp'])
        #    #    #print(p,q,norm(q-p))
        #    #    ax.quiver(*p, *(.02*(q - p)/norm(q-p)),  color=color)
        #    # if l[-1]['v'] != 0:
        #    #    color = 'green' if l[0]['v'] > 0 else 'red'
        #    #    p, q = self.S(*l[-1]['fp']), self.S(*l[-2]['fp'])
        #    #    #print(p,q,norm(q-p))
        #    #    ax.quiver(*p, *(.02*(q - p)/norm(q-p)),  color=color)

        #    line = []
        #    for j in range(len(l) - 1):
        #        p, q = l[j]["fp"], l[j + 1]["fp"]
        #        line.append(p)
        #        for s, v in self.line_bks[i][j]:
        #            line.append((1 - s) * p + s * q)
        #            # ax.scatter(*self.S(*line[-1]))
        #            if vis == 0:
        #                visible_lines.append(line)
        #            else:
        #                # if vis > 0:
        #                #    plt.close(fig)
        #                #    return 'fail'
        #                plot_uv_line(
        #                    ax, type, line, vis
        #                )  # , 0 if type !='?' else 'navajowhite')
        #            vis = vis + v
        #            line = [line[-1]]
        #    line.append(l[-1]["fp"])
        #    if vis == 0:
        #        visible_lines.append(line)
        #    else:
        #        # if vis > 0:
        #        #    plt.close(fig)
        #        #    return 'fail'
        #        plot_uv_line(ax, type, line, vis)
        #for line in visible_lines:
        #    plot_uv_line(ax, type, line, 0)

        # ax.quiver(*self.S(self.u_grid[0,0], self.v_grid[0,0]), *(axis*10), color='black') # flèche vers l'observateur

        # ax.plot_surface(self.S_grid[0,0,:, :], self.S_grid[1,0,:, :], self.S_grid[2,0,:, :]) # la surface en facettes

        # for i, p in enumerate(self.c_pts): # des flèches vers la surface
        #    ax.quiver(*self.S(*p['fp'])[:,0], *(.03 * self.XYZ(p.d)), color='black')
        # for i, e in enumerate(self.beds):
        #    ax.quiver(*self.S(*e['fp'])[:,0], *(.02 * self.XYZ(self.dirint(e))),  color='black')

        # for l in self.b_lines: # des gros points là où il y a changement de visibilité au bord
        #    for i, idx in enumerate(l):
        #        index = idx if idx >= 0 else len(self.beds) + idx
        #        if index in self.breaks['b']:
        #            for s, v in self.breaks['b'][index]:
        #                point = (1 - s) * self.dom_pt(self.beds[idx].p) + s * self.dom_pt(self.beds[idx].q)
        #                ax.scatter(*self.S(*point))

        # for l in self.c_lines: # des gros points là où il y a changement de visibilité sur les contours
        #    for i in range(len(l)-1):
        #        if l[i] in self.breaks['c']:
        #            for s, v in self.breaks['c'][l[i]]:
        #                point = (1 - s) * self.c_pts[l[i]]['fp'] + s * self.c_pts[l[i+1]]['fp']
        #                ax.scatter(*self.S(*point))



        limits = np.array([ax.get_xlim(), ax.get_ylim(), ax.get_zlim()])
        ax.set_box_aspect(aspect=np.ptp(limits, axis=1))
        self.print("[%0.3fs] %s" % (time.perf_counter() - t0, "préparation du dessin"))
        self.print("[%0.3fs] %s" % (time.perf_counter() - T0, "Total"))
        # plt.gca()
        # plt.close(fig)
        # return 'pass'
        plt.show()

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
        t = self.XY(self.dS(*fp, *e.fdp))
        t = t / norm(t)
        vec = vec - np.inner(vec, t) * t
        return vec / norm(vec)

    def dirvec(self, p):  # direction de la surface dans l'image au point de contour p
        return self.kerdS(*p)[1]

    def bvis_chge(self, e, s):  # e est une arête de bord
        p, dp, Nb = e["fp"], e["fdp"], e["dir"]  # segment dans le domaine

        p0 = p + s * dp  # le point où la visibilité peut changer
        # p0 = (1-s) * p + s * q # le point où la visibilité peut changer

        ker = self.kerdS(*p0)[0]
        # on le dirige vers le pli supérieur, càd les z croissants
        if self.Z(self.dS(*p0, *ker)) < 0:
            ker = -ker

        Tb = dp / norm(dp)  # tangent au bord
        # Tb = (q-p)/norm(q-p) # tangent au bord
        # Nb = Tb.perp().dir(r-p) # normale pointant vers la surface

        Np = np.array(self.SD_jac(*p0, *self.axis))  # normale au pli
        if np.inner(Np, ker) < 0:
            Np = -Np  # dirigée vers le pli supérieur

        if (
            np.inner(ker, Nb) > 0
        ):  # il y a changement de visibilité, la surface est vers le pli supérieur
            v = 1 if np.inner(Tb, Np) > 0 else -1
            return v
        return 0

    def vis_chge(
        self, p, dir
    ):  # changement de visi quand part d'un pli en p dans la direction dir
        ker = self.kerdS(*p)[0]
        # on le dirige vers le pli supérieur, càd les z croissants
        if np.inner(self.dS(*p, *ker), self.axis) < 0:
            ker = -ker
        Np = np.array(self.SD_jac(*p, *self.axis))  # normale au pli
        if np.inner(Np, ker) < 0:
            Np = -Np  # dirigée vers le pli supérieur


        return 0 if np.inner(Np, dir) >= 0 else -1

    def kerdS(self, u, v):  # u et v peuvent être des float ou des vecteurs de floats
        # matrice(s) de la différentielle de la projection de S sur le plan
        du, dv = vect_prod(self.Su(u, v), self.axis), vect_prod(
            self.Sv(u, v), self.axis
        )
        A = np.stack((du.T, dv.T), axis=-1)
        # singular value decomp. vh est une matrice dont les lignes sont les vecteurs correspondant
        # aux valeurs singulières en ordre décroissant
        ux, s, vh = svd(A)
        # vecteur correspondant à la dernière val sing, on transpose pour avoir le format [u, v]
        ker = vh[..., -1, :].T
        diff2 = self.XY(self.d2S(u, v, ker[0], ker[1]))
        im = self.XY(self.dS(u, v, -ker[1], ker[0]))
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
                if zp + s * dzp > zq + t * dzq:
                    v = -e["n"] if np.inner(dq, e["d"]) > 0 else e["n"]
                    if v != 0:  # v= 0 pour les arêtes fictives
                        self.add_line_bk(
                            f["l_idx"], f["e_idx"], t, v
                        )  # e est audessus de f

                elif zp + s * dzp < zq + t * dzq:
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

    def propagate(self, l_idx, idx, visib): 
        # remplacé par "propagation"
        # propage la visibilité sur une ligne l_idx à partir du point idx de visibilité visib
        # attention, visib est la visibilité quand on est déjà dans la ligne
        l = self.lines[l_idx]
        visa = 0
        for bks in self.line_bks[l_idx][
            idx:
        ]:  # attention, idx peut-être négatif, y faut pas utiliser range
            for s, v in bks:
                visa = visa + v
        # on soustrait le changement de visi à l'extrêmité
        self.visibilities[tuple(l[-1]["ixfp"])] = visib + visa - l[-1]["v"]
        visb = 0
        for bks in self.line_bks[l_idx][0:idx]:
            for s, v in bks:
                visb = visb + v
        self.visibilities[tuple(l[0]["ixfp"])] = visib - visb - l[0]["v"]
        # print('propagated: line %i, length %i, point %i with visibility %i. Result: start visibility %i, end visibility %i'%
        #      (l_idx, len(self.lines[l_idx]), idx, visib, visib - visb - l[0]['v'], visib + visa - l[-1]['v']))
        # string = ' visibility changes : ' + str(l[0]['v'])+' '
        # for bks in self.line_bks[l_idx]:# attention, idx peut-être négatif, y faut pas utiliser range
        #    for s,v in bks:
        #        string += str(v) + ' '
        # string += str(-l[-1]['v'])
        # print(string)
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

    def plot_lines(self, ax):
        visible_lines = []  # tracées en dernier
        for i, l in enumerate(self.lines):
            type = "?"
            if l[0]["type"] == "??":
                continue
            vis = l[0]["v"]
            if tuple(l[0]["ixfp"]) in self.visibilities:
                type = l[0]["type"][1]
                vis += self.visibilities[tuple(l[0]["ixfp"])]
            line = []
            for j in range(len(l) - 1):
                p, q = l[j]["fp"], l[j + 1]["fp"]
                line.append(p)
                for s, v in self.line_bks[i][j]:
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

