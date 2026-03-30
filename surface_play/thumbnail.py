import numpy as np
from .silhouette import Surface

THUMBNAIL_RES = 100


def _ij_from_angles(elev, azim, inplane=0.0):
    """Compute camera right (I) and up (J) unit vectors from Euler angles in degrees."""
    phi   = elev   * np.pi / 180
    theta = azim   * np.pi / 180
    psi   = inplane * np.pi / 180
    I0 = np.array([-np.sin(theta), np.cos(theta), 0.0])
    J0 = np.array([-np.sin(phi)*np.cos(theta), -np.sin(phi)*np.sin(theta), np.cos(phi)])
    I  = (I0 * np.cos(psi) + J0 * np.sin(psi)).tolist()
    J  = (-I0 * np.sin(psi) + J0 * np.cos(psi)).tolist()
    return I, J


def compute_thumbnail(rec, I=None, J=None, eye=None, elev=35.0, azim=45.0, inplane=0.0):
    """Generate an SVG thumbnail.

    If I and J are provided (camera right/up world vectors), they are used directly.
    Otherwise, I/J are derived from elev/azim/inplane angles.
    eye is the camera position for perspective; None means orthographic.
    """
    if I is None or J is None:
        I, J = _ij_from_angles(elev, azim, inplane)
        if eye is None:
            # build a default perspective eye from angles for backwards-compat
            pass  # stays orthographic unless caller supplies eye

    surf = Surface(
        rec.X, rec.Y, rec.Z, rec.parameter_names,
        bounds=(rec.u_min, rec.u_max, rec.v_min, rec.v_max),
        quotient=(rec.u_identify, rec.v_identify),
        domain_type=rec.domain_type, coord_type=rec.coord_type,
        r_min=rec.r_min, r_max=rec.r_max,
        output_type=rec.output_type,
    )
    surf.triangulate(THUMBNAIL_RES)
    surf.set_axis(I, J, eye=eye)
    surf.traitement()

    lines = surf.plot_for_browser()
    xy_grid = surf.XY(surf.S_grid)  # (2, N) — projected surface points
    x_min, x_max = float(xy_grid[0].min()), float(xy_grid[0].max())
    y_min, y_max = float(xy_grid[1].min()), float(xy_grid[1].max())
    origin = np.array([(x_min + x_max) / 2, (y_min + y_max) / 2])
    r = max(x_max - x_min, y_max - y_min) / 2

    parts = []
    vis_keys = sorted(lines)
    for idx, vis in enumerate(vis_keys):
        visible = idx == len(vis_keys) - 1
        opacity = '1' if visible else '0.35'
        dash = ' stroke-dasharray="0.05 0.05"' if not visible else ''
        for l in lines[vis]:
            pts = ' '.join(
                f'{(p[0] - origin[0]) / r:.4f},{-(p[1] - origin[1]) / r:.4f}'
                for p in l
            )
            if pts:
                parts.append(
                    f'<polyline points="{pts}" fill="none" stroke="steelblue"'
                    f' stroke-width="0.04" stroke-opacity="{opacity}"{dash}/>'
                )

    return (
        '<svg viewBox="-1.1 -1.1 2.2 2.2" xmlns="http://www.w3.org/2000/svg"'
        ' style="width:100%;height:100%;">'
        + ''.join(parts)
        + '</svg>'
    )
