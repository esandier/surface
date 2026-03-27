import numpy as np
from .silhouette import Surface

THUMBNAIL_RES = 100

# Default viewpoint: elev=35°, azim=45°, inPlane=0°
# I = camera right, J = camera up (same convention as set_axis)
_phi   = 35 * np.pi / 180
_theta = 45 * np.pi / 180
DEFAULT_I = [-np.sin(_theta),                    np.cos(_theta),                 0]
DEFAULT_J = [-np.sin(_phi) * np.cos(_theta), -np.sin(_phi) * np.sin(_theta), np.cos(_phi)]


def compute_thumbnail(rec):
    surf = Surface(
        rec.X, rec.Y, rec.Z, rec.parameter_names,
        bounds=(rec.u_min, rec.u_max, rec.v_min, rec.v_max),
        quotient=(rec.u_identify, rec.v_identify),
    )
    surf.triangulate(THUMBNAIL_RES)
    surf.set_axis(DEFAULT_I, DEFAULT_J)
    surf.traitement()

    lines = surf.plot_for_browser()
    origin = surf.XY(np.array(surf.center))
    r = surf.radius

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
