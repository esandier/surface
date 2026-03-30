# pages/views.py
from django.views.generic import TemplateView, ListView
from django.urls import reverse_lazy
from django.views.generic.edit import CreateView, DeleteView, UpdateView
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, redirect, render
from django.conf import settings

from collections import defaultdict

from surface_play.models import SurfaceRecord
from surface_play.forms import SurfaceRecordForm
from surface_play.thumbnail import compute_thumbnail
from .silhouette import *
import json

_surface_cache = {}  # pk -> (cache_key, triangulated Surface)

class SurfacePlayView(TemplateView):
    template_name = "play.html"

    # def dispatch(self, request, *args, **kwargs):
    #     self.pk = self.kwargs['pk']
    #     return super().dispatch(request, *args, **kwargs)

    def get(self, request, pk):
        rec = get_object_or_404(SurfaceRecord, pk=pk)
        surf = Surface(rec.X, rec.Y, rec.Z, rec.parameter_names,
                       bounds=(rec.u_min, rec.u_max, rec.v_min, rec.v_max),
                       quotient=(rec.u_identify, rec.v_identify),
                       domain_type=rec.domain_type, coord_type=rec.coord_type,
                       r_min=rec.r_min, r_max=rec.r_max,
                       output_type=rec.output_type)
        positions, normals, faces, center, radius = surf.for_3js()
        context = {}
        context['positions'] = positions
        context['normals'] = normals
        context['faces'] = faces
        context['center'] = center
        context['radius'] = radius
        context['surface'] = rec
        context['pk'] = pk
        context['initial_elev']    = rec.initial_elev
        context['initial_azim']    = rec.initial_azim
        context['initial_inplane'] = rec.initial_inplane
        context['initial_zoom']    = rec.initial_zoom
        context['debug_ui']        = settings.DEBUG

        return render(request, self.template_name, context)
        
    def post(self, request, pk):
        rec = get_object_or_404(SurfaceRecord, pk=pk)
        data = json.loads(request.body.decode())
        debug = data.get('debug', {})
        grid_res = int(debug.get('grid_res', settings.RESOLUTION))
        cache_key = (rec.X, rec.Y, rec.Z, rec.parameter_names,
                     rec.u_min, rec.u_max, rec.v_min, rec.v_max,
                     rec.u_identify, rec.v_identify,
                     rec.domain_type, rec.coord_type, rec.r_min, rec.r_max,
                     rec.output_type, grid_res)

        if _surface_cache.get(pk, (None,))[0] != cache_key:
            surf = Surface(rec.X, rec.Y, rec.Z, rec.parameter_names,
                           bounds=(rec.u_min, rec.u_max, rec.v_min, rec.v_max),
                           quotient=(rec.u_identify, rec.v_identify),
                           domain_type=rec.domain_type, coord_type=rec.coord_type,
                           r_min=rec.r_min, r_max=rec.r_max,
                           output_type=rec.output_type)
            surf.triangulate(grid_res)
            _surface_cache[pk] = (cache_key, surf)
        surf = _surface_cache[pk][1]

        I   = data['I']
        J   = data['J']
        O   = data['O']
        eye = data.get('eye')  # None for orthographic, [x,y,z] for perspective

        surf.set_axis(I, J, eye=eye)
        surf.traitement(
            use_lp=debug.get('use_lp', True),
            use_newton=debug.get('use_newton', True),
            simplify_pts=debug.get('simplify_pts', None),
        )

        lines = surf.plot_for_browser()
        lines_by_visibility = []
        for vis in sorted(lines):
            lines_out = []
            for i, l in enumerate(lines[vis]):
                string = ""
                for p in l:
                    string += " %f,%f " % (p[0], p[1])
                lines_out.append(string)
            lines_by_visibility.append(lines_out)

        response = {}
        response['lines_by_visibility'] = lines_by_visibility
        if eye is None:
            response['origin'] = surf.XY(np.array(O)).tolist()
        else:
            response['origin'] = surf.XY(np.array(surf.center)).tolist()
        # response['domain'] = surf.domain_data()

        return HttpResponse(json.dumps(response), content_type='application/json')        



def save_view(request, pk):
    if request.method != 'POST':
        return HttpResponse(status=405)
    rec = get_object_or_404(SurfaceRecord, pk=pk)
    data = json.loads(request.body.decode())
    rec.initial_elev    = float(data['elev'])
    rec.initial_azim    = float(data['azim'])
    rec.initial_inplane = float(data['inplane'])
    rec.initial_zoom    = float(data['zoom'])
    try:
        rec.thumbnail = compute_thumbnail(
            rec, elev=rec.initial_elev, azim=rec.initial_azim, inplane=rec.initial_inplane
        )
        rec.save(update_fields=['initial_elev', 'initial_azim', 'initial_inplane', 'initial_zoom', 'thumbnail'])
    except Exception:
        rec.save(update_fields=['initial_elev', 'initial_azim', 'initial_inplane', 'initial_zoom'])
    return HttpResponse('{}', content_type='application/json')


class SurfaceRecordListView(ListView):
    model = SurfaceRecord
    context_object_name = 'surfaces'

class SurfaceRecordCreateView(CreateView):
    model = SurfaceRecord
    form_class = SurfaceRecordForm
    context_object_name = 'surface'

    def form_valid(self, form):
        response = super().form_valid(form)
        try:
            self.object.thumbnail = compute_thumbnail(self.object)
            self.object.save(update_fields=['thumbnail'])
        except Exception:
            pass
        return response

class SurfaceRecordUpdateView(UpdateView):
    model = SurfaceRecord
    form_class = SurfaceRecordForm
    context_object_name = 'surface'

    def form_valid(self, form):
        response = super().form_valid(form)
        try:
            self.object.thumbnail = compute_thumbnail(self.object)
            self.object.save(update_fields=['thumbnail'])
        except Exception:
            pass
        return response

class SurfaceRecordDeleteView(DeleteView):
    model = SurfaceRecord
    success_url = reverse_lazy('surfaces')
    context_object_name = 'surface'

#surf = Surface('u','v', 'log(cos(v)/cos(u))', 'u v')
#surf.triangulate(bounds = (-1.5,1.5,-1.5,1.5), res = 53)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('(1.2+ cos(u))*cos(v)','(1.2 + cos(u))*sin(v) ','sin(u)', 'u v')
#surf.triangulate(bounds = (0,6,0, 5.5), res = 147, quotient=('no','no'))
#surf.set_axis(elev = 30, azim = -30)
#surf.traitement()
#surf.plot(T0)

#surf = Surface(
#"u*cos(v)", "u*sin(v)", "v", "u v"
#)  # demande une perturbation pas trop grande à azimuth 0.
#surf.triangulate(bounds=(-5, 5, 0, 25), res=200)
#surf.set_axis(elev=10, azim=8)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u', 'v', 'u*u*u + v*u', 'u v')
#surf.triangulate(bounds = (-1,1,-1,1), res = 32)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','u*v', 'u v' )
#surf.triangulate(bounds = (-2,2,-2,2), res = 61)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','2*sin(u*v)', 'u v' ) # fonctionne avec res = 600, res ligne = 200, bornes -5,5,-5,5
#surf.triangulate(bounds = (-5,5,-5,5), res = 700)
#surf.set_axis(elev = 30, azim = 38)
#surf.traitement()
#surf.plot(T0)
## mauvais résultat pour certains points de vue

#surf = Surface('u','v','3*sin(u)*sin(v)', 'u v' ) # marche bien à 200
#surf.triangulate(bounds = (0,15,0,15), res = 203)
#surf.set_axis(elev = 9, azim = 67)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','(1 - (u*u + v*v))**2', 'u v' )  # marche à 100, elev 20, azim 28
#surf.triangulate(bounds = (-1.3,1.3,-1.3,1.3), res = 33)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','3*(sin(u)+sin(v))**2', 'u v' )
#surf.triangulate(bounds = (0,20,0, 20), res = 506) # à (0,20,0,20) res 500, elev 25, azim 35 ça marche pas, même à 1500
#surf.set_axis(elev = 25, azim = 35)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','u*u', 'u v' )
#surf.triangulate(bounds = (-2,2,-2,2), res = 51)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','u*u - v*v', 'u v' )
#surf.triangulate(bounds = (-2,2,-2,2), res = 29)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','4*sin(u)', 'u v' )
#surf.triangulate(bounds = (0,10,-2,2), res = 43)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u*(4+ sin(10*exp(-5*(u*u+v*v))))/sqrt(u*u+v*v+.0001)','3*(u*u+v*v)','v*(4+ sin(10*exp(-5*(u*u+v*v))))/sqrt(u*u+v*v+.0001)', 'u v' )
#surf.triangulate(bounds = (-1,1,-1,1), res = 367)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','2*sin(u*v)', 'u v' )
#surf.triangulate(bounds = (-4,4,-4,4), res = 176)
#surf.set_axis(elev = 31, azim = 39)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('(3+ u*cos(3*v))*cos(v)','(3+ u*cos(3*v))*sin(v) ','u*sin(3*v)', 'u v')
#surf.triangulate(bounds = (-1,1,0, 6.28318530717958647), res = 143, quotient=('no','cy'))
#surf.set_axis(elev = 30, azim = -30)
#surf.traitement()
#surf.plot(T0)

# surf = Surface('u','v','sin(u*u + v*v)', 'u v' ) # temps à battre : entre 10 et 12,5 secondes
# surf.triangulate(bounds = (-5,5,-5,5), res = 500) # 600 OK avec x*x/(1+x), même pour -6,6,-6,6. 800 pour x^{3/2}/(1+sqrt(x))
# surf.set_axis(elev = 20, azim = 48)
# surf.traitement()
# surf.plot(T0)


#################################################################################################################
#################################################################################################################
#################################################################################################################


#surf = Surface('(3+ cos(u))*cos(v)','(3 + cos(u))*sin(v) ','sin(u)', 'u v')
#surf.triangulate(bounds = (.3,6, 0,6.2831853071795864769252), res = 100, quotient=('no','cy'))
#surf.set_axis(elev = 30, azim = -30)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('(2.5+ cos(u))*cos(v)','(2.5 + cos(u))*sin(v) ','sin(u)', 'u v')
#surf.triangulate(bounds = (0,6.2831853071795864769252, 0,6.2831853071795864769252), res = 200, quotient=('cy','cy'))
#surf.set_axis(elev = 26.6, azim = 14.88)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('(2+ sin(v))*cos(u)','(2+sin(v))*sin(u) ','v', 'u v')
#surf.triangulate(bounds = (0,6.2831853071795864769252, 1,10), res = 200, quotient=('cy','no'))
#surf.set_axis(elev = 10, azim = 40)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('(1.2+ sin(2*v))*cos(u)','(1.2+sin(2*v))*sin(u) ','v', 'u v') 
#surf.triangulate(bounds = (0,6.2831853071795864769252, 0,10), res = 100, quotient=('cy','no'))
#surf.set_axis(elev = 32.667, azim = -45.333)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('cos(u)','sin(u) ','v', 'u v')
#surf.triangulate(bounds = (0,6.2831853071795864769252, 3,9), res = 200, quotient=('cy','no'))
#surf.set_axis(elev = 10, azim = 40)
#surf.traitement()
#surf.plot(T0)

#surf = Surface("(2+ u*cos(v/2))*cos(v)", "(2+ u*cos(v/2))*sin(v) ", "u*sin(v/2)", "v u") # deconne elev=40, azim=40 si res = 600, uniquement si raccord "mo, no"
#surf.triangulate(bounds=(0, 6.28318, -1, 1), res=100, quotient=("mo", "no"))
#surf.set_axis(elev=40, azim=40)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('(2+ u*cos(3*v/2))*cos(v)','(2+ u*cos(3*v/2))*sin(v) ','u*sin(3*v/2)', 'u v')
#surf.triangulate(bounds = (-1.4,1.4,0, 6.2831), res = 100, quotient=('no','mo'))
#surf.set_axis(elev = 8, azim = 50)
#surf.traitement()
#surf.plot(T0)