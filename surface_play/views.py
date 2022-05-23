# pages/views.py
from django.views.generic import TemplateView
from django.http import HttpResponse, HttpResponseRedirect
from .silhouette import *
import json



class HomePageView(TemplateView):
    template_name = "home.html"
    surf = Surface('(2+ u*cos(3*v/2))*cos(v)','(2+ u*cos(3*v/2))*sin(v) ','u*sin(3*v/2)', 'u v', bounds = (-1.4,1.4,0, 6.2831), quotient=('no','mo'))

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        positions, normals, faces, center, radius = self.surf.for_3js() 
        context['positions'] = positions
        context['normals'] = normals
        context['faces'] = faces
        context['center'] = center
        context['radius'] = radius
        return context
        
    def post(self, request):
        data = json.loads(request.body.decode())
        I = data['I']
        J = data['J']
        zoom = data['zoom']

        self.surf.triangulate(70)
        self.surf.set_axis(I,J)
        self.surf.traitement()
        
        lines, center, radius = self.surf.plot_for_browser()
        for vis in lines:
            for i, l in enumerate(lines[vis]):
                string = ""
                for p in l:
                    string += " %f,%f " % (p[0], p[1])
                lines[vis][i] = string

        response = {}
        response['lines'] = lines[0]
        response['center'] = center
        response['radius'] = radius
        response['zoom'] = zoom

        return HttpResponse(json.dumps(response), content_type='application/json')        


class AboutPageView(TemplateView):
    template_name = "about.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # surf = Surface('u','v','2*sin(u)', 'u v', bounds = (-3,3,-3,3))
        # surf.triangulate(70)
        surf = Surface('(2+ u*cos(3*v/2))*cos(v)','(2+ u*cos(3*v/2))*sin(v) ','u*sin(3*v/2)', 'u v', bounds = (-1.4,1.4,0, 6.2831), quotient=('no','mo'))
        surf.triangulate(100)        
        
        surf.set_axis(elev=40, azim=40)
        surf.traitement()
        lines, center, radius = surf.plot_for_browser()
        for vis in lines:
            for i, l in enumerate(lines[vis]):
                string = ""
                for p in l:
                    string += " %f,%f " % (p[0], p[1])
                lines[vis][i] = string
        context['lines'] = lines[0]
        context['center'] = center
        context['radius'] = radius

        return context    
    