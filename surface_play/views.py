# pages/views.py
from django.views.generic import TemplateView
from .silhouette import *


class HomePageView(TemplateView):
    template_name = "home.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        surf = Surface('u','v','sin(u*u + v*v)', 'u v' ) 
        positions, normals, faces, center, radius = surf.for_3js(bounds = (-5,5,-5,5), res = 100) 
        context['positions'] = positions
        context['normals'] = normals
        context['faces'] = faces
        context['center'] = center
        context['radius'] = radius

        return context    

class AboutPageView(TemplateView):
    template_name = "about.html"