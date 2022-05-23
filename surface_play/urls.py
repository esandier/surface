from django.urls import path
from .views import HomePageView, AboutPageView
from django.views.decorators.csrf import csrf_exempt

urlpatterns = [
    path("", csrf_exempt(HomePageView.as_view()), name="home"),
    path("about/", AboutPageView.as_view(), name="about"),
]