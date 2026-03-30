from django.urls import path
from .views import SurfaceRecordListView, SurfaceRecordCreateView, SurfaceRecordUpdateView, SurfaceRecordDeleteView, SurfacePlayView, save_view
from django.views.decorators.csrf import csrf_exempt

urlpatterns = [
    path("", SurfaceRecordListView.as_view(), name='surfaces'),
    path('surface/add/', SurfaceRecordCreateView.as_view(), name='surface-add'),
    path('surface/<int:pk>/', SurfaceRecordUpdateView.as_view(), name='surface-update'),
    path('surfaceplay/<int:pk>/', csrf_exempt(SurfacePlayView.as_view()), name='surface-play'),
    path('surface/<int:pk>/delete/', SurfaceRecordDeleteView.as_view(), name='surface-delete'),
    path('surfaceplay/<int:pk>/saveview/', csrf_exempt(save_view), name='surface-save-view'),
]
