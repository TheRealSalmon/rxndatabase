from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('transformation/', views.TransformationListView.as_view(), name='transformation'),
    path('transformation/add', views.addtfm, name='addtransformation'),
    path('transformation/<int:pk>', views.conditions_by_transformation, name='cndbytfm'),
]
