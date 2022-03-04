from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('home/', views.home, name='home'),
    path('transformation/', views.TransformationListView.as_view(), name='transformation'),
    path('transformation/<int:pk>/', views.condition_by_transformation, name='condition_by_transformation'),
    path('condition/<int:pk>/', views.condition_by_index, name='condition_by_index'),
    path('substrate/<int:pk>/', views.substrate_by_index, name='substrate_by_index'),

    # path('add/', views.add, name='add')
    path('add/transformation/', views.add_transformation, name='add_transformation'),
    path('add/chemical/', views.add_chemical, name='add_chemical'),
    path('add/condition/', views.add_condition, name='add_condition'),
    path('add/substrate/', views.add_substrate, name='add_substrate'),
    path('add/reaction/', views.add_reaction, name='add_reaction'),
]