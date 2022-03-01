from django.shortcuts import render, get_object_or_404
from django.views import generic
from django.http import HttpResponseRedirect
from django.urls import reverse
from .models import Transformation, Condition, Reaction

from rxndatapp.forms import AddTransformation

# Create your views here.
def home(request):
    return render(request, 'home.html', None)

class TransformationListView(generic.ListView):
    model = Transformation

def conditions_by_transformation(request, pk):
    conditions_list = Condition.objects.filter(transformation__id = pk)
    context = {
        'conditions_by_transformation_list': conditions_list,
    }
    return render(request, 'rxndatapp/conditions_by_transformation.html', context)

def addtfm(request):

    if request.method == 'POST':
        form = AddTransformation(request.POST)

        if form.is_valid():
            new_tfm = Transformation(**form.cleaned_data)
            new_tfm.save()

        return HttpResponseRedirect(reverse('addtransformation'))

    else:
        initial_form = {
            'name': 'the reaction',
            'reactant_smarts': 'C',
            'product_smarts': 'C',
        }
        form = AddTransformation(initial=initial_form)

    context = {
        'form': form,
    }

    return render(request, 'rxndatapp/addtfm.html', context)
