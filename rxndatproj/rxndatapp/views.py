from django.shortcuts import render
from django.views import generic
from django.http import HttpResponseRedirect
from django.urls import reverse
from django.contrib.admin.views.decorators import staff_member_required
from .models import Transformation, Chemical, Condition, Substrate, Reaction
from .forms import TransformationForm, ChemicalForm, ConditionForm, SubstrateForm, ReactionForm

from rdkit import Chem
from rdkit.Chem import Draw

import base64
# import io

# Create your views here.
def home(request):
    latest_substrates = reversed(Substrate.objects.order_by('id')[:10])
    context = {
        'substrate_list': latest_substrates,
    }
    return render(request, 'home.html', context)

def get_base64_image_from_smiles(smi):
    # taken from https://iwatobipen.wordpress.com/2020/01/17/draw-rdkit-mol-reaction-object-on-html-without-static-png-image-rdkit-memo/
    m = Chem.MolFromSmiles(smi)
    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(150,150)
    drawer.DrawMolecule(m)
    drawer.FinishDrawing()

    # bytes_io = io.BytesIO()
    text = drawer.GetDrawingText()

    return base64.b64encode(text).decode('utf8')

class TransformationListView(generic.ListView):
    model = Transformation

    def get_context_data(self, *args, **kwargs):
        transformation_list = list(Transformation.objects.order_by('name'))
        context = super(TransformationListView, self).get_context_data(*args, **kwargs)
        transformation_details = []
        for transformation in transformation_list:
            related_conditions = Condition.objects.filter(transformation__id = transformation.id)
            num_related_substrates = 0
            for condition in related_conditions:
                related_substrates = Substrate.objects.filter(condition__id = condition.id)
                num_related_substrates += len(related_substrates)
            transformation_details.append({
                                           'records': f'{len(related_conditions)} condition(s), {num_related_substrates} substrate(s)',
                                           'reactants': get_base64_image_from_smiles(transformation.reactant_smiles),
                                           'products': get_base64_image_from_smiles(transformation.product_smiles),
                                         })
        context['transformation_list_with_details'] = zip(transformation_list, transformation_details)
        return context

def condition_by_transformation(request, pk):
    condition_list = Condition.objects.filter(transformation__id = pk)
    substrate_list = {}
    for condition in condition_list:
        substrate_list[condition.id] = Substrate.objects.filter(condition__id = condition.id)
    context = {
        'condition_and_substrate_list': zip(condition_list, substrate_list)
    }
    return render(request, 'rxndatapp/condition_by_transformation.html', context)

def condition_by_index(request, pk):
    condition = Condition.objects.filter(id = pk)
    context = {
        'cd': condition[0],
    }
    return render(request, 'rxndatapp/condition_by_index.html', context)

def substrate_by_index(request, pk):
    substrate = Substrate.objects.filter(id = pk)
    context = {
        'sb': substrate[0],
        'cd': substrate[0].condition,
    }
    return render(request, 'rxndatapp/substrate_by_index.html', context)

@staff_member_required
def add_transformation(request):
    if request.method == 'POST':
        form = TransformationForm(request.POST)

        if form.is_valid():
            transformation = Transformation(**form.cleaned_data)
            transformation.save()

        return HttpResponseRedirect(reverse('add_transformation'))

    else:
        initial = {}
        form = TransformationForm(initial = initial)

    context = {
        'form': form,
    }

    return render(request, 'rxndatapp/transformation_form.html', context)

@staff_member_required
def add_chemical(request):
    if request.method == 'POST':
        form = ChemicalForm(request.POST)

        if form.is_valid():
            chemical = Chemical(**form.cleaned_data)
            chemical.save()

        return HttpResponseRedirect(reverse('add_chemical'))

    else:
        initial = {}
        form = ChemicalForm(initial = initial)

    context = {
        'form': form,
    }

    return render(request, 'rxndatapp/chemical_form.html', context)

@staff_member_required
def add_condition(request):
    if request.method == 'POST':
        form = ConditionForm(request.POST)

        if form.is_valid():
            condition = Condition(**form.cleaned_data)
            condition.save()

        return HttpResponseRedirect(reverse('add_condition'))

    else:
        initial = {}
        form = ConditionForm(initial = initial)

    context = {
        'form': form,
    }

    return render(request, 'rxndatapp/condition_form.html', context)

@staff_member_required
def add_substrate(request):
    if request.method == 'POST':
        form = SubstrateForm(request.POST)

        if form.is_valid():
            substrate = Substrate(**form.cleaned_data)
            substrate.save()

        return HttpResponseRedirect(reverse('add_substrate'))

    else:
        initial = {
            'eln_page': 'EXP-##-AA####',
        }
        form = SubstrateForm(initial = initial)

    context = {
        'form': form,
    }

    return render(request, 'rxndatapp/substrate_form.html', context)

def add_reaction(request):
    if request.method == 'POST':
        form = ReactionForm(request.POST)

        if form.is_valid():
            reaction = Reaction(**form.cleaned_data)
            reaction.save()

        return HttpResponseRedirect(reverse('add_reaction'))

    else:
        initial = {
            'eln_page': 'EXP-##-AA####',
        }
        form = ReactionForm(initial = initial)

    context = {
        'form': form,
    }

    return render(request, 'rxndatapp/reaction_form.html', context)