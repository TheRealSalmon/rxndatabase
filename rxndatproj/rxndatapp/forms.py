from django import forms
from .models import Transformation, Chemical, Condition

class TransformationForm(forms.Form):
    name = forms.CharField(max_length=100)
    reactant_smarts1 = forms.CharField(max_length=50)
    reactant_smarts2 = forms.CharField(max_length=50, required=False)
    reactant_smarts3 = forms.CharField(max_length=50, required=False)
    reactant_smarts4 = forms.CharField(max_length=50, required=False)
    reactant_smarts5 = forms.CharField(max_length=50, required=False)
    product_smarts = forms.CharField(max_length=50)

class ChemicalForm(forms.Form):
    common_name = forms.CharField(max_length=300, required=False)
    iupac_name = forms.CharField(max_length=300)
    smiles = forms.CharField(max_length=300)
    cas_number = forms.CharField(max_length=20)

class ConditionForm(forms.Form):
    AIR_AMBIENT = 0
    AIR_PURGE = 1
    AIR_BACKFILL = 2
    air_free_choices = [
        (AIR_AMBIENT, 'ambient air'),
        (AIR_PURGE, 'argon purge'),
        (AIR_BACKFILL, 'nitrogen/argon backfill'),
    ]
    WATER_AMBIENT = 0
    WATER_OVENDRY = 1
    water_free_choices = [
        (WATER_AMBIENT, 'no dry'),
        (WATER_OVENDRY, 'oven-dried'),
    ]
    
    transformation = forms.ModelChoiceField(
        queryset=Transformation.objects.all(),
    )
    reactant_smarts1 = forms.CharField(max_length=50)
    reactant_smarts2 = forms.CharField(max_length=50, required=False)
    reactant_smarts3 = forms.CharField(max_length=50, required=False)
    reactant_smarts4 = forms.CharField(max_length=50, required=False)
    reactant_smarts5 = forms.CharField(max_length=50, required=False)
    catalyst1 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    catalyst2 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    catalyst3 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    reagent1 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    reagent2 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    reagent3 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    reagent4 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    reagent5 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    solvent1 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    solvent2 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    solvent3 = forms.ModelChoiceField(
        queryset=Chemical.objects.all(),
        required=False
    )
    product_smarts = forms.CharField(max_length=50)
    air_free = forms.ChoiceField(choices=air_free_choices)
    water_free = forms.ChoiceField(choices=water_free_choices)
    comment = forms.CharField(max_length=300, required=False)

class SubstrateForm(forms.Form):
    AIR_AMBIENT = 0
    AIR_PURGE = 1
    AIR_BACKFILL = 2
    air_free_choices = [
        (AIR_AMBIENT, 'ambient air'),
        (AIR_PURGE, 'argon purge'),
        (AIR_BACKFILL, 'nitrogen/argon backfill'),
    ]
    WATER_AMBIENT = 0
    WATER_OVENDRY = 1
    water_free_choices = [
        (WATER_AMBIENT, 'no dry'),
        (WATER_OVENDRY, 'oven-dried'),
    ]

    condition = forms.ModelChoiceField(
        queryset=Condition.objects.all()
    )
    reactant_smiles1 = forms.CharField(max_length=300)
    reactant_smiles2 = forms.CharField(max_length=300, required=False)
    reactant_smiles3 = forms.CharField(max_length=300, required=False)
    reactant_smiles4 = forms.CharField(max_length=300, required=False)
    reactant_smiles5 = forms.CharField(max_length=300, required=False)
    cat1_eq = forms.FloatField(required=False)
    cat2_eq = forms.FloatField(required=False)
    cat3_eq = forms.FloatField(required=False)
    rgt1_eq = forms.FloatField(required=False)
    rgt2_eq = forms.FloatField(required=False)
    rgt3_eq = forms.FloatField(required=False)
    rgt4_eq = forms.FloatField(required=False)
    rgt5_eq = forms.FloatField(required=False)
    svt1_vol = forms.FloatField(required=False)
    svt2_vol = forms.FloatField(required=False)
    svt3_vol = forms.FloatField(required=False)
    product_smiles = forms.CharField(max_length=300)
    air_free = forms.ChoiceField(choices=air_free_choices)
    water_free = forms.ChoiceField(choices=water_free_choices)
    rxn_yield = forms.FloatField()
    comment = forms.CharField(max_length=300, required=False)
    eln_page = forms.CharField(max_length=20)

class ReactionForm(forms.Form):
    name = forms.CharField(max_length=100)
    reactant_smiles1 = forms.CharField(max_length=300)
    reactant_smiles2 = forms.CharField(max_length=300, required=False)
    reactant_smiles3 = forms.CharField(max_length=300, required=False)
    reactant_smiles4 = forms.CharField(max_length=300, required=False)
    reactant_smiles5 = forms.CharField(max_length=300, required=False)
    product_smiles = forms.CharField(max_length=300)
    reagents = forms.CharField(max_length=300)
    solvents = forms.CharField(max_length=50)
    procedure = forms.CharField(max_length=500)
    rxn_yield = forms.FloatField()
    eln_page = forms.CharField(max_length=20)