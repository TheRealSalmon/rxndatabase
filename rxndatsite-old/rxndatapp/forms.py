from django import forms

class AddTransformation(forms.Form):
    name = forms.CharField(max_length=200)
    reactant_smarts = forms.CharField(max_length=50)
    product_smarts = forms.CharField(max_length=50)
