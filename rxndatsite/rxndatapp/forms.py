from django import forms

# from django.core.exceptions import ValidationError
# from django.utils.translation import ugettext_lazy as _

class AddTransformation(forms.Form):
    name = forms.CharField(max_length=200)
    reactant_smarts = forms.CharField(max_length=50)
    product_smarts = forms.CharField(max_length=50)
