from django.db import models

# Create your models here.

class Transformation(models.Model):
    reactant_smarts = models.CharField(max_length=50)
    product_smarts = models.CharField(max_length=50)
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name

class Condition(models.Model):
    transformation = models.ForeignKey(
        Transformation,
        on_delete=models.CASCADE,
    )
    reagents = models.CharField(max_length=200)
    solvent = models.CharField(max_length=20)
    degree_of_caution = models.CharField(max_length=50)
    comments = models.CharField(max_length=200)

    def __str__(self):
        return self.comments

class Reaction(models.Model):
    condition = models.ForeignKey(
        Condition,
        on_delete=models.CASCADE,
    )
    reactant1_smiles = models.CharField(max_length=200)
    reactant2_smiles = models.CharField(max_length=200)
    reactant3_smiles = models.CharField(max_length=200)
    reactant4_smiles = models.CharField(max_length=200)
    reactant5_smiles = models.CharField(max_length=200)
    product_smiles = models.CharField(max_length=200)
    rxn_yield = models.IntegerField()

    def __str__(self):
        return self.product_smiles
