from django.db import models

# Create your models here.
class Transformation(models.Model):
    name = models.CharField(max_length=100)
    alias = models.CharField(max_length=300, blank=True)
    reactant_smarts1 = models.CharField(max_length=50)
    reactant_smarts2 = models.CharField(max_length=50, blank=True)
    reactant_smarts3 = models.CharField(max_length=50, blank=True)
    reactant_smarts4 = models.CharField(max_length=50, blank=True)
    reactant_smarts5 = models.CharField(max_length=50, blank=True)
    product_smarts = models.CharField(max_length=50)
    reactant_smiles = models.CharField(max_length=300, blank=True)
    product_smiles = models.CharField(max_length=300, blank=True)

    def __str__(self):
        return self.name

class Chemical(models.Model):
    common_name = models.CharField(max_length=300, blank=True)
    iupac_name = models.CharField(max_length=300)
    smiles = models.CharField(max_length=300)
    cas_number = models.CharField(max_length=20)

    def __str__(self):
        if self.common_name != '':
            return self.common_name
        else:
            return self.iupac_name

    def get_name_with_cas_number(self):
        if self.common_name != '':
            return f'{self.common_name} - {self.cas_number}'
        else:
            return f'{self.iupac_name} - {self.cas_number}'

class Condition(models.Model):
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
    
    transformation = models.ForeignKey(
        Transformation, 
        on_delete=models.CASCADE
    )
    reactant_smarts1 = models.CharField(max_length=50)
    reactant_smarts2 = models.CharField(max_length=50, blank=True)
    reactant_smarts3 = models.CharField(max_length=50, blank=True)
    reactant_smarts4 = models.CharField(max_length=50, blank=True)
    reactant_smarts5 = models.CharField(max_length=50, blank=True)
    catalyst1 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='cat1',
    )
    catalyst2 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='cat2',
    )
    catalyst3 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='cat3',
    )
    reagent1 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='rgt1',
    )
    reagent2 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='rgt2',
    )
    reagent3 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='rgt3',
    )
    reagent4 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='rgt4',
    )
    reagent5 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='rgt5',
    )
    solvent1 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='svt1',
    )
    solvent2 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='svt2',
    )
    solvent3 = models.ForeignKey(
        Chemical,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='svt3',
    )
    product_smarts = models.CharField(max_length=50)
    air_free = models.IntegerField(choices=air_free_choices, default=0)
    water_free = models.IntegerField(choices=water_free_choices, default=0)
    comment = models.CharField(max_length=300, blank=True)

    def __str__(self):
        return f'{self.transformation.name} ({self.comment})'

class Substrate(models.Model):
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

    condition = models.ForeignKey(
        Condition, 
        on_delete=models.CASCADE,
    )
    reactant_smiles1 = models.CharField(max_length=300)
    reactant_smiles2 = models.CharField(max_length=300, blank=True)
    reactant_smiles3 = models.CharField(max_length=300, blank=True)
    reactant_smiles4 = models.CharField(max_length=300, blank=True)
    reactant_smiles5 = models.CharField(max_length=300, blank=True)
    cat1_eq = models.FloatField(blank=True, null=True)
    cat2_eq = models.FloatField(blank=True, null=True)
    cat3_eq = models.FloatField(blank=True, null=True)
    rgt1_eq = models.FloatField(blank=True, null=True)
    rgt2_eq = models.FloatField(blank=True, null=True)
    rgt3_eq = models.FloatField(blank=True, null=True)
    rgt4_eq = models.FloatField(blank=True, null=True)
    rgt5_eq = models.FloatField(blank=True, null=True)
    svt1_vol = models.FloatField(blank=True, null=True)
    svt2_vol = models.FloatField(blank=True, null=True)
    svt3_vol = models.FloatField(blank=True, null=True)
    product_smiles = models.CharField(max_length=300)
    air_free = models.IntegerField(choices=air_free_choices, default=0)
    water_free = models.IntegerField(choices=water_free_choices, default=0)
    rxn_yield = models.FloatField()
    comment = models.CharField(max_length=300, blank=True)
    eln_page = models.CharField(max_length=20)

    def __str__(self):
        return f'{self.condition} - {self.eln_page}'

class Reaction(models.Model):
    name = models.CharField(max_length=100)
    suggested = models.BooleanField(default=True)
    reactant_smiles1 = models.CharField(max_length=300)
    reactant_smiles2 = models.CharField(max_length=300, blank=True)
    reactant_smiles3 = models.CharField(max_length=300, blank=True)
    reactant_smiles4 = models.CharField(max_length=300, blank=True)
    reactant_smiles5 = models.CharField(max_length=300, blank=True)
    reagents = models.CharField(max_length=300)
    solvents = models.CharField(max_length=50)
    procedure = models.CharField(max_length=500)
    rxn_yield = models.FloatField()
    eln_page = models.CharField(max_length=20)