# Generated by Django 3.1.2 on 2022-02-20 07:15

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('rxndatapp', '0002_auto_20220220_0656'),
    ]

    operations = [
        migrations.CreateModel(
            name='Condition',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('reagents', models.CharField(max_length=200)),
                ('solvent', models.CharField(max_length=20)),
                ('degree_of_caution', models.CharField(max_length=50)),
                ('comments', models.CharField(max_length=200)),
                ('transformation', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rxndatapp.transformation')),
            ],
        ),
        migrations.CreateModel(
            name='Reaction',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('reactant1_smiles', models.CharField(max_length=200)),
                ('reactant2_smiles', models.CharField(max_length=200)),
                ('reactant3_smiles', models.CharField(max_length=200)),
                ('reactant4_smiles', models.CharField(max_length=200)),
                ('reactant5_smiles', models.CharField(max_length=200)),
                ('product_smiles', models.CharField(max_length=200)),
                ('rxn_yield', models.IntegerField()),
                ('condition', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rxndatapp.condition')),
            ],
        ),
    ]