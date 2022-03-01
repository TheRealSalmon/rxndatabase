# Generated by Django 3.1.2 on 2022-02-27 00:48

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('rxndatapp', '0002_auto_20220225_2101'),
    ]

    operations = [
        migrations.CreateModel(
            name='Reaction',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('suggested', models.BooleanField(default=True)),
                ('reactant_smiles1', models.CharField(max_length=300)),
                ('reactant_smiles2', models.CharField(blank=True, max_length=300)),
                ('reactant_smiles3', models.CharField(blank=True, max_length=300)),
                ('reactant_smiles4', models.CharField(blank=True, max_length=300)),
                ('reactant_smiles5', models.CharField(blank=True, max_length=300)),
                ('reagents', models.CharField(max_length=300)),
                ('solvents', models.CharField(max_length=50)),
                ('procedure', models.CharField(max_length=500)),
                ('rxn_yield', models.FloatField()),
                ('eln_page', models.CharField(max_length=20)),
            ],
        ),
        migrations.RemoveField(
            model_name='condition',
            name='suggested',
        ),
        migrations.RemoveField(
            model_name='substrate',
            name='suggested',
        ),
        migrations.RemoveField(
            model_name='transformation',
            name='suggested',
        ),
    ]
