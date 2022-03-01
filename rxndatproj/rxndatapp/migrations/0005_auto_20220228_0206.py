# Generated by Django 3.1.2 on 2022-02-28 02:06

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('rxndatapp', '0004_auto_20220227_0059'),
    ]

    operations = [
        migrations.AlterField(
            model_name='substrate',
            name='cat1_eq',
            field=models.FloatField(blank=True, default=0.0),
        ),
        migrations.AlterField(
            model_name='substrate',
            name='cat2_eq',
            field=models.FloatField(blank=True, default=0.0),
        ),
        migrations.AlterField(
            model_name='substrate',
            name='cat3_eq',
            field=models.FloatField(blank=True, default=0.0),
        ),
        migrations.AlterField(
            model_name='substrate',
            name='rgt1_eq',
            field=models.FloatField(blank=True, default=0.0),
        ),
        migrations.AlterField(
            model_name='substrate',
            name='rgt2_eq',
            field=models.FloatField(blank=True, default=0.0),
        ),
        migrations.AlterField(
            model_name='substrate',
            name='rgt3_eq',
            field=models.FloatField(blank=True, default=0.0),
        ),
        migrations.AlterField(
            model_name='substrate',
            name='rgt4_eq',
            field=models.FloatField(blank=True, default=0.0),
        ),
        migrations.AlterField(
            model_name='substrate',
            name='rgt5_eq',
            field=models.FloatField(blank=True, default=0.0),
        ),
        migrations.AlterField(
            model_name='substrate',
            name='svt1_vol',
            field=models.FloatField(blank=True, default=0.0),
        ),
        migrations.AlterField(
            model_name='substrate',
            name='svt2_vol',
            field=models.FloatField(blank=True, default=0.0),
        ),
        migrations.AlterField(
            model_name='substrate',
            name='svt3_vol',
            field=models.FloatField(blank=True, default=0.0),
        ),
    ]