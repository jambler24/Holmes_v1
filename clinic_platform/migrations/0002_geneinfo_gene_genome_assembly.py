# Generated by Django 3.1.1 on 2021-09-20 11:18

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('clinic_platform', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='geneinfo',
            name='gene_genome_assembly',
            field=models.CharField(default='NA', max_length=100),
        ),
    ]
