# Generated by Django 2.2.4 on 2019-09-12 12:15

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('clinic_platform', '0008_auto_20190912_1156'),
    ]

    operations = [
        migrations.AddField(
            model_name='cdsinfo',
            name='cds_coverage',
            field=models.TextField(default='None', max_length=2000),
        ),
        migrations.AddField(
            model_name='exoninfo',
            name='exon_coverage',
            field=models.TextField(default='None', max_length=2000),
        ),
        migrations.AlterField(
            model_name='transcriptinfo',
            name='gene_info',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='clinic_platform.GeneInfo'),
        ),
    ]