# Generated by Django 2.2.4 on 2019-08-30 09:20

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('clinic_platform', '0003_auto_20190830_0840'),
    ]

    operations = [
        migrations.AddField(
            model_name='geneinfo',
            name='gene_description',
            field=models.CharField(default='NA', max_length=200),
        ),
        migrations.AlterField(
            model_name='transcriptinfo',
            name='gene_info',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='clinic_platform.GeneInfo'),
        ),
    ]