from django import forms
from .models import PanelGeneList
from django.forms import ModelForm


class UploadFileForm(forms.Form):
    title = forms.CharField(max_length=50)
    file = forms.FileField()


class GeneQueryForm(forms.Form):
    gene_name = forms.CharField(max_length=50)


class PanelListForm(ModelForm):
    class Meta:
        model = PanelGeneList
        fields = ['panel_id', 'phenotype', 'gene_list', 'reference_genome_used']


