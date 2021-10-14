from __future__ import unicode_literals

from django.db import models

# Create your models here.


class Experiment(models.Model):
	"""This describes the experimental setup of the gene expression data"""

	# Fields

	experiment_shortname = models.CharField(max_length=30, help_text='Enter field documentation')
	comparisons = models.TextField(help_text='Enter the conditions separated by a comma.')
	time_series = models.BooleanField()

	def __str__(self):
		return self.experiment_shortname


class ExpressionSet(models.Model):
	comparison = models.CharField(max_length=160, help_text='TBD')
	source_file = models.CharField(max_length=160, help_text='TBD')
	experiment = models.ForeignKey(Experiment, on_delete=models.CASCADE)

	def __str__(self):
		return self.comparison


class GeneExpression(models.Model):
	expression_set = models.ForeignKey(ExpressionSet, on_delete=models.CASCADE)
	gene_id = models.CharField(max_length=60, help_text='')
	p_value = models.FloatField(null=True)
	q_value = models.FloatField(null=True)
	fold_change = models.FloatField(null=True)

	def __str__(self):
		return self.gene_id


class CurrentSettings(models.Model):
	"""The settings for the current user and session"""

	# Fields

	current_experimental_setup = models.CharField(max_length=30, help_text='Currently selected experimental conditions')


class PanelGeneList(models.Model):
	"""A list of gene related to a panel"""

	panel_id = models.CharField(max_length=60, help_text='Descriptor of panel')
	phenotype = models.CharField(max_length=60, help_text='Phenotypic information')
	gene_list = models.TextField(max_length=10000, help_text='Comma separated list of gene IDs')

	bam_files = models.TextField(max_length=2000, default='None')

	vcf_files = models.TextField(max_length=2000, default='NA')

	comparison = models.CharField(max_length=160, default='None')

	reference_genome_used = models.CharField(max_length=160, default='None')

	last_threshold = models.IntegerField(default=-1)

	def __str__(self):
		return self.panel_id


class GeneInfo(models.Model):
	"""Info for genes loaded into the system"""

	gene_id = models.CharField(max_length=60, help_text='')
	GO_annotation = models.CharField(max_length=200, help_text='', default='NA')
	gene_name = models.CharField(max_length=200, help_text='', default='NA')
	gene_description = models.CharField(max_length=200, help_text='', default='NA')
	gene_chrom = models.CharField(max_length=200, help_text='', default='NA')

	gene_ENS_protein = models.CharField(max_length=60, help_text='', default='NA')
	gene_ENS_transcript = models.CharField(max_length=60, help_text='', default='NA')

	gene_biotype = models.CharField(max_length=60, help_text='', default='NA')
	gene_strand = models.CharField(max_length=10, help_text='', default='NA')

	gene_genome_assembly = models.CharField(max_length=100, help_text='', default='NA')

	def __str__(self):
		return self.gene_id


class TranscriptInfo(models.Model):

	transcript_id = models.CharField(max_length=60, help_text='')

	#exon_positions = models.CharField(max_length=10000, help_text='', default='NA')

	gene_info = models.ForeignKey(GeneInfo, on_delete=models.CASCADE)

	def __str__(self):
		return self.transcript_id


class ExonInfo(models.Model):

	exon_id = models.CharField(max_length=60, help_text='')

	exon_start = models.IntegerField(default=-1)

	exon_stop = models.IntegerField(default=-1)

	transcript_info = models.ForeignKey(TranscriptInfo, on_delete=models.CASCADE)

	coverage = models.TextField(max_length=2000, default='None')

	covOverThreshold = models.TextField(max_length=2000, default='None')

	pb_coverage = models.TextField(max_length=1000000, default='None')

	selected_cov = models.IntegerField(default=-1)

	# Format is "thisbam:value,otherbam:value"

	def __str__(self):
		return self.exon_id


class CDSInfo(models.Model):

	cds_id = models.CharField(max_length=60, help_text='')

	protein_id = models.CharField(max_length=60, help_text='')

	cds_start = models.IntegerField(default=-1)

	cds_stop = models.IntegerField(default=-1)

	cds_frame = models.IntegerField(default=-1)

	coverage = models.TextField(max_length=2000, default='None')

	covOverThreshold = models.TextField(max_length=2000, default='None')

	# Per pase coverage
	pb_coverage = models.TextField(max_length=1000000, default='None')

	selected_cov = models.IntegerField(default=-1)

	transcript_info = models.ForeignKey(TranscriptInfo, on_delete=models.CASCADE)

	def __str__(self):
		return self.cds_id


class SampleInfo(models.Model):

	experiment = models.ForeignKey(PanelGeneList, on_delete=models.CASCADE)

	sample_id = models.TextField(max_length=200, default='None')

	phenotype_info = models.TextField(max_length=1000000, default='None')

	def __str__(self):
		return self.sample_id




