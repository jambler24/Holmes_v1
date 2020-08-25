from django.http import HttpResponse
from django.shortcuts import render
from django.conf import settings
from django.shortcuts import redirect
from django.core.files.storage import FileSystemStorage
import json
import os
import networkx as nx
from networkx.readwrite import json_graph
#import coverage_tools as cov_tools
import clinic_platform.coverage_tools as cov_tools
from .processing import *
from .forms import UploadFileForm, PanelListForm
from clinic_platform.models import GeneExpression, ExpressionSet, Experiment, CurrentSettings, PanelGeneList, TranscriptInfo, GeneInfo, CDSInfo
from os import listdir
from os.path import isfile, join
from bs4 import BeautifulSoup
import requests
import pandas as pd

# Container Paths

anno_folder = settings.ANNO_FOLDER
bam_files_dir = settings.BAM_FILES_DIR
ref_genome_folder = settings.REF_GENOME_FOLDER
vcf_folder = settings.VARIANT_FOLDER
custom_anno_folder = anno_folder


def home(request):

	loaded_genomes = check_loaded_genomes(ref_genome_folder)

	loaded_bam_files = check_loaded_bam(bam_files_dir)

	panel_count = len(list(PanelGeneList.objects.all()))

	return render(request, 'index.html', {
		'ref_files': loaded_genomes,
		'loaded_panel_count': panel_count,
		'loaded_bam': len(loaded_bam_files),
	})


def uploads(request):
	if request.method == 'POST' and request.FILES['myfile']:
		#myfile = request.FILES['myfile']
		#fs = FileSystemStorage()
		#filename = fs.save(myfile.name, myfile)
		#uploaded_file_url = fs.url(filename)
		uploaded_file_url = 'Temp disabled'
		file_type = 'gff'

		if file_type == 'gff':

			gff2network('/Users/panix/iCloud/programs/Holmes/holmes_core/uploads/sequence.gff3',
						'/Users/panix/iCloud/programs/Holmes/holmes_core/processed')

		file_type = 'expression'

		if file_type == 'expression':

			genome_graph = nx.read_graphml('processed/test_result_network.xml')

			expression2network('uploads/gene_exp.diff', genome_graph, 'processed')

			print("processing expression data: complete")

		file_type = 'variants'

		if file_type == 'variants':
			print("processing variant data")

			genome_graph = nx.read_graphml('processed/test_exp_network.xml')


			variants2network('uploads/S5527_aligned_to_h37rv_RG_sorted_filtered_ann_unique.vcf', genome_graph, 'processed')
			variants2network('uploads/S507_aligned_to_h37rv_RG_sorted_filtered_ann_unique.vcf', genome_graph, 'processed')

			print("processing variant data: complete")


		# Convert gff to GML

		return render(request, 'uploads.html', {'uploaded_file_url': uploaded_file_url})

	return render(request, 'uploads.html')


def gene_list_input(request):

	if request.method == 'POST':
		# Extract info

		form = PanelListForm(request.POST)

		form_panel_id = request.POST['panel_id']

		selected_bam = request.POST.getlist('bamFiles')

		ref_genome_gff = request.POST['anno_selection']

		# Add to DB
		if form.is_valid():
			u = form.save()
			new_panel_obj = PanelGeneList.objects.get(panel_id=form_panel_id)

			new_panel_obj.bam_files = ','.join(selected_bam)
			new_panel_obj.save()

			# Add genes to the database

			panel_val = form.cleaned_data['panel_id']

			in_file_path = ref_genome_folder + ref_genome_gff

			import_result = subset_gff_to_db(in_file_path, panel_val)

			print(import_result)

			return redirect('/')

	else:

		form_class = PanelListForm

		all_bam_files = os.listdir(bam_files_dir)
		clean_bam_list = []

		for a_bam in all_bam_files:
			if a_bam[-3:].lower() == 'bam':
				clean_bam_list.append(a_bam)

		# Get ref anno files

		all_ref_gff_files = os.listdir(ref_genome_folder)
		clean_gff_list = []

		for a_file in all_ref_gff_files:
			if a_file[-3:].lower() == 'gff':
				clean_gff_list.append(a_file)

		return render(request, 'gene_list_input.html', {'form': form_class, 'bam_list': clean_bam_list, 'gff_list': clean_gff_list})


def load_data_to_network(request):
	if request.method == 'POST' and request.FILES['myfile']:
		#myfile = request.FILES['myfile']
		#fs = FileSystemStorage()
		#filename = fs.save(myfile.name, myfile)
		#uploaded_file_url = fs.url(filename)
		uploaded_file_url = 'Temp disabled'
		file_type = 'gff'

		if file_type == 'gff':

			gff2network('/Users/panix/iCloud/programs/Holmes/holmes_core/uploads/sequence.gff3',
						'/Users/panix/iCloud/programs/Holmes/holmes_core/processed')

		file_type = 'expression'

		if file_type == 'expression':

			genome_graph = nx.read_graphml('processed/test_result_network.xml')

			expression2network('uploads/gene_exp.diff', genome_graph, 'processed')

			print("processing expression data: complete")

		file_type = 'variants'

		if file_type == 'variants':
			print("processing variant data")

			genome_graph = nx.read_graphml('processed/test_exp_network.xml')


			variants2network('uploads/S5527_aligned_to_h37rv_RG_sorted_filtered_ann_unique.vcf', genome_graph, 'processed')
			variants2network('uploads/S507_aligned_to_h37rv_RG_sorted_filtered_ann_unique.vcf', genome_graph, 'processed')

			print("processing variant data: complete")


		# Convert gff to GML

		return render(request, 'uploads.html', {'uploaded_file_url': uploaded_file_url})

	return render(request, 'uploads.html')


def gene_search(request):

	if request.method == 'POST':
		# Extract gene subnet
		a_gene = request.POST['gene_name']
		genome_graph = nx.read_graphml('processed/test_var_network.xml')
		graph_net = extract_subview_genome(genome_graph, a_gene)

		# extract expression of all genes in subnet
		res_info = subnet_to_json_dict(graph_net)
		json_res_info = json.dumps(res_info)
		print(json_res_info)

		# Dump to JSON for vis
		json_dict = subnet2json(graph_net)

		# data = json.dumps({"nodes": json_data["nodes"], "links": json_data["links"]})
		data = json.dumps(json_dict)

		#print(json_data)
		#outfile = open('processed/tempGraph.json', 'w')
		#outfile.write(data)
		#outfile.close()

		# Get expression data
		expression_table_data = extract_and_format_gene_expression_table(graph_net, a_gene)

		# Format

		return redirect('/')

	else:

		return render(request, 'search_genes.html')


def gene_view(request, gene_id):
	# Extract gene subnet
	genome_graph = nx.read_graphml('processed/test_var_network.xml')
	graph_net = extract_subview_genome(genome_graph, gene_id)

	# extract expression of all genes in subnet
	res_info = subnet_to_json_dict(graph_net)
	json_res_info = json.dumps(res_info)
	print(json_res_info)

	# Dump to JSON for vis
	json_dict = subnet2json(graph_net)

	# data = json.dumps({"nodes": json_data["nodes"], "links": json_data["links"]})
	data = json.dumps(json_dict)

	# print(json_data)
	# outfile = open('processed/tempGraph.json', 'w')
	# outfile.write(data)
	# outfile.close()

	# Get expression data
	expression_table_data = extract_and_format_gene_expression_table(graph_net, gene_id)

	# Format

	return render(request, 'gene_view.html', {'JSON_data': data, 'table': expression_table_data, 'gene_info': json_res_info, 'gene_name': gene_id})


def sub_graphs(request):

	genome_graph = nx.read_graphml('processed/test_var_network.xml')

	subnet_list = extract_subnets(genome_graph, condition='507midlogtreated_5527midlogtreated')

	request.session['subnet_data'] = subnet_list

	return render(request, 'subnet_list.html', {'table_data': subnet_list})


def variant_network_overview(request):

	#gene_graph = nx.read_graphml('processed/test_var_network.xml')
	#genome_graph = nx.read_graphml('processed/test3genome.xml')

	variant_info = get_variant_info(gene_graph, genome_graph)

	#subnet_list = extract_subnets(genome_graph, condition='507midlogtreated_5527midlogtreated')

	#request.session['subnet_data'] = subnet_list

	return render(request, 'variant_info.html', {'table_data': 'this'})


def variant_overview(request, variant='default', panel='default'):

	#anno_set_selection = "Exon"
	#anno_set_selection = "CDS"
	#anno_set_selection = 'test_1_bed'

	vcf_headder_list = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

	# Get available panels

	panel_obj = list(PanelGeneList.objects.all())
	panel_list = []

	for gene_panel in panel_obj:
		panel_list.append(str(gene_panel))

	#gene_list = panel_obj.gene_list.replace('\r\n', '').split(',')


	if request.method == 'POST':

		panel_obj = PanelGeneList.objects.get(panel_id=request.POST['panel_selection'])

		# process appropriate vcf

		vcf_file_list = panel_obj.vcf_files.split(',')

		vcf_info_list = []

		gene_list = []
		var_list = []
		sample_list = []
		var_mutations = {}

		for a_vcf_file in vcf_file_list:

			file_name = '_'.join(a_vcf_file.split('.')[:-1])

			vcf_info = input_parser(vcf_folder + a_vcf_file)

			vcf_dict = {
				'file_name': file_name,
				'vcf_info': vcf_info
			}

			vcf_info_list.append(vcf_dict)

		# Find overlaps where vars are in Exons for the panel

		panel_gene_list = panel_obj.gene_list.split(',')

		result_dict = {}

		for a_gene in panel_gene_list:

			# Get transcripts

			gene_obj = GeneInfo.objects.get(gene_id=a_gene)

			gene_chrom = gene_obj.gene_chrom.replace('chr','')

			transcripts = list(TranscriptInfo.objects.filter(gene_info=gene_obj))

			for trans in transcripts:
				# To the exon level

				trans_exons = list(ExonInfo.objects.filter(transcript_info=trans))

				trans_exons.sort(key=lambda x: x.exon_start)

				for a_exon in trans_exons:

					for a_vcf in vcf_info_list:

						for a_var in a_vcf['vcf_info']:

							# Per variant, per gene

							if str(a_var['CHROM']) == str(gene_chrom):

								if int(a_exon.exon_start) <= int(a_var['POS']) <= int(a_exon.exon_stop):

									# POLD1_AMPL7154402731, chr19, 50 920 347 - 50 920 673

									var_pos = 'pos_' + str(a_var['POS'] + '_chr_' + str(gene_chrom))

									if var_pos not in var_list:
										var_list.append(var_pos)
										var_mutations[var_pos] = []

									if a_gene not in gene_list:
										gene_list.append(a_gene)

									# Also account for multiple samples in a vcf

									for a_sample in a_var.keys():

										# Per variant in the parsed vcf

										if a_sample not in vcf_headder_list:

											# Per sample in the vcf file

											if a_sample not in sample_list:
												sample_list.append(a_sample)

											# Also check alleles

											ref_allele = a_var['REF']
											alt_alleles = a_var['ALT'].split(',')
											all_alleles = [ref_allele] + alt_alleles
											try:
												sample_allele = a_var[a_sample]
											except KeyError:
												sample_allele = '0/0'

											if sample_allele != '0/0':

												mutation_string = ref_allele + '->' + all_alleles[int(sample_allele.split('/')[0])] + '/' + all_alleles[int(sample_allele.split('/')[1])]

												mut_simple_string_1 = ref_allele + '>' + all_alleles[int(sample_allele.split('/')[0])]
												mut_simple_string_2 = ref_allele + '>' + all_alleles[int(sample_allele.split('/')[1])]

												if mut_simple_string_1 not in var_mutations[var_pos]:
													var_mutations[var_pos].append(mut_simple_string_1)
												if mut_simple_string_2 not in var_mutations[var_pos]:
													var_mutations[var_pos].append(mut_simple_string_2)

												sample_var_dict = {
													'quality': a_var['QUAL'],
													'genotype': sample_allele,
													'mutation': mutation_string,
													'gene': a_gene,
													'info': 'info string'
												}

												if a_gene not in result_dict.keys():

														# No gene, fresh start
														result_dict[a_gene] = {var_pos: {'samples': {a_sample: sample_var_dict}, 'effect': a_var['INFO']}}
												else:
													if var_pos not in result_dict[a_gene].keys():
														# Existing gene, new variant position
														result_dict[a_gene][var_pos] = {'samples': {a_sample: sample_var_dict}, 'effect': a_var['INFO']}

													else:
														# Existing gene, existing var position, new sample
														result_dict[a_gene][var_pos]['samples'][a_sample] = sample_var_dict


		# Show sample

		df = pd.DataFrame(index=sample_list, columns=var_list)

		for a_gene in result_dict.keys():
			for a_var in result_dict[a_gene].keys():
				for a_sample in result_dict[a_gene][a_var]['samples']:
					df[a_var][a_sample] = result_dict[a_gene][a_var]['samples'][a_sample]['mutation']

		render_table = var_results_to_html_table(result_dict, sample_list)

		variant_info_dict = {
			'selected': variant,
			'panel': panel
		}

		# variation view
		# <iframe src="https://www.ncbi.nlm.nih.gov/variation/view/?chr=19&q=&assm=GCF_000001405.38&from=87863348&to=87863348" title="Var details" width="100%" height="600"></iframe>

		var_info_block_dict = {}

		for a_var in var_mutations.keys():
			var_info_block_dict[a_var] = myvariant_html(a_var, var_mutations[a_var])

		return render(request, 'variant_info.html', {
			'panels': panel_list,
			'var_table': render_table,
			'var_detail': variant_info_dict,
			'var_list': var_info_block_dict
		}
			)

	else:
		# Select which sample to look at


		return render(request, 'variant_info.html', {
			'selected_panel': 'None',
			'panels': panel_list,
			'var_detail': variant
		}
			)


def sub_graph_detail(request):

	subnet_id = request.path.rsplit('/', 1)[-1]
	subnet_list = request.session.get('subnet_data', None)

	for subnet in subnet_list:
		if subnet_id == subnet['label']:
			displayed_subnet = subnet

	return render(request, 'subnet_view.html', {'subnet_data': displayed_subnet})


def coverage_summary_old(request):

	if request.method == 'POST':

		loaded_custom_tracks = check_loaded_annotations('/annotations/')

		anno_path = anno_folder + loaded_custom_tracks['custom_annotation']

		q_gene = request.POST['a_gene_selection']
		cov_threshold = int(request.POST['coverage_threshold'])

		# Getting gene info
		target_page = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' + q_gene
		r = requests.get(target_page)
		data = r.text
		soup = BeautifulSoup(data, 'html.parser')

		mydivs = soup.find("div")
		print(mydivs)

		sample_cov_result = cov_tools.calc_gene_coverage_info(anno_path, q_gene, bam_files_dir, threshold=cov_threshold)

		exon_cot_row = []
		cot_matrix = []

		# This needs changing
		for a_sample in sample_cov_result.keys():
			rep_sample_name = a_sample

		for an_exon in sample_cov_result[rep_sample_name].keys():
			exon_cot_row.append(int(an_exon))

		exon_sample_row = sorted(exon_cot_row)

		for a_sample in sample_cov_result.keys():
			exon_cot_row = [{'cot': a_sample, 'avc': a_sample, 'cot_perc': 0}]
			for exon in exon_sample_row:
				#print(exon)
				#print(sample_cov_result[a_sample][str(exon_sample_row)]['coverage_over_threshold'])
				#print(sample_cov_result[a_sample][str(exon_sample_row)]['average_coverage'])
				cot = sample_cov_result[a_sample][str(exon)]['coverage_over_threshold']
				avc = sample_cov_result[a_sample][str(exon)]['average_coverage']
				cot_perc = cot / 100
				res_dict = {'cot': cot, 'avc': avc, 'cot_perc': cot_perc}
				exon_cot_row.append(res_dict)
			cot_matrix.append(exon_cot_row)

		exon_sample_row = ['Sample'] + exon_sample_row

		print(exon_sample_row)

		view_info = {'gene': q_gene, 'cov_threshold': cov_threshold}

		return render(request, 'coverage_summary.html', {'JSON_data': q_gene, 'exon_cot_row': exon_sample_row, 'cot_matrix': cot_matrix, 'view_info': view_info})

	else:

		loaded_custom_tracks = check_loaded_annotations('/annotations/')

		anno_obj = cov_tools.parse_annotation_file(anno_folder + loaded_custom_tracks['custom_annotation'])

		gene_list = []

		for a_gene in anno_obj.keys():
			gene_list.append(a_gene)

		return render(request, 'coverage_summary_select.html', {'gene_list': gene_list, 'track_name': loaded_custom_tracks['custom_annotation']})


def coverage_summary_gene(request):

	"""
	This is the view for a single gene
	:param request:
	:return:
	"""

	if request.method == 'POST':

		q_gene = request.POST['a_gene_selection']
		selected_panel = request.POST['panel_selection']

		# Get transcripts related to the gene

		gene_obj = GeneInfo.objects.get(gene_id=q_gene)
		panel_obj = PanelGeneList.objects.get(panel_id=selected_panel)
		cov_threshold = panel_obj.last_threshold

		panel_bam_files = panel_obj.bam_files.replace('\r\n', '').split(',')

		# Then calculate the per gene info

		## First the exons

		# Get transcript info
		transcripts = list(TranscriptInfo.objects.filter(gene_info=gene_obj))

		trans_cov_matrix = {}

		for trans in transcripts:

			# To the exon level

			trans_exons = list(ExonInfo.objects.filter(transcript_info=trans))
			trans_cds = list(CDSInfo.objects.filter(transcript_info=trans))
			trans_cov_matrix[trans.transcript_id] = {'header': [], 'samples': []}

			trans_exons.sort(key=lambda x: x.exon_start)

			for a_sample_bam in panel_bam_files:

				#a_sample_bam = a_sample_bam.split('.')[0]

				exon_sample_row = ['Sample']

				exon_cot_row = [{'cot': a_sample_bam.split('.')[0], 'avc': a_sample_bam.split('.')[0], 'cot_perc': 0}]

				exon_count = 1

				for an_exon in trans_exons:

					exon_sample_row.append(str(exon_count))

					exon_sample_cot_dict = {}

					bam_cot_list = an_exon.covOverThreshold.split(',')[:-1]

					for a_sample in bam_cot_list:
						split_ele = a_sample.split(':')
						exon_sample_cot_dict[split_ele[0]] = split_ele[1]

					print(exon_sample_cot_dict)

					exon_sample_cov_dict = {}
					bam_cov_list = an_exon.coverage.split(',')[:-1]

					for a_sample in bam_cov_list:
						split_ele = a_sample.split(':')
						exon_sample_cov_dict[split_ele[0]] = split_ele[1]

					print(exon_sample_cov_dict)

					cot_perc = float(exon_sample_cot_dict[a_sample_bam]) / 100

					res_dict = {
						'cot': float(exon_sample_cot_dict[a_sample_bam]),
						'avc': float(exon_sample_cov_dict[a_sample_bam]),
						'cot_perc': cot_perc,
						'exon_no': exon_count
					}

					exon_cot_row.append(res_dict)

					exon_count += 1

				trans_cov_matrix[trans.transcript_id]['samples'].append(exon_cot_row)

				trans_cov_matrix[trans.transcript_id]['header'] = exon_sample_row

		print(trans_cov_matrix)

		view_info = {'gene': q_gene, 'cov_threshold': cov_threshold}

		return render(
			request,
			'coverage_summary.html',
			{
				'JSON_data': q_gene,
				'trans_cov_matrix': trans_cov_matrix,
				'view_info': view_info
			}
		)

	else:

		gene_panels = PanelGeneList.objects.values_list('panel_id', flat=True)

		gene_list = list(GeneInfo.objects.all())

		return render(request, 'coverage_summary_select.html', {'gene_list': gene_list, 'gene_panels': gene_panels})


def coverage_summary_sample(request):

	"""
	This is the view for looking at the average coverage for all samples.
	:param request:
	:return:
	"""

	if request.method == 'POST':

		# q_gene = request.POST['a_gene_selection']
		#cov_threshold = int(request.POST['coverage_threshold'])

		anno_set_selection = request.POST['set_selection']

		selected_level = request.POST['set_level']

		panel_obj = PanelGeneList.objects.get(panel_id=anno_set_selection)

		gene_list = panel_obj.gene_list.replace('\r\n', '').split(',')

		table_matrix = {}

		# Get a list of samples to iterate over

		bam_file_list = panel_obj.bam_files.replace('\r\n', '').split(',')

		sample_list = []

		for filename in bam_file_list:
			sample_list.append(filename.split('.')[0])

		for a_gene in gene_list:

			gene_cov_result = cov_tools.calc_average_gene_coverage_info(a_gene, bam_file_list, level=selected_level)

			print(gene_cov_result)

			table_matrix[a_gene] = []

			for a_trans in gene_cov_result.keys():

				# This is the first column of the row
				gene_cot_row = [{'cot': a_trans, 'avc': a_trans, 'cot_perc': 0}]

				print(a_trans)

				sample_list = []

				for filename in bam_file_list:
					sample_list.append(filename.split('.')[0])

				for a_sample in sample_list:

					print(a_sample)

					sample_cot_row = [{'cot': gene_cov_result[a_trans][a_sample]['percentage_over_cov_thres'], 'avc': gene_cov_result[a_trans][a_sample]['average_coverage'], 'cot_perc': gene_cov_result[a_trans][a_sample]['percentage_over_cov_thres_perc']}]

					gene_cot_row = gene_cot_row + sample_cot_row

				table_matrix[a_gene].append(gene_cot_row)

		cov_threshold = panel_obj.last_threshold

		sample_list = ['Transcript'] + sample_list

		view_info = {'cov_threshold': cov_threshold, 'anno_file': anno_set_selection}

		print(table_matrix)

		return render(request, 'coverage_summary_samples.html', {'gene_list': gene_list, 'sample_list': sample_list, 'cot_matrix': table_matrix, 'view_info': view_info})

	else:

		#anno_file_list = [f for f in listdir(anno_folder) if isfile(join(anno_folder, f)) and f[-3:] == 'bed']
		#anno_file_list = anno_file_list + [f for f in listdir(anno_folder) if isfile(join(anno_folder, f)) and f[-3:] == 'gff']

		gene_set_list = list(PanelGeneList.objects.all())

		gene_set_objs = PanelGeneList.objects.all()

		passed_sets = []

		for a_set in gene_set_objs:
			if a_set.last_threshold != -1:
				passed_sets.append(a_set)

		return render(request, 'coverage_summary_select_samples.html', {'gene_set_list': passed_sets})


def view_region(request):

	return render(request, 'view_read_align_region.html', {'gene_list': ''})


def create_subset_file(request):

	loaded_genomes = check_loaded_genomes(ref_genome_folder)
	loaded_custom_tracks = check_loaded_annotations(anno_folder)

	anno_file_path = loaded_genomes['annotation']
	tar_feat_path = loaded_custom_tracks['gene_list']

	extract_gff_subset(ref_genome_folder + anno_file_path, anno_folder + tar_feat_path, anno_folder)

	return redirect('/')


def gene_list_to_db(request):

	if request.method == 'POST':

		panel_val = request.POST['panel_val']

		in_file_path = ref_genome_folder + ref_genome_gff

		import_result = subset_gff_to_db(in_file_path, panel_val)

		return redirect('/')

	else:

		gene_panels = PanelGeneList.objects.values_list('panel_id', flat=True)

		return render(request, 'gene_list_2_db.html', {'gene_panels': gene_panels})


def process_sample_coverage(request):

	if request.method == 'POST':
		# Extract info

		panel_val = request.POST['panel_val']

		cutoff_val = int(request.POST['cutoff_val'])

		print('Processing ' + panel_val)

		# First get the list of genes

		panel_info = PanelGeneList.objects.get(panel_id=panel_val)

		panel_info.last_threshold = cutoff_val
		panel_info.save()

		gene_list = panel_info.gene_list.replace('\r\n', '').split(',')
		bam_list = panel_info.bam_files.replace('\r\n', '').split(',')

		gene_obj_list = []

		for a_gene in gene_list:
			gene_obj_list += list(GeneInfo.objects.filter(gene_id=a_gene))

		# Then calculate the per gene info

		## First the exons

		for a_gene in gene_obj_list:
			gene_chrom = a_gene.gene_chrom
			# Get transcript info
			transcripts = list(TranscriptInfo.objects.filter(gene_info=a_gene))
			#print('-------')
			#print(transcripts)
			#print('-------')

			for trans in transcripts:

				#print(trans)

				# To the exon level

				trans_exons = list(ExonInfo.objects.filter(transcript_info=trans))
				trans_cds = list(CDSInfo.objects.filter(transcript_info=trans))

				#print(trans_exons)
				## First the exons
				print('Processing Exons')

				for an_exon in trans_exons:

					cov_list = ''
					cov_OverThres = ''
					base_cov = ''

					for a_bam in bam_list:

						bam_file = bam_files_dir + a_bam

						cov_res = cov_tools.get_region_coverage_info(
							bam_file,
							an_exon.exon_start,
							an_exon.exon_stop,
							gene_chrom,
							threshold=cutoff_val
						)

						cov_list += a_bam + ':' + str(cov_res['average_coverage']) + ','
						cov_OverThres += a_bam + ':' + str(cov_res['coverage_over_threshold']) + ','
						base_cov = base_cov + a_bam + ':'
						for a_base in cov_res['coverage']:
							base_cov = base_cov + str(a_base) + ','
						base_cov = base_cov[:-1]
						base_cov += ';'

					an_exon.coverage = cov_list
					an_exon.covOverThreshold = cov_OverThres
					an_exon.pb_coverage = base_cov[:-1]
					an_exon.selected_cov = cutoff_val
					an_exon.save()

					## Then the CDS
				print('Processing CDS')
				for a_cds in trans_cds:

					cov_list = ''
					cov_OverThres = ''
					base_cov = ''

					for a_bam in bam_list:

						bam_file = bam_files_dir + a_bam

						cov_res = cov_tools.get_region_coverage_info(
							bam_file,
							a_cds.cds_start,
							a_cds.cds_stop,
							gene_chrom,
							threshold=cutoff_val
						)

						cov_list += a_bam + ':' + str(cov_res['average_coverage']) + ','
						cov_OverThres += a_bam + ':' + str(cov_res['coverage_over_threshold']) + ','
						for a_base in cov_res['coverage']:
							base_cov = base_cov + str(a_base) + ','
						base_cov = base_cov[:-1]
						base_cov += ';'

					#print(cov_list)
					#print(cov_OverThres)
					a_cds.coverage = cov_list
					a_cds.covOverThreshold = cov_OverThres
					a_cds.pb_coverage = base_cov[:-1]
					a_cds.selected_cov = cutoff_val
					a_cds.save()

			# Then the per sample info

		return redirect('/')

	else:

		gene_panels = PanelGeneList.objects.values_list('panel_id', flat=True)

		return render(request, 'process_bam.html', {'gene_panels': gene_panels})


def gene_file_to_db(request):
	if request.method == 'POST':

		# Process genes to db

		selected_file = request.POST['anno_file_selection']
		selected_bam = request.POST.getlist('bamFiles')
		selected_vcf = request.POST.getlist('vcfFiles')

		if selected_file.split('.')[-1].lower() == 'bed':

			file_obj = open(custom_anno_folder + selected_file)

			gene_list = []

			for line in file_obj:

				info = line.split('\t')

				chrom = info[0]
				start = info[1]
				stop = info[2]
				feature = info[3]
				exon = info[4]
				strand = info[5]

				gene_list.append(feature)

				if not GeneInfo.objects.filter(gene_id=feature).exists():

					gene_obj = GeneInfo(
						gene_id=feature,
						gene_name=feature,
						gene_description='custom bed annotation',
						gene_chrom=chrom,
						gene_strand=strand
					)
					gene_obj.save()

				gene_obj = GeneInfo.objects.get(gene_id=feature)

				if not TranscriptInfo.objects.filter(transcript_id=feature).exists():

					transcript_obj = TranscriptInfo(
						gene_info=gene_obj,
						transcript_id=feature
					)

					transcript_obj.save()

				transcript_obj = TranscriptInfo.objects.get(transcript_id=feature)

				if not ExonInfo.objects.filter(exon_id=feature + '_' + exon).exists():

					new_exon = ExonInfo(
						exon_id=feature + '_' + exon,
						exon_start=start,
						exon_stop=stop,
						transcript_info=transcript_obj
					)
					print('exon creates')
					new_exon.save()


		# Create the PanelGeneList entry

		new_panel_obj = PanelGeneList(
			panel_id=request.POST['panelName'],
			phenotype=request.POST['phenotype'],
			gene_list=','.join(gene_list),
			bam_files=','.join(selected_bam),
			vcf_files=','.join(selected_vcf),
		)
		new_panel_obj.save()

		return redirect('/')

	else:

		custom_anno_files = []

		all_custom_anno_files = os.listdir(custom_anno_folder)

		all_bam_files = os.listdir(bam_files_dir)
		clean_bam_list = []

		all_vcf_files = os.listdir(vcf_folder)
		clean_vcf_list = []

		allowed_files = ['gff', 'bed', 'gff3']

		for file in all_custom_anno_files:
			if file.split('.')[-1].lower() in allowed_files:
				custom_anno_files.append(file)

		for a_bam in all_bam_files:
			if a_bam[-3:].lower() == 'bam':
				clean_bam_list.append(a_bam)

		for a_vcf in all_vcf_files:
			if a_vcf[-3:].lower() == 'vcf':
				clean_vcf_list.append(a_vcf)

		return render(request, 'custom_gene_set_upload.html', {'anno_file_list': custom_anno_files, 'bam_list': clean_bam_list, 'vcf_list': clean_vcf_list})


def panel_detail_view(request):

	panel_list = list(PanelGeneList.objects.all())

	if request.method == 'POST':

		selected_panel = request.POST['panel_selection']

		panel_info = PanelGeneList.objects.get(panel_id=selected_panel)

		gene_list = panel_info.gene_list.replace('\r\n', '').split(',')
		bam_list = panel_info.bam_files.replace('\r\n', '').split(',')

		gene_obj_list = []

		for a_gene in gene_list:
			gene_obj_list += list(GeneInfo.objects.filter(gene_id=a_gene))

		print(gene_list)
		print(bam_list)

		gene_status_list = []

		for a_gene in gene_list:

			res_dict = {'gene_id': a_gene}

			# look in DB

			try:
				gene_obj = GeneInfo.objects.get(gene_id=a_gene)
			except:
				res_dict['db'] = 'Not found'
			else:
				res_dict['db'] = 'Found'

			if res_dict['db'] == 'Found':

				transcripts = list(TranscriptInfo.objects.filter(gene_info=gene_obj))

				res_dict['trans_no'] = len(transcripts)

			else:
				res_dict['trans_no'] = 0

			gene_status_list.append(res_dict)

		return render(request, 'panel_details.html', {'panel_list': panel_list, 'panel_name': selected_panel, 'gene_info': gene_status_list, 'bam_list': bam_list})

	else:

		selected_panel = 'none'

		return render(request, 'panel_details.html', {'panel_list': panel_list, 'panel_name': selected_panel})


def add_experiment(request):

	if request.method == 'POST':

		if request.POST['timeSeries'] == 'True':
			time_series = True
		else:
			time_series = False

		if not Experiment.objects.filter(experiment_shortname=request.POST['experimentName']).exists():

			experiment_obj = Experiment(
				experiment_shortname=request.POST['experimentName'],
				comparisons=request.POST['comparison'],
				time_series=time_series,
			)
			experiment_obj.save()

		return redirect('/')

	else:

		return render(request, 'create_experiment.html')


def add_expression_data(request):

	if request.method == 'POST':

		experiment_name = request.POST['experimentName']
		expression_comparison = request.POST['comparison']
		expression_file = request.POST['sourceFile']

		experiment_obj = Experiment.objects.get(experiment_shortname=experiment_name)

		exp_set_obj = ExpressionSet(
			comparison=expression_comparison,
			source_file=expression_file,
			experiment=experiment_obj
		)

		exp_set_obj.save()

		path_to_expression_file = bam_files_dir + expression_file

		expression_info = expression_file_parser(path_to_expression_file)

		expression_info['exp_set_obj'] = exp_set_obj

		expression_info.apply(save_expression_data, axis=1)



		'''

		experiment_obj = Experiment(
			experiment_shortname=request.POST['experimentName'],
			comparisons=request.POST['comparison'],
			time_series=time_series,
		)
		experiment_obj.save()
		'''

		return redirect('/')

	else:

		# Get list of existing experiments to link

		experiment_obj_list = list(Experiment.objects.all())
		experiment_list = []

		for thing in experiment_obj_list:
			experiment_list.append(thing)

		# Get a list of the available expression files

		all_expression_files = os.listdir(bam_files_dir)
		clean_expression_list = []

		allowed_files = ['csv']

		for file in all_expression_files:
			if file.split('.')[-1].lower() in allowed_files:
				clean_expression_list.append(file)

		print(clean_expression_list)

		return render(request, 'add_expression_data.html', {
			'expression_file_list': clean_expression_list,
			'experiments': experiment_list
			}
		)


def help_gene_input(request):
	return render(request, "helppages/gene_input_help.html")

