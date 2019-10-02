# extract_sub_gtf

tar_feat_file = open('target_genes.txt', 'r')

search_gene_list =  []


for line in tar_feat_file:
	if len(line) > 1:
		search_gene_list.append(line.strip())


anno_file = open('/Volumes/External/Genomes/human/Homo_sapiens.GRCh37.87.chr.gff3', 'r')

out_anno_file = open("result_refined.gtf", 'w')


for line in anno_file:

	if line[0] == '#' and line[2] != '#':
		out_anno_file.write(line)

for a_feature in search_gene_list:

	for line in anno_file:


		line_list = line.split('\t')

		feature = line_list[2]

		info_list = line_list[8].split(';')

		if feature == 'gene':

		for a_feature in search_gene_list:

			if a_feature in line:
				out_anno_file.write(line)

