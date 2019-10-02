

#in_bed_path = '/Users/panix/Library/Mobile Documents/com~apple~CloudDocs/programs/Holmes/annotations/test_1_renamed.bed'
in_bed_path = '/Users/panix/Library/Mobile Documents/com~apple~CloudDocs/programs/Holmes/annotations/CRC_panel_designed.bed'


new_file = open(in_bed_path[:-4] + '_clean.bed', 'w')

chrom_rename = True

case = 2

bed_file = open(in_bed_path, 'r')


exon_count_dict = {}

for line in bed_file:

	info = line.split('\t')

	tab_val = '\t'

	if case == 1:
		chrom = info[0]
		start = info[1]
		stop = info[2]
		strand = info[5]

		feature_name = info[3].split('_')[0]

		exon = info[3].split('_')[3]

		out_line = chrom + tab_val + start + tab_val + stop + tab_val + feature_name + tab_val + exon + tab_val + strand

		new_file.write(out_line)

	if case == 2:

		if chrom_rename:
			chrom = info[0][3:]
			start = info[1]
			stop = info[2]
			strand = info[5]

			feature_name = info[3].split('_')[0]

			if feature_name not in exon_count_dict.keys():
				exon_count_dict[feature_name] = 1

			else:
				exon_count_dict[feature_name] = exon_count_dict[feature_name] + 1

			exon = str(exon_count_dict[feature_name])

			out_line = chrom + tab_val + start + tab_val + stop + tab_val + feature_name + tab_val + exon + tab_val + strand

			new_file.write(out_line)

