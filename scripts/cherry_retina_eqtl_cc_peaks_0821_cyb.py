#!/usr/bin/env python
import subprocess
import os
import glob
'''
info:
qsub -Iq cdbrmq -l mem=20gb,ncpus=1 -P 19833a08-f6fb-4bea-8526-8a79069da878
'''

##parameters
delim = '\t'
##programs
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
##files etc
bt_genome = '/home/atimms/ngs_data/references/10x/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.chrom_sizes'


##methods
def make_eqtl_bed(q_infile, all_infile, out_bed):
	qvalue_dict, eqtl_dict = {}, {}
	with open(q_infile, "r") as qin_fh:
		lc = 0
		for qline in qin_fh:
			# print(line)
			qline = qline.rstrip().split(delim)
			rs_number = qline[0]
			pval = qline[2]
			gene = qline[11]
			biotype = qline[14]
			if rs_number in qvalue_dict:
				qvalue_dict[rs_number][0].append(gene)
				qvalue_dict[rs_number][1].append(pval)
				qvalue_dict[rs_number][2].append(biotype)
			else:
				qvalue_dict[rs_number] = [[gene], [pval], [biotype]]


	with open(all_infile, "r") as ain_fh:
		for line in ain_fh:
			# print(line)
			line = line.rstrip().split(delim)
			biotype = line[14]
			rs_number = line[0]
			pos = '_'.join(line[6:8])
			##remove header lines
			if biotype != 'biotype':
				if rs_number not in eqtl_dict:
					if rs_number in qvalue_dict:
						eqtl_dict[rs_number] = [pos] + qvalue_dict[rs_number]
					else:
						eqtl_dict[rs_number] = [pos, [], [], []]
	vc, evc = 0, 0
	with open(out_bed, "w") as out_fh:
		for rs in eqtl_dict:
			##get counts and values
			vc += 1
			cse = ['chr' + eqtl_dict[rs][0].split('_')[0], str(int(eqtl_dict[rs][0].split('_')[1]) - 1), eqtl_dict[rs][0].split('_')[1]]
			p_counts = len(eqtl_dict[rs][2])

			if p_counts > 0:
				evc += 1
				genes = '_'.join(eqtl_dict[rs][1])
				pvals = '_'.join(eqtl_dict[rs][2])
				biotypes = '_'.join(eqtl_dict[rs][3])
			else:
				genes = '.'
				pvals = '.'
				biotypes = '.'
			out_fh.write(delim.join(cse + [rs, str(p_counts), genes, pvals, biotypes]) + '\n')
	print(vc, evc)
		



def bt_intersect_snps_bed(snp_bed, peak_beds, out_bed):
	##merge sort eqtl bed
	# temp_bed = 'temp1.bed'
	# with open(temp_bed, "w") as temp_fh:
	# 	sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', snp_bed], stdout=subprocess.PIPE)
	# 	bt_merge = subprocess.Popen([bedtools, 'merge', '-i', '-'], stdin=sort_file.stdout, stdout=temp_fh)
	# 	bt_merge.wait()
	##bedtools intersect 
	with open(out_bed, "w") as out_fh: 
		# hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', region_bed, '-b'] + summit_peak_beds + ['-wa', '-wb', '-filenames'], stdout=out_fh)
		# hom_bt_intersect.wait()
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', snp_bed, '-b'] + peak_beds + ['-C', '-filenames'], stdout=out_fh)
		hom_bt_intersect.wait()


def cat_merge_narrowpeaks(in_beds, out_bed):
	t1_bed = 'temp.a.a.a.bed'
	#combine beds
	with open(t1_bed, 'w') as t1_fh:
		cat_beds = subprocess.Popen(['cat'] + in_beds, stdout=t1_fh)
		cat_beds.wait()
	##sort and merge
	with open(out_bed, "w") as out_fh:
		sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', t1_bed], stdout=subprocess.PIPE)
		bt_merge = subprocess.Popen([bedtools, 'merge', '-i', '-'], stdin=sort_file.stdout, stdout=out_fh)
		bt_merge.wait()


def get_counts_per_cell_class(infile, outfile):
	count_dict = {'all': [0,0]}
	with open(infile, "r") as in_fh:
		for line in in_fh:
			# print(line)
			line = line.rstrip().split(delim)
			cell_class = line[8].split('.')[1] + '_' + line[8].split('.')[4]
			intersection = int(line[9])
			evar = line[4]
			##add data for all
			if cell_class == 'Mature_combined_narrowPeak':
				count_dict['all'][0] += 1
				if evar != '0':
					count_dict['all'][1] += 1
			if intersection > 0:
				if cell_class in count_dict:
					count_dict[cell_class][0] += 1
					if evar != '0':
						count_dict[cell_class][1] += 1
				else:
					if evar == '0':
						count_dict[cell_class] = [1,0]
					else:
						count_dict[cell_class] = [1,1]
	for cc in count_dict:
		print(cc, count_dict[cc])
	with open(outfile, "w") as out_fh:
		out_fh.write(delim.join(['Cell Class', 'Variants', 'eVariants', 'Percentage']) + '\n')
		for cc in count_dict:
			percentage = (float(count_dict[cc][1]) / float(count_dict[cc][0])) * 100
			out_fh.write(delim.join([cc, str(count_dict[cc][0]), str(count_dict[cc][1]), str(round(percentage, 3))]) + '\n')


def bt_slop_summit_beds(in_beds, slop_sizes):
	for slop_size in slop_sizes:
		for in_bed in in_beds:
			##bedtools slop
			out_bed = in_bed.rsplit('.', 1)[0] + '.ext' + slop_size + '.bed'
			with open(out_bed, "w") as out_fh: 
				hom_bt_intersect = subprocess.Popen([bedtools, 'slop', '-i', in_bed, '-g', bt_genome, '-b', slop_size], stdout=out_fh)
				hom_bt_intersect.wait()


def make_p2g_bed(in_file, out_file):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			# print(line)
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				cse = line[6].replace('-', ':').split(':')
				gene = line[7]
				out_fh.write(delim.join(cse + [gene]) + '\n')


def make_evars_bed(in_file, out_file):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			# print(line)
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				chrom = 'chr' + line[6]
				start = str(int(line[7]) - 1)
				end = str(int(line[7]))
				gene = line[11]
				out_fh.write(delim.join([chrom, start, end, gene]) + '\n')

def bt_intersect_snp_p2g_beds(snp_bed, peak_bed, out_bed):
	##bedtools intersect 
	with open(out_bed, "w") as out_fh: 
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', peak_bed, '-b', snp_bed, '-wa', '-wb', '-filenames'], stdout=out_fh)
		hom_bt_intersect.wait()

##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_retina_eqtl_cc_peaks_0821'
os.chdir(working_dir)

##params
##peaks copied from /active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/human_cell_class_macs2_0621
cell_class_narrowpeaks = ['human.Mature_Amacrines.macs2_q0.000001_peaks.narrowPeak', 'human.Mature_Bipolars.macs2_q0.000001_peaks.narrowPeak', 
		'human.Mature_Cones.macs2_q0.000001_peaks.narrowPeak', 'human.Mature_Horizontals.macs2_q0.000001_peaks.narrowPeak', 
		'human.Mature_Mullers.macs2_q0.000001_peaks.narrowPeak', 'human.Mature_Rods.macs2_q0.000001_peaks.narrowPeak']
cell_class_combined_narrowpeak = 'human.Mature_combined.macs2_q0.000001_peaks.narrowPeak'
all_cell_class_narrowpeaks = cell_class_narrowpeaks + [cell_class_combined_narrowpeak]
cell_class_summit_beds = ['human.Mature_Amacrines.macs2_q0.000001_summits.bed', 'human.Mature_Bipolars.macs2_q0.000001_summits.bed', 
		'human.Mature_Cones.macs2_q0.000001_summits.bed', 'human.Mature_Horizontals.macs2_q0.000001_summits.bed', 
		'human.Mature_Mullers.macs2_q0.000001_summits.bed', 'human.Mature_Rods.macs2_q0.000001_summits.bed'] 
# slop_sizes = ['25','50','75','100','125','150','175','200','225','250']
slop_sizes = ['100', '200', '300', '400', '500', '600', '700', '800', '900', '1000']
all_cell_class_files = all_cell_class_narrowpeaks + glob.glob('human*ext*bed')
eqtl_q_txt = 'Retina_merged3_hg38_FastQTL_eVariants_chrALL_FDR_0.05.txt'
eqtl_all_txt = 'Retina_merged3_hg38_FastQTL_nominal_combined.txt'
eqtl_bed = 'Retina_merged3_hg38_FastQTL_nominal_combined.bed'
eqtl_cc_peaks_bed = 'Retina_merged3_hg38_FastQTL.human_mature_peaks.bed'
eqtl_cc_counts_file = 'retina_eqtl_human_mature_peaks_counts.xls'


##1. make bed from eqtl data
# make_eqtl_bed(eqtl_q_txt, eqtl_all_txt, eqtl_bed)

##2. combine peak data, make diffent peak sizes and intersect with eqtl bed
# cat_merge_narrowpeaks(cell_class_narrowpeaks,cell_class_combined_narrowpeak)
# bt_slop_summit_beds(cell_class_summit_beds, slop_sizes)
# bt_intersect_snps_bed(eqtl_bed, all_cell_class_files, eqtl_cc_peaks_bed)

##3. get numbers of vars and evars
# get_counts_per_cell_class(eqtl_cc_peaks_bed, eqtl_cc_counts_file)


##4. compare evars to peak2gene data
human_peak2gene_file = 'human_0621.peak2gene.all.xls'
human_peak2gene_bed = 'human_0621.peak2gene.bed'
eqtl_p2g_bed = 'retina_eqtl_evars.peak2gene.intersect.bed'
evars_bed = 'retina_eqtl_evars.bed'
##make bed from peak2 gene data
# make_p2g_bed(human_peak2gene_file, human_peak2gene_bed)
##make bed from evars
# make_evars_bed(eqtl_q_txt, evars_bed)
##intersect with all eqtl snps
bt_intersect_snp_p2g_beds(evars_bed, human_peak2gene_bed, eqtl_p2g_bed)
