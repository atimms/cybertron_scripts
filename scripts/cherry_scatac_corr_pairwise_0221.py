#!/usr/bin/env python
import subprocess
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

'''
info....

##for macs2.. not used in this script
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load local_python/3.7.6 
source activate macs2

##homer
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da87
module load biobuilds
module load homer/4.9.1

##for graphing
qsub -Iq cdbrmq -l mem=40gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
conda activate pd_np_plt_etc


'''

##parameters
delim = '\t'
thread_number = '20'
##programs
sinto = '/home/atimms/programs/sinto/scripts/sinto'
genrich = '/home/atimms/programs/Genrich/Genrich'
bedtools = '/home/atimms/programs/bedtools2.28/bin/bedtools'
bg_to_bw = '/home/atimms/programs/bedGraphToBigWig'


##methods
def merge_peaks(out_prefix, d_value, peak_files, out_suffix):
	out_file = out_prefix + d_value + 'd' + out_suffix
	print(len(peak_files), out_file)
	with open(out_file, 'w') as out_fh:
		# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
		run_merge_peaks = subprocess.Popen(['mergePeaks', '-d', d_value] + peak_files, stdout=out_fh)
		run_merge_peaks.wait()

def annotate_peaks(samples, peak_file, d_value, tag_suffix, out_prefix, sizes_req, peak_suffix):
	tag_dirs = []
	##get list of tag_dirs
	for sample in samples:
		tag_dir = sample + tag_suffix
		tag_dirs.append(tag_dir)
	print(len(tag_dirs), peak_file)
	##annotate
	for size_req in sizes_req:
		out_file_using_size = out_prefix + d_value + 'd.' + size_req + 'size' + peak_suffix
		with open(out_file_using_size, 'w') as out_fh:
			mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-size', size_req, '-d'] + tag_dirs, stdout=out_fh)
			mk_tag_dir.wait()

def get_correlation_info_homer(cc_dict, d_values, size_values, bed_suffix, tag_dir_suffix):
	for cell_class in cc_dict:
		samples = cc_dict[cell_class]
		##combine bed files
		merge_prefix = cell_class + '_merged.'
		annotate_prefix = cell_class + '_annotated.'
		annotate_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
		corr_prefix = cell_class + '_correlation.'
		'''
		comb_beds = []
		for s in samples:
			bed = s + bed_suffix
			comb_beds.append(bed)
		# for d in d_values:
		# 	merge_peaks(merge_prefix, d, comb_beds, bed_suffix)
		
		##annotate files
		tag_dir_suffix = '.tag_dir'
		# correlation_prefix = 'correlation.'
		for d in d_values:
			merged_bed = merge_prefix + d + 'd' + bed_suffix
			
			##annotate those peaks, and then make a correlation file
			# annotate_peaks(samples, merged_bed, d, tag_dir_suffix, annotate_prefix, size_values, annotate_suffix)
		'''
		##format annotation files
		for size_req in size_values:
			for d_value in d_values:
				infile = annotate_prefix + d_value + 'd.' + size_req + 'size' + annotate_suffix
				out_file = corr_prefix + infile.split('.', 1)[1]
				# print(infile, out_file)
				with open(out_file, 'w') as out_fh, open(infile, 'r') as in_fh:
					lc = 0
					for line in in_fh:
						lc += 1
						line = line.strip('\n').split(delim)
						if lc == 1:
							sample_info = line[19:]
							sample_info = [s.split('.')[0] + '.' + s.split('.')[1] for s in sample_info]
							out_fh.write(delim.join(['peak_id'] + sample_info) + '\n')
							# print(sample_info)
						else:
							chrom = line[1]
							if '_' not in chrom and chrom != 'chrM' and chrom != 'chrY':
								human_count = float(line[19])
								org_count = float(line[20])
								# if min([human_count, org_count]) > 1:
								line_out = [line[0]] + line[19:]
								out_fh.write(delim.join(line_out) + '\n')

def graph_scatterplot(cc_dict, d_values, size_values, bed_suffix):
	for cell_class in cc_dict:
		corr_prefix = cell_class + '_correlation.'
		annotate_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
		for size_req in size_values:
			for d_value in d_values:
				in_file = corr_prefix + d_value + 'd.' + size_req + 'size' + annotate_suffix
				pdf_name = in_file.rsplit('.', 1)[0] + '.png'
				# pdf_name = in_file.rsplit('.', 1)[0] + '.svg'
				int_data = pd.read_table(in_file, index_col=0 )
				ccs = list(int_data.columns)
				correlation = round(int_data[ccs[0]].corr(int_data[ccs[1]]),3)
				##without reg line
				# ag = sns.regplot(x = int_data[ccs[0]], y = int_data[ccs[1]], fit_reg=False, scatter_kws={'s':2})
				ag = sns.regplot(x = int_data[ccs[0]], y = int_data[ccs[1]], ci=None, 
					scatter_kws={'s':2}, line_kws={"color": "red"})
				ag.set_xlim(1,)
				ag.set_ylim(1,)
				# ag.set(xscale="log", yscale="log")
				plt.title("Pearson's r = " + str(correlation))
				plt.savefig(pdf_name)
				# plt.show()
				plt.close()

def filter_split_peaks_for_great(cc_dict, d_values, size_values, bed_suffix, tag_count_req, fold_change_req):
	for cell_class in cc_dict:
		samples = cc_dict[cell_class]
		##combine bed files
		annotate_prefix = cell_class + '_annotated.'
		annotate_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
		for size_req in size_values:
			for d_value in d_values:
				cs = []
				human_peak_count, org_peak_count, shared_peak_count = 0,0,0
				annotate_file = annotate_prefix + d_value + 'd.' + size_req + 'size' + annotate_suffix
				great_peak_file = annotate_prefix + d_value + 'd.' + size_req + 'size.' + str(tag_count_req) + 'min.' + str(fold_change_req) + 'fc.peaks.txt'
				great_peak_human_bed = great_peak_file.rsplit('.',2)[0] + '.human_peaks.bed'
				great_peak_org_bed = great_peak_file.rsplit('.',2)[0] + '.organoid_peaks.bed'
				great_peak_shared_bed = great_peak_file.rsplit('.',2)[0] + '.shared_peaks.bed'
				all_peak_bed = great_peak_file.rsplit('.',2)[0] + '.all_peaks.bed'
				with open(great_peak_file, 'w') as out_fh, open(annotate_file, 'r') as in_fh, open(great_peak_human_bed, 'w') as gphb_fh, open(great_peak_org_bed, 'w') as gpob_fh, open(great_peak_shared_bed, 'w') as gpsb_fh, open(all_peak_bed, 'w') as apb_fh :
					lc = 0
					for line in in_fh:
						lc += 1
						line = line.strip('\n').split(delim)
						if lc == 1:
							out_fh.write(delim.join(line + ['group']) + '\n')
						else:
							human_count = float(line[19])
							org_count = float(line[20])
							chrom = line[1]
							cse = line[1:4] + [line[0]]
							if '_' not in chrom and chrom != 'chrM' and chrom != 'chrY':
								# if chrom not in cs:
								# 	cs.append(chrom)
								##sufficent coverage
								if max([human_count, org_count]) > tag_count_req:
									##write to background file
									# print(cse, org_count, human_count)
									apb_fh.write(delim.join(cse) + '\n')
									##if enrichred in a type
									if org_count == 0.0 or (human_count / org_count > fold_change_req):
										human_peak_count += 1
										gphb_fh.write(delim.join(cse) + '\n')
										out_fh.write(delim.join(line + ['human']) + '\n')
									elif human_count == 0.0 or (org_count / human_count > fold_change_req):
										org_peak_count += 1
										gpob_fh.write(delim.join(cse) + '\n')
										out_fh.write(delim.join(line + ['organoid']) + '\n')
									elif (org_count / human_count < 1.1) and (human_count / org_count < 1.1):
										shared_peak_count += 1
										gpsb_fh.write(delim.join(cse) + '\n')
										out_fh.write(delim.join(line + ['shared']) + '\n')

			print(annotate_file, lc, human_peak_count, org_peak_count, shared_peak_count)
			# print(cs)

def homer_find_motifs_from_beds(cc_dict, d_values, size_values, bed_types):
	for cell_class in cc_dict:
		samples = cc_dict[cell_class]
		##combine bed files
		annotate_prefix = cell_class + '_annotated.'
		annotate_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
		for size_req in size_values:
			for d_value in d_values:
				for bed_type in bed_types:
					in_bed_file = annotate_prefix + d_value + 'd.' + size_req + 'size.10min.5fc.' + bed_type + '_peaks.bed'
					out_dir = in_bed_file.rsplit('_', 1)[0] + '.find_motifs'
					print(in_bed_file)
					#findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
					find_motifs = subprocess.Popen(['findMotifsGenome.pl', in_bed_file, 'hg38', out_dir, '-size', 'given', '-preparsedDir', './temp4'])
					find_motifs.wait()



##run methods

##params etc
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scatac_corr_pairwise_0221'
os.chdir(working_dir)

cell_class_dict = {'cones':['human.Mature_Cones', 'organoid.Cones'], 'rods':['human.Mature_Rods', 'organoid.Rods'],
		'bipolars':['human.Mature_Bipolars', 'organoid.Bipolar_Cells'], 'mullers':['human.Mature_Mullers', 'organoid.Muller_Glia'],
		'horizontals':['human.Mature_Horizontals', 'organoid.Horizontal_Cells'], 'amacrines':['human.Mature_Amacrines', 'organoid.Amacrine_Cells'],
		'early_progenitor':['human.Early_Progenitors', 'organoid.Early_Progenitors'], 'late_progenitor':['human.Late_Progenitors', 'organoid.Late_Progenitors']}

# cell_class_dict = {'cones':['human.Mature_Cones', 'organoid.Cones']}

d_values_wanted = ['100', '500']
size_values_wanted = ['500', '2000'] ##does none work for homer
bed_suffix = '.macs2_q0.01_summits.bed'
tag_dir_suffix = '.tag_dir'
# get_correlation_info_homer(cell_class_dict, d_values_wanted, size_values_wanted, bed_suffix, tag_dir_suffix)

##and graph
# graph_scatterplot(cell_class_dict, d_values_wanted, size_values_wanted, bed_suffix)

##great analysis
tag_count_wanted = 10
fold_change_wanted = 5
# d_values_wanted = ['500']
# size_values_wanted = ['2000'] ##does none work for homer
# cell_class_dict = {'rods':['human.Mature_Rods', 'organoid.Rods']}
# filter_split_peaks_for_great(cell_class_dict, d_values_wanted, size_values_wanted, bed_suffix, tag_count_wanted, fold_change_wanted)

##motif enruchment

size_values_motif = ['500']
bed_types = ['human', 'organoid', 'shared']
# cell_class_dict = {'cones':['human.Mature_Cones', 'organoid.Cones'], 'rods':['human.Mature_Rods', 'organoid.Rods']}
# cell_class_dict = {'bipolars':['human.Mature_Bipolars', 'organoid.Bipolar_Cells'], 'mullers':['human.Mature_Mullers', 'organoid.Muller_Glia']}
# cell_class_dict = {'horizontals':['human.Mature_Horizontals', 'organoid.Horizontal_Cells'], 'amacrines':['human.Mature_Amacrines', 'organoid.Amacrine_Cells']}
cell_class_dict = {'early_progenitor':['human.Early_Progenitors', 'organoid.Early_Progenitors'], 'late_progenitor':['human.Late_Progenitors', 'organoid.Late_Progenitors']}



homer_find_motifs_from_beds(cell_class_dict, d_values_wanted, size_values_motif, bed_types)
