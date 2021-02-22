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
				int_data = pd.read_table(in_file, index_col=0 )
				ccs = list(int_data.columns)
				correlation = round(int_data[ccs[0]].corr(int_data[ccs[1]]),3)
				ag = sns.regplot(x = int_data[ccs[0]], y = int_data[ccs[1]], fit_reg=False)
				ag.set_xlim(1,)
				ag.set_ylim(1,)
				ag.set(xscale="log", yscale="log")
				plt.title("Pearson's r = " + str(correlation))
				plt.savefig(pdf_name)
				# plt.show()
				plt.close()

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
graph_scatterplot(cell_class_dict, d_values_wanted, size_values_wanted, bed_suffix)



