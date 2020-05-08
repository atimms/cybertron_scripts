#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import pandas as pd
import scipy.stats as stats


##note
'''

load modules:
module load homer/4.9.1
module load local_python/3.6.4

'''


##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/misc/cherry_atac_correlation/cherry_atac_correlation_1218'
os.chdir(working_dir)
threads = '15'

def make_tag_dirs(name, bam):
	outdir = name + '.tag_dir'
	mk_tag_dir = subprocess.Popen(['makeTagDirectory', outdir, bam])
	mk_tag_dir.wait()


def find_peaks(name, tag_suffix, peak_suffix):
	tag_dir = name + tag_suffix
	peak_file = name + peak_suffix
	mk_tag_dir = subprocess.Popen(['findPeaks', tag_dir, '-o', peak_file])
	mk_tag_dir.wait()

def merge_peaks(file_dict, out_prefix, d_value):
	peak_files = [file_dict[i][0] for i in file_dict]
	out_file = out_prefix + d_value + 'd.txt'
	with open(out_file, 'w') as out_fh:
		mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
		mk_tag_dir.wait()

def annotate_peaks(names, peak_prefix, d_value, tag_suffix, out_prefix, size_req):
	peak_file = peak_prefix + d_value + 'd.txt'
	tag_dirs = [i + tag_suffix for i in names]
	out_file = out_prefix + d_value + 'd.' + size_req + 'size.txt'
	with open(out_file, 'w') as out_fh:
		#annotatePeaks.pl pu1peaks.txt mm8 -size 400 -d Macrophage-PU.1/ Bcell-PU.1/ > output.txt
		##this way works
		mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg19', '-size', size_req, '-d'] + tag_dirs, stdout=out_fh)
		# mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg19', '-size', size_req], stdout=out_fh)
		# mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg19', '-noann', '-size', size_req, '-norm', '10000000', '-d'] + tag_dirs, stdout=out_fh)
		mk_tag_dir.wait()

def compute_correlation(annotate_prefix, d_value, size_req, names, out_prefix):
	in_file = annotate_prefix + d_value + 'd.' + size_req + 'size.txt'
	temp_file = 'temp.' + d_value + 'd.' + size_req + 'size.txt'
	out_file = out_prefix + d_value + 'd.' + size_req + 'size.txt'
	print(in_file)
	with open(temp_file, 'w') as out_fh, open(in_file, 'r') as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.strip('\n').split(delim)
			if lc == 1:
				sample_info = line[19:]
				sample_info = [s.split('.')[0] for s in sample_info]
				out_fh.write(delim.join(['peak_id'] + sample_info + ['\n']))
				print(sample_info)
			else:
				line_out = [line[0]] + line[19:]
				out_fh.write(delim.join(line_out + ['\n']))
	##imprt to dataframe in pandas
	df = pd.read_csv(temp_file,sep='\t',header=(0))
	with open(out_file, 'w') as out_fh2:
		out_fh2.write(delim.join(['sample1', 'sample2', 'spearmanr', '\n']))
		for n1 in names:
			for n2 in names:
				print('analyzing:', n1, n2)
				# print(df[n1].corr(df[n2]))
				sp_r = df[n1].corr(df[n2], method= 'spearman')
				# sp_r = df[n1].corr(df[n2])
				print(sp_r)
				out_fh2.write(delim.join([n1, n2, str(sp_r), '\n']))
				# print(stats.pearsonr(df[n1], df[n2]))











##from summit 
d_values = ['100', '200', '500', '1000']
size_values = ['200', '400', '1000']
sample_dict = {'103dayhuret':['103dayhuret_withBackground_fixchr_summits.bed', '103dayhuret_fixchr_ext200_bedsort.bed'], 
		'11a_combined':['11α_OrgCombined_postidr_summits_withBackground.bed', '11a_combined_fixchr_ext200_bedsort.bed'],
		'125dayhuret': ['125dayhuret_withBackground_fixchr_summits.bed', '125dayhuret_fixchr_ext200_bedsort.bed'],
		'74dayhuret': ['74dayhuret_withBackground_fixchr_summits.bed', '74dayhuret_fixchr_ext200_bedsort.bed'],
		'85dayhuret': ['85dayhuret_withBackground_fixchr_summits.bed', '85dayhuret_fixchr_ext200_bedsort.bed'],
		'89dayhuret': ['89dayhuret_withBackground_fixchr_summits.bed', '89dayhuret_fixchr_ext200_bedsort.bed'],
		'bcells_combined': ['bcells_combined_postidr_summits.bed', 'bcells_combined_fixchr_ext200_bedsort.bed'],
		'adult_huret_combined': ['adult_huret_combined_postidr_summits.bed', 'huret_combined_fixchr_ext200_bedsort.bed'],
		'adult_RPE_combined': ['adult_RPE_combined_postidr_summits.bed', 'RPE_combined_fixchr_ext200_bedsort.bed'],
		'KE_org_4_5_combined': ['KE_org_4_5_combined_postidr_summits_withBackground.bed', 'KE45_combined_fixchr_ext200_bedsort.bed'],
		'KE_org_6_7_combined': ['KE_org_6_7_combined_postidr_summits_withBackground.bed', 'KE67_combined_fixchr_ext200_bedsort.bed'],
		'muscle_cells_combined': ['muscle_cells_combined_postidr_summits.bed', 'muscle_cells_combined_fixchr_ext200_bedsort.bed']}
tag_dir_suffix = '.tag_dir'
merge_prefix = 'merged.'
annotate_prefix = 'annotated.'
correlation_prefix = 'correlation.'

##make tag dirs
# for sample in sample_dict:
# 	make_tag_dirs(sample, sample_dict[sample][1])


##merge files
# '''
for d in d_values:
	# merge_peaks(sample_dict, merge_prefix, d)
	for size in size_values:
		annotate_peaks(sample_dict, merge_prefix, d, tag_dir_suffix,annotate_prefix, size)
		compute_correlation(annotate_prefix, d, size, sample_dict, correlation_prefix)
# '''
