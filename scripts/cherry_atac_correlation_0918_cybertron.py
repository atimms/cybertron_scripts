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
working_dir = '/home/atimms/ngs_data/misc/cherry_atac_correlation_0918'
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

def merge_peaks(names, peak_file_suffix, out_prefix, d_value):
	peak_files = [i + peak_file_suffix for i in names]
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
		mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg19', '-size', size_req, '-d'] + tag_dirs, stdout=out_fh)
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










##original analysis
samples = ['125dayhuret', 'bcells_combined', 'KE_org_6_7_combined']
bed_suffix = '_fixchr_bedsort.bed'
tag_dir_suffix = '.tag_dir'
peak_file_suffix = '.peaks'
d_values = ['100', '200', '500', '1000']
size_values = ['200', '400', '1000']
merge_prefix = 'merged.'
annotate_prefix = 'annotated.'
correlation_prefix = 'correlation.'
##make tag dir and find peaks for all samples
'''
for sample in samples:
	# make_tag_dirs(sample, sample + bed_suffix)
	find_peaks(sample, tag_dir_suffix, peak_file_suffix)
'''
##merge files
'''
for d in d_values:
	# merge_peaks(samples, peak_file_suffix, merge_prefix, d)
	for size in size_values:
		# annotate_peaks(samples, merge_prefix, d, tag_dir_suffix,annotate_prefix, size)
		compute_correlation(annotate_prefix, d, size, samples, correlation_prefix)
'''

##testing
# compute_correlation(annotate_prefix, '1000', '1000', samples, 'testing.')

# annotate_peaks(samples, merge_prefix, '1000', tag_dir_suffix,'testing.', '1000')



##from summit 
d_values = ['100', '200', '500', '1000']
size_values = ['200', '400', '1000']
samples = ['bcells_combined', 'KE_org_6_7_combined', '125dayhuret']
tag_dir_suffix = '.tag_dir'
peak_file_suffix = '_summits.bed'
merge_prefix = 'merged_sum.'
annotate_prefix = 'annotated_sum.'
correlation_prefix = 'correlation_sum.'
##merge files
# '''
for d in d_values:
	merge_peaks(samples, peak_file_suffix, merge_prefix, d)
	for size in size_values:
		annotate_peaks(samples, merge_prefix, d, tag_dir_suffix,annotate_prefix, size)
		compute_correlation(annotate_prefix, d, size, samples, correlation_prefix)
# '''
