#!/usr/bin/env python
import subprocess
import os
import glob

'''
info....

##homer
load modules:
qsub -Iq cdbrmq -l mem=100gb,ncpus=5,walltime=10:00:00 -P 19833a08-f6fb-4bea-8526-8a79069da87
module load biobuilds
module load homer/4.9.1

'''

##parameters
delim = '\t'
thread_number = '20'

##files etc
tag_suffix = '.tag_dir'

##methods

def annotate_peaks(tag_directories, peak_file, sizes_req, out_suffix):
	out_prefix = peak_file.rsplit('.',1)[0]
	for size_req in sizes_req:
		##not needed as none == size=given
		# if size_req == 'none':
		# 	out_file_no_size = out_prefix + out_suffix
		# 	with open(out_file_no_size, 'w') as out_ns_fh:
		# 		ann_peaks = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-d'] + tag_directories, stdout=out_ns_fh)
		# 		ann_peaks.wait()
		out_file_using_size = out_prefix + '.' + size_req + 'size' + out_suffix
		with open(out_file_using_size, 'w') as out_fh:
			#annotatePeaks.pl pu1peaks.txt mm8 -size 400 -d Macrophage-PU.1/ Bcell-PU.1/ > output.txt
			##alternatives
			# mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg19', '-size', size_req], stdout=out_fh)
			# mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg19', '-noann', '-size', size_req, '-norm', '10000000', '-d'] + tag_dirs, stdout=out_fh)
			ann_peaks = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-size', size_req, '-d'] + tag_directories, stdout=out_fh)
			ann_peaks.wait()

def find_genes_for_set_of_peaks(samples, cell_types, sizes, peak_files):
	##define tag dirs to use
	tag_dirs = []
	for cell_type in cell_types:
		for sample in samples:
			tag_dir = sample + '_' + cell_type + tag_suffix
			tag_dirs.append(tag_dir)
	##annoatate peaks
	ann_peak_suffix = '.ann_peak.xls'
	for peak_file in peak_files:
		annotate_peaks(tag_dirs, peak_file, sizes, ann_peak_suffix)




##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_define_enhancer_genes_0820'
os.chdir(working_dir)

###testing on a few peak file and human_adult set
##samples and tissue
sample_names = ['hu5', 'hu7', 'hu8']
cell_type_names = ['Amacrines', 'Bipolars', 'Cones', 'Ganglions', 'Horizontals', 'Mullers', 'Rods']
##so for summit.bed files best to use size = 500 or 2000, for IDR so actual peaks use size = given
# size_values = ['500', '2000', 'given']
size_values = ['given', '2000']
##try just using real peaks make bed from narrowpeak i.e.
##cut -f1-4 Rods_hu5.macs2.bampe_q0.01_keepdups_peaks.narrowPeak  > Rods_hu5.macs2.bampe_q0.01_keepdups_peaks.bed
peak_files = ['Rods_hu5.macs2.bampe_q0.01_keepdups_peaks.bed', 'Rods.macs2.bampe_p0.01_keepdups.idr_combined.bed']

##master method
find_genes_for_set_of_peaks(sample_names, cell_type_names, size_values, peak_files)



