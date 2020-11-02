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
ann_peaks = '/home/atimms/programs/homer/bin/annotatePeaks.pl'

def annotate_peaks(peak_file):
	outfile = peak_file.split('.')[0] + '.homer.txt'
	with open(outfile, 'w') as out_fh:
		mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'mm9'], stdout=out_fh)
		# mk_tag_dir = subprocess.Popen([ann_peaks, peak_file, 'mm9'], stdout=out_fh)
		mk_tag_dir.wait()


##run methods
working_dir = '/home/atimms/ngs_data/misc/lisa_pbx_chip_1020'
os.chdir(working_dir)
pbx_peak_file = 'penkov_pbx_peaks.bed'
annotate_peaks(pbx_peak_file)