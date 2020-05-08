#!/usr/bin/env python
import sys
import subprocess
import os
import glob



'''
##load modules required for analysis
module load biobuilds/2017.11
module load R/3.5.0
'''


##parameters
delim = '\t'
adtex = '/home/atimms/programs/ADTEx.v.2.0/ADTEx.py'



def run_adtex_from_bam(sample_bam, control_bam, bed, outdir):
	# python ADTEx.py --normal normal_sample.BAM --tumor tumor_sample.BAM --bed target.bed --out output_folder
	run_adtex = subprocess.Popen(['python', adtex, '-n', control_bam, '-t', sample_bam, '-o', outdir, '-b', bed ])
	run_adtex.wait()



##run methos
ss_bed = 'SureSelectXT_hg19_V5_UTRs_padded_no_chr.bed'

proband_bam = 'LR15-377a1.bwa_gatk.bam'
parent_bam = 'LR15-377m.bwa_gatk.bam'
out_folder = proband_bam.split('.')[0] + '.adtex_results'

run_adtex_from_bam(proband_bam, parent_bam, ss_bed,out_folder)


proband_bam = 'LR15-377a2.bwa_gatk.bam'
parent_bam = 'LR15-377r.bwa_gatk.bam'
out_folder = proband_bam.split('.')[0] + '.adtex_results'

run_adtex_from_bam(proband_bam, parent_bam, ss_bed,out_folder)



