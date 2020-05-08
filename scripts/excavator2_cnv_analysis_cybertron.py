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
TargetPerla = '/home/atimms/programs/EXCAVATOR2_Package_v1.1.2/TargetPerla.pl'
SourceTarget = '/home/atimms/programs/EXCAVATOR2_Package_v1.1.2/SourceTarget.txt'
refgene_x_bed = '/home/atimms/programs/EXCAVATOR2_Package_v1.1.2/data/targets/hg19/hg19_refgene_chrX_nc.bed'

def run_excavator(bed, outdir):
	# python ADTEx.py --normal normal_sample.BAM --tumor tumor_sample.BAM --bed target.bed --out output_folder
	run_adtex = subprocess.Popen(['perl', TargetPerla, SourceTarget, bed, outdir, '50000', 'hg19'])
	run_adtex.wait()



##run methods


run_excavator(refgene_x_bed, 'test')




