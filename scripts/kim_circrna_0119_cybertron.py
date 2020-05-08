#!/usr/bin/env python
import sys
import subprocess
import os
import glob


##note
'''
analysis of circrna in rnaseq data
load modules:
module load biobuilds/2017.11

'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/rnaseq/kim_circrna_1218'
os.chdir(working_dir)

##programs and files
ciri2 = '/home/atimms/programs/CIRI_v2.0.6/CIRI2.pl'
fasta = '/home/atimms/ngs_data/references/igenomes/hg19/genome.fa'
# gtf = '/home/atimms/ngs_data/references/igenomes/hg19/genes.gtf'
##igenome gtf gave an error message, so used the UCSC table browser to get refgene genes
gtf = 'hg19_refgene.gtf'




def map_using_bwa(samples, r1_suffix, r2_suffix, sam_suffix):
	for sample in samples:
		r1_fq = sample + r1_suffix
		r2_fq = sample + r2_suffix
		sam_file = sample + sam_suffix
		##bwa and convert to bam
		bwa_pe = subprocess.Popen(['bwa', 'mem', '-T', '19', '-t', '18', '-o', sam_file, fasta, r1_fq, r2_fq])
		bwa_pe.wait()

def run_ciri2(samples, sam_suffix, out_suffix):
	for sample in samples:
		sam_file = sample + sam_suffix
		out_file = sample + out_suffix
		##bwa and convert to bam
		bwa_pe = subprocess.Popen(['perl', ciri2, '-I', sam_file, '-O', out_file, '-F', fasta, '-A', gtf])
		bwa_pe.wait()






##run methods

# sample_names = ['H26857-EGL1-tube13']
# sample_names = ['H26857-PK1-1-tube11', 'H26857_RL_combined', 'H26857-RL-tube14', 'H26857-VZ-tube16', 'H26857-Whm-tube17', 'H26857-Wver-tube15']

sample_names = ['H26857-EGL1-tube13', 'H26857-PK1-1-tube11', 'H26857_RL_combined', 'H26857-RL-tube14', 'H26857-VZ-tube16', 'H26857-Whm-tube17', 'H26857-Wver-tube15']
fq_r1_suffix = '.r1.fq.gz'
fq_r2_suffix = '.r2.fq.gz'
sam_file_suffix = '.sam'
out_file_suffix = '.out'
##map files
# map_using_bwa(sample_names, fq_r1_suffix, fq_r2_suffix, sam_file_suffix)
run_ciri2(sample_names, sam_file_suffix, out_file_suffix)