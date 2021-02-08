#!/usr/bin/env python
import os
import subprocess
import glob
import shutil


'''
load these:

'''

##set input variables and parameters
delim = '\t'


##program
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
ann_var = '/home/atimms/programs/annovar_1019/annotate_variation.pl'
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'

##annovar parameters --- add cosmic
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene']
av_operation = ['-operation', 'g']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']

##methods
def annotate_filtered_pisces_vcf(samples):
	for sample in samples:
		in_vcf = sample + '_mutect2.filtered.vcf.gz'
		##run_table_annovar
		command = [table_annovar] + av_buildver + [in_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', sample]
		annovar = subprocess.Popen(command)
		annovar.wait()

##run methods
work_dir = '/home/atimms/ngs_data/genomes/ghayda_ucb_somatic_comparison_0221'
os.chdir(work_dir)
all_samples = ['UC1401004', 'UC1401005', 'UC1401007', 'UC1401013_2', 'UC1401017', 'UC1401039', 'UC1402089', 
	'UC1405044', 'UC1701038-2', 'UC1708005', 'UC1711041', 'UC1812014']

annotate_filtered_pisces_vcf(all_samples)