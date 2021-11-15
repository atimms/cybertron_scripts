#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
requires:

'''
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'


##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,gnomad211_exome,gnomad211_genome,dbnsfp41a,clinvar_20210123,generic', '-genericdbfile', 'Geno2MP.avinput']
av_operation = ['-operation', 'g,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.']


def run_annovar(avinput_file, av_out):
	##annotate vcfs with annovar i.e. run_table_annovar
	command = [table_annovar] + av_buildver + [avinput_file] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_out]
	annovar = subprocess.Popen(command)
	annovar.wait()


##for SMARCA1 1121
# working_dir = '/archive/mirzaa_g/misc/SMARCA1_var_annotation_1121'
# os.chdir(working_dir)
# run_annovar('smarca1_vars_1121.avinput', 'smarca1_vars_1121')



##for ANKLE2 1121
working_dir = '/archive/mirzaa_g/misc/ankle2_var_annotation_1121'
os.chdir(working_dir)
run_annovar('ankle2_vars_1121.avinput', 'ankle2_vars_1121')


