#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
requires:
qsub -Iq cdbrmq -l mem=30gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
#java 1.8 for snpeff
module load java/1.8.0_202 
#biobuilds for tabix etc
module load biobuilds/2017.11
'''
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'


av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,knownGene,ensGene,gnomad211_exome']
av_operation = ['-operation', 'g,g,g,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-arg', '-splicing 10 ,,,']


def run_annovar(avinput_file, av_out):
	##annotate vcfs with annovar i.e. run_table_annovar
	command = [table_annovar] + av_buildver + [avinput_file] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_out]
	annovar = subprocess.Popen(command)
	annovar.wait()



working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalyze_exomes_0121/var_check_0821'
os.chdir(working_dir)


run_annovar('var_check_0821.avinput', 'var_check_0821')