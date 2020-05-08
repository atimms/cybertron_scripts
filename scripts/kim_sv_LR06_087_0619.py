#!/usr/bin/python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/exomes/working/kim_sv_LR06-087_0619'
os.chdir(working_dir)

##file names etc
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
exome_capture_bed = ref_dir + 'dobyns_exome.in_any_padded_target.1015.bed'
svaba = '/home/atimms/programs/svaba/bin/svaba'

def run_svaba(case_bams, control_bams, out_prefix):
	## eg targets.bed is a set of exome capture regions
	#svaba run -t $BAM -k targets.bed -a exome_cap -G $REF
	svaba_run = subprocess.Popen([svaba, 'run', '-t'] + case_bams + ['-n'] + control_bams + ['-k', exome_capture_bed,'-a', out_prefix, '-G', fasta, '-p', '10'])
	svaba_run.wait()

##run methods
case_files = ['LR06-087.bwa_mkdup.bam']
control_files = ['LR06-087f.bwa_mkdup.bam', 'LR06-087m.bwa_mkdup.bam']
# bam_files = ['LR06-087.bwa_mkdup.bam']
out_file_prefix = 'LR06-087'
run_svaba(case_files, control_files, out_file_prefix)