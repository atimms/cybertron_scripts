#!/usr/bin/env python
import sys
import subprocess
import os
import glob




##file names etc
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'


##annovar parameters
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
# av_genome = 'hg19'
# av_buildver = ['-buildver', av_genome]
# av_ref_dir = ['/Users/atimms/Desktop/ngs/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,knownGene,ensGene,rmsk,genomicSuperDups,ljb26_all,esp6500si_all,esp6500si_aa,esp6500si_ea,exac03,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_amr,1000g2014oct_eas,1000g2014oct_eur,1000g2014oct_sas,avsnp147']
# av_operation = ['-operation', 'g,g,g,r,f,f,f,f,f,f,f,f,f,f,f,f,f']
# av_options = ['-remove', '-otherinfo']

av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147,gnomad211_genome,gnomad211_exome']
# av_protocol_pisces = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147,vcf', '-vcfdbfile']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f']
# av_operation_pisces = ['-operation', 'g,r,r,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.','-arg', '-splicing 10 ,,,,,,,']
# av_options_vcf = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput','-arg', '-splicing 10 ,,,,,']



##methods
def annotate_variants(avinput, working_d):
	os.chdir(working_d)
	out_prefix = avinput.rsplit('.',1)[0]
	command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()



##run methods


##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/misc/add_annotation/ghayda_znf292_0519'
##manually make avinput file from ghayda's file
avinput_file = 'ZNF292.avinput'
annotate_variants(avinput_file, working_dir)