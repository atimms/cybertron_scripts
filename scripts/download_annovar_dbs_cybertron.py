#!/usr/bin/env python
import sys
import subprocess
import os


##parameters
delim = '\t'
ann_var = '/home/atimms/programs/annovar/annotate_variation.pl'
av_genome = 'hg19'
av_ref_dir = '/home/atimms/ngs_data/references/annovar/' + av_genome

def download_annovar_dbs(db_name):
	#$ann_var -buildver hg19 -downdb clinvar_20170130 $av_ref_hg19 -webfrom annovar
	dl_ann = subprocess.Popen([ann_var, '-buildver', av_genome,  '-downdb', db_name, av_ref_dir, '-webfrom', 'annovar'])
	dl_ann.wait()

##download
db = 'clinvar_20170905'
download_annovar_dbs(db)