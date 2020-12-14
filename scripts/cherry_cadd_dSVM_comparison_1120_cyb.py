#!/usr/bin/env python
import sys
import subprocess
import os
import glob

'''
set up:
interactive session i.e.
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load biobuilds

get cadd data, and place in references:
wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi


'''

##parameters
delim = '\t'
cadd_data = '/home/atimms/ngs_data/references/cadd/GRCh38_v1.6/whole_genome_SNVs.tsv.gz'

##methods
def get_cadd_score_for_bed(bed, outfile):
	with open(outfile, 'w') as out_fh:
		run_tabix = subprocess.Popen(['tabix', '-R', bed, cadd_data], stdout=out_fh)
		run_tabix.wait()	

##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_cadd_dSVM_comparison_1120'
os.chdir(working_dir)
bed_file = 'ALL_disease_summits_ext300_sorted_chr_removed.bed' #sorted with bedtools sort and chr removed from start of line
cadd_scores_file = 'ALL_disease_summits_ext300.CADD_score.tsv'

get_cadd_score_for_bed(bed_file, cadd_scores_file)

