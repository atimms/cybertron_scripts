#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


'''
set up:
interactive session i.e.
qsub -Iq cdbrmq -l mem=20gb,ncpus=1 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load biobuilds
conda activate pd_np_plt_etc

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

def combine_dsvm_cadd(cadd_file, dvsm_file, outfile):
	##make dict from cadd file
	cadd_dict = {}
	with open(cadd_file, "r") as cadd_fh:
		for line in cadd_fh:
			line = line.rstrip().split(delim)
			mut = '_'.join(line[:4])
			cadd = line[4:6]
			cadd_dict[mut] = cadd
	with open(dvsm_file, "r") as dsvm_fh, open(outfile, "w") as out_fh:
		out_fh.write(delim.join(['chr', 'pos', 'ref', 'alt', 'dvsm', 'cadd_raw', 'cadd_phred']) + '\n')
		for line in dsvm_fh:
			line = line.rstrip().split(delim)
			##need to add 1 to dsvm position values
			mut_list = [line[0], str(int(line[1]) + 1)] + line[3:5]
			mut2 = '_'.join(mut_list)
			dvsm = line[7]
			if mut2 in cadd_dict:
				line_out = mut2.split('_') + [dvsm] + cadd_dict[mut2]
				out_fh.write(delim.join(line_out) + '\n')

def make_scatter_plot(infile, pdf):
	int_data = pd.read_table(infile)
	randsample = int_data.sample(n=10000)
	plt.figure(figsize=(16,9))
	sns.relplot(x="dvsm", y="cadd_raw", data=randsample);
	plt.savefig(pdf)

##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_cadd_dSVM_comparison_1120'
os.chdir(working_dir)
bed_file = 'ALL_disease_summits_ext300_sorted_chr_removed.bed' #sorted with bedtools sort and chr removed from start of line
cadd_scores_file = 'ALL_disease_summits_ext300.CADD_score.tsv'

##get cadd score for bed
# get_cadd_score_for_bed(bed_file, cadd_scores_file)

##combine cadd and dsvm scores and grapg
dsvm_leah = 'dsvm_simplified.txt'
cadd_dsvm_leah = 'leah_dsvm_simplified.cadd.txt'
cadd_dsvm_leah_scatter = 'leah_dsvm_simplified.cadd.pdf'
##had to format with... 
##paste dsvm_simplified.txt ALL_diseaseLoci_deltaSVM_predict_ATAC_ChIP.txt | cut -f1-7,10 > ALL_diseaseLoci_deltaSVM_predict_ATAC_ChIP.formatted.txt
dsvm_scrambled = 'ALL_diseaseLoci_deltaSVM_predict_ATAC_ChIP.formatted.txt'
cadd_dsvm_scrambled = 'scrambled_dsvm.cadd.txt'
cadd_dsvm_scrambled_scatter = 'scrambled_dsvm.cadd.pdf'

##combine
# combine_dsvm_cadd(cadd_scores_file, dsvm_leah, cadd_dsvm_leah)
# combine_dsvm_cadd(cadd_scores_file, dsvm_scrambled, cadd_dsvm_scrambled)


##scatter plot
make_scatter_plot(cadd_dsvm_leah, cadd_dsvm_leah_scatter)
make_scatter_plot(cadd_dsvm_scrambled, cadd_dsvm_scrambled_scatter)


