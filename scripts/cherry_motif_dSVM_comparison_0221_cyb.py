#!/usr/bin/env python
import sys
import subprocess
import os
# import glob
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns


'''
set up:
interactive session i.e.
qsub -Iq cdbrmq -l mem=40gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load biobuilds
conda activate pd_np_plt_etc

get motif data from homer and decompress
get http://homer.ucsd.edu/homer/data/motifs/homer.KnownMotifs.hg38.191020.bed.gz
gunzip -d homer.KnownMotifs.hg38.191020.bed.gz

'''

##parameters
delim = '\t'
motif_data = 'homer.KnownMotifs.hg38.191020.bed'

def convert_dvsm_to_bed(infile, all_bed, neg_bed):
	score_dict = {}
	with open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			chrom = line[0]
			start = line[1]
			end = str(int(line[1]) + 1)
			cse = '_'.join([chrom,start,end])
			dvsm = float(line[7])
			ref = line[3]
			mut = line[4]
			# print(cse, ref, mut, dvsm, line)
			##
			if ref != mut:
				if cse in score_dict:
					score_dict[cse][0].append(dvsm)
					if dvsm < 0:
						score_dict[cse][1].append(dvsm)
				else:
					if dvsm < 0:
						score_dict[cse] = [[dvsm], [dvsm], ref]
					else:
						score_dict[cse] = [[dvsm], [], ref]
	# print(len(score_dict))
	with open(all_bed, "w") as all_fh, open(neg_bed, "w") as neg_fh:
		for s in score_dict:
			# print(s, score_dict[s])
			all_ave = sum(score_dict[s][0]) / len(score_dict[s][0])
			if len(score_dict[s][1]) == 0:
				neg_ave = 0
			else:
				neg_ave = sum(score_dict[s][1]) / len(score_dict[s][1])
			all_fh.write(delim.join(s.split('_') + [str(all_ave), ref]) + '\n')
			neg_fh.write(delim.join(s.split('_') + [str(neg_ave), ref]) + '\n')















##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_motif_dSVM_comparison_0221'
os.chdir(working_dir)
## motif bed from scanMotifGenomeWide.pl or from homer.KnownMotifs.hg38.191020.bed 
dsvm_leah = 'dsvm_simplified.txt'
# dsvm_leah = 'temp.txt'
dsvm_random_all_bed = 'dsvm.random.all_scores.bed'
dsvm_random_neg_bed = 'dsvm.random.neg_scores.bed'
dsvm_scrambled = 'ALL_diseaseLoci_deltaSVM_predict_ATAC_ChIP.formatted.txt'
dsvm_scrambled_all_bed = 'dsvm.scrambled.all_scores.bed'
dsvm_scrambled_neg_bed = 'dsvm.scrambled.neg_scores.bed'
homer_bed = 'homer.KnownMotifs.hg38.191020.bed'

##make dvsm bed files
convert_dvsm_to_bed(dsvm_leah, dsvm_random_all_bed, dsvm_random_neg_bed)
convert_dvsm_to_bed(dsvm_scrambled, dsvm_scrambled_all_bed, dsvm_scrambled_neg_bed)

##get motif specific bed
##CRX, OTX2, NRL, MEF2D, RORB, CREB, CTCF

##bedtools intersect with homer beds - complete or split?


##split by motif and graph in different ways
##all motif, pos in motif
##negative only, pos only and all




