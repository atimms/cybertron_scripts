#!/usr/bin/env python
import sys
import subprocess
import os
import glob
from itertools import combinations
# from Bio import SeqIO
import random

'''
set up:
interactive session i.e.
qsub -Iq cdbrmq -l mem=40gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
for lsgkm:
qsub -Iq cdbrmq -l mem=80gb,ncpus=4 -P 19833a08-f6fb-4bea-8526-8a79069da878 
for splitting the fasta file:
conda e
##add biopython to pyfasta env
conda install -c anaconda biopython
'''

##parameters
delim = '\t'
##programs
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
gkmtrain = '/home/atimms/programs/lsgkm/src/gkmtrain'
gkmpredict = '/home/atimms/programs/lsgkm/src/gkmpredict'
gkmmatrix = '/home/atimms/programs/lsgkm/src/gkmmatrix'

##files
hg38_fasta = '/home/atimms/ngs_data/references/hg38/hg38.fa'


###methods
def get_unique_and_combined_peaks(beds, bed_names):
	##get all combinations of 11 beds out of the 12
	bed_combs = combinations(beds, 11)
	##make dict for names
	bed_dict = dict(zip(beds, bed_names))
	for bed_comb in list(bed_combs):
		print(list(bed_comb))
		missing_bed = list(set(beds) - set(bed_comb))
		print(missing_bed)
		unique_bed = 'Huret_' + bed_dict[missing_bed[0]] + '.unique.bed'
		nonunique_bed = 'Huret_' + bed_dict[missing_bed[0]] + '.non_unique.bed'
		with open(nonunique_bed, "w") as nub_fh:
			bt_int = subprocess.Popen([bedtools, 'intersect', '-wa', '-a', missing_bed[0], '-b'] + list(bed_comb), stdout=nub_fh)
			bt_int.wait()
		with open(unique_bed, "w") as ub_fh:
			bt_int = subprocess.Popen([bedtools, 'intersect', '-wa', '-v', '-a', missing_bed[0], '-b'] + list(bed_comb), stdout=ub_fh)
			bt_int.wait()
		# ##returns final list with all beds (weird)
		# final_bed = 'Huret.minus_' + bed_dict[missing_bed[0]] + '.bed'
		# t1_bed = 'Huret.minus_' + bed_dict[missing_bed[0]] + 'temp1.bed'
		# # t2_bed = 'Huret.minus_' + bed_dict[missing_bed[0]] + 'temp2.bed'
		# print(missing_bed, final_bed)
		# ##combine beds
		# with open(t1_bed, 'w') as t1_fh:
		# 	cat_beds = subprocess.Popen(['cat'] + list(bed_comb), stdout=t1_fh)
		# 	cat_beds.wait()
		# ##sort and merge
		# with open(final_bed, "w") as out_fh:
		# 	sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', t1_bed], stdout=subprocess.PIPE)
		# 	bt_merge = subprocess.Popen([bedtools, 'merge', '-i', '-'], stdin=sort_file.stdout, stdout=out_fh)
		# 	bt_merge.wait()

##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_lsgkm_subsets_0421'
os.chdir(working_dir)
##filenames
tenk_beds = ['adult_huret_ATAC_combined_postidr_summits_t10k_300bp.bed', 'Huret_ChIP_CREB_combined_postidr_summits_t10k_300bp.bed', 
	'Huret_ChIP_CRX_B11X_combined_postidr_summits_t10k_300bp.bed', 'Huret_ChIP_CRX_H120X_combined_postidr_summits_t10k_300bp.bed', 
	'Huret_ChIP_CTCF_combined_postidr_summits_t10k_300bp.bed', 'Huret_ChIP_H3K27ac_combined_postidr_summits_t10k_300bp.bed', 
	'Huret_ChIP_H3K4me2_combined_postidr_summits_t10k_300bp.bed', 'Huret_ChIP_MEF2D_combined_postidr_summits_t10k_300bp.bed', 
	'Huret_ChIP_NRL_combined_postidr_summits_t10k_300bp.bed', 'Huret_ChIP_OTX2_abcam_combined_postidr_summits_t10k_300bp.bed', 
	'Huret_ChIP_OTX2_thermo_combined_postidr_summits_t10k_300bp.bed', 'Huret_ChIP_RORB_combined_postidr_summits_t10k_300bp.bed']
tenk_bed_names = ['ATAC', 'ChIP_CREB', 'ChIP_CRX_B11X', 'ChIP_CRX_H120X', 'ChIP_CTCF', 'ChIP_H3K27ac', 
	'ChIP_H3K4me2', 'ChIP_MEF2D', 'ChIP_NRL', 'ChIP_OTX2_abcam', 'ChIP_OTX2_thermo', 'ChIP_RORB']

###
get_unique_and_combined_peaks(tenk_beds, tenk_bed_names)




