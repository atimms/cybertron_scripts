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

##methods
def make_hold_out_beds(beds, bed_names):
	bed_combs = combinations(beds, 11)
	bed_dict = dict(zip(beds, bed_names))
	for bed_comb in list(bed_combs):
		# print(list(bed_comb))
		missing_bed = list(set(beds) - set(bed_comb))
		# print(missing_bed)
		##returns finla list with all beds (weird)
		final_bed = 'Huret.minus_' + bed_dict[missing_bed[0]] + '.bed'
		t1_bed = 'Huret.minus_' + bed_dict[missing_bed[0]] + 'temp1.bed'
		# t2_bed = 'Huret.minus_' + bed_dict[missing_bed[0]] + 'temp2.bed'
		print(missing_bed, final_bed)
		##combine beds
		with open(t1_bed, 'w') as t1_fh:
			cat_beds = subprocess.Popen(['cat'] + list(bed_comb), stdout=t1_fh)
			cat_beds.wait()
		##sort and merge
		with open(final_bed, "w") as out_fh:
			sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', t1_bed], stdout=subprocess.PIPE)
			bt_merge = subprocess.Popen([bedtools, 'merge', '-i', '-'], stdin=sort_file.stdout, stdout=out_fh)
			bt_merge.wait()

def split_neg_fasta(fasta):
	bt_merge = subprocess.Popen(['pyfasta', 'split', '-n', '9', fasta])
	bt_merge.wait()

def make_fastas_from_bed(beds):
	for bed in beds:
		out_fasta = bed.rsplit('.',1)[0] + '.fasta'
		#bedtools getfasta -fi test.fa -bed test.bed -fo test.fa.out
		bt_gf = subprocess.Popen([bedtools, 'getfasta', '-fi', hg38_fasta, '-bed', bed, '-fo', out_fasta])
		bt_gf.wait()

def run_lsgkm_crossvalidation(pos_set, neg_set, out_prefix):
	#gkmtrain -l 11 -k 7 -d 3 -t 2 -e 0.005 -x 5 (also added 20gb cache and 4 threads)
	lsgkm_cv = subprocess.Popen([gkmtrain, '-l', '11', '-k', '7', '-d', '3', '-t', '2', '-e', '0.005', '-x', '5', '-m', '20000', '-T', '4', pos_set, neg_set, out_prefix])
	lsgkm_cv.wait()

def run_lsgkm_without_crossvalidation(pos_set, neg_set, out_prefix):
	##gkmtrain
	#gkmtrain -l 11 -k 7 -d 3 -t 2 -e 0.005 -x 5 (also added 20gb cache and 4 threads)
	lsgkm_cv = subprocess.Popen([gkmtrain, '-l', '11', '-k', '7', '-d', '3', '-t', '2', '-e', '0.005', '-m', '20000', '-T', '4', pos_set, neg_set, out_prefix])
	lsgkm_cv.wait()
	##gkmpredict on pos and neg sets
	#gkmpredict [options] <test_seqfile> <model_file> <output_file>
	model_file = out_prefix + '.model.txt'
	pos_predictions = out_prefix + '.positive_prediction.txt'
	neg_predictions = out_prefix + '.negative_prediction.txt'
	lsgkm_predict_a = subprocess.Popen([gkmpredict, '-T', '4', pos_set, model_file, pos_predictions])
	lsgkm_predict_a.wait()
	lsgkm_predict_b = subprocess.Popen([gkmpredict, '-T', '4', neg_set, model_file, neg_predictions])
	lsgkm_predict_b.wait()
	##gkmmatrix
	out_matrix = out_prefix + '.matrix.txt'
	#gkmmatrix -l 11 -k 7 -d 3 -t 2 -T 4 (also 4 threads)
	# lsgkm_matrix = subprocess.Popen([gkmmatrix, '-l', '11', '-k', '7', '-d', '3', '-t', '2', pos_set, neg_set, out_matrix])
	# lsgkm_matrix.wait()

def run_lsgkm_without_crossvalidation_random_10k(pos_set, neg_fasta, out_prefix):
	random_10k_fasta = out_prefix + '.fasta'
	'''
	##make dict from 86k fasta
	neg_dict = {rec.id : rec.seq for rec in SeqIO.parse(neg_fasta, "fasta")}
	##get random 10k sequence
	tenk_rand_dict = dict(random.sample(neg_dict.items(), 10000))
	# for i in tenk_rand_dict:
	# 	print(i, tenk_rand_dict[i])
	##write 10k to fasta
	with open(random_10k_fasta, "w") as out_fh:
		for name in tenk_rand_dict:
			out_fh.write('>' + name + '\n')
			out_fh.write(str(tenk_rand_dict[name]) + '\n')
	#gkmtrain -l 11 -k 7 -d 3 -t 2 -e 0.005 -x 5 (also added 20gb cache and 4 threads)
	lsgkm_cv = subprocess.Popen([gkmtrain, '-l', '11', '-k', '7', '-d', '3', '-t', '2', '-e', '0.005', '-m', '20000', '-T', '4', pos_set, random_10k_fasta, out_prefix])
	lsgkm_cv.wait()
	'''
	##gkmpredict on pos and neg sets
	#gkmpredict [options] <test_seqfile> <model_file> <output_file>
	model_file = out_prefix + '.model.txt'
	pos_predictions = out_prefix + '.positive_prediction.txt'
	neg_predictions = out_prefix + '.negative_prediction.txt'
	# lsgkm_predict_a = subprocess.Popen([gkmpredict, '-T', '4', pos_set, model_file, pos_predictions])
	# lsgkm_predict_a.wait()
	# lsgkm_predict_b = subprocess.Popen([gkmpredict, '-T', '4', random_10k_fasta, model_file, neg_predictions])
	# lsgkm_predict_b.wait()
	##gkmmatrix
	out_matrix = out_prefix + '.matrix.txt'
	#gkmmatrix -l 11 -k 7 -d 3 -t 2 -T 4 (also 4 threads)
	lsgkm_matrix = subprocess.Popen([gkmmatrix, '-l', '11', '-k', '7', '-d', '3', '-t', '2', '-T', '4', pos_set, random_10k_fasta, out_matrix])
	lsgkm_matrix.wait()


##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_lsgkm_subsets_0321'
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
random_neg_86k = 'filterN_rand_100k_filt_cleaned.fa'
random_neg_9k = 'filterN_rand_100k_filt_cleaned.0.fa'
tenk_fastas = [i.rsplit('.',1)[0] + '.fasta' for i in tenk_beds]
hold_out_beds = ['Huret.minus_ATAC.bed', 'Huret.minus_ChIP_CREB.bed', 'Huret.minus_ChIP_CRX_B11X.bed', 
		'Huret.minus_ChIP_CRX_H120X.bed', 'Huret.minus_ChIP_CTCF.bed', 'Huret.minus_ChIP_H3K27ac.bed', 
		'Huret.minus_ChIP_H3K4me2.bed', 'Huret.minus_ChIP_MEF2D.bed', 'Huret.minus_ChIP_NRL.bed', 
		'Huret.minus_ChIP_OTX2_abcam.bed', 'Huret.minus_ChIP_OTX2_thermo.bed', 'Huret.minus_ChIP_RORB.bed']
hold_out_fastas = [i.rsplit('.',1)[0] + '.fasta' for i in hold_out_beds]
##combine beds, leaving out one set at a time
# make_hold_out_beds(tenk_beds, tenk_bed_names)

##makes fastas from all bed files
# make_fastas_from_bed(tenk_beds)
# make_fastas_from_bed(hold_out_beds)

##make 10k negative training set for individual positive sets
# split_neg_fasta(random_neg_86k)

##run lsgkm cross validation on single and leave out data
##for 10k set
# for tenk_fasta in tenk_fastas:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0]
# 	run_lsgkm_crossvalidation(tenk_fasta, random_neg_9k, outfile_prefix)

##for 86k set etc
# ho_fastas = hold_out_fastas[:3]
# ho_fastas = hold_out_fastas[3:6]
# ho_fastas = hold_out_fastas[6:9]
# ho_fastas = hold_out_fastas[9:]
# for ho_fasta in ho_fastas:
# 	outfile_prefix = ho_fasta.rsplit('.',1)[0]
# 	run_lsgkm_crossvalidation(ho_fasta, random_neg_86k, outfile_prefix)


##for 10k set against 86k set without cross validation
# for tenk_fasta in tenk_fastas:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.86k_neg'
# 	run_lsgkm_without_crossvalidation(tenk_fasta, random_neg_86k, outfile_prefix)


##for 10k set against random 10 from 86k set 
# for tenk_fasta in tenk_fastas:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.10k_rand_neg'
# 	run_lsgkm_without_crossvalidation_random_10k(tenk_fasta, random_neg_86k, outfile_prefix)


##run lsgkm on complete
merged_pos_peaks = 'merged_huret_ATAC_ChIP_postidr_summits_ext300_noNs.fa'
outfile_prefix = merged_pos_peaks.rsplit('.',1)[0] + '.86k_neg'
run_lsgkm_without_crossvalidation(merged_pos_peaks, random_neg_86k, outfile_prefix)



##temp stuffff....

tenk_fastas1 = tenk_fastas[:4]
tenk_fastas2 = tenk_fastas[4:8]
tenk_fastas3 = tenk_fastas[8:]

# for tenk_fasta in tenk_fastas1:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.10k_rand_neg'
# 	run_lsgkm_without_crossvalidation_random_10k(tenk_fasta, random_neg_86k, outfile_prefix)
# for tenk_fasta in tenk_fastas2:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.10k_rand_neg'
# 	run_lsgkm_without_crossvalidation_random_10k(tenk_fasta, random_neg_86k, outfile_prefix)
# for tenk_fasta in tenk_fastas3:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.10k_rand_neg'
# 	run_lsgkm_without_crossvalidation_random_10k(tenk_fasta, random_neg_86k, outfile_prefix)

tenk_fastas1 = tenk_fastas[:2]
tenk_fastas2 = tenk_fastas[2:4]
tenk_fastas3 = tenk_fastas[4:6]
tenk_fastas4 = tenk_fastas[6:8]
tenk_fastas5 = tenk_fastas[8:10]
tenk_fastas6 = tenk_fastas[10:]

# for tenk_fasta in tenk_fastas1:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.86k_neg'
# 	run_lsgkm_without_crossvalidation(tenk_fasta, random_neg_86k, outfile_prefix)
# for tenk_fasta in tenk_fastas2:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.86k_neg'
# 	run_lsgkm_without_crossvalidation(tenk_fasta, random_neg_86k, outfile_prefix)
# for tenk_fasta in tenk_fastas3:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.86k_neg'
# 	run_lsgkm_without_crossvalidation(tenk_fasta, random_neg_86k, outfile_prefix)
# for tenk_fasta in tenk_fastas4:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.86k_neg'
# 	run_lsgkm_without_crossvalidation(tenk_fasta, random_neg_86k, outfile_prefix)
# for tenk_fasta in tenk_fastas5:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.86k_neg'
# 	run_lsgkm_without_crossvalidation(tenk_fasta, random_neg_86k, outfile_prefix)
# for tenk_fasta in tenk_fastas6:
# 	outfile_prefix = tenk_fasta.rsplit('.',1)[0] + '.86k_neg'
# 	run_lsgkm_without_crossvalidation(tenk_fasta, random_neg_86k, outfile_prefix)





