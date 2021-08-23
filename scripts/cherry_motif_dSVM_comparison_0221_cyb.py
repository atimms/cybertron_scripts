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
qsub -Iq cdbrmq -l mem=40gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
conda activate pd_np_plt_etc

get motif data from homer and decompress
get http://homer.ucsd.edu/homer/data/motifs/homer.KnownMotifs.hg38.191020.bed.gz
gunzip -d homer.KnownMotifs.hg38.191020.bed.gz

'''

##parameters
delim = '\t'
##programs
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
##files etc
bt_genome = '/home/atimms/ngs_data/references/hg38/hg38.nochr.chrom.sizes'


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
	all_bed_temp = all_bed.rsplit('.', 1)[0] + 'temp.bed'
	neg_bed_temp = neg_bed.rsplit('.', 1)[0] + 'temp.bed'
	with open(all_bed_temp, "w") as allt_fh, open(neg_bed_temp, "w") as negt_fh:
		for s in score_dict:
			# print(s, score_dict[s])
			all_ave = sum(score_dict[s][0]) / len(score_dict[s][0])
			if len(score_dict[s][1]) == 0:
				neg_ave = 'na'
			else:
				neg_ave = sum(score_dict[s][1]) / len(score_dict[s][1])
			allt_fh.write(delim.join(s.split('_') + [str(all_ave), score_dict[s][2]]) + '\n')
			negt_fh.write(delim.join(s.split('_') + [str(neg_ave), score_dict[s][2]]) + '\n')
	##writw final files
	with open(all_bed, "w") as all_fh:
		sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', all_bed_temp], stdout=all_fh)
		sort_file.wait()
	with open(neg_bed, "w") as neg_fh:
		sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', neg_bed_temp], stdout=neg_fh)
		sort_file.wait()

def get_specific_motifs_from_homer_bed(infile, nomerge_out_prefix, merge_out_prefix, motif_dict):
	for motif in motif_dict:
		alt_motif = motif_dict[motif]
		no_merge_bed = nomerge_out_prefix + alt_motif + '.bed'
		final_bed = merge_out_prefix + alt_motif + '.bed'
		with open(infile, "r") as in_fh, open(no_merge_bed, "w") as tb_fh:
			for line in in_fh:
				line = line.rstrip().split(delim)
				chrom = line[0].replace('chr','')
				start = str(int(line[1]) - 1)
				end = line[2]
				homer_motif = line[3]
				if homer_motif == motif:
					tb_fh.write(delim.join([chrom, start, end] + [motif]) + '\n')
		with open(final_bed, "w") as out_fh:
			sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', no_merge_bed], stdout=subprocess.PIPE)
			bt_merge = subprocess.Popen([bedtools, 'merge', '-i', '-'], stdin=sort_file.stdout, stdout=out_fh)
			bt_merge.wait()


def bedtools_intersect_dsvm_motif(dsvm_beds, motifs_beds, out_suffix):
	for dsvm_bed in dsvm_beds:
		for motifs_bed in motifs_beds:
			motif = motifs_bed.split('.')[3]
			out_bed = dsvm_bed.rsplit('.',1)[0] + '.' + motif + out_suffix
			with open(out_bed, "w") as naf_fh: 
				hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', motifs_bed, '-b', dsvm_bed, '-wa', '-wb'], stdout=naf_fh)
				hom_bt_intersect.wait()




def get_counts_per_dvsm_and_motif(dsvm_beds, motifs_beds, out_file):
	dsvm_bed_count = 0
	header = ['dsvm_type']
	with open(out_file, "w") as out_fh:
		for dsvm_bed in dsvm_beds:
			dsvm_type = dsvm_bed.rsplit('.',1)[0]
			dsvm_bed_count += 1
			dsvm_scores = []
			for motifs_bed in motifs_beds:
				motif = motifs_bed.split('.')[3]
				header.append(motif)
				dvsm_motif_bed = dsvm_type + '.' + motif + '.bed'
				dvsm_list = []
				with open(dvsm_motif_bed, "r") as dmb_fh:
					for line in dmb_fh:
						line = line.rstrip().split(delim)
						score = line[6]
						if score != 'na':
							dvsm_list.append(float(score))
				average_dvsm = (sum(dvsm_list)/len(dvsm_list))
				dsvm_scores.append(str(round(average_dvsm,3)))
			all_bed = dsvm_type + '.bed'
			header.append('all')
			dvsm_list = []
			with open(all_bed, "r") as ab_fh:
				for line in ab_fh:
					line = line.rstrip().split(delim)
					score = line[3]
					if score != 'na':
						dvsm_list.append(float(score))
			average_dvsm = (sum(dvsm_list)/len(dvsm_list))
			dsvm_scores.append(str(round(average_dvsm,3)))

			if dsvm_bed_count == 1:
				out_fh.write(delim.join(header) + '\n')
			line_out = [dsvm_type] + dsvm_scores
			out_fh.write(delim.join(line_out) + '\n')


def bedtools_slop_motif_beds(in_beds):
	for in_bed in in_beds:
		##bedtools slop
		out_bed = in_bed.rsplit('.', 1)[0] + '.ext25.bed'
		with open(out_bed, "w") as out_fh: 
			hom_bt_intersect = subprocess.Popen([bedtools, 'slop', '-i', in_bed, '-g', bt_genome, '-b', '25'], stdout=out_fh)
			hom_bt_intersect.wait()

def get_ave_dvsm_per_bp_around_motif(dvsm_motif_beds):
	for dvsm_motif_bed in dvsm_motif_beds:
		##make dict with all score per motif regions
		motif_region_svsm_dict = {}
		with open(dvsm_motif_bed, "r") as dmb_fh:
			for line in dmb_fh:
				line = line.rstrip().split(delim)
				motif_region = '_'.join(line[:3])
				motif_region_len = int(line[2]) - int(line[1])
				score = line[7]
				if motif_region in motif_region_svsm_dict:
					motif_region_svsm_dict[motif_region][0].append(score)
				else:
					motif_region_svsm_dict[motif_region] = [[score], motif_region_len]
		##make new dict with just the regions that have dsvm score for all the reion
		new_motif_dvsm_dict = {}
		for mr in motif_region_svsm_dict:
			# print(dvsm_motif_bed, mr, motif_region_svsm_dict[mr][1], len(motif_region_svsm_dict[mr][0]))
			if motif_region_svsm_dict[mr][1] == len(motif_region_svsm_dict[mr][0]):
				new_motif_dvsm_dict[mr] = motif_region_svsm_dict[mr][0]
		##make dict into a panda's dataframe and get 
		motif_df = pd.DataFrame.from_dict(new_motif_dvsm_dict, orient='index')
		##write file
		dvsm_motif_csv = dvsm_motif_bed.rsplit('.', 1)[0] + '.csv'
		motif_df.to_csv(dvsm_motif_csv)
		##convert columns to floats
		cols = motif_df.columns
		motif_df[cols] = motif_df[cols].apply(pd.to_numeric, errors='coerce')
		##get averages
		motif_df_means = motif_df.mean(axis = 0, skipna = True)
		# print(motif_df_means)
		# plt.figure(figsize=(6,4))
		sns.lineplot(data=motif_df_means)
		pdf_name = dvsm_motif_bed.rsplit('.', 1)[0] + '.pdf' 
		plt.savefig(pdf_name)
		plt.close()


def get_ave_dvsm_per_bp_only_motif(dvsm_motif_beds):
	for dvsm_motif_bed in dvsm_motif_beds:
		##make dict with all score per motif regions
		motif_region_svsm_dict = {}
		with open(dvsm_motif_bed, "r") as dmb_fh:
			for line in dmb_fh:
				line = line.rstrip().split(delim)
				motif_region = '_'.join(line[:3])
				motif_region_len = int(line[2]) - int(line[1])
				score = line[7]
				if motif_region in motif_region_svsm_dict:
					motif_region_svsm_dict[motif_region][0].append(score)
				else:
					motif_region_svsm_dict[motif_region] = [[score], motif_region_len]
		##make new dict with just the regions that have dsvm score for all the reion
		new_motif_dvsm_dict = {}
		for mr in motif_region_svsm_dict:
			# print(dvsm_motif_bed, mr, motif_region_svsm_dict[mr][1], len(motif_region_svsm_dict[mr][0]))
			if motif_region_svsm_dict[mr][1] == len(motif_region_svsm_dict[mr][0]):
				new_motif_dvsm_dict[mr] = motif_region_svsm_dict[mr][0]
		##make dict into a panda's dataframe and get 
		motif_df = pd.DataFrame.from_dict(new_motif_dvsm_dict, orient='index')
		##write file
		dvsm_motif_csv = dvsm_motif_bed.rsplit('.', 1)[0] + '.csv'
		motif_df.to_csv(dvsm_motif_csv)
		##convert columns to floats
		cols = motif_df.columns
		motif_df[cols] = motif_df[cols].apply(pd.to_numeric, errors='coerce')
		##get averages
		# motif_df_means = motif_df.mean(axis = 0, skipna = True)
		# print(motif_df_means)
		# plt.figure(figsize=(6,4))
		sns.barplot(data=motif_df, color = 'blue', ci=None)
		pdf_name = dvsm_motif_bed.rsplit('.', 1)[0] + '.pdf' 
		plt.savefig(pdf_name)
		plt.close()



##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_motif_dSVM_comparison_0221'
os.chdir(working_dir)
## motif bed from scanMotifGenomeWide.pl or from homer.KnownMotifs.hg38.191020.bed 
dsvm_leah = 'dsvm_simplified.txt'
# dsvm_leah = 'temp.txt'
dsvm_random_all_bed = 'dsvm.random.all_scores.bed'
dsvm_random_neg_bed = 'dsvm.random.neg_scores.bed'
dsvm_random_motif_bed = 'dsvm.random.neg_scores.motifs.bed'
dsvm_scrambled = 'ALL_diseaseLoci_deltaSVM_predict_ATAC_ChIP.formatted.txt'
dsvm_scrambled_all_bed = 'dsvm.scrambled.all_scores.bed'
dsvm_scrambled_neg_bed = 'dsvm.scrambled.neg_scores.bed'
dsvm_bed_files = [dsvm_random_all_bed, dsvm_random_neg_bed, dsvm_scrambled_all_bed, dsvm_scrambled_neg_bed]
homer_bed = 'homer.KnownMotifs.hg38.191020.bed'
homer_motif_dict = {'CRX(Homeobox)': 'CRX', 'Otx2(Homeobox)': 'OTX2', 'MafF(bZIP)': 'NRL', 
		'Mef2d(MADS)':'MEF2D', 'RORgt(NR)': 'RORB', 'CREB5(bZIP)':'CREB', 
		'CTCF(Zf)': 'CTCF'}
homer_motif_prefix = 'homer.hg38.191020.'
homer_motifs_beds = [homer_motif_prefix + i + '.bed' for i in homer_motif_dict.values()]
motif_counts_file = 'dvsm_motifs.counts.txt'
dvsm_motifs_ext25_beds = glob.glob('dsvm*ext25.bed')
dvsm_motif_only_beds = glob.glob('dsvm*.motif_only.bed')
homer_motif_nomerge_prefix = 'homer.hg38.no_merge.'
homer_motifs_nomerge_beds = [homer_motif_nomerge_prefix + i + '.bed' for i in homer_motif_dict.values()]
homer_motifs_ext25_beds = [homer_motif_nomerge_prefix + i + '.ext25.bed' for i in homer_motif_dict.values()]


##make dvsm bed files
# convert_dvsm_to_bed(dsvm_leah, dsvm_random_all_bed, dsvm_random_neg_bed)
# convert_dvsm_to_bed(dsvm_scrambled, dsvm_scrambled_all_bed, dsvm_scrambled_neg_bed)

##get motif specific bed files, with and without merge
##CRX, OTX2, NRL, MEF2D, RORB, CREB, CTCF
# get_specific_motifs_from_homer_bed(homer_bed, homer_motif_nomerge_prefix, homer_motif_prefix, homer_motif_dict)

##bedtools intersect dsvm scores with homer beds
# print(homer_motifs_beds)
# bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_beds, '.bed')

##get average dvsm per motif for the different score
# get_counts_per_dvsm_and_motif(dsvm_bed_files, homer_motifs_beds, motif_counts_file)


##look at specific postions within motif

##extend motif beds +- 25
# bedtools_slop_motif_beds(homer_motifs_nomerge_beds)
##bedtools intersect dsvm scores with homer beds
# print(homer_motifs_beds)
# bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_ext25_beds, '.ext25.bed')
# bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_nomerge_beds, '.motif_only.bed')


##then count data and make graphs

##using +- 25bps bed files
# print(dvsm_motifs_ext25_beds)
# get_ave_dvsm_per_bp_around_motif(dvsm_motifs_ext25_beds)
# get_ave_dvsm_per_bp_only_motif(dvsm_motif_only_beds)


##repeat analysis 0721 with new dvsm 
dsvm_0721 = 'ATAC_dsvm20pout_simple.txt'
dsvm_0721_all_bed = 'ATAC_dsvm20pout.all_scores.bed'
dsvm_0721_neg_bed = 'ATAC_dsvm20pout.neg_scores.bed'
dsvm_bed_files = [dsvm_0721_all_bed, dsvm_0721_neg_bed]
motif_counts_0721_file = 'ATAC_dsvm20pout_motifs_0721.counts.txt'
atac_dvsm_motifs_ext25_beds = glob.glob('ATAC_dsvm20pout*ext25.bed')
atac_ddvsm_motif_only_beds = glob.glob('ATAC_dsvm20pout*.motif_only.bed')


##make dvsm bed files
# convert_dvsm_to_bed(dsvm_0721, dsvm_0721_all_bed, dsvm_0721_neg_bed)

##bedtools intersect dsvm scores with homer beds
# print(homer_motifs_beds)
# bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_beds, '.bed')

##get average dvsm per motif for the different score
# get_counts_per_dvsm_and_motif(dsvm_bed_files, homer_motifs_beds, motif_counts_0721_file)


##bedtools intersect dsvm scores with homer beds
# print(homer_motifs_beds)
# bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_ext25_beds, '.ext25.bed')
# bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_nomerge_beds, '.motif_only.bed')

##then count data and make graphs

##using +- 25bps bed files
# print(atac_dvsm_motifs_ext25_beds)
# get_ave_dvsm_per_bp_around_motif(atac_dvsm_motifs_ext25_beds)
# get_ave_dvsm_per_bp_only_motif(atac_ddvsm_motif_only_beds)


##repeat analysis 0821 with new dvsm 
dsvm_0821 = 'ATAC_shuffled_dsvm_simplified.txt'
dsvm_0821_all_bed = 'ATAC_shuffled_dsvm.all_scores.bed'
dsvm_0821_neg_bed = 'ATAC_shuffled_dsvm.neg_scores.bed'
dsvm_bed_files = [dsvm_0821_all_bed, dsvm_0821_neg_bed]
motif_counts_0821_file = 'ATAC_shuffled_dsvm_motifs_0821.counts.txt'
atac_dvsm_motifs_ext25_beds = glob.glob('ATAC_shuffled_dsvm*ext25.bed')
atac_ddvsm_motif_only_beds = glob.glob('ATAC_shuffled_dsvm*.motif_only.bed')


##make dvsm bed files
# convert_dvsm_to_bed(dsvm_0821, dsvm_0821_all_bed, dsvm_0821_neg_bed)

##bedtools intersect dsvm scores with homer beds
# print(homer_motifs_beds)
# bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_beds, '.bed')

##get average dvsm per motif for the different score
# get_counts_per_dvsm_and_motif(dsvm_bed_files, homer_motifs_beds, motif_counts_0821_file)


##bedtools intersect dsvm scores with homer beds
# print(homer_motifs_beds)
# bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_ext25_beds, '.ext25.bed')
# bedtools_intersect_dsvm_motif(dsvm_bed_files, homer_motifs_nomerge_beds, '.motif_only.bed')

##then count data and make graphs

##using +- 25bps bed files and just motifs
print(atac_dvsm_motifs_ext25_beds)
get_ave_dvsm_per_bp_around_motif(atac_dvsm_motifs_ext25_beds)
get_ave_dvsm_per_bp_only_motif(atac_ddvsm_motif_only_beds)

