#!/usr/bin/env python
import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from collections import Counter


'''
info:
##need to use correct env
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878
conda activate pd_np_plt_etc
'''
##parameters
delim = '\t'
##programs
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
##files etc
bt_genome = '/home/atimms/ngs_data/references/10x/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.chrom_sizes'




##methods
def make_combined_region_bed_summarize_data(in_files, region_bed):
	region_dict = {}
	for in_file in in_files:
		with open(in_file, "r") as in_fh:
			for line in in_fh:
				# print(line)
				line = line.rstrip().split(delim)
				chrom = line[0]	
				start = int(line[1])
				end = int(line[2])
				region = line[3]
				if region in region_dict:
					if chrom == region_dict[region][0]:
						region_dict[region][3] += 1
						if start < region_dict[region][1]:
							region_dict[region][1] = start
						if end > region_dict[region][2]:
							region_dict[region][2] = end
					else:
						print('issue with:', line)						

				else:
					region_dict[region] = [chrom, start, end, 1]
	with open(region_bed, "w") as out_fh:
		for r in region_dict:
			print(region_dict[r] + [r])
			line_out = [str(i) for i in region_dict[r] + [r]]
			out_fh.write(delim.join(line_out) + '\n')

def make_combined_snp_bed(in_files, out_file):
	with open(out_file, "w") as out_fh:
		for in_file in in_files:
			with open(in_file, "r") as in_fh:
				for line in in_fh:
					line = line.rstrip().split(delim)
					line_out = line[:3] + [line[3].split('_')[0]]
					out_fh.write(delim.join(line_out) + '\n')



def pandas_to_summarize_comb_bed(c_bed):
	##load data and format
	c_bed = 'combined_regions.bed'
	c_data = pd.read_table(c_bed, header=None)
	##add header
	c_data.columns=["chr", "start", "end", "snp_count", "region"]
	##add 2 new columns
	c_data['size'] = c_data.end - c_data.start
	c_data['study'] = c_data['region'].str.split("_").str[0]
	# c_data.head()
	##graph using seaborn
	fig, ax=plt.subplots(1,2)
	sns.boxplot(x="study", y="snp_count", data=c_data, ax=ax[0])
	sns.boxplot(x="study", y="size", data=c_data, ax=ax[1])
	plt.savefig("study_regions.boxplots.pdf")

def bt_regions_vs_cc_summit_peaks(region_bed, summit_peak_beds, out_bed):
	##bedtools intersect 
	with open(out_bed, "w") as out_fh: 
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', region_bed, '-b'] + summit_peak_beds + ['-C', '-filenames'], stdout=out_fh)
		hom_bt_intersect.wait()

def bt_snps_vs_cc_extended_peaks(region_bed, summit_peak_beds, out_bed):
	##bedtools intersect 
	with open(out_bed, "w") as out_fh: 
		# hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', region_bed, '-b'] + summit_peak_beds + ['-wa', '-wb', '-filenames'], stdout=out_fh)
		# hom_bt_intersect.wait()
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', region_bed, '-b'] + summit_peak_beds + ['-C', '-filenames'], stdout=out_fh)
		hom_bt_intersect.wait()

def graph_heatmaps_regions_cc_peaks(int_file):
	int_data = pd.read_table(int_file, header=None)
	##add header
	int_data.columns=["chr", "start", "end", "snp_count", "region", "filename", "peak_count"]
	##add 3 new columns
	##split filename to just peak names
	int_data['cc_peaks'] = int_data['filename'].str.split("/").str[-1].str.rsplit(".",3).str[0]

	##normalize peak count
	int_data['cc_peaks_snp_count'] = int_data.peak_count / int_data.snp_count
	int_data['cc_peaks_length'] = int_data.peak_count / (int_data.end - int_data.start)

	##split peak count column every 54th col i.e. how many peak beds
	row_names = int_data.region.unique()
	col_names = int_data.cc_peaks.unique()
	df_for_hm = pd.DataFrame(int_data.peak_count.values.reshape(-1, 54))
	df_for_hm.columns=col_names
	df_for_hm['region'] = row_names
	df_for_hm = df_for_hm.set_index('region')
	df_for_hm.head()
	##ussing normalized data
	df_for_hm2 = pd.DataFrame(int_data.cc_peaks_snp_count.values.reshape(-1, 54))
	df_for_hm2.columns=col_names
	df_for_hm2['region'] = row_names
	df_for_hm2 = df_for_hm2.set_index('region')
	df_for_hm2.head()
	df_for_hm3 = pd.DataFrame(int_data.cc_peaks_length.values.reshape(-1, 54))
	df_for_hm3.columns=col_names
	df_for_hm3['region'] = row_names
	df_for_hm3 = df_for_hm3.set_index('region')
	df_for_hm3.head()
	plt.figure(figsize=(16,9))
	sns.heatmap(df_for_hm)
	plt.savefig("regions.summit_peaks.heatmap.pdf")
	plt.figure(figsize=(16,9))
	sns.heatmap(df_for_hm2)
	plt.savefig("regions.summit_peaks.heatmap.norm_snpcount.pdf")
	plt.figure(figsize=(16,9))
	sns.heatmap(df_for_hm3)
	plt.savefig("regions.summit_peaks.heatmap.norm_size.pdf")


def extend_summit_peak_files(in_beds):
	for in_bed in in_beds:
		##bedtools slop
		out_bed = in_bed.split('/')[-1].rsplit('.', 1)[0] + '.ext150.bed'
		with open(out_bed, "w") as out_fh: 
			hom_bt_intersect = subprocess.Popen([bedtools, 'slop', '-i', in_bed, '-g', bt_genome, '-b', '150'], stdout=out_fh)
			hom_bt_intersect.wait()


def compare_snps_vs_cc_extended_peaks(in_file, out_file):
	snp_count_dict = {}
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			# print(line)
			line = line.rstrip().split(delim)
			snp = '_'.join(line[:3])
			study = line[3]
			peak = line[4].split('.')[0] + '.' + line[4].split('.')[1]
			count = int(line[5])
			
			if study in snp_count_dict:
				##if snp seen before
				if snp in snp_count_dict[study]:
					if count == 0:
						snp_count_dict[study][snp][1].append(peak)
					else:
						snp_count_dict[study][snp][0] += 1
						snp_count_dict[study][snp][2].append(peak)
				else:
					if count == 0:
						snp_count_dict[study][snp] = [0, [peak], []]
					else:
						snp_count_dict[study][snp] = [1, [], [peak]]	
			##if study not seen before 
			else:
				if count == 0:
					snp_count_dict[study] = {snp:[0, [peak], []]}
				else:
					snp_count_dict[study] = {snp:[1, [], [peak]]}
	# print(snp_count_dict)
	cc_study_count_dict = {}
	header = ['cell class']
	for st in snp_count_dict:
		header.append(st)
		print('total snps in', st, "=", len(snp_count_dict[st]))
		total_snps = len(snp_count_dict[st])
		all_study_peaks_with_snps = []
		for s in snp_count_dict[st]:
			# print(st, s, snp_count_dict[st][s][0], len(snp_count_dict[st][s][1]), len(snp_count_dict[st][s][2]))
			all_study_peaks_with_snps.extend(snp_count_dict[st][s][2])
		# print(st, len(Counter(all_study_peaks_with_snps)), Counter(all_study_peaks_with_snps))
		##get counts of all snp in cc peaks
		snp_peak_counter = Counter(all_study_peaks_with_snps)
		# print(snp_peak_counter)
		for peak in snp_peak_counter:
			# print(peak, snp_peak_counter[peak])
			intersect_count = snp_peak_counter[peak]
			intersect_fraction = (intersect_count/total_snps) * 100
			# print(peak, intersect_count, intersect_fraction, total_snps, st)
			if peak in cc_study_count_dict:
				cc_study_count_dict[peak].append(intersect_fraction)
			else:
				cc_study_count_dict[peak] = [intersect_fraction]
	# print(header)
	# print(cc_study_count_dict)
	with open(out_file, "w") as out_fh:
		out_fh.write(delim.join(header) + '\n')
		for s in cc_study_count_dict:
			print(s, cc_study_count_dict[s])
			line_out = [s] + [str(i) for i in cc_study_count_dict[s]]
			out_fh.write(delim.join(line_out) + '\n')

def graph_heatmaps_snps_cc_peaks(int_file, plot_name):
	int_data = pd.read_table(int_file, index_col = 0)
	sns.clustermap(int_data, yticklabels = 1)
	plt.savefig(plot_name)

##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_gwas_peaks_1220'
os.chdir(working_dir)

##bed files with gwas snps
bed_files = ['fritsche_2016.bed', 'scerri_2017.bed', 'jansen_2019.bed', 'kunkle_2019.bed']
##removed region from janson that squed the data and repeat
adjusted_bed_files = ['fritsche_2016.bed', 'scerri_2017.bed', 'jansen_2019_alt.bed', 'kunkle_2019.bed']

##region/peak files
combined_snp_bed = 'combined_snps.bed'
combined_region_bed = 'combined_regions.bed'

peak_beds_dir = '/active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/archr_macs2_peak_beds_1020/'
summit_peak_beds = glob.glob(peak_beds_dir + 'human*') + glob.glob(peak_beds_dir + 'organoid*')
regions_summit_bt_bed = 'combined_regions.summit_peaks.bt_int.bed'
extended_peak_beds = glob.glob('human*bed') + glob.glob('organoid*bed')
snps_extended_bt_bed = 'combined_snps.extended_peaks.bt_int.bed'
snps_extended_sum = 'combined_snps.extended_peaks.counts_for_heatmap.txt'
snps_extended_heatmap = "snps.cc_peaks.heatmap.pdf"
snps_extended_sum_norm = 'combined_snps.extended_peaks.counts_for_heatmap.norm_peak_count.txt'
snps_extended_norm_heatmap = "snps.cc_peaks.heatmap.norm_peak_count.pdf"

##region analysis
## 1. make a combined region bed file from all snp beds
# make_combined_region_bed_summarize_data(adjusted_bed_files, combined_region_bed)

## 2. summarize data from different studies
# pandas_to_summarize_comb_bed(combined_region_bed)

## 3. bedtools on regions against peak files
# print(summit_peak_beds)
# bt_regions_vs_cc_summit_peaks(combined_region_bed, summit_peak_beds, regions_summit_bt_bed)

## 4. make heatmaps from region/peak data
# graph_heatmaps_regions_cc_peaks(regions_summit_bt_bed)

##snp analysis
## 1. extend summit peak bed files
# extend_summit_peak_files(summit_peak_beds)

## 2. combine all snp files
# make_combined_snp_bed(bed_files, combined_snp_bed)

## 3. intersect/compare summit peaks with gwas SNPs
# print(extended_peak_beds)
##use bt intersect
# bt_snps_vs_cc_extended_peaks(combined_snp_bed, extended_peak_beds, snps_extended_bt_bed)
##compare studies
# compare_snps_vs_cc_extended_peaks(snps_extended_bt_bed, snps_extended_sum)
##draw heatmap
graph_heatmaps_snps_cc_peaks(snps_extended_sum, snps_extended_heatmap)
graph_heatmaps_snps_cc_peaks(snps_extended_sum_norm, snps_extended_norm_heatmap)

## 3. motifs in summit peaks: intersect motifs with summit peaks
## motif bed from scanMotifGenomeWide.pl or from homer.KnownMotifs.hg38.191020.bed 
## OTX2, CRX, MEF2C, MEF2D, MEF2A, PAX6, RORB, NRL, NR2E3, VSX2, POU4F2, POU4F1, TFAP2A/B/C, LHX1, ONECUT1/2/3, ESRRB, SOX2, LHX2, LHX4, RAX




