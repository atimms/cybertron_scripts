#!/usr/bin/python
import os
import subprocess
import glob
import shutil
import math

##set input variables and parameters
delim = '\t'

def filter_get_ratios_from_norm_counts(norm_file, counts_req):
	ratio_dict = {}
	ratio_file = norm_file.rsplit('.', 1)[0] + '.ratios.txt'
	with open(norm_file, "r") as nf_fh, open(ratio_file, "w") as rf_fh:
		lc = 0
		for line in nf_fh:
			line = line.replace('"', '').rstrip().split(',')
			lc += 1
			if lc == 1:
				ratio_samples = line[1::2]
				# ratio_sample_list = [i.split('_')[0] + '_' + sample_dict[i.split('_')[0]] for i in ratio_samples]
				ratio_sample_list = [i.replace('_Strained', '')for i in ratio_samples]
				print(ratio_samples, ratio_sample_list)
				rf_fh.write(delim.join(['gene'] + ratio_sample_list) + '\n')
			else:
				gene = line[0]
				first_counts = [float(i) for i in line[1:]]
				##covert all 
				counts = [0.01 if x==0 else x for x in first_counts]
				##if any cellline has norm counts >50
				# if max(counts) >= 50 and min(counts) > 0:
				if max(counts) >= counts_req:
					# print(gene, max(counts))
					ratios = []
					##make list 0-x for each sample and then split into lists of 2 values
					total_sample_numbers = list(range(0, len(line[1:])))
					sample_groups = [total_sample_numbers[i:i + 2] for i in range(0, len(total_sample_numbers), 2)]
					# print(total_sample_numbers, sample_groups)
					for sample_group in sample_groups:
						# print(line)
						if counts[sample_group[0]] >= counts[sample_group[1]]:
							ratio1 = counts[sample_group[0]]/counts[sample_group[1]]
						else:
							ratio1 = -(counts[sample_group[1]]/counts[sample_group[0]])
						ratios.append(ratio1)
					ratios_out = [str(i) for i in ratios]
					rf_fh.write(delim.join([gene] + ratios_out) + '\n')


def filter_get_ratios_from_norm_counts_matrix(norm_file, counts_req):
	ratio_dict = {}
	ratio_file = norm_file.rsplit('.', 1)[0] + '.ratios.txt'
	with open(norm_file, "r") as nf_fh, open(ratio_file, "w") as rf_fh:
		lc = 0
		for line in nf_fh:
			line = line.replace('"', '').rstrip().split(',')
			lc += 1
			if lc == 1:
				ratio_samples = line[1::2]
				# ratio_sample_list = [i.split('_')[0] + '_' + sample_dict[i.split('_')[0]] for i in ratio_samples]
				ratio_sample_list = [i.replace('_100_', '_')for i in ratio_samples]
				ratio_sample_list = [i.replace('_0_', '_')for i in ratio_sample_list]
				print(ratio_samples, ratio_sample_list)
				rf_fh.write(delim.join(['gene'] + ratio_sample_list) + '\n')
			else:
				gene = line[0]
				first_counts = [float(i) for i in line[1:]]
				##covert all 
				counts = [0.01 if x==0 else x for x in first_counts]
				##if any cellline has norm counts >50
				# if max(counts) >= 50 and min(counts) > 0:
				if max(counts) >= counts_req:
					# print(gene, max(counts))
					ratios = []
					##make list 0-x for each sample and then split into lists of 2 values
					total_sample_numbers = list(range(0, len(line[1:])))
					sample_groups = [total_sample_numbers[i:i + 2] for i in range(0, len(total_sample_numbers), 2)]
					# print(total_sample_numbers, sample_groups)
					for sample_group in sample_groups:
						# print(line)
						if counts[sample_group[0]] >= counts[sample_group[1]]:
							ratio1 = counts[sample_group[0]]/counts[sample_group[1]]
						else:
							ratio1 = -(counts[sample_group[1]]/counts[sample_group[0]])
						ratios.append(ratio1)
					ratios_out = [str(i) for i in ratios]
					rf_fh.write(delim.join([gene] + ratios_out) + '\n')


def get_genenames_from_file(infile):
	genes = []
	with open(infile, "r") as in_fh:
		for line in in_fh:
			gene = line.rstrip()
			genes.append(gene)
	return(genes)

def filter_ratio_file_genelist(ratio_file, genelists):
	for genelist in genelists:
		gl_genes = get_genenames_from_file(genelist)
		out_file = ratio_file.split('.', 1)[0] + '.' + genelist.split('.', 1)[0]  + '.all_samples.txt'
		with open(ratio_file, "r") as in_fh, open(out_file, "w") as out_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				if lc == 1:
					out_fh.write(line)
				else:
					gene = line.split(delim)[0]
					if gene in gl_genes:
						out_fh.write(line)

def filter_order_ratio_files_by_sample(ratio_files, samplelist_dict):
	for ratio_file in ratio_files:
		for samplelist in samplelist_dict:
			samples_wanted = samplelist_dict[samplelist]
			out_file = ratio_file.rsplit('.', 2)[0] + '.' + samplelist + '.heatmap.txt'
			with open(ratio_file, "r") as in_fh, open(out_file, "w") as out_fh:
				out_fh.write(delim.join(['gene'] + samples_wanted) + '\n')
				lc = 0
				for line in in_fh:
					line = line.rstrip().split(delim)
					lc += 1
					if lc == 1:
						header = line
					else:
						gene = line[0]
						ratios = []
						for sample_wanted in samples_wanted:
							ratio = line[header.index(sample_wanted)]
							ratios.append(ratio)
						out_fh.write(delim.join([gene] + ratios) + '\n')


##run methods

##flexcell data
working_dir = '/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0621'
os.chdir(working_dir)

norm_counts_file = 'cunn_flexcell_rnaseq_0621.norm_counts.csv'
norm_count_ratio_file = 'cunn_flexcell_rnaseq_0621.norm_counts.ratios.txt'
genelist_files = ['GOBP_osteoblast_DIFFERENTIATION.txt', 'GOBP_CELL_matrix_ADHESION.txt',
		'GOBP_chondrocyte_DIFFERENTIATION.txt', 'GOBP_FOCAL_ADHESION_ASSEMBLY.txt']
genelist_ratio_files = ['cunn_flexcell_rnaseq_0621.GOBP_osteoblast_DIFFERENTIATION.all_samples.txt', 'cunn_flexcell_rnaseq_0621.GOBP_CELL_matrix_ADHESION.all_samples.txt',
		'cunn_flexcell_rnaseq_0621.GOBP_chondrocyte_DIFFERENTIATION.all_samples.txt', 'cunn_flexcell_rnaseq_0621.GOBP_FOCAL_ADHESION_ASSEMBLY.all_samples.txt']
sample_dict = {'all': ['OST08', 'OST16', 'OST22', 'OST31', 'OST35', 'C2064', 'C1722', 'C3049_1', 'C3049_2', 'C2012_1', 'C2012_2', 'C4019_1', 'C4019_2', 'C3066_1', 'C3066_2', 'C2013_1', 'C2013_2', 'C1001', 'C1653_1', 'C1653_2', 'C1957'],
		'original': ['OST08', 'OST16', 'OST22', 'C2064', 'C1722', 'C3049_1', 'C2012_1', 'C4019_1', 'C3066_1', 'C2013_1', 'C1001', 'C1653_1', 'C1957'],
		'new': ['OST31', 'OST35', 'C3049_2', 'C2012_2', 'C4019_2', 'C3066_2', 'C2013_2', 'C1653_2']}


##filter that at least one cell line has counts >50 and get ratios for each sample
# filter_get_ratios_from_norm_counts(norm_counts_file, 50)

##filter by genelist
# filter_ratio_file_genelist(norm_count_ratio_file, genelist_files)

##filter/order by sample
# filter_order_ratio_files_by_sample(genelist_ratio_files, sample_dict)




##matrix data
working_dir = '/home/atimms/ngs_data/rnaseq/cunn_rnaseq_0220'
os.chdir(working_dir)

#renamed and removed 15 data using cut
norm_counts_file = 'matrix_rnaseq_0220.norm_counts.csv' 
norm_count_ratio_file = 'matrix_rnaseq_0220.norm_counts.ratios.txt'
genelist_files = ['GOBP_osteoblast_DIFFERENTIATION.txt', 'GOBP_CELL_matrix_ADHESION.txt',
		'GOBP_chondrocyte_DIFFERENTIATION.txt', 'GOBP_FOCAL_ADHESION_ASSEMBLY.txt']
genelist_ratio_files = ['matrix_rnaseq_0220.GOBP_osteoblast_DIFFERENTIATION.all_samples.txt', 'matrix_rnaseq_0220.GOBP_CELL_matrix_ADHESION.all_samples.txt',
		'matrix_rnaseq_0220.GOBP_chondrocyte_DIFFERENTIATION.all_samples.txt', 'matrix_rnaseq_0220.GOBP_FOCAL_ADHESION_ASSEMBLY.all_samples.txt']
sample_dict = {'all': ['S_AUT25_1', 'S_AUT25_2', 'S_AUT25_3', 'S_OST85_1', 'S_OST85_2', 'S_OST85_3', 'S_STL19_1', 'S_STL19_2', 
		'S_STL19_3', 'S_C1625_1', 'S_C1625_2', 'S_C1625_3', 'S_C2084_1', 'S_C2084_2', 'S_C2084_3', 'S_C3049_1', 'S_C3049_2', 
		'S_C3049_3', 'S_1032_1', 'S_1032_2', 'S_1032_3', 'S_4032_1', 'S_4032_2', 'S_4032_3', 'S_C3066_1', 'S_C3066_2', 
		'S_C3066_3', 'S_C1671_1', 'S_C1671_2', 'S_C1671_3', 'S_C1859_1', 'S_C1859_2', 'S_C1859_3', 'S_C1926_1', 'S_C1926_2', 
		'S_C1926_3', 'S_C1957_1', 'S_C1957_2', 'S_C1957_3']}


##filter that at least one cell line has counts >50 and get ratios for each sample
# filter_get_ratios_from_norm_counts_matrix(norm_counts_file, 50)

##filter by genelist
# filter_ratio_file_genelist(norm_count_ratio_file, genelist_files)

##filter/order by sample
# filter_order_ratio_files_by_sample(genelist_ratio_files, sample_dict)


##matrix data average the duplicates
norm_counts_file = 'matrix_rnaseq_0220.norm_counts.csv' 
norm_counts_combined_file = 'matrix_rnaseq_combined_0220.norm_counts.csv' 
norm_count_ratio_file = 'matrix_rnaseq_combined_0220.norm_counts.ratios.txt'
genelist_files = ['GOBP_osteoblast_DIFFERENTIATION.txt', 'GOBP_CELL_matrix_ADHESION.txt',
		'GOBP_chondrocyte_DIFFERENTIATION.txt', 'GOBP_FOCAL_ADHESION_ASSEMBLY.txt']
genelist_ratio_files = ['matrix_rnaseq_combined_0220.GOBP_osteoblast_DIFFERENTIATION.all_samples.txt', 'matrix_rnaseq_combined_0220.GOBP_CELL_matrix_ADHESION.all_samples.txt',
		'matrix_rnaseq_combined_0220.GOBP_chondrocyte_DIFFERENTIATION.all_samples.txt', 'matrix_rnaseq_combined_0220.GOBP_FOCAL_ADHESION_ASSEMBLY.all_samples.txt']
sample_dict = {'averaged': ['S_AUT25', 'S_OST85', 'S_STL19', 'S_C1625', 'S_C2084', 'S_C3049', 'S_1032', 'S_4032', 'S_C3066', 'S_C1671', 'S_C1859', 'S_C1926', 'S_C1957']}


def average_replicates(in_file, out_file):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			line = line.replace('"', '').rstrip().split(',')
			lc += 1
			if lc == 1:
				ratio_samples = line[1::3]
				# ratio_sample_list = [i.split('_')[0] + '_' + sample_dict[i.split('_')[0]] for i in ratio_samples]
				ratio_sample_list = [i.replace('_0_1', '_0')for i in ratio_samples]
				ratio_sample_list = [i.replace('_100_2', '_100')for i in ratio_sample_list]
				print(ratio_samples, ratio_sample_list)
				out_fh.write(','.join(['gene'] + ratio_sample_list) + '\n')
			else:
				gene = line[0]
				counts = [float(i) for i in line[1:]]
				averages = []
				##make list 0-x for each sample and then split into lists of 6 values
				total_sample_numbers = list(range(0, len(line[1:])))
				sample_groups = [total_sample_numbers[i:i + 6] for i in range(0, len(total_sample_numbers), 6)]
				if lc ==2:
					print(total_sample_numbers, sample_groups)
				for sample_group in sample_groups:
					# print(line)
					zeros = [counts[sample_group[0]], counts[sample_group[2]], counts[sample_group[4]]]
					hundereds = [counts[sample_group[1]], counts[sample_group[3]], counts[sample_group[5]]]
					average_0 = sum(zeros) / len(zeros)
					average_100 = sum(hundereds) / len(hundereds)
					averages.extend([average_0, average_100])
				averages_out = [str(i) for i in averages]
				out_fh.write(','.join([gene] + averages_out) + '\n')


def filter_get_ratios_from_norm_counts_matrix_averages(norm_file, counts_req):
	ratio_dict = {}
	ratio_file = norm_file.rsplit('.', 1)[0] + '.ratios.txt'
	with open(norm_file, "r") as nf_fh, open(ratio_file, "w") as rf_fh:
		lc = 0
		for line in nf_fh:
			line = line.replace('"', '').rstrip().split(',')
			lc += 1
			if lc == 1:
				ratio_samples = line[1::2]
				# ratio_sample_list = [i.split('_')[0] + '_' + sample_dict[i.split('_')[0]] for i in ratio_samples]
				ratio_sample_list = [i.replace('_100', '')for i in ratio_samples]
				ratio_sample_list = [i.replace('_0', '')for i in ratio_sample_list]
				print(ratio_samples, ratio_sample_list)
				rf_fh.write(delim.join(['gene'] + ratio_sample_list) + '\n')
			else:
				gene = line[0]
				first_counts = [float(i) for i in line[1:]]
				##covert all 
				counts = [0.01 if x==0 else x for x in first_counts]
				##if any cellline has norm counts >50
				# if max(counts) >= 50 and min(counts) > 0:
				if max(counts) >= counts_req:
					# print(gene, max(counts))
					ratios = []
					##make list 0-x for each sample and then split into lists of 2 values
					total_sample_numbers = list(range(0, len(line[1:])))
					sample_groups = [total_sample_numbers[i:i + 2] for i in range(0, len(total_sample_numbers), 2)]
					# print(total_sample_numbers, sample_groups)
					for sample_group in sample_groups:
						# print(line)
						if counts[sample_group[0]] >= counts[sample_group[1]]:
							ratio1 = counts[sample_group[0]]/counts[sample_group[1]]
						else:
							ratio1 = -(counts[sample_group[1]]/counts[sample_group[0]])
						ratios.append(ratio1)
					ratios_out = [str(i) for i in ratios]
					rf_fh.write(delim.join([gene] + ratios_out) + '\n')

##merge the replicates
# average_replicates(norm_counts_file, norm_counts_combined_file)

##filter that at least one cell line has counts >50 and get ratios for each sample
# filter_get_ratios_from_norm_counts_matrix_averages(norm_counts_combined_file, 50)

##filter by genelist
filter_ratio_file_genelist(norm_count_ratio_file, genelist_files)

##filter/order by sample
filter_order_ratio_files_by_sample(genelist_ratio_files, sample_dict)
