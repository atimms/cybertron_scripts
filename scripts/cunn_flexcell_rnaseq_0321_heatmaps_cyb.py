#!/usr/bin/python
import os
import subprocess
import glob
import shutil
import math

##set input variables and parameters
delim = '\t'
sample_dict = {'OST08': 'CTRL', 'OST16': 'CTRL', 'OST22': 'CTRL', 'C2064': 'FLNA', 'C1722': 'FLNA', 'C3049': 'FLNA', 
		'C2012': 'FLNB', 'C4019': 'FLNB', 'C3066': 'FLNB', 'C1589': 'FLNC', 'C1677': 'FLNC', 'C1192': 'FLNC', 
		'C2013': 'PIEZO1', 'C1957': 'PIEZO1', 'C1001': 'PIEZO1', 'C1653': 'PIEZO1'}

def filter_get_ratios_from_norm_counts(norm_files, counts_req):
	for norm_file in norm_files:
		ratio_dict, logfc_dict = {}, {}
		new_norm_file = norm_file.rsplit('.', 1)[0] + '.counts.txt'
		ratio_file = norm_file.rsplit('.', 1)[0] + '.ratios.txt'
		log2fc_file = norm_file.rsplit('.', 1)[0] + '.log2fc.txt'
		with open(norm_file, "r") as nf_fh, open(new_norm_file, "w") as nnf_fh:
			lc = 0
			for line in nf_fh:
				line = line.replace('"', '').rstrip().split(',')
				lc += 1
				if lc == 1:
					ratio_samples = line[1::2]
					ratio_sample_list = [i.split('_')[0] + '_' + sample_dict[i.split('_')[0]] for i in ratio_samples]
					all_sample_list = [i + '_' + sample_dict[i.split('_')[0]] for i in line[1:]]
					nnf_fh.write(delim.join(['gene'] + all_sample_list) + '\n')
					print(ratio_samples, ratio_sample_list, all_sample_list)
				else:
					gene = line[0]
					first_counts = [float(i) for i in line[1:]]
					##covert all 
					counts = [0.01 if x==0 else x for x in first_counts]
					##if any cellline has norm counts >50
					# if max(counts) >= 50 and min(counts) > 0:
					if max(counts) >= counts_req:
						nnf_fh.write(delim.join(line) + '\n')
						# print(line)
						if counts[0] >= counts[1]:
							ratio1 = counts[0]/counts[1]
						else:
							ratio1 = -(counts[1]/counts[0])
						if counts[2] >= counts[3]:
							ratio2 = counts[2]/counts[3]
						else:
							ratio2 = -(counts[3]/counts[2])
						if counts[4] >= counts[5]:
							ratio3 = counts[4]/counts[5]
						else:
							ratio3 = -(counts[5]/counts[4])
						if counts[6] >= counts[7]:
							ratio4 = counts[6]/counts[7]
						else:
							ratio4 = -(counts[7]/counts[6])
						if counts[8] >= counts[9]:
							ratio5 = counts[8]/counts[9]
						else:
							ratio5 = -(counts[9]/counts[8])
						if counts[10] >= counts[11]:
							ratio6 = counts[10]/counts[11]
						else:
							ratio6 = -(counts[11]/counts[10])
						
						if len(ratio_sample_list) == 7:
							if counts[12] >= counts[13]:
								ratio7 = counts[12]/counts[13]
							else:
								ratio7 = -(counts[13]/counts[12])
							ratios = [str(ratio1), str(ratio2), str(ratio3), str(ratio4), str(ratio5), str(ratio6), str(ratio7)]
						else:
							ratios = [str(ratio1), str(ratio2), str(ratio3), str(ratio4), str(ratio5), str(ratio6)]
						ratio_dict[gene] = ratios
						log2fc1 = math.log(counts[0]/counts[1], 2)
						log2fc2 = math.log(counts[2]/counts[3], 2)
						log2fc3 = math.log(counts[4]/counts[5], 2)
						log2fc4 = math.log(counts[6]/counts[7], 2)
						log2fc5 = math.log(counts[8]/counts[9], 2)
						log2fc6 = math.log(counts[10]/counts[11], 2)
						
						if len(ratio_sample_list) == 7:
							log2fc7 = math.log(counts[12]/counts[13], 2)
							logfcs = [str(log2fc1), str(log2fc2), str(log2fc3), str(log2fc4), str(log2fc5), str(log2fc6), str(log2fc7)]
						else:
							logfcs = [str(log2fc1), str(log2fc2), str(log2fc3), str(log2fc4), str(log2fc5), str(log2fc6)]

						logfc_dict[gene] = logfcs
						# print(ratios)
		with open(ratio_file, "w") as rf_fh, open(log2fc_file, "w") as lfc_fh:
			rf_fh.write(delim.join(['gene'] + ratio_sample_list) + '\n')
			lfc_fh.write(delim.join(['gene'] + ratio_sample_list) + '\n')
			for gene in ratio_dict:
				rf_fh.write(delim.join([gene] + ratio_dict[gene]) + '\n')
			for gene in logfc_dict:
				lfc_fh.write(delim.join([gene] + logfc_dict[gene]) + '\n')


def make_genelist_dict(in_file):
	gl_dict = {}
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			gl_name = line[0]
			genes = line[1:]
			gl_dict[gl_name] = genes
	# print(gl_dict)
	return(gl_dict)

def filter_by_genelists(in_files, gl_dict):
	for in_file in in_files:
		for gl in gl_dict:
			out_file = in_file.rsplit('.', 1)[0] + '.' + gl + '.heatmap.txt'
			with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
				lc = 0
				for line in in_fh:
					lc += 1
					if lc == 1:
						out_fh.write(line)
					else:
						gene = line.split(delim)[0]
						if gene in gl_dict[gl]:
							out_fh.write(line)


def combine_ratio_heatmap_files(filename_parts, gene_group_names, genelist_names, sample_names):
	for genelist_name in genelist_names:
		gene_dict = {}
		out_file = 'all_samples' + filename_parts[1] + genelist_name + filename_parts[2]
		print(out_file)
		for gene_group_name in gene_group_names:
			in_file = filename_parts[0] + gene_group_name + filename_parts[1] + genelist_name + filename_parts[2]
			print(in_file)
			lc = 0
			with open(in_file, "r") as in_fh:
				for line in in_fh:
					line = line.rstrip().split(delim)
					lc += 1
					if lc == 1:
						samples_from_header = line[1:]
					else:
						##set up dict
						gene = line[0]
						if gene not in gene_dict:
							gene_dict[gene] = ['.'] * len(sample_names)
						counts = line[1:]
						##add counts to dict
						for sample in samples_from_header:
							count = counts[samples_from_header.index(sample)]
							gene_dict[gene][sample_names.index(sample)] = count

		with open(out_file, "w") as out_fh:
			out_fh.write(delim.join(['gene'] + sample_names) + '\n')
			for g in gene_dict:
				if '.' not in gene_dict[g]:
					# print(g, gene_dict[g])
					out_fh.write(delim.join([g] + gene_dict[g]) + '\n')








##run methods
working_dir = '/home/atimms/ngs_data/rnaseq/cunn_flexcell_rnaseq_0321'
os.chdir(working_dir)

norm_counts_files = ['cunn_flexcell_rnaseq_0321_ctl_flna.norm_counts.csv', 'cunn_flexcell_rnaseq_0321_ctl_flnb.norm_counts.csv', 
	'cunn_flexcell_rnaseq_0321_ctl_flnc.norm_counts.csv', 'cunn_flexcell_rnaseq_0321_ctl_piezo1.norm_counts.csv']

vst_counts_files = ['cunn_flexcell_rnaseq_0321_ctl_flna.vst_counts.csv', 'cunn_flexcell_rnaseq_0321_ctl_flnb.vst_counts.csv', 
	'cunn_flexcell_rnaseq_0321_ctl_flnc.vst_counts.csv', 'cunn_flexcell_rnaseq_0321_ctl_piezo1.vst_counts.csv']

genelist_file = 'genelist_040821.txt'

##test
# vst_counts_files = ['cunn_flexcell_rnaseq_0321_ctl_flna.vst_counts.csv']
# norm_counts_files = ['cunn_flexcell_rnaseq_0321_ctl_flna.norm_counts.csv']
# genelist_file = 'gl_temp.txt'



##filter that at least one cell line has counts >50 and get ratios/counts for each sample
# filter_get_ratios_from_norm_counts(norm_counts_files, 50)
# filter_get_ratios_from_norm_counts(vst_counts_files, 5)

##make dict of all genelist
genelist_dict = make_genelist_dict(genelist_file) 

##filter counts/ratio files for each gene list
counts_and_ratio_files = glob.glob('*.norm_counts.*.txt')
# counts_and_ratio_files = glob.glob('*.vst_counts.*.txt')
# print(counts_and_ratio_files)
# filter_by_genelists(counts_and_ratio_files, genelist_dict)

##combine ratio files in order
ratio_hm_filenames = ['cunn_flexcell_rnaseq_0321_ctl_', '.norm_counts.ratios.', '.heatmap.txt']
gene_groups = ['flna', 'flnb', 'flnc', 'piezo1']
genelists = genelist_dict.keys()
samples = ['OST08_CTRL', 'OST16_CTRL', 'OST22_CTRL', 'C1001_PIEZO1', 'C2064_FLNA', 'C1722_FLNA', 'C3066_FLNB',
		'C1589_FLNC', 'C1677_FLNC', 'C1192_FLNC', 'C2013_PIEZO1', 'C1957_PIEZO1', 'C1653_PIEZO1', 'C3049_FLNA',
		'C2012_FLNB', 'C4019_FLNB']
combine_ratio_heatmap_files(ratio_hm_filenames, gene_groups, genelists, samples)


###rm *.norm_counts.*.txt




