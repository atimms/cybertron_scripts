#!/usr/bin/env python
import sys
import subprocess
import os
import glob


##note
'''
analysis of f1

load modules:
module load biobuilds/2017.11
module load local_python/3.6.4

'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_pisces_1218/working'
os.chdir(working_dir)


def filter_pisces_files(in_files, mtor_genes, summary_file):
	with open(summary_file, "w") as sum_fh:
		sum_fh.write(delim.join(['sample', 'total', 'mtor', 'clinvar_pathogenic', 'clinvar_pathogenic_or_uncertain', 'cosmic83_5', 'cosmic83_10']) + '\n')
		for in_file in in_files:
			line_count, mtor_lists, cosmic5_list, cosmic10_list, clinvar_path_lists,clinvar_p_and_un_lists  = 0,[],[],[],[],[]
			sample = in_file.split('.')[0]
			mtor_out_file = in_file.rsplit('.', 1)[0] + '.mtor.xls'
			cosmic5_out_file = in_file.rsplit('.', 1)[0] + '.cosmic83_5.xls'
			cosmic10_out_file = in_file.rsplit('.', 1)[0] + '.cosmic83_10.xls'
			clin_p_out_file = in_file.rsplit('.', 1)[0] + '.clinvar_pathogenic.xls'
			clin_p_un_out_file = in_file.rsplit('.', 1)[0] + '.clinvar_pathogenic_and_uncertain.xls'
			print(in_file)
			with open(in_file, "r") as in_fh:
				for line in in_fh:
					line_count +=1
					if line_count == 1:
						header = line
					else:
						line = line.split(delim)
						gene = line[6]
						clinsig = line[66]
						cosmic_coding = line[72]
						rmsk = line[10]
						supdup = line[11]
						##remove duplicates
						# if rmsk == '.' and supdup == '.':
						##mtor gene
						if gene in mtor_genes:
							mtor_lists.append(line)
						##if pathogenic in clinsigdb
						if 'Pathogenic' in clinsig:
							clinvar_path_lists.append(line)
						##if pathogenic or uncertain
						if 'Pathogenic' in clinsig or 'Uncertain' in clinsig:
							clinvar_p_and_un_lists.append(line)
						##count occurances in cosmic
						if cosmic_coding != '.':
							cosmic_coding = cosmic_coding.split(';')[1]
							# print(cosmic_coding)	
							occurances = [int(i) for i in cosmic_coding if i.isdigit()]
							# print(occurances)
							occurances =sum(occurances)
							# print(occurances)
							if occurances >= 5:
								cosmic5_list.append(line)
							if occurances >= 10:
								cosmic10_list.append(line)
			info = [sample, line_count -1,len(mtor_lists), len(clinvar_path_lists), len(clinvar_p_and_un_lists), len(cosmic5_list), len(cosmic10_list)]
			info = [str(i) for i in info]
			sum_fh.write(delim.join(info) +'\n')
			print(info)
			##write out files 
			if len(mtor_lists) >= 1:
				with open(mtor_out_file, "w") as out_fh:
					out_fh.write(header)
					for line_out in mtor_lists:
						out_fh.write(delim.join(line_out))
			if len(clinvar_path_lists) >= 1:
				with open(clin_p_out_file, "w") as out_fh:
					out_fh.write(header)
					for line_out in clinvar_path_lists:
						out_fh.write(delim.join(line_out))
			if len(clinvar_p_and_un_lists) >= 1:
				with open(clin_p_un_out_file, "w") as out_fh:
					out_fh.write(header)
					for line_out in clinvar_p_and_un_lists:
						out_fh.write(delim.join(line_out))
			if len(cosmic5_list) >= 1:
				with open(cosmic5_out_file, "w") as out_fh:
					out_fh.write(header)
					for line_out in cosmic5_list:
						out_fh.write(delim.join(line_out))
			if len(cosmic10_list) >= 1:
				with open(cosmic10_out_file, "w") as out_fh:
					out_fh.write(header)
					for line_out in cosmic10_list:
						out_fh.write(delim.join(line_out))





##run methods
kegg_mtor_genes = ['AKT1', 'AKT2', 'AKT3', 'BRAF', 'CAB39', 'CAB39L', 'DDIT4', 'EIF4B', 'EIF4E', 'EIF4E1B', 'EIF4E2', 'EIF4EBP1', 'FIGF', 
		'HIF1A', 'IGF1', 'INS', 'MAPK1', 'MAPK3', 'MLST8', 'MTOR', 'PDPK1', 'PGF', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 
		'PIK3R3', 'PIK3R5', 'PRKAA1', 'PRKAA2', 'RHEB', 'RICTOR', 'RPS6', 'RPS6KA1', 'RPS6KA2', 'RPS6KA3', 'RPS6KA6', 'RPS6KB1', 'RPS6KB2', 
		'RPTOR', 'STK11', 'STRADA', 'TSC1', 'TSC2', 'ULK1', 'ULK2', 'ULK3', 'VEGFA', 'VEGFB', 'VEGFC']

pisces_files = glob.glob('*.pisces.xls')
# pisces_files = ['LR18-492.pisces.xls', 'LR17-490_1801044.pisces.xls']
count_file = 'pisces_filtering_summary_0119.xls'

filter_pisces_files(pisces_files, kegg_mtor_genes, count_file)