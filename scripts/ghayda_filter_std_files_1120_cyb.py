#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'

##methods
def filter_cadd_maf(in_files, out_prefix, sd_peds):
	trio_outfile = out_prefix + '.trio.xls'
	single_duo_outfile = out_prefix + '.single_duo.xls'
	with open(trio_outfile, "w") as tout_fh, open(single_duo_outfile, "w") as sdout_fh:
		fc = 0
		for in_file in in_files:
			ped = in_file.split('/')[-1].split('.')[0]
			with open(in_file, "r") as in_fh:
				lc = 0
				fc += 1
				for line in in_fh:
					lc += 1
					line = line.strip('\n').split(delim)
					if lc == 1:
						if fc == 1:
							tout_fh.write(delim.join(['ped'] + line[:47] + line[-7:]) + '\n')
							sdout_fh.write(delim.join(['ped'] + line[:47] + line[-7:]) + '\n')
					else:
						anal = line[-2]
						cadd = line[11]
						gnomad_ex_af = line[33]
						annovar = line[13]
						print(ped, anal, cadd, gnomad_ex_af)
						if cadd == 'None':
							cadd = 25
						else:
							cadd = float(cadd)
						if gnomad_ex_af == 'na' or gnomad_ex_af == '.':
							gnomad_ex_af = 0
						else:
							gnomad_ex_af = float(gnomad_ex_af)
						if cadd >= 20 and gnomad_ex_af <= 0.01 and annovar != '.':
							line_out = [ped] + line[:47] + line[-7:]
							if ped in sd_peds:
								sdout_fh.write(delim.join(line_out) + '\n')
							else:
								tout_fh.write(delim.join(line_out) + '\n')

def filter_by_gene(in_files, out_prefix, gene_file, sd_peds):
	candidiate_genes = []
	with open(gene_file, "r") as gf_fh:
		for line in gf_fh:
			gene = line.rstrip().split('\t')[0]
			if gene not in candidiate_genes:
				candidiate_genes.append(gene)
	print(candidiate_genes, len(candidiate_genes))
	trio_outfile = out_prefix + '.trio.xls'
	single_duo_outfile = out_prefix + '.single_duo.xls'
	with open(trio_outfile, "w") as tout_fh, open(single_duo_outfile, "w") as sdout_fh:
		fc = 0
		for in_file in in_files:
			ped = in_file.split('/')[-1].split('.')[0]
			with open(in_file, "r") as in_fh:
				lc = 0
				fc += 1
				for line in in_fh:
					lc += 1
					line = line.strip('\n').split(delim)
					if lc == 1:
						if fc == 1:
							tout_fh.write(delim.join(['ped'] + line[:47] + line[-7:]) + '\n')
							sdout_fh.write(delim.join(['ped'] + line[:47] + line[-7:]) + '\n')
					else:
						gene = line[5]
						# anal = line[-2]
						# cadd = line[11]
						# gnomad_ex_af = line[33]
						# annovar = line[13]
						# print(ped, anal, cadd, gnomad_ex_af)
						if gene in candidiate_genes:
							line_out = [ped] + line[:47] + line[-7:]
							if ped in sd_peds:
								sdout_fh.write(delim.join(line_out) + '\n')
							else:
								tout_fh.write(delim.join(line_out) + '\n')


##run methods
working_dir = '/archive/mirzaa_g/exomes/result_files_at/exome_project_1120/filtered_analysis_files_1120'
os.chdir(working_dir)
rss_std_files =  glob.glob('/archive/mirzaa_g/exomes/result_files_at/exome_project_1120/results_files_1120/files_from_rss/*std_analysis.xls')
new_std_files =  glob.glob('/archive/mirzaa_g/exomes/result_files_at/exome_project_1120/results_files_1120/0920_files/*std_analysis.xls')
std_files = rss_std_files + new_std_files
incomplete_peds = ['LR02-263', 'LR03-214', 'LR03-384', 'LR04-022', 'LR04-185', 'LR04-199', 'LR04-414', 'LR05-162', 'LR06-105', 'LR07-082', 
		'LR08-418', 'LR10-270', 'LR11-019', 'LR11-173', 'LR12-018', 'LR12-204', 'LR12-214', 'LR12-224', 'LR12-346', 'LR12-401', 'LR14-098', 
		'LR15-287', 'LR16-029', 'LR16-406', 'LR16-503', 'LR17-059', 'LR17-260', 'LR17-439', 'LR18-030', 'LR18-416', 'LR18-427', 'LR18-340',
		'LR19-354', 'LR15-215', 'LR16-004', 'LR16-157', 'LR17-176', 'LR05-375', 'LP99-100', 'LR06-018', 'LR09-416', 'LR11-465', 'LP97-141', 
		'LP98-078', 'LR00-016', 'LR00-225', 'LR01-173', 'LR01-242', 'LR01-279', 'LR01-306', 'LR01-314', 'LR02-085', 'LR03-039', 'LR03-304', 
		'LR03-340', 'LR04-222', 'LR04-315', 'LR04-341', 'LR05-054', 'LR05-265', 'LR07-037', 'LR08-056', 'LR10-260', 'LR11-151', 'LR11-429', 
		'LR12-109', 'LR12-475', 'LR13-076', 'LR13-084', 'LR14-075', 'LR15-137', 'LR16-306', 'LR16-441', 'LR16-461', 'LR16-462', 'LR16-463', 
		'LR16-510', 'LR16-511', 'LR16-512', 'LR16-514', 'LR16-517', 'LR17-100', 'LR17-101', 'LR17-102', 'LR17-103', 'LR17-478', 'LR18-339', 
		'LR18-370', 'LR18-529', 'LR19-273', 'LR07-054', 'LR14-255', 'LR16-053', 'LR19-515', '3C-10', '3C-11', '3C-12', '3C-14', '3C-3', 
		'3C-7', 'DWM10', 'DWM13', 'DWM3']


##look at files
# print('rss, 0320:',len(rss_std_files))
# print('new, 0920:',len(new_std_files))


##filter by cadd/maf
cadd_maf_results_prefix = 'exome_data.cadd20_maf0.01.1120'
# filter_cadd_maf(std_files, cadd_maf_results_prefix, incomplete_peds)


##filter by gene
kims_genes = 'genelist_kim_1120.txt'
kim_results_prefix = 'exome_data.kims_genes.1120'
ghayda_genes = 'genelist_ghayda_1120.txt'
ghayda_results_prefix = 'exome_data.ghayda_genes.1120'
centrosome_list_1220 = 'centrosome_list_1220.txt'
centrosome_results_prefix = 'exome_data.centrosome_genes_1220'
# filter_by_gene(std_files, kim_results_prefix, kims_genes, incomplete_peds)
# filter_by_gene(std_files, ghayda_results_prefix, ghayda_genes, incomplete_peds)
filter_by_gene(std_files, centrosome_results_prefix, centrosome_list_1220, incomplete_peds)


