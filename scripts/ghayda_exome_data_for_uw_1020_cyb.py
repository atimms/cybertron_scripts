#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil
import re

##parameters
delim = '\t'


##methods
def filter_solved_for_mygene2(in_file, out_file):
	##get all denovo vars in one file
	with open(in_file, "U") as in_fh, open(out_file, "w") as out_fh:
		##add header
		header = ['Ped', 'Gene', 'Ped Type', 'DxGroup1', 'Variant Info', 'Inheritance',	'Genomic Coordinates', 'Interpretation']
		out_fh.write(delim.join(header) + '\n')
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				# print(line)
				# line = line.replace("\r","")
				# line = line.replace("\n","")
				line = line.rstrip('\n').rstrip('\r').split(delim)
				ped_id = line[0]
				# print(ped_id, len(line))
				solved_status = line[35].split(' ')[0]
				solved_gene = line[36].rstrip()
				ped_type = line[12].rstrip()
				dx = line[23]
				# print(ped_id, solved_status, solved_gene, ped_type)
				variant_info = line[37:40]
				if solved_status == 'candidate' or solved_status == 'solved':
					line_out = [ped_id, solved_gene, ped_type, dx] + variant_info + [solved_status]
					out_fh.write(delim.join(line_out) + '\n')

def filter_missing_vars_project(in_file, out_file, ad_dir):
	##get list of all unsolved trios
	unsolved_peds = []
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip('\n').rstrip('\r').split(delim)
				ped_id = line[0]
				# print(ped_id, len(line))
				solved_status = line[35].split(' ')[0]
				solved_gene = line[36].rstrip()
				ped_type = line[12].rstrip()
				dx = line[23]
				# print(ped_id, solved_status, solved_gene, ped_type)
				variant_info = line[37:40]
				if solved_status != 'candidate' and solved_status != 'solved':
					if ped_type == 'Trio' or ped_type == 'Trio*' or ped_type == 'Quad':
						unsolved_peds.append(ped_id)
	print(len(unsolved_peds))
	with open(out_file, "w") as out_fh:
		fc = 0
		for ped in unsolved_peds:
			ad_file = ad_dir + ped + '.auto_dominant_analysis.xls'
			fc += 1
			with open(ad_file, "r") as ad_fh:
				lc2 = 0
				for line in ad_fh:
					lc2 += 1
					line = line.strip('\n').split(delim)
					if lc2 == 1:
						if fc ==1:
							out_fh.write(delim.join(['pedigree'] + line) + '\n')
					else:
						# clinsig = line[18].split('/')[0].split(',')[0]
						clinsig = re.split("/|,", line[18])[0]
						print(line[18], clinsig)

						gnomad_ex_af = line[33]
						if gnomad_ex_af == 'na' or gnomad_ex_af == '.':
							gnomad_ex_af = '0'
						if clinsig == 'Pathogenic' or clinsig == 'Likely_pathogenic': 
							if float(gnomad_ex_af) <= 0.01:
								out_fh.write(delim.join([ped] + line) + '\n')




##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/exome_project_20'
os.chdir(working_dir)
ped_info_file = 'master_ped_tracking_101520.txt'
mygene2_file = 'mirzaa_mygene2_1020.xls'
missing_vars_file = 'mirzaa_missing_vars_1020.xls'
auto_dom_dir = '/home/atimms/ngs_data/exomes/working/exome_project_20/src_files/auto_dom/'


##step 0. get the data
##step1. manually format and get files
## copy src_files.zip to cybertron
## use clean() function to copy data from master_ped_tracking_10320 to new excel file and then copy to text file


##step1. filter data for 'solved' exomes in format for mygene2
##want: Ped, Gene, Nucleotide change, Amino acid change, Transcript, Zygosity, Inheritance, Interpretation
# filter_solved_for_mygene2(ped_info_file, mygene2_file)

##step2. get data for missing variant project from auto dom files
filter_missing_vars_project(ped_info_file, missing_vars_file, auto_dom_dir)






