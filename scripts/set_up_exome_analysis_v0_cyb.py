#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parms
delim = '\t'
##pipeline versions 
germline_version = 'dobyns_exome_pipeline_cybertron_v11'
mosaic_version = 'dobyns_exome_mosaic_pipeline_cybertron_v5'



##methods
def make_ped_files(input_dict):
	for ped in input_dict:
		# print ped, input_dict[ped]
		outfile = ped + '.ped'
		header = ['#family_id', 'name', 'paternal_id', 'maternal_id', 'sex', 'phenotype', '\n']
		with open(outfile, "w") as out_fh:
			out_fh.write(delim.join(header))
			for outline in input_dict[ped]:
				out_fh.write(delim.join(outline) + '\n')

def make_analysis_dicts(ped_info_file):
	ped_file_dict, analysis_dict = {}, {}
	with open(ped_info_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				ped_name = line[0]
				sample_name = line[4]
				ped_type = line[1]
				anal_type = line[2]
				input_file_type = line[3]
				seq = line[11]
				##get input file depending on type
				if input_file_type == 'bam' or input_file_type == 'messy_bam':
					input_files = line[5]
				elif input_file_type == 'fastq':
					input_files = [line[5], line[6]]
				else:
					print('input_file_type not recognized:', input_file_type)
				##general ped info
				anal_info = [ped_type, anal_type, input_file_type]
				# print ped_name, line
				##add data to ped file dict
				ped_file_info = [line[0], line[4]] + line[7:11]
				# print ped_file_info
				if ped_name in ped_file_dict:
					ped_file_dict[ped_name].append(ped_file_info)
				else:
					ped_file_dict[ped_name] = [ped_file_info]
				##add data to analysis_dict
				if seq == 'yes':
					if ped_name in analysis_dict:
						analysis_dict[ped_name][1][sample_name] = input_files
					else:
						analysis_dict[ped_name] = [anal_info, {sample_name:input_files}]
	return ped_file_dict, analysis_dict	


def make_analysis_info_file(work_dir, in_dict, outfile):
	with open(outfile, "w") as out_fh:
		for ped in in_dict:
			print(ped, in_dict[ped])
			ped_type = in_dict[ped][0][0]
			anal_type = in_dict[ped][0][1]
			file_type = in_dict[ped][0][2]
			ped_file_info_dict = in_dict[ped][1]
			##fill in 2 lists for the bams to test and the parental bams
			samples = ped_file_info_dict.keys()
			parent_bams, proband_bams = [], []
			for s in samples:
				if s.endswith('f') or s.endswith('m'):
					parent_bams.append(s + '.bwa_gatk.bam')
				else:
					proband_bams.append(s + '.bwa_gatk.bam')
			if anal_type == 'germline' or anal_type == 'both':
				gl_line_out = germline_version + '.call_all_exome_methods_inc_gemini("' + work_dir + '", "' + ped + '", "' + file_type + '", ' + str(ped_file_info_dict) + ', "' + ped_type + '")' + '\n'
				out_fh.write(gl_line_out)
			if anal_type == 'mosaic' or anal_type == 'both':
				mosaic_line_out = mosaic_version + '.run_mosaic_variant_calling("' + work_dir + '", "' + ped + '", ' + str(proband_bams) + ', ' + str(parent_bams) + ', "", "", "pisces_v2")' + '\n'
				out_fh.write(mosaic_line_out)			

def make_files_master(working_dir, analysis_file):
	os.chdir(working_dir)
	project = analysis_file.split('.')[0]
	##read in analysis file and make 2 dicts
	ped_file_dict, analysis_dict = make_analysis_dicts(analysis_file)
	make_ped_files(ped_file_dict)
	# for ped in analysis_dict:
	# 	print(ped, analysis_dict[ped])
	anal_cmds_file = project + '.cmds.txt'
	make_analysis_info_file(working_dir, analysis_dict, anal_cmds_file)


##run methods

##ghayda_broad_exomes_0920
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_broad_exomes_0920'
# input_file = 'ghayda_test_0920.txt'
# input_file = 'ghayda_broad_0920.txt'
input_file = 'ghayda_broad_2_0920.txt'
make_files_master(work_dir, input_file)





