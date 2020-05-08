#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import daa_run_analysis_v01


'''
##load modules required for analysis
module load java/1.8.0_121 
module load biobuilds/2016.11
module load vt/0.5772-60f436c3 
module load dotnet/1.0.5 
'''



##parameters
delim = '\t'

def make_analysis_dicts(input_file):
	project = input_file.split('.')[0]
	ped_file_dict, analysis_dict = {}, {}
	with open(input_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				ped_name = line[0]
				# print ped_name, line
				##add data to ped file dict
				ped_file_info = [line[0], line[4]] + line[7:11]
				# print ped_file_info
				if ped_name in ped_file_dict:
					ped_file_dict[ped_name].append(ped_file_info)
				else:
					ped_file_dict[ped_name] = [ped_file_info]
				##add data to script command dict
				anal_info = line[1:7] + [line[11]]
				# print anal_info
				ped_input_type = line[3]
				if ped_input_type == 'fastq':
					if ped_name in analysis_dict:
						analysis_dict[ped_name][1][anal_info[3]] = anal_info[4:6]
					else:
						analysis_dict[ped_name] = [anal_info[:3] + [anal_info[-1]], {anal_info[3]: anal_info[4:6]}] 
				elif ped_input_type == 'bam':
					if ped_name in analysis_dict:
						analysis_dict[ped_name][1][anal_info[3]] = anal_info[4]
					else:
						analysis_dict[ped_name] = [anal_info[:3] + [anal_info[-1]], {anal_info[3]: anal_info[4]}] 	
				else:
					print 'ped_input_type %s not recognized'%ped_input_type	
	return ped_file_dict, analysis_dict		

def make_ped_files(input_dict):
	for ped in input_dict:
		# print ped, input_dict[ped]
		outfile = ped + '.ped'
		header = ['#family_id', 'name', 'paternal_id', 'maternal_id', 'sex', 'phenotype', '\n']
		with open(outfile, "w") as out_fh:
			out_fh.write(delim.join(header))
			for outline in input_dict[ped]:
				out_fh.write(delim.join(outline) + '\n')




##master method
def run_amplicon_analysis(working_dir, analysis_file):
	os.chdir(working_dir)
	project = analysis_file.split('.')[0]
	##read in analysis file and make 2 dicts
	ped_file_dict, analysis_dict = make_analysis_dicts(analysis_file)
	##make ped files from first
	make_ped_files(ped_file_dict)
	##run analysis from 2nd dict
	daa_run_analysis_v01.run_analysis(working_dir, analysis_dict, project)
	
def just_make_peds(working_dir, analysis_file):
	os.chdir(working_dir)
	project = analysis_file.split('.')[0]
	print(project)
	##read in analysis file and make 2 dicts
	ped_file_dict, analysis_dict = make_analysis_dicts(analysis_file)
	##make ped files from first
	make_ped_files(ped_file_dict)

def just_run_gemini(working_dir, analysis_file):
	os.chdir(working_dir)
	project = analysis_file.split('.')[0]
	##read in analysis file and make 2 dicts
	ped_file_dict, analysis_dict = make_analysis_dicts(analysis_file)
	##make ped files from first
	make_ped_files(ped_file_dict)
	##run analysis from 2nd dict
	daa_run_analysis_v01.run_gemini(working_dir, analysis_dict, project)

##run methods
##comes from pbs script
# work_dir = sys.argv[1]
# input_file = sys.argv[2]


##analaysis for zach 0718 -- test from MiSeq data
##for testing
# work_dir = '/home/atimms/ngs_data/targetted/zach_miseq_test_0718'
# input_file = 'zach_miseq_test_0718.txt'
# input_file = 'test.txt'

##plate1
work_dir = '/home/atimms/ngs_data/targetted/kim_cbl_plate1_0918'
input_file = 'cbl_plate1.txt'
# run_amplicon_analysis(work_dir, input_file)
##repeat using intersected vars (0619)
# just_run_gemini(work_dir, input_file)

##plate2
work_dir = '/home/atimms/ngs_data/targetted/kim_cbl_plate2_1018'
input_file = 'cbl_plate2.txt'
# # input_file = 'test1.txt'
# run_amplicon_analysis(work_dir, input_file)
##repeat using intersected vars (0619)
# just_run_gemini(work_dir, input_file)

##plate3
work_dir = '/home/atimms/ngs_data/targetted/kim_cbl_plate3_1218'
# work_dir = '/home/atimms/ngs_data/targetted/kim_cbl_plate3_1218_temp'
# input_file = 'cbl_plate3.txt'
# input_file = 'cbl_plate3_2.txt'
input_file = 'cbl_plate3_3.txt'
##just make ped files, so can modify
# just_make_peds(work_dir, input_file)
##then run the analysis
# run_amplicon_analysis(work_dir, input_file)
##repeat using intersected vars (0619)
# just_run_gemini(work_dir, input_file)

##plate4
work_dir = '/home/atimms/ngs_data/targetted/kim_cbl_plate4_0319'
input_file = 'cbl_plate4.txt'
# run_amplicon_analysis(work_dir, input_file)
##repeat using intersected vars (0619)
# just_run_gemini(work_dir, input_file)


##plate5
work_dir = '/home/atimms/ngs_data/targetted/kim_cbl_plate5_1019'
input_file = 'cbl_plate5.txt'
run_amplicon_analysis(work_dir, input_file)

