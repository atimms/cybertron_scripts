#!/usr/bin/python
import os
import subprocess
import glob
import shutil

##set input variables and parameters
delim = '\t'


def get_list_of_peds_from_filename(file_suffix):
	filenames = glob.glob('*' + file_suffix)
	peds = [p.split('.')[0] for p in filenames]
	return(peds)


def mkdir_and_move_files(peds):
	temp_dir_name = 'temp_dir'
	for ped in peds:
		os.mkdir(temp_dir_name)
		files_to_move = glob.glob('*' + ped + '*')
		print(files_to_move)
		for ftm in files_to_move:
			shutil.move(ftm, temp_dir_name + '/' + ftm)
		shutil.move(temp_dir_name, ped)



















##run methods
work_dir = '/archive/dobyns_w/exome_data/all_exome_files'
# work_dir = '/archive/dobyns_w/exome_data/temp_bams'
os.chdir(work_dir)
std_anal_suffix = '.std_analysis.xls'

##get names of all ped with a std analysis file
peds_to_combine = get_list_of_peds_from_filename(std_anal_suffix)
print(peds_to_combine)

##put into directory with ped name
mkdir_and_move_files(peds_to_combine)