#!/usr/bin/env python
import sys
import subprocess
import os.path

##notes
'''
move to folder where data is i.e. rawlings_d/Rawlings_lab/Stem_Cell_Therapies/Globin_bbb/Sickle/20181016_MiSeq/181008_M04866_0180_000000000-C4N59_SCL
'''

##parameters to add:
##the prefix of any summary files
# project_name = 'HBB'
project_name = 'HBD'
##filenames... no need to change
delim = '\t'
results_file_location = 'Analysis/'
mapping_stats_file = 'Mapping_statistics.txt'
combined_mapping_stat = project_name + '.mapping_statistics.xls'
quant_file = 'Quantification_of_editing_frequency.txt'
combined_quant_file = project_name + '.quantification_of_editing_frequency.xls'
indel_hist_file = 'indel_histogram.txt'
combined_indel_hist = project_name + '.indel_histogram.xls'
framshift_file = 'Frameshift_analysis.txt'
combined_framshift = project_name + '.frameshift_analysis.xls'

##methods
def get_sample_names_from_analysis_files(file_location):
	sample_name_dict = {}
	results_dirs = os.listdir(file_location)
	for r in results_dirs:
		sample_name = r.split('_')[2]
		s2 = r.split('_')[2] + '_' + r.split('_')[3]
		if sample_name in sample_name_dict:
			sample_name_dict[s2] = file_location + r
		else:
			sample_name_dict[sample_name] = file_location + r
	return(sample_name_dict)

def combine_mapping_stats(sample_dir_dict, infile_name, outfile):
	out_dict = {}
	for sample in sample_dir_dict:
		infile = sample_dir_dict[sample] + '/' + infile_name
		if os.path.exists(infile):
			with open(infile, "r") as infh:
				for line in infh:
					info = line.rstrip().split(':')[1]
					if sample in out_dict:
						out_dict[sample].append(info)
					else:
						out_dict[sample] = [info]
		else:
			out_dict[sample] = ['-', '-', '-']
	with open(outfile, "w") as outfh:
		outfh.write(delim.join(['sample', 'reads in inputs', 'reads after preprocessing', 'reads aligned']) + '\n')
		for s in out_dict:
			outfh.write(delim.join([s] + out_dict[s]) + '\n')

def combine_quant_stats(sample_dir_dict, infile_name, outfile):
	out_dict = {}
	for sample in sample_dir_dict:
		infile = sample_dir_dict[sample] + '/' + infile_name
		if os.path.exists(infile):
			with open(infile, "r") as infh:
				lc = 0
				info_out = []
				for line in infh:
					lc += 1
					if lc < 6:
						info = line.rstrip().split(':')[1]
						if lc == 2:
							info_out.append(info.split(' ')[0])
							print info_out
						if lc == 3 or lc == 4 or lc == 5:
							info_out.append(info.split(' ')[0])
							info_out.append(info.split(' ')[2].split('(')[1])
							info_out.append(info.split(' ')[6])
							info_out.append(info.split(' ')[10])
							print info_out
			##write to dict
			out_dict[sample] = info_out
		else:
			out_dict[sample] = ['-', '-', '-', '-','-', '-', '-', '-','-', '-', '-', '-', '-']
	with open(outfile, "w") as outfh:
		outfh.write(delim.join(['sample', 'Unmodified reads', 'NHEJ total', 'NHEJ with insertions', 'NHEJ with deletions', 
			'NHEJ with substitutions', 'HDR total', 'HDR with insertions', 'HDR with deletions', 'HDR with substitutions', 
			'Mixed total', 'Mixed with insertions', 'Mixed with deletions', 'Mixed with substitutions']) + '\n')
		for s in out_dict:
			outfh.write(delim.join([s] + out_dict[s]) + '\n')

def combine_indel_hist(sample_dir_dict, infile_name, outfile):
	out_dict = {}
	header_list  = []
	for sample in sample_dir_dict:
		infile = sample_dir_dict[sample] + '/' + infile_name
		if os.path.exists(infile):
			header_list.append(sample)
			with open(infile, "r") as infh:

				lc = 0
				for line in infh:
					lc += 1
					if lc > 1:
						line = line.rstrip().split(delim)
						i_size = line[0]
						count = line[1]
						if i_size in out_dict:
							out_dict[i_size].append(count)
						else:
							out_dict[i_size] = [count]
	with open(outfile, "w") as outfh:
		outfh.write(delim.join(['indel size'] + header_list) + '\n')
		for isize in out_dict:
			outfh.write(delim.join([isize] + out_dict[isize]) + '\n')

def combine_framshift_files(sample_dir_dict, infile_name, outfile):
	out_dict = {}
	for sample in sample_dir_dict:
		infile = sample_dir_dict[sample] + '/' + infile_name
		if os.path.exists(infile):
			with open(infile, "r") as infh:
				lc = 0
				for line in infh:
					lc += 1
					if lc >1:
						info = line.rstrip().split(':')[1].split(' ')[0]
						if sample in out_dict:
							out_dict[sample].append(info)
						else:
							out_dict[sample] = [info]
		else:
			out_dict[sample] = ['-', '-', '-']
	with open(outfile, "w") as outfh:
		outfh.write(delim.join(['sample', 'Noncoding mutation', 'In-frame mutation', 'Frameshift mutation']) + '\n')
		for s in out_dict:
			outfh.write(delim.join([s] + out_dict[s]) + '\n')

##run methods
##get sample name and folder name
sample_dict = get_sample_names_from_analysis_files(results_file_location)
##check dict
# for s in sample_dict:
# 	print(s, sample_dict[s])

##get mapping stats files
combine_mapping_stats(sample_dict, mapping_stats_file, combined_mapping_stat)
combine_quant_stats(sample_dict, quant_file, combined_quant_file)
combine_indel_hist(sample_dict, indel_hist_file, combined_indel_hist)
combine_framshift_files(sample_dict, framshift_file, combined_framshift)


