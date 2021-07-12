#!/usr/bin/env python
import subprocess
import os


'''
info:
qsub -Iq cdbrmq -l mem=20gb,ncpus=1 -P 19833a08-f6fb-4bea-8526-8a79069da878
'''

##parameters
delim = '\t'


##methods
def reformat_db_data(in_file, out_file):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			##reformat some errors
			line = line.replace(';', ',')
			line = line.replace('"', '')
			line = line.replace('2104038. s2104039', '2104038,s2104039')
			line = line.replace('1612059. 1906002', '1612059,1906002')

			##continue
			line = line.rstrip().split(delim)
			lc += 1
			if lc >1:
				##if we have an id 
				if len(line) == 2:
					lr = line[0].strip()
					##get ids and strip whitespaces
					ids = [i.replace(" ", "") for i in line[1].split(',')]
					##remove empty strings
					id_names = [x for x in ids if x]
					print(line, lr, ids, id_names)
					for id_name in id_names:
						##remove anything after a space
						id_name = id_name.split('(')[0]
						##remove/format weird ids
						if id_name.startswith(('xq')):
							id_name2 = id_name[2:]
						elif id_name.endswith(('lb')):
							id_name2 = id_name[1:-2]
						elif id_name.endswith(('pp')):
							id_name2 = id_name[:-2]
						elif id_name.startswith(('sch')):
							id_name2 = 'na'
						elif id_name == '13091m24':
							id_name2 = '1309124'
						elif id_name.startswith(('x', 'q', 's', 'd')):
							id_name2 = id_name[1:]
						elif id_name.startswith(('U', 'c', 'p', 'D', 'f', 'n', 'o')):
							id_name2 = 'na'
						elif id_name.endswith(('t', 'n', 'm', 'd', 'e', 's', 'o', 'd')):
							id_name2 = 'na'
						elif id_name == '0':
							id_name2 = 'na'
						else:
							id_name2 = id_name
						if id_name2 != 'na':
							out_fh.write(delim.join([lr, id_name2]) + '\n')

def make_sample_location_dict(in_file1, in_file2):
	s_dict = {}
	with open(in_file1, "r") as in1_fh, open(in_file2, "r") as in2_fh:
		for line in in1_fh:
			line = line.rstrip().split(delim)
			box_name = line[0]
			samples = line[1:]
			for sample in samples:
				sample = sample.split(' ')[0].split('-')[0]
				if sample.startswith(('x', 'q', 's', 'd')):
					sample = sample[1:]
				if sample in s_dict:
					s_dict[sample].append(box_name)
				else:
					s_dict[sample] = [box_name]
		lc = 0
		for line in in2_fh:
			lc += 1
			if lc >1:
				line = line.rstrip().split(delim)
				box_name = line[1]
				sample = line[0]
				sample = sample.split(' ')[0].split('-')[0]
				if sample.startswith(('x', 'q', 's', 'd')):
					sample = sample[1:]
				if sample in s_dict:
					s_dict[sample].append(box_name)
				else:
					s_dict[sample] = [box_name]
	return(s_dict)


def compare_samples(in_file, sample_dict, out_file):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		header = ['LR number', 'ID', 'box location']
		out_fh.write(delim.join(header) + '\n')
		for line in in_fh:
			line = line.rstrip().split(delim)
			sample = line[1]
			if sample in sample_dict:
				out_fh.write(delim.join(line + [', '.join(sample_dict[sample])]) + '\n')
			else:
				out_fh.write(delim.join(line + ['']) + '\n')


##run methods
working_dir = '/home/atimms/ngs_data/misc/ghayda_sample_xc_0621'
os.chdir(working_dir)

db_data_file = 'db_info.txt'
db_data_file_reformatted = 'db_info_formatted.xls'
boxes_4c_file = '4c_boxes.txt'
new_samples_file = 'new_samples_pull.txt'
final_file = 'sample_xc_062421.xls'


##make new file with tidy names
# reformat_db_data(db_data_file, db_data_file_reformatted)

##make dict from sample names in 4c boxes
box_dict = make_sample_location_dict(boxes_4c_file, new_samples_file)

##make final file
compare_samples(db_data_file_reformatted, box_dict, final_file)

