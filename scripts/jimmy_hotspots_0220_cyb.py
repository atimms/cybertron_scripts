#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'


def make_tables_all_hotspots(infiles, outfile):
	hs_dict = {}
	samples = [i.split('.')[0].split('_',2)[2] for i in infiles]
	print(samples)
	for infile in infiles:
		lc = 0
		sample_name = infile.split('.')[0].split('_',2)[2]
		print(sample_name)
		with open(infile, 'r') as in_fh:
			for line in in_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc == 1:
					header = line[:3]
				else:
					chr_pos = '_'.join(line[:2])
					info = line[:3]
					allele_info = '_'.join(line[10:13])
					if chr_pos in hs_dict:
						hs_dict[chr_pos][1][samples.index(sample_name)] = allele_info
					else:
						hs_dict[chr_pos] = [info, ['na_na_na'] * len(samples)]
						hs_dict[chr_pos][1][samples.index(sample_name)] = allele_info
	with open(outfile, 'w') as out_fh:
		extra_header = []
		for s in samples:
			extra_header.extend([s + ' RefCount', s + ' AltCount', s + ' AltAlleleFraction'])
		# print(header, extra_header)
		out_fh.write(delim.join(header + extra_header) + '\n')
		for cp in hs_dict:
			# print(hs_dict[cp])
			hs_info = hs_dict[cp][0]
			counts = []
			for c in hs_dict[cp][1]:
				cc = c.split('_')
				counts.extend(cc)
			# print(counts)
			out_fh.write(delim.join(hs_info + counts) + '\n')

def make_tables_all_hotspots_take2(infiles, outfile):
	hs_dict = {}
	samples = [i.split('.')[0].split('_',2)[2] for i in infiles]
	print(samples)
	for infile in infiles:
		lc = 0
		sample_name = infile.split('.')[0].split('_',2)[2]
		print(sample_name)
		with open(infile, 'r') as in_fh:
			for line in in_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc == 1:
					header = line[:3]
				else:
					chr_pos = '_'.join(line[:2])
					info = line[:3]
					allele_info = line[12]
					if chr_pos in hs_dict:
						hs_dict[chr_pos][1][samples.index(sample_name)] = allele_info
					else:
						hs_dict[chr_pos] = [info, ['na'] * len(samples)]
						hs_dict[chr_pos][1][samples.index(sample_name)] = allele_info
	with open(outfile, 'w') as out_fh:
		extra_header = []
		for s in samples:
			extra_header.append(s)
		# print(header, extra_header)
		out_fh.write(delim.join(header + extra_header) + '\n')
		for cp in hs_dict:
			# print(hs_dict[cp])
			hs_info = hs_dict[cp][0]
			counts = hs_dict[cp][1]
			# print(counts)
			out_fh.write(delim.join(hs_info + counts) + '\n')








##run methods
working_dir = '/home/atimms/ngs_data/misc/jimmy_hotspots_0220'
os.chdir(working_dir)
hotspot_files = glob.glob('Hotspot_report_*')


##make table of all hotspots
hs_table_outfile = 'all_hotspots_0220.table1.xls'
make_tables_all_hotspots(hotspot_files, hs_table_outfile)
hs_table_outfile_aaf = 'all_hotspots_0220.table_aaf.xls'
make_tables_all_hotspots_take2(hotspot_files, hs_table_outfile_aaf)