#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'

##methods
def extract_hotspot_data(infiles, regions_file, outfile):
	region_list = []
	with open(regions_file, 'r') as reg_fh:
		for line in reg_fh:
			line = line.split(delim)
			c_pos = '_'.join(line[:2])
			region_list.append(c_pos)
	print(region_list)
	fc = 0
	with open(outfile, 'w') as out_fh:
		for infile in infiles:
			fc += 1
			lc = 0
			with open(infile, 'r') as in_fh:
				for line in in_fh:
					line = line.rstrip().split(delim)
					lc += 1
					if lc == 1:
						if fc == 1:
							out_fh.write(delim.join(line) + '\n')
					else:
						chr_pos = '_'.join(line[:2])
						if chr_pos in region_list:
							out_fh.write(delim.join(line) + '\n')

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


def count_all_alleles_per_hotspot(infiles, outfile):
	hs_dict = {}
	for infile in infiles:
		lc = 0
		sample_name = infile.split('.')[0].split('_',2)[2]
		print(sample_name)
		with open(infile, 'r') as in_fh:
			for line in in_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc == 1:
					header = line[:3] + [line[5]]
				else:
					chr_pos = '_'.join(line[:2])
					info = line[:3] + [line[5]]
					allele_counts = [int(i) for i in line[6:10]]
					if chr_pos in hs_dict:
						if sample_name not in hs_dict[chr_pos][1]:
							hs_dict[chr_pos][1].append(sample_name)
							hs_dict[chr_pos][2][0].append(allele_counts[0])
							hs_dict[chr_pos][2][1].append(allele_counts[1])
							hs_dict[chr_pos][2][2].append(allele_counts[2])
							hs_dict[chr_pos][2][3].append(allele_counts[3])
					else:
						hs_dict[chr_pos] = [info, [sample_name], [[],[],[],[]]]
						hs_dict[chr_pos][2][0].append(allele_counts[0])
						hs_dict[chr_pos][2][1].append(allele_counts[1])
						hs_dict[chr_pos][2][2].append(allele_counts[2])
						hs_dict[chr_pos][2][3].append(allele_counts[3])
	with open(outfile, 'w') as out_fh:
		extra_header = ['total samples','A count', 'C count', 'T count', 'G count', 'A samples>4', 'C samples>4', 'T samples>4', 'G samples>4']
		# print(header, extra_header)
		out_fh.write(delim.join(header + extra_header) + '\n')
		for cp in hs_dict:
			# print(hs_dict[cp])
			hs_info = hs_dict[cp][0]
			counts, samples_more_zero = [], []
			for c in hs_dict[cp][2]:
				sumc = sum(c)
				more_than_zero = len([i for i in c if i > 4])
				counts.append(str(sumc))
				samples_more_zero.append(str(more_than_zero))
				print(cp, c, more_than_zero)
			total_samples = str(len(c))
			# print(counts)
			out_fh.write(delim.join(hs_info + [total_samples] + counts + samples_more_zero) + '\n')



##run methods

working_dir = '/home/atimms/ngs_data/misc/jimmy_hotspot_format_1120'
os.chdir(working_dir)
hotspot_files = glob.glob('Hotspot_report_*')


##just extract hotspots
hotspots_wanted = 'hotspots_wanted_112320.txt'
hotspot_outfile = 'hotspots_112320.combined.xls'
# extract_hotspot_data(hotspot_files, hotspots_wanted, hotspot_outfile)

##make table of all hotspots
hs_table_outfile = 'all_hotspots.table1.xls'
# make_tables_all_hotspots(hotspot_files, hs_table_outfile)
hs_table_outfile_aaf = 'all_hotspots.table_aaf.xls'
make_tables_all_hotspots_take2(hotspot_files, hs_table_outfile_aaf)

##counts all alleles per hotspot
hs_allele_counts_outfile = 'all_hotspots.allele_counts.xls'
# count_all_alleles_per_hotspot(hotspot_files, hs_allele_counts_outfile)




