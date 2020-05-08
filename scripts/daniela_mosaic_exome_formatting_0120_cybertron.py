#!/usr/bin/env python
import sys
import subprocess
import os
import glob



'''
##load modules required for analysis
'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/exomes/working/daniela_mosaic_exomes_redo_1219'
os.chdir(working_dir)


def combine_in_both_tissue_only_filter_parent_vaf(samples, in_both_outfile, in_tissue_outfile):
	in_both_dict, in_tissue_dict = {}, {}
	sc = 0
	for sample in samples:
		sc += 1
		var_dict = {}
		tissue_infile = sample + '-T.pisces.xls'
		blood_infile = sample + '-B.pisces.xls'
		with open(blood_infile, "r") as bin_fh:
			line_count = 0
			for line in bin_fh:
				line = line.rstrip().split(delim)
				line_count += 1
				if line_count == 1:
					print(sample, len(line))
					if sc ==1:
						header = line
				else:
					if len(line) == 83:
						max_parents_alt_reads = int(line[82])
					else:
						max_parents_alt_reads == 0
					var = '_'.join(line[:5])
					if max_parents_alt_reads < 2:
						var_dict[var] = line
		with open(tissue_infile, "r") as tin_fh:
			line_count = 0
			for line in tin_fh:
				line = line.rstrip().split(delim)
				line_count += 1
				if line_count == 1:
					print(sample, len(line))
				else:
					if len(line) == 83:
						max_parents_alt_reads = int(line[82])
					else:
						max_parents_alt_reads == 0
					var = '_'.join(line[:5])
					if max_parents_alt_reads < 2:
						if var in var_dict:
							if sample in in_both_dict:
								in_both_dict[sample][var] = [var_dict[var], line]
							else:
								in_both_dict[sample] = {var:[var_dict[var], line]}
						else:
							if sample in in_tissue_dict:
								in_tissue_dict[sample][var] = line
							else:
								in_tissue_dict[sample] = {var:line}
	with open(in_tissue_outfile, "w") as tout_fh:
		tout_fh.write(delim.join(['sample'] + header) + '\n')
		for samp in in_tissue_dict:
			print(samp, len(in_tissue_dict[samp]))
			for var in in_tissue_dict[samp]:
				# print(samp, var)
				tout_fh.write(delim.join([samp] + in_tissue_dict[samp][var]) + '\n')
	with open(in_both_outfile, "w") as bout_fh:
		extra_header = ['blood_q', 'blood_coverage', 'blood_filter', 'blood_format', 'blood_info', 'blood_var_combined', 'blood_sample_alt_reads', 
				'blood_parent1', 'blood_parent2', 'blood_max_parents_alt_reads', 'tissue_q', 'tissue_coverage', 'tissue_filter', 'tissue_format', 
				'tissue_info', 'tissue_var_combined', 'tissue_sample_alt_reads', 'tissue_parent1', 'tissue_parent2', 'tissue_max_parents_alt_reads']
		bout_fh.write(delim.join(['sample'] + header[:73] + extra_header) + '\n')
		for samp in in_both_dict:
			print(samp, len(in_both_dict[samp]))
			for var in in_both_dict[samp]:
				# print(samp, var)
				if len(in_both_dict[samp][var][0]) == 83:
					bout_fh.write(delim.join([samp] + in_both_dict[samp][var][0] + in_both_dict[samp][var][1][73:]) + '\n')
				else:
					bout_fh.write(delim.join([samp] + in_both_dict[samp][var][0] + ['.','.','.'] + in_both_dict[samp][var][1][73:] + ['.','.','.']) + '\n')

def combine_filter_by_genelists(in_files, out_file, genes):
	with open(out_file, "w") as out_fh:
		fc = 0
		for in_file in in_files:
			sample = in_file.split('.')[0]
			fc += 1
			with open(in_file, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line_count += 1
					if line_count == 1:
						if fc == 1:
							out_fh.write('sample' + delim + line)
					else:
						line = line.rstrip().split(delim)
						gene = line[6]
						if gene in genes:
							out_fh.write(delim.join([sample] + line) + '\n')

def filter_trio_by_coverage(samples, outfile_10x, outfile_20x, genes):
	var_dict = {}
	sc = 0
	for sample in samples:
		sc += 1
		tissue_infile = sample + '-T.pisces.xls'
		blood_infile = sample + '-B.pisces.xls'
		with open(blood_infile, "r") as bin_fh:
			line_count = 0
			for line in bin_fh:
				line = line.rstrip().split(delim)
				line_count += 1
				if line_count == 1:
					print(sample, len(line))
					if sc ==1:
						header = line[:73] + []
				else:
					var = '_'.join([sample] + line[:5])
					s_coverage = line[74]
					s_alt = line[79]
					s_aaf = str(int(s_alt) / float(s_coverage))
					p1_alleles = line[80].split(':')[1]
					# print(p1_alleles)
					if p1_alleles == '.' or line[80].split(':')[0] == './.':
						p1_cov, p1_alt, p1_aaf  = '0', '0', '0'
					else:
						p1_cov = str(int(p1_alleles.split(',')[0]) + int(p1_alleles.split(',')[1]))
						p1_alt = p1_alleles.split(',')[1]
						p1_aaf = str(int(p1_alt) / float(p1_cov))
					p2_alleles = line[81].split(':')[1]
					# print(p2_alleles)
					if p2_alleles == '.' or line[81].split(':')[0] == './.':
						p2_cov, p2_alt, p2_aaf  = '0', '0', '0'
					else:					
						p2_cov = str(int(p2_alleles.split(',')[0]) + int(p2_alleles.split(',')[1]))
						p2_alt = p2_alleles.split(',')[1]
						p2_aaf = str(int(p2_alt) / float(p2_cov))
					var_dict[var] = [sample] + line[:73] + [s_coverage, s_alt, s_aaf, p1_cov, p1_alt, p1_aaf, p2_cov, p2_alt, p2_aaf]

		with open(tissue_infile, "r") as tin_fh:
			line_count = 0
			for line in tin_fh:
				line = line.rstrip().split(delim)
				line_count += 1
				if line_count == 1:
					print(sample, len(line))
				else:
					var = '_'.join([sample] + line[:5])
					s2_coverage = line[74]
					s2_alt = line[79]
					s2_aaf = str(int(s2_alt) / float(s2_coverage))
					p1_alleles = line[80].split(':')[1]
					# print(p1_alleles, line[80])
					if p1_alleles == '.' or line[80].split(':')[0] == './.':
						p1_cov, p1_alt, p1_aaf  = '0', '0', '0'
					else:
						p1_cov = str(int(p1_alleles.split(',')[0]) + int(p1_alleles.split(',')[1]))
						p1_alt = p1_alleles.split(',')[1]
						p1_aaf = str(int(p1_alt) / float(p1_cov))
					p2_alleles = line[81].split(':')[1]
					# print(p2_alleles)
					if p2_alleles == '.' or line[81].split(':')[0] == './.':
						p2_cov, p2_alt, p2_aaf  = '0', '0', '0'
					else:					
						p2_cov = str(int(p2_alleles.split(',')[0]) + int(p2_alleles.split(',')[1]))
						p2_alt = p2_alleles.split(',')[1]
						p2_aaf = str(int(p2_alt) / float(p2_cov))
					if var in var_dict:
						var_dict[var].extend([s2_coverage, s2_alt, s2_aaf])
					else:
						var_dict[var] = [sample] + line[:73] + ['.', '.', '.', p1_cov, p1_alt, p1_aaf, p2_cov, p2_alt, p2_aaf, s2_coverage, s2_alt, s2_aaf]
	with open('temp_vars.xls', "w") as temp_outfh:
		extra_header = ['blood_coverage', 'blood_alt', 'blood_aaf', 'parent1_coverage', 'parent1_alt', 'parent1_aaf',
				'parent2_coverage', 'parent2_alt', 'parent2_aaf', 'tissue_coverage', 'tissue_alt', 'tissue_aaf', 'in_genelist']
		temp_outfh.write(delim.join(['sample'] + header + extra_header) + '\n')
		for var in var_dict:
			print(len(var_dict[var]))
			if len(var_dict[var]) == 86:
				temp_outfh.write(delim.join(var_dict[var]) + '\n')
			else:
				temp_outfh.write(delim.join(var_dict[var] + ['.', '.', '.']) + '\n')

	with open(outfile_10x, "w") as ten_outfh, open(outfile_20x, "w") as twen_outfh, open('temp_vars.xls', "r") as temp_infh:
		line_count = 0
		for line in temp_infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				ten_outfh.write(delim.join(line) + '\n')
				twen_outfh.write(delim.join(line) + '\n')
			else:
				gene = line[7]
				rmsk = line[11]
				supdup = line[12]
				p1_alt = int(line[78])
				p2_alt = int(line[81])
				print(len(line))
				coverages = [line[74],line[77],line[80],line[83]]
				coverages2 = ['200' if x=='.' else x for x in coverages]
				coverages3 = [int(x) for x in coverages2]
				print(coverages, coverages2, coverages3)
				if gene in genes:
					gl = 'yes'
				else:
					gl = 'no'
				if rmsk == '.' and supdup == '.' and p1_alt <2 and p2_alt <2:
					if min(coverages3) >= 10:
						ten_outfh.write(delim.join(line + [gl]) + '\n')
					if min(coverages3) >= 20:
						twen_outfh.write(delim.join(line + [gl]) + '\n')





##run methods
project_name = 'daniela_mosaic_exomes_1219'
pisces_files = ['CFM-MOS-01-01-B.pisces.xls', 'CFM-MOS-01-01-T.pisces.xls', 'CFM-MOS-06-01-B.pisces.xls', 'CFM-MOS-06-01-T.pisces.xls', 'CFM-MOS-09-01-B.pisces.xls', 
		'CFM-MOS-09-01-T.pisces.xls', 'CFM-MOS-13-01-B.pisces.xls', 'CFM-MOS-13-01-T.pisces.xls', 'CFM-MOS-14-01-B.pisces.xls', 'CFM-MOS-14-01-T.pisces.xls', 
		'CFM-MOS-15-01-B.pisces.xls', 'CFM-MOS-15-01-T.pisces.xls', 'CFM-MOS-16-01-B.pisces.xls', 'CFM-MOS-16-01-T.pisces.xls', 'CFM-MOS-18-01-B.pisces.xls', 
		'CFM-MOS-18-01-T.pisces.xls']
samples = ['CFM-MOS-01-01', 'CFM-MOS-06-01', 'CFM-MOS-09-01', 'CFM-MOS-13-01', 'CFM-MOS-14-01', 'CFM-MOS-15-01', 'CFM-MOS-16-01']
trio_samples = ['CFM-MOS-01-01', 'CFM-MOS-06-01', 'CFM-MOS-13-01', 'CFM-MOS-14-01', 'CFM-MOS-15-01']
in_tissue_and_blood_outfile = project_name + '.in_tissue_and_blood.xls' 
in_tissue_only_outfile = project_name + '.in_tissue_not_blood.xls' 
in_genelists_outfile = project_name + '.in_genelists.xls'
new_project_name = 'daniela_mosaic_exomes_0220'
trio_10x_file = project_name + '.trios_10x_in_all.xls'
trio_20x_file = project_name + '.trios_20x_in_all.xls'
genelist = ['TWIST2', 'POLR1A', 'PLCB4', 'GNAI3', 'TSR2', 'RPS28', 'RPS26', 'HOXA1', 'TFAP2A', 'EYA', 'SIX1', 'EYA1', 'SIX5', 'CHD7', 'SEMA3E', 'HOXA2', 'FGF3', 
		'FANCL', 'FRAS1', 'FREM2', 'GRIP1', 'KMT2D', 'GDF6', 'FGFR2', 'FGFR3', 'FGF10', 'EFTUD2', 'EDNRA', 'ORC1', 'ORC4', 'ORC6', 'CDT1', 'CDC6', 'GMNN', 'DHODH', 
		'SF3B4', 'HMX1', 'SALL4', 'GLI3', 'SALL1', 'TCOF1', 'POLR1D', 'POLR1C', 'DCHS1', 'FAT4', 'POMT1', 'FGF13', 'COG1', 'GLI2', 'FIG4', 'WNT7A', 'A2M', 'ABR', 
		'AKAP13', 'ARAP1', 'ARAP2', 'ARAP3', 'ARHGAP1', 'ARHGAP10', 'ARHGAP11A', 'ARHGAP11B', 'ARHGAP12', 'ARHGAP15', 'ARHGAP17', 'ARHGAP18', 'ARHGAP19', 'ARHGAP20', 
		'ARHGAP21', 'ARHGAP22', 'ARHGAP23', 'ARHGAP24', 'ARHGAP25', 'ARHGAP26', 'ARHGAP27', 'ARHGAP28', 'ARHGAP29', 'ARHGAP30', 'ARHGAP31', 'ARHGAP32', 'ARHGAP33', 
		'ARHGAP35', 'ARHGAP36', 'ARHGAP39', 'ARHGAP4', 'ARHGAP40', 'ARHGAP42', 'ARHGAP44', 'ARHGAP5', 'ARHGAP6', 'ARHGAP8', 'ARHGAP9', 'ARHGDIA', 'ARHGDIB', 'ARHGDIG', 
		'ARHGEF1', 'ARHGEF10', 'ARHGEF10L', 'ARHGEF11', 'ARHGEF12', 'ARHGEF15', 'ARHGEF16', 'ARHGEF17', 'ARHGEF18', 'ARHGEF19', 'ARHGEF2', 'ARHGEF26', 'ARHGEF3', 'ARHGEF33', 
		'ARHGEF35', 'ARHGEF37', 'ARHGEF38', 'ARHGEF39', 'ARHGEF4', 'ARHGEF40', 'ARHGEF5', 'ARHGEF6', 'ARHGEF7', 'ARHGEF9', 'BCR', 'CDC42', 'CHN1', 'CHN2', 'DEPDC1B', 'DEPDC7', 
		'DLC1', 'ECT2', 'FAM13A', 'FAM13B', 'FGD1', 'FGD2', 'FGD3', 'FGD4', 'GDI1', 'GDI2', 'GDP', 'GMIP', 'GNA13', 'GTP', 'HMHA1', 'ITSN1', 'KALRN', 'MCF2L', 'MYO9A', 'MYO9B', 
		'NET1', 'NGEF', 'OBSCN', 'OCRL', 'OPHN1', 'PIK3R2', 'PLEKHG2', 'PLEKHG5', 'PREX1', 'RAC1', 'RAC2', 'RAC3', 'RACGAP1', 'RALBP1', 'RASGRF2', 'RHOA', 'RHOB', 'RhoBTB', 
		'RHOC', 'RHOD', 'RHOF', 'RHOG', 'RHOH', 'RHOJ', 'RHOQ', 'RHOT1', 'RHOT2', 'RHOU', 'RHOV', 'SOS1', 'SOS2', 'SRGAP1', 'SRGAP2', 'SRGAP3', 'STARD13', 'STARD8', 'SYDE1', 
		'SYDE2', 'TAGAP', 'TIAM1', 'TIAM2', 'TRIO', 'TRIP10', 'VAV1', 'VAV2', 'VAV3', 'ARHGEF28', 'ARHGEF25', 'ARHGEF8', 'ARHGEF29', 'DOCK2', 'DOCK3', 'DOCK180', 'DOCK10', 
		'ALS2', 'ARHGEF34P', 'DNMBP', 'ECT2L', 'FARP1', 'FARP2', 'FGD5', 'FGD6', 'ITSN2', 'MCF2', 'MCF2L2', 'PLEKHG1', 'PLEKHG3', 'PLEKHG4', 'PREX2', 'RASGRF1', 'SPATA13', 
		'ARHGAP16P', 'ARHGAP45', 'SH3BP1', 'RHOBTB1', 'RHOBTB2', 'RND1', 'RND2', 'RND3']


##get filtered version
# combine_in_both_tissue_only_filter_parent_vaf(samples, in_tissue_and_blood_outfile, in_tissue_only_outfile)
# combine_filter_by_genelists(pisces_files, in_genelists_outfile, genelist)

##filter trio sample by x coverage in all, and reformat file
filter_trio_by_coverage(trio_samples, trio_10x_file, trio_20x_file, genelist)


