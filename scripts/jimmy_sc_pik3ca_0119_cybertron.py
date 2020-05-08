#!/usr/bin/env python
import sys
import subprocess
import os
import glob


##note
'''
analysis of f1

load modules:
module load biobuilds/2017.11
module load local_python/3.6.4

'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/misc/jimmy_sc_pik3ca_0119'
os.chdir(working_dir)
bedtools = 'bedtools'

##get exome coverage for bam files
def calculate_coverage(samples, bam_suffix, exome_bed, prefix):
	##parameters and header
	coverage_required = [1,5,10,20,50]
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
	output_file = prefix + '.coverage.xls'
	print(output_file)
	with open(output_file, "w") as outfh:
		outfh.write(header)
		##go through bam files and get coverage histogram
		for sample in samples:
			bam = sample + bam_suffix
			print('calculating coverge for bam file', bam)
			##get temp coverage files
			coverage_histogram = sample + '.hist.temp'
			hist_fh = open(coverage_histogram, 'w')
			# bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist'], stdout=subprocess.PIPE)
			grep_cmd = subprocess.Popen(['grep', '^all'], stdin=bedtools_cov.stdout, stdout=hist_fh)
			grep_cmd.wait()
			hist_fh.close()
			##calculate numbers from histogram
			with open(coverage_histogram, "r") as cov_temp:
				seq_total = 0
				for line in cov_temp:
					line = line.strip('\n').split(delim)
					target_size = float(line[3])
					if line[1] != '0':
						seq_total += int(line[1]) *  int(line[2])
			cov_stats = []
			for cr in coverage_required:
				coverage_bp = 0
				with open(coverage_histogram, "r") as cov_temp:
					for line in cov_temp:
						line = line.strip('\n').split(delim)
						target_size = float(line[3])
						if int(line[1]) >= cr:
							coverage_bp += int(line[2])
					cov_stats.extend([str(coverage_bp), str((coverage_bp / target_size) * 100)])
			line_out = delim.join([bam, str(target_size), str(seq_total), str(seq_total / target_size)] + cov_stats + ['\n'])
			outfh.write(line_out)






##run methods
samples = ['BC_0076', 'BC_0077', 'BC_0078', 'BC_0079', 'BC_0080', 'BC_0081', 'BC_0082', 'BC_0083', 'BC_0246', 'BC_0247', 'BC_0248', 
		'BC_0249', 'BC_0250', 'BC_0251', 'BC_0252', 'BC_0253', 'BC_0254', 'BC_0255', 'BC_0256', 'BC_0257']
# samples = ['BC_0076', 'BC_0077', 'BC_0078']
bam_file_suffix = '.sorted.bam'
pik3ca_bed = 'hg19_pik3ca.bed'
pik3ca_hot_spot_bed = 'hg19_pik3ca_hs.bed'
whole_gene = 'pik3ca_gene_cov_0119'
hs_list = ['E542K', 'E545K', 'H1047']
##calculate coverage over whole gene
# calculate_coverage(samples, bam_file_suffix, pik3ca_bed, whole_gene)

##calulate for individual hot spots, make bed file for each hot spot
for hs in hs_list:
	calculate_coverage(samples, bam_file_suffix, hs + '.bed', 'pik3ca_' + hs + '_cov_0119')







