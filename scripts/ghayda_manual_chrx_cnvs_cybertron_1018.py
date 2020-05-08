#!/usr/bin/env python
import os
import subprocess
import sys
import pybedtools as pbt

working_dir = '/home/atimms/ngs_data/exomes/working/cnv_analysis/ghayda_1018'
os.chdir(working_dir)
delim = '\t'

'''
load modules:
module load local_python/3.6.4
and install pybedtools
conda install -c bioconda pybedtools
'''


def get_coverage_data(bam_list, output_file, bed_file):
	b = pbt.BedTool(bed_file)
	b_multi = b.multi_bam_coverage(bams=bam_list, output='coverage.temp')
	with open('coverage.temp', "r") as cov_temp, open(output_file, 'w') as outfile:
		header_list = ['chr', 'start', 'end', 'exon', '', 'strand']
		for bam in bam_list:
			bam = bam.split('.')[0]
			header_list.append(bam)
		header = delim.join(header_list + ['\n'])
		outfile.write(header)
		for line in cov_temp:
			outfile.write(line)

def filter_coverage_data(infile, outfile):
	with open(infile, "r") as in_fh, open(outfile, 'w') as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc == 1:
				out_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				pro_covs = [int(i) for i in line[6:8]]
				par_covs = [int(i) for i in line[8:]]
				# print(pro_covs, par_covs, line)
				if sum(pro_covs) < 20 and sum(par_covs) > 50:
					out_fh.write(delim.join(line) + '\n')


bams = ['LR15-377a1.bwa_gatk.bam', 'LR15-377a2.bwa_gatk.bam', 'LR15-377m.bwa_gatk.bam', 'LR15-377r.bwa_gatk.bam']
bed = 'hg19_refgene_chrX_nc.bed'
all_exon_coverage = 'LR15-377.all_x_exon_coverage.txt'
filtered_exon_coverage = 'LR15-377.filtered_x_exon_coverage.txt'
##test
# bed = 'temp.bed'
# all_exon_coverage = 'test_cov.txt'

##get coverage for all chrx exons 
# get_coverage_data(bams, all_exon_coverage, bed )
##filter those
filter_coverage_data(all_exon_coverage, filtered_exon_coverage)
