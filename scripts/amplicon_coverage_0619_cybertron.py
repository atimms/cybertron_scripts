#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'

##program
bedtools = '/home/atimms/programs/bedtools2.28/bin/bedtools'

##methods
def calculate_exome_coverage(samples, bam_suffix, exome_bed, prefix):
	##parameters and header
	coverage_required = [1,5,10,20,50,100,200]
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
	output_file = prefix + '.coverage_data.txt'
	print output_file
	with open(output_file, "w") as outfh:
		outfh.write(header)
		##go through bam files and get coverage histogram
		for sample in samples:
			bam = sample + bam_suffix
			print 'calculating coverge for bam file', bam
			##get temp coverage files
			coverage_histogram = sample + '.hist.temp'
			hist_fh = open(coverage_histogram, 'w')
			# bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=subprocess.PIPE)
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


def calculate_coverage(wd, bed, output_file):
	os.chdir(wd)
	bams = glob.glob('*.bwa.bam')
	header = delim.join(['sample', 'gene', 'coverage', 'percentage_covered_10x']) + '\n'
	with open(output_file, "w") as outfh:
		outfh.write(header)
		for bam in bams:
			sample = bam.split('.')[0] 
			print 'calculating coverge for bam files', bam
			##get temp coverage files
			coverage_histogram = sample + 'coverage.hist.temp'
			with open(coverage_histogram, 'w') as hist_fh:
				bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', bed, '-b', bam, '-hist'], stdout=hist_fh)
				bedtools_cov.wait()
			##make cov dict for each gene
			gene_cov_dict = {}
			with open(coverage_histogram, "r") as cov_temp:
				for line in cov_temp:

					line = line.strip('\n').split(delim)
					##don't use all
					if line[0] != 'all':				
						amplicon = '_'.join(line[:3])
						gene = line[3]
						target_size = float(line[6])
						depth = int(line[4])
						bases_at_depth = int(line[5])
						##get total bps sequenced
						if depth == 0:
							bps = 0
						else:
							bps = depth * bases_at_depth
						##get bases covered >=10
						bases_ge_10 = 0
						if depth >= 10:
							bases_ge_10 = bases_at_depth
						##add info to dict
						if gene in gene_cov_dict:
							gene_cov_dict[gene][2] += bps
							gene_cov_dict[gene][3] += bases_ge_10
							##only add to target size if a new amplicon
							if amplicon not in gene_cov_dict[gene][0]:
								gene_cov_dict[gene][0].append(amplicon)
								gene_cov_dict[gene][1] += target_size

						#start dict, add
						else:
							gene_cov_dict[gene] = [[amplicon], target_size, bps, bases_ge_10]
			for g in gene_cov_dict:
				# print(g, gene_cov_dict[g])
				coverage = gene_cov_dict[g][2] / gene_cov_dict[g][1]
				percentage_covered_10x = gene_cov_dict[g][3] / gene_cov_dict[g][1]
				line_out = delim.join([sample, g, str(coverage), str(percentage_covered_10x)]) + '\n'
				outfh.write(line_out)



##run methods
bedfile = '/home/atimms/ngs_data/references/hg19/TruSeq_CAT_0718_genes.nochr.bed'
# bamfiles = glob.glob('L*.bwa.bam')
# bamfiles = ['LR02-148.bwa.bam', 'LR02-148f.bwa.bam', 'LR02-148m.bwa.bam']
##plate1
working_dir = '/home/atimms/ngs_data/targetted/kim_cbl_plate1_0918'
results_file = 'kim_cbl_plate1_0918.coverage.xls'
# calculate_coverage(working_dir, bedfile, results_file)

##plate2
working_dir = '/home/atimms/ngs_data/targetted/kim_cbl_plate2_1018'
results_file = 'kim_cbl_plate2_1018.coverage.xls'
# calculate_coverage(working_dir, bedfile, results_file)

##plate3
working_dir = '/home/atimms/ngs_data/targetted/kim_cbl_plate3_1218'
results_file = 'kim_cbl_plate3_1218.coverage.xls'
# calculate_coverage(working_dir, bedfile, results_file)

##plate4
working_dir = '/home/atimms/ngs_data/targetted/kim_cbl_plate4_0319'
results_file = 'kim_cbl_plate4_0319.coverage.xls'
# calculate_coverage(working_dir, bedfile, results_file)

##plate5
working_dir = '/home/atimms/ngs_data/targetted/backedup/kim_cbl_plate5_1019'
results_file = 'kim_cbl_plate5_1019.coverage.xls'
calculate_coverage(working_dir, bedfile, results_file)



