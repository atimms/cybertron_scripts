#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##info
'''
##for daniela data
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 9aa67182-09d6-405c-b7f7-93b0320828c1
module load biobuilds


'''


##parameters
delim = '\t'



##get exome coverage for bam files
def calculate_exome_coverage(samples, bam_dir, bam_suffix, exome_bed, prefix):
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
			bam = bam_dir + sample + bam_suffix
			print 'calculating coverge for bam file', bam
			##get temp coverage files
			coverage_histogram = sample + '.hist.temp'
			hist_fh = open(coverage_histogram, 'w')
			# bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen(['bedtools', 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=subprocess.PIPE)
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


##daniela exomes 1221 
working_dir = '/archive/heike_c/CFM/exomes_1221'
os.chdir(working_dir)
sample_names = ['101-0046-01.642210', '101-0046-11.642211', '101-0046-21.642212', '101-0048-01.642213', '101-0048-11.642214', '101-0048-21.642215', '101-0049-01.642216', '101-0049-11.642217', '101-0049-21.642218', '101-0052-01.642219', '101-0052-11.642220', '101-0052-21.642543', '101-0053-01.642221', '101-0053-11.642222', '101-0053-21.642223', '101-0054-01.642224', '101-0054-11.642225', '101-0054-21.642226', '101-0057-01.642227', '101-0057-11.642228', '101-0057-21.642229', '101-0067-01.642233', '101-0067-11.642234', '101-0067-21.642235', '101-0068-01.642236', '101-0068-11.642237', '101-0068-21.642238', '101-0069-01.642239', '101-0069-11.642240', '101-0069-21.642241', '101-0070-01.642242', '101-0070-02.642243', '101-0070-11.642244', '101-0070-21.642245', '101-0071-01.642246', '101-0071-21.642248', '101-0073-01.642249', '101-0073-11.642250', '101-0073-21.642251', '101-0074-01.642252', '101-0074-11.642253', '101-0074-21.642254', '103-0035-01.642255', '103-0035-11.642256', '103-0035-21.642257', '103-0037-01.642258', '103-0037-11.642259', '103-0037-21.642260', '103-0038-01.642261', '103-0038-11.642262', '103-0038-21.642263', '103-0042-01.642264', '103-0042-11.642265', '103-0042-21.642266', '103-0046-01.642267', '103-0046-11.642268', '103-0046-21.642269', '103-0048-01.642270', '103-0048-11.642271', '103-0048-21.642272', '103-0049-01.642273', '103-0049-11.642274', '103-0049-21.642275', '103-0052-01.642276', '103-0052-11.642277', '103-0052-21.642278', '103-0053-01.642279', '103-0053-11.642280', '103-0053-21.642281', '103-0054-01.642282', '103-0054-11.642283', '103-0054-21.642284', '103-0055-01.642285', '103-0055-11.642286', '103-0055-21.642287', '103-0056-01.642288', '103-0056-11.642289', '103-0056-21.642290', '103-0057-01.642291', '103-0057-11.642292', '103-0057-21.642293', '103-0059-01.642294', '103-0059-11.642295', '103-0059-21.642296', '103-0061-01.642297', '103-0061-11.642298', '103-0061-21.642299', '103-0063-01.642300', '103-0063-11.642301', '103-0063-21.642302', '103-0064-01.642303', '103-0064-11.642304', '103-0064-21.642305', '103-0066-01.642306', '103-0066-11.642307', '103-0066-21.642308', '103-0067-01.642309', '103-0067-11.642310', '103-0067-21.642311', '103-0068-01.642312', '103-0068-11.642313', '103-0068-21.642314', '103-0069-01.642315', '103-0069-11.642316', '103-0069-21.642317', '103-0073-01.642318', '103-0073-11.642319', '103-0073-21.642320', '103903.642323', '104003.718155', '104502.642331', '104503.642332', '106-0035-01.642336', '106-0035-11.642337', '106-0035-21.642338', '106-0040-01.642345', '106-0040-11.642346', '106-0040-21.642347', '106-0050-01.642354', '106-0050-11.642355', '106-0050-21.642356', '107-0023-01.642360', '107-0023-11.642361', '107-0023-21.642362', '107-0025-01.642366', '107-0025-11.642367', '107-0025-21.642368', '107-0026-01.642369', '107-0026-11.642370', '107-0026-21.642371', '107-0027-01.642372', '107-0027-11.642373', '107-0027-21.642374', '107-0028-01.642375', '107-0028-11.642376', '107-0028-21.642377', '107-0040-01.642384', '107-0040-11.642385', '107-0040-21.642386', '107-0043-01.642390', '107-0043-11.642391', '107-0043-21.642392', '107-0044-01.642393', '107-0044-11.642394', '107-0044-21.642395', '309303.718157', '401-0017-01.642402', '401-0017-11.642403', '401-0017-21.642404', '401-0019-11.642406', '401-0019-21.642407', '402202.642409', '402203.718159', 'CA-01-002-01.642411', 'CA-01-002-11.642412', 'CA-01-002-21.642413', 'CA-01-003-01.642414', 'CA-01-003-11.642415', 'CA-01-003-21.642416', 'CA-01-008-01.642417', 'CA-01-008-11.642418', 'CA-01-008-21.642419', 'CA-01-013-01.642420', 'CA-01-013-11.642421', 'CA-01-013-21.642422', 'CA-01-015-01.642423', 'CA-01-015-11.642424', 'CA-01-015-21.642425', 'CA-01-023-01.642426', 'CA-01-023-11.642427', 'CA-01-023-21.642428', 'CA-01-026-01.642429', 'CA-01-026-11.642430', 'CA-01-026-21.642431', 'CA-02-006-01.642432', 'CA-02-006-11.642433', 'CA-02-006-21.642434', 'CA-03-004-01.642438', 'CA-03-004-11.642439', 'CA-03-004-21.642440', 'CA-03-005-21.642443', 'CA-03-009-11.642451', 'CA-03-009-21.642452', 'CA-03-010-01.642453', 'CA-03-010-11.642454', 'CA-03-010-21.642455', 'CA-03-017-01.642456', 'CA-03-017-11.642457', 'CA-03-017-21.642458', 'CA-03-022-01.642459', 'CA-03-022-11.642460', 'CA-03-022-21.642461', 'CA-03-025-21.642464', 'CA-03-026-11.642466', 'CA-03-026-21.642467', 'CA-03-028-11.642469', 'CA-03-028-21.642470', 'CA-03-029-01.642471', 'CA-03-029-11.642472', 'CA-03-032-01.642474', 'CA-03-032-21.642476', 'CA-03-041-21.642485', 'CA-03-044-11.642490', 'CA-03-045-11.642493', 'CA-03-050-01.642498', 'CA-03-050-11.642499', 'CA-03-050-21.642500', 'CA-03-051-01.642501', 'CA-03-052-01.718177', 'CA-03-052-21.642506', 'CA-03-054-11.642511', 'CA-03-054-21.642512', 'CA-03-057-01.642513', 'CA-03-057-21.642515', 'CA-04-001-01.642516', 'CA-04-001-11.642517', 'CA-04-001-21.642518', 'CA-04-003-21.642521', 'CA-04-005-01.718181', 'CA-04-005-11.642523', 'CA-04-005-21.642524', 'CA-04-006-01.642525', 'CA-04-006-11.642526', 'CA-04-006-21.718182', 'CA-04-009-11.642529', 'CA-04-011-01.642531', 'CA-04-011-11.642532', 'CA-04-011-21.642533', 'CA-04-012-01.642534', 'CA-04-012-11.642535', 'CA-04-012-21.642536', 'CA-04-017-01.642540', 'CA-04-017-11.642541', 'CA-04-017-21.642542']
bamfile_dir = '/archive/heike_c/CFM/exomes_1221/bamFiles/'
twist_bed = '/home/atimms/ngs_data/references/hg19/Twist_Exome_RefSeq_targets_hg19_0.no_chr.bed'
bamfile_suffix = '.bam'
project_prefix = 'exomes_1221'
calculate_exome_coverage(sample_names, bamfile_dir, bamfile_suffix, twist_bed, project_prefix)




