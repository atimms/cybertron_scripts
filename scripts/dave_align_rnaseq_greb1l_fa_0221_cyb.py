#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'
star_threads = '18'

##programs
star = '/home/atimms/programs/STAR-2.5.3a/bin/Linux_x86_64_static/STAR'

##ref files and file names
genome_name = 'mm10_greb1l_0221'
star_index_dir = '/home/atimms/ngs_data/references/star/' + genome_name
fa_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name +'/genes.gtf'

##methods
def make_star_index_files_no_gtf(star_genome_dir, genome_fas, threads_to_use):
	star_index = subprocess.Popen([star, '--runMode', 'genomeGenerate', '--genomeDir', star_genome_dir, '--genomeFastaFiles', genome_fas, '--runThreadN', threads_to_use])
	star_index.wait()

def star_align_paired_end(sample_name, fq_files, star_genome_dir):
	r1_fq = fq_files[0]
	r2_fq = fq_files[1]
	print(sample_name, r1_fq, r2_fq, star_genome_dir)
	# star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align.wait()

def call_rnaseq_methods_for_deseq_analysis(work_dir, info_file, genome):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	print('working on project:', project_name)
	##get info on project
	with open(info_file, "U") as infh:
		metadata_dict = {}
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print(conditions)
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				metadata_dict[sample] = condition_by_sample
				sample_list.append(sample)
				fq1_info = line[1].split(',')
				fq2_info = line[2].split(',')
				# print fq1_info, len(fq1_info)
				# print fq2_info, len(fq2_info)
				##do we need to combine fqs
				if len(fq1_info) == 1:
					fq1 = fq1_info[0]
				elif len(fq1_info) >= 1:
					fq1 = sample + '.r1.fq.gz'
					combine_fqs(fq1_info, fq1)
				##single or paired sequencing
				if fq2_info == ['']:
					seq_type = 'single_end'
				else:
					seq_type = 'paired_end'
					##if same fqs for r1 and r2
					if len(fq1_info) == len(fq2_info):
						##combine files
						if len(fq2_info) == 1:
							fq2 = fq2_info[0]
						elif len(fq2_info) > 1:
							fq2 = sample + '.r2.fq.gz'
							combine_fqs(fq2_info, fq2)
					else:
						print('issue with sample %s, there are different number of fq files for read 1 and 2'%sample)
				##run star
				print('seq_type = ', seq_type)
				if seq_type == 'single_end':
					star_align_single_end(sample, fq1, star_index_dir)
				elif seq_type == 'paired_end':
					star_align_paired_end(sample, [fq1, fq2], star_index_dir)


##run methods
##make star index files
# make_star_index_files_no_gtf(star_index_dir, fa_file, star_threads)

##then align reads
working_dir = '/home/atimms/ngs_data/rnaseq/dave_greb1l_es_cell_rnaseq_0121'
call_rnaseq_methods_for_deseq_analysis(working_dir, 'dave_greb1l_esc_0121.txt', 'mm10')


##get count number for each bam using samtools
#samtools idxstats aln.bam | cut -f 1,3
