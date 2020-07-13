#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'
star_threads = '16'

"""
to run:
module load java/1.8.0_121 ##? needed
module load biobuilds/2017.11 ##needed for star/picard, but can't be used for rnaseqc
"""

##programs
star = '/home/atimms/programs/STAR-2.5.3a/bin/Linux_x86_64_static/STAR'
picard = '/home/atimms/programs/picard.jar'
qorts = '/home/atimms/programs/hartleys-QoRTs-39cd1fc/QoRTs.jar'
rnaseqc = '/home/atimms/programs/RNA-SeQC_v1.1.8.jar'

##ref files etc
genome_name = 'hg19'
star_ref_dir = '/home/atimms/ngs_data/references/star/' + genome_name
fa_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name +  '/genes.gtf'
rrna_gtf = 'hg19_rRNA.gtf'

##methods
def star_align_paired_end(sample_name, fq_files, star_genome_dir):
	r1_fq = fq_files[0]
	r2_fq = fq_files[1]
	print sample_name, r1_fq, r2_fq, star_genome_dir
	# star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align.wait()

def star_align_paired_end_not_gz(sample_name, fq_files, star_genome_dir):
	r1_fq = fq_files[0]
	r2_fq = fq_files[1]
	print sample_name, r1_fq, r2_fq, star_genome_dir
	# star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align.wait()

def add_rg_index_bam(sample_name, in_bam_suffix, out_bam_suffix):
	in_bam = sample_name + in_bam_suffix
	out_bam = sample_name + out_bam_suffix
	picard_arg = subprocess.Popen(['picard', 'AddOrReplaceReadGroups', 'INPUT=' + in_bam, 'OUTPUT=' + out_bam, 'SO=coordinate', 'RGID=' + sample_name, 'RGLB=' + sample_name, 'RGPL=Illumina', 'RGPU=machine1', 'RGSM=' + sample_name ])
	picard_arg.wait()
	samtools_index = subprocess.Popen(['samtools', 'index', out_bam])
	samtools_index.wait()
	# picard_rs = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'ReorderSam', 'INPUT=' + in_bam, 'OUTPUT=' + out_bam, 'CREATE_INDEX=true', 'REFERENCE=' + fa_file])
	# picard_rs = subprocess.Popen(['picard', 'SortSam', 'INPUT=' + arg_bam, 'OUTPUT=' + out_bam, 'CREATE_INDEX=true', 'SORT_ORDER=coordinate'])
	# picard_rs.wait()

def run_rnaseqc(sample_file, out_dir):
	# module_rm = subprocess.Popen(['module', 'rm', 'java/1.8.0_121'])
	# module_rm.wait()
	rnaseqc_cmd = subprocess.Popen(['java', '-jar', rnaseqc, '-s', sample_file, '-t', gtf_file, '-r', fa_file, '-o', out_dir, '-rRNA', rrna_gtf])
	rnaseqc_cmd.wait()

def combine_rqc_rpkm_files(sample_names, outfile, results_dir):
	gene_super_dict = {}
	for sample in sample_names:
		ind_gene_dict = {}
		infile = results_dir + sample + '/' + sample + '.metrics.tmp.txt.rpkm.gct'
		line_count = 0
		with open(infile, "r") as in_fh:
			for line in in_fh:
				line_count += 1
				if line_count > 3:
					line = line.rstrip().split(delim)
					gene = line[1]
					rpkm = float(line[2])
					if gene in ind_gene_dict:
						ind_gene_dict[gene].append(rpkm)
					else:
						ind_gene_dict[gene] = [rpkm]
		for g in ind_gene_dict:
			# print g, ind_gene_dict[g]
			ave_rpkm = sum(ind_gene_dict[g]) / len(ind_gene_dict[g])
			if g in gene_super_dict:
				gene_super_dict[g].append(str(ave_rpkm))
			else:
				gene_super_dict[g] = [str(ave_rpkm)]

	print 'dict contains %s genes'%len(gene_super_dict)
	with open(outfile, "w") as out_fh:
		out_fh.write(delim.join(['gene'] + sample_names + ['\n']))
		for g in gene_super_dict:
			if g != "":
				print g, gene_super_dict[g]
				out_fh.write(delim.join([g] + gene_super_dict[g] + ['\n']))

def get_samplenames(infile):
	samples = []
	with open(infile, "r") as infh:
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count > 1:
				sample = line[0]
				samples.append(sample)
	return(samples)


def run_rnaseqc_master(sample_filename, result_dir, rpkm_out_file):
	run_rnaseqc(sample_filename, result_dir)
	samples = get_samplenames(sample_filename)
	combine_rqc_rpkm_files(samples, rpkm_out_file, result_dir)


def analyze_strain_fqs(info_file):
	##get fq files etc
	with open(info_file, "U") as infh:
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count != 1:
				sample = line[0]
				fq1 = line[1]
				fq2 = line[2]
				##run star
				star_align_paired_end(sample, [fq1, fq2], star_ref_dir)
				add_rg_index_bam(sample, '.Aligned.sortedByCoord.out.bam', '.hg19.bam')

def convert_bam_fastq_picard(bamfile, r1_fastq, r2_fastq):
	# picard_sam_fq = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
	picard_sam_fq = subprocess.Popen(['picard', 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
	picard_sam_fq.wait()

def analyze_original_bams(info_file):
	##get fq files etc
	with open(info_file, "U") as infh:
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count != 1:
				sample = line[0]
				bam = line[1]
				fq1 = sample + '.temp.r1_fq'
				fq2 = sample + '.temp.r2_fq'
				##convert bam to fq
				convert_bam_fastq_picard(bam, fq1, fq2)
				##run star then add read group and index
				star_align_paired_end_not_gz(sample, [fq1, fq2], star_ref_dir)
				add_rg_index_bam(sample, '.Aligned.sortedByCoord.out.bam', '.hg19.bam')
##run methods
##setup working directory where results will be
working_dir = '/home/atimms/ngs_data/misc/cunningham_rpkm_calculation_0720'
os.chdir(working_dir)

##remap strain rnaseq to hg19
##just needs file with sample name, fq1 and fq2 (with header)
# strain_fq_info = 'strain_rnaseq_test_fqs.txt'
strain_fq_info = 'strain_rnaseq_fqs.txt'
# analyze_strain_fqs(strain_fq_info)

##make fq and remap original rnaseq to hg19
##just needs file with sample name, bam (no header)
# original_bam_info = 'original_rnaseq_bam_test.txt'
# original_bam_info = 'original_rnaseq_bam_1.txt'
# original_bam_info = 'original_rnaseq_bam_2.txt'
# original_bam_info = 'original_rnaseq_bam_3.txt'
# analyze_original_bams(original_bam_info)

##run rnaseqc on all samples(can't use java 8)
##manually make rna-seqc sample file, 3 column sample name, bam, notes (header)
# sample_info = 'test.txt'
# outdir = 'rnaseqc_results_test/'
# rnaseqc_rpkm_file = 'test.rpkm.txt'
##need to make this !!!!
sample_info = 'rnaseqc_info_all.txt'
outdir = 'rnaseqc_results_all/'
rnaseqc_rpkm_file = 'cunn_all_rnaseq.rpkm.0720.txt'
run_rnaseqc_master(sample_info, outdir, rnaseqc_rpkm_file)


