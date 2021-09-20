#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
requires:
qsub -Iq cdbrmq -l mem=120gb,ncpus=10 -P 19833a08-f6fb-4bea-8526-8a79069da878
#java 1.8 for gatk
module load java/1.8.0_202 


'''

##set input variables and parameters
delim = '\t'
thread_number = '10'

##programs
# star = '/home/atimms/programs/STAR-2.5.3a/bin/Linux_x86_64_static/STAR'
star = '/home/atimms/programs/STAR-2.7.9a/bin/Linux_x86_64_static/STAR'
picard = '/home/atimms/programs/picard_2.19/picard.jar'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
samtools = '/home/atimms/programs/samtools-1.11/bin/samtools'


##ref files and file names
genome_name = 'GRCh38'
star_index_dir = '/home/atimms/ngs_data/references/star21/' + genome_name
fa_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name +'/genes.gtf'
star_bam_suffix = '.Aligned.sortedByCoord.out.bam'
proccessed_bam_suffix = '.star.gatk.bam'
bamlist = 'bams.list'
ref_dir = '/home/atimms/ngs_data/references/hg38_gatk/'
# fasta = ref_dir + 'Homo_sapiens_assembly38.fasta'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
indels_1000g = ref_dir + '1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf'
dbsnp = ref_dir + 'Homo_sapiens_assembly38.dbsnp138.vcf'


##methods
def make_star_index_files(star_genome_dir, genome_fas, genome_gtf):
	star_index = subprocess.Popen([star, '--runMode', 'genomeGenerate', '--genomeDir', star_genome_dir, '--genomeFastaFiles', genome_fas, '--sjdbGTFfile', genome_gtf, '--runThreadN', thread_number])
	star_index.wait()

def star_align_paired_end_2_pass(sample_name, r1_fq, r2_fq):
	print(sample_name, r1_fq, r2_fq)
	# star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--outSAMmapqUnique', '60', '--readFilesCommand', 'zcat', '--twopassMode', 'Basic', '--outFileNamePrefix', sample_name + '.', '--runThreadN', threads_to_use])
	star_align = subprocess.Popen([star,'--genomeDir', star_index_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMmapqUnique', '60', '--readFilesCommand', 'zcat', '--twopassMode', 'Basic', '--outFileNamePrefix', sample_name + '.', '--runThreadN', thread_number])
	star_align.wait()
	bam_file = sample_name + star_bam_suffix
	samtools_index = subprocess.Popen([samtools, 'index', bam_file])
	samtools_index.wait()


def proccess_bam_files(sample_name):
	in_bam = sample_name + star_bam_suffix
	arg_bam = sample_name + '.with_rg.temp.bam'
	mkdup_bam = sample_name + '.mkdup.temp.bam'
	reorder_bam = sample_name + '.reorder.temp.bam'
	split_bam = sample_name + '.split.temp.bam'
	final_bam = sample_name + proccessed_bam_suffix
	picard_arg = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'AddOrReplaceReadGroups', 'INPUT=' + in_bam, 'OUTPUT=' + arg_bam, 'SO=coordinate', 'RGID=' + sample_name, 'RGLB=' + sample_name, 'RGPL=Illumina', 'RGPU=machine1', 'RGSM=' + sample_name ])
	picard_arg.wait()	
	picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'INPUT=' + arg_bam, 'OUTPUT=' + mkdup_bam,'CREATE_INDEX=true', 'METRICS_FILE=' + sample_name + '.metrics'])
	picard_md.wait()
	picard_rs = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'ReorderSam', 'INPUT=' + mkdup_bam, 'OUTPUT=' + reorder_bam, 'CREATE_INDEX=true', 'REFERENCE=' + fa_file])
	picard_rs.wait()
	gatk_snt = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'SplitNCigarReads',  '-R', fa_file, '-I', reorder_bam, '-o', split_bam, '-U', 'ALLOW_N_CIGAR_READS'])
	gatk_snt.wait()
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '6', '-R', fa_file, '-I', split_bam, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample_name + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '6', '-R', fa_file, '-I', split_bam, '-BQSR', sample_name + '.recal_data.table', '-o', final_bam])
	gatk_pr.wait()
	##remove intermeidate files
	# files_to_go = glob.glob('*with_rg.bam') + glob.glob('*mkdup.bam') + glob.glob('*mkdup.bai') + glob.glob('*metrics') + glob.glob('*reorder.bam')+ glob.glob('*reorder.bai') + glob.glob('*split.bam') + glob.glob('*split.bai')
	# print('removing files:')
	# for f in files_to_go:
	# 	os.remove(f)
	# 	print(f)

def make_sample_fq_dict(in_file):
	fq_dict = {}
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc > 1:
				sample = line[0]
				fq1 = line[1]
				fq2 = line[2]
				fq_dict[sample] = [fq1, fq2]
	return(fq_dict)

def genotyping_pipeline(info_file):
	sample_fq_dict = make_sample_fq_dict(info_file)
	project = info_file.split('_')[0]
	for sample in sample_fq_dict:
		print(sample, sample_fq_dict[sample])
		fq1 = sample_fq_dict[sample][0]
		fq2 = sample_fq_dict[sample][1]
		##align with star
		star_align_paired_end_2_pass(sample, fq1, fq2)
		##process files 
		proccess_bam_files(sample)
	##call variants
	# make_list_of_bams(sample_fq_dict, sample_fq_dict, bamlist)
	# get_vars_around_rs3184504(bamlist, project)

##run methods
working_dir = '/home/atimms/ngs_data/rnaseq/kim_rnaseq_varcalling_0921'
os.chdir(working_dir)


##make new reference for latest star version
# make_star_index_files(star_index_dir, fa_file, gtf_file)

##analyze the files
sample_fq_file = 'kim_rnaseq_varcalling_0921.txt'
# sample_fq_file = 'kim_rnaseq_varcalling_0921_temp.txt' ##testing 2 sample
# sample_fq_file = 'kim_rnaseq_varcalling_0921a.txt' ##split into two for mapping
sample_fq_file = 'kim_rnaseq_varcalling_0921b.txt' ##split into two for mapping
genotyping_pipeline(sample_fq_file)

