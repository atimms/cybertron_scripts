#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##info
'''
##convert_bam_fastq_picard
module load biobuilds

##xengsort - already made index
qsub -Iq longq -l mem=200gb,ncpus=20
module load local_python/3.7.6
source activate xengsort


'''


##parameters
delim = '\t'

##dirs and file
xengsort_index = '/home/atimms/ngs_data/references/xengsort/mouse_human_ref_0720.h5'


##methods

def convert_bam_fastq_picard(s_dict):
	for sample in s_dict:
		bamfile = s_dict[sample]
		r1_fastq = sample + '.r1.fastq'
		r2_fastq = sample + '.r2.fastq'
		# picard_sam_fq = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
		picard_sam_fq = subprocess.Popen(['picard', 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
		picard_sam_fq.wait()
		# gzip_one = subprocess.Popen(['gzip', r1_fastq])
		# gzip_one.wait()
		# gzip_two = subprocess.Popen(['gzip', r2_fastq])
		# gzip_two.wait()

def run_xengsort(s_dict):
	for sample in s_dict:
		r1_fastq = sample + '.r1.fastq'
		r2_fastq = sample + '.r2.fastq'
		outdir = './' + sample + '_xe/'
		out_prefix = outdir + sample
		os.mkdir(outdir) 
		#xengsort classify --index <index> --classification new -T <threads> --fastq <fq1> --pairs <fq2> -o <out_dir>
		xengsort_classify = subprocess.Popen(['xengsort', 'classify', '--index', xengsort_index, '--classification', 'new', '-T', '8', '--fastq', r1_fastq, '--pairs', r2_fastq, '--prefix', out_prefix])
		xengsort_classify.wait()

def align_with_bwa_one_at_time(sample, r1_fq, r2_fq):
	rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
	post_bwa_bam = sample + '.bwa.bam'
	sort_bam = sample + '.bwa_sort.bam'
	mkdup_bam = sample + '.bwa_mkdup.bam'
	realigned_bam = sample + '.bwa_religned.bam'
	gatk_bam = sample + final_bam_suffix
	# mkdup_bam = sample + '.bwa_mkdup.bam'
	##bwa and convert to bam
	bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
	st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
	st_sam_bam_pe.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam, '-@', '10', '-m', '10G', post_bwa_bam])
	st_sort_pe.wait()
	##mark duplicates
	# picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md = subprocess.Popen([picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md.wait()
	##realign around indels
	gatk_ir = subprocess.Popen(['java', '-Xmx80g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam, '-L', exome_capture_bed, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
	gatk_ir.wait()
	##bqsr
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-L', exome_capture_bed, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-L', exome_capture_bed, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
	gatk_pr.wait()

def exome_analysis_pipeline(s_dict):
	for sample in s_dict:
		read1_fastq = sample + '.r1.fastq'
		read2_fastq = sample + '.r2.fastq'
		align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)

##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/rich_exomes_0720'
os.chdir(working_dir)

sample_dict = {'TZ001': 'TZ001.368077.bam', 'TZ002': 'TZ002.368078.bam', 'TZ003': 'TZ003.368079.bam', 'TZ004': 'TZ004.368080.bam', 
		'TZ005': 'TZ005.368081.bam', 'TZ006': 'TZ006.368082.bam', 'TZ007': 'TZ007.370762.bam', 'TZ008': 'TZ008.368084.bam', 
		'TZ009': 'TZ009.368085.bam'}

##make fastq files from bams
# convert_bam_fastq_picard(sample_dict)

##run xengsort
# sample_dict = ['TZ001']
run_xengsort(sample_dict)

##align, call vars and annotate
# exome_analysis_pipeline(sample_dict)



