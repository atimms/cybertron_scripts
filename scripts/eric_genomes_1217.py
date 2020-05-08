#!/usr/bin/python
import sys
import subprocess
import os
import glob
import dobyns_gemini_pipeline_cybertron_v7

##parameters
delim = '\t'


##file names etc
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.b37.vcf'
indels_1000g = ref_dir + '1000G_phase1.indels.b37.vcf'
rtc_intervals = ref_dir + 'Mills_1000G.intervals'
dbsnp = ref_dir + 'dbsnp_138.b37.vcf'
# exome_capture_bed = ref_dir + 'dobyns_exome.in_any_padded_target.1015.bed'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
exome_bed_for_coverage = ref_dir + 'dobyns_exome.in_all_targets.1015.bed'


##programs
samtools = 'samtools'
bedtools = 'bedtools'
picard = 'picard'
bwa = 'bwa'
gatk = '/cm/shared/apps/GATK/GenomeAnalysisTK.jar'
bcftools = 'bcftools'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
bgzip = 'bgzip'

plink = '/home/atimms/programs/plink'


def relign_target_creater_knowns():
	gatk_rtc = subprocess.Popen(['java', '-Xmx3g', '-jar', gatk, '-T', 'RealignerTargetCreator', '-R', fasta, '-nt', '15', '-o', rtc_intervals, '-known', indels_mills, '-known', indels_1000g])
	gatk_rtc.wait()

def convert_bam_fastq_bedtools(bamfile, r1_fastq, r2_fastq):
	read_sorted_bam = bamfile[:-4] + 'n_sorted.bam'
	st_n_sort = subprocess.Popen([samtools, 'sort', '-nO', 'bam', '-o', read_sorted_bam, '-@', '10', '-m', '10G', '-T', 'tempy', bamfile])
	st_n_sort.wait()
	bam_fq = subprocess.Popen([bedtools, 'bamtofastq', '-i', read_sorted_bam, '-fq', r1_fastq, '-fq2', r2_fastq])
	bam_fq.wait()

def convert_bam_fastq_picard(bamfile, r1_fastq, r2_fastq):
	# picard_sam_fq = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
	picard_sam_fq = subprocess.Popen([picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])

	picard_sam_fq.wait()

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
	gatk_ir = subprocess.Popen(['java', '-Xmx80g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
	gatk_ir.wait()
	##bqsr
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
	gatk_pr.wait()


def call_all_methods(working_dir, sample_list):
	os.chdir(working_dir)
	for sample in sample_list:
		bam_file = sample + '.bam'
		read1_fastq = sample + '.r1.fastq'
		read2_fastq = sample + '.r2.fastq'
		convert_bam_fastq_picard(bam_file, read1_fastq, read2_fastq)
		align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)







working_dir = '/home/atimms/ngs_data/genomes/eric_genomes_1217'
samples = ['SEA001_Patient_I.a', 'SEA001_Patient_II.a', 'SEA001_Patient_II.d', 'SEA001_Patient_I.b', 'SEA001_Patient_II.b']
##complete in one shot
call_all_methods(working_dir, samples)
