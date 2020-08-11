#!/usr/bin/env python
import subprocess
import os
import glob

'''
info....


'''

##parameters
delim = '\t'
thread_number = '18'

##ref files and file names
genome_name = 'GRCh38'
star_index_dir = '/home/atimms/ngs_data/references/star/' + genome_name
fa_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name +'/genes.gtf'
# star_bam_suffix = '.Aligned.out.bam'
star_bam_suffix = '.Aligned.sortedByCoord.out.bam'
proccessed_bam_suffix = '.star.gatk.bam'
bamlist = 'bams.list'
rs3184504_bed = 'rs3184504.bed'

##programs
picard = '/home/atimms/programs/picard_2.19/picard.jar'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
star = '/home/atimms/programs/STAR-2.5.3a/bin/Linux_x86_64_static/STAR'

##methods
def make_star_index_files(star_genome_dir, genome_fas, genome_gtf, threads_to_use):
	star_index = subprocess.Popen([star, '--runMode', 'genomeGenerate', '--genomeDir', star_genome_dir, '--genomeFastaFiles', genome_fas, '--sjdbGTFfile', genome_gtf, '--runThreadN', threads_to_use])
	star_index.wait()

def star_align_paired_end_2_pass(sample_name, star_genome_dir, threads_to_use, r1_fq, r2_fq):
	print sample_name, r1_fq, r2_fq
	# star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--outSAMmapqUnique', '60', '--readFilesCommand', 'zcat', '--twopassMode', 'Basic', '--outFileNamePrefix', sample_name + '.', '--runThreadN', threads_to_use])
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMmapqUnique', '60', '--readFilesCommand', 'zcat', '--twopassMode', 'Basic', '--outFileNamePrefix', sample_name + '.', '--runThreadN', threads_to_use])
	star_align.wait()
	bam_file = sample_name + star_bam_suffix
	samtools_index = subprocess.Popen(['samtools', 'index', bam_file])
	samtools_index.wait()


def proccess_bam_files(sample_name):
	in_bam = sample_name + star_bam_suffix
	arg_bam = sample_name + '.with_rg.bam'
	mkdup_bam = sample_name + '.mkdup.bam'
	reorder_bam = sample_name + '.reorder.bam'
	split_bam = sample_name + '.split.bam'
	final_bam = sample_name + proccessed_bam_suffix
	picard_arg = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'AddOrReplaceReadGroups', 'INPUT=' + in_bam, 'OUTPUT=' + arg_bam, 'SO=coordinate', 'RGID=' + sample_name, 'RGLB=' + sample_name, 'RGPL=Illumina', 'RGPU=machine1', 'RGSM=' + sample_name ])
	picard_arg.wait()	
	picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'INPUT=' + arg_bam, 'OUTPUT=' + mkdup_bam,'CREATE_INDEX=true', 'METRICS_FILE=' + sample_name + '.metrics'])
	picard_md.wait()
	picard_rs = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'ReorderSam', 'INPUT=' + mkdup_bam, 'OUTPUT=' + reorder_bam, 'CREATE_INDEX=true', 'REFERENCE=' + fa_file])
	picard_rs.wait()
	gatk_snt = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'SplitNCigarReads',  '-R', fa_file, '-I', reorder_bam, '-o', split_bam, '-U', 'ALLOW_N_CIGAR_READS'])
	gatk_snt.wait()
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '6', '-R', fa_file, '-I', split_bam, '-knownSites', indels_mills, '-knownSites', dbsnp, '-o', sample_name + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '6', '-R', fa_file, '-I', split_bam, '-BQSR', sample_name + '.recal_data.table', '-o', final_bam])
	gatk_pr.wait()
	##remove intermeidate files
	files_to_go = glob.glob('*with_rg.bam') + glob.glob('*mkdup.bam') + glob.glob('*mkdup.bai') + glob.glob('*metrics') + glob.glob('*reorder.bam')+ glob.glob('*reorder.bai') + glob.glob('*split.bam') + glob.glob('*split.bai')
	print 'removing files:'
	for f in files_to_go:
		os.remove(f)
		print f


def add_rg_to_bam(sample_name):
	in_bam = sample_name + star_bam_suffix
	arg_bam = sample_name + '.with_rg.bam'
	picard_arg = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'AddOrReplaceReadGroups', 'INPUT=' + in_bam, 'OUTPUT=' + arg_bam, 'SO=coordinate','CREATE_INDEX=true', 'RGID=' + sample_name, 'RGLB=' + sample_name, 'RGPL=Illumina', 'RGPU=machine1', 'RGSM=' + sample_name ])
	picard_arg.wait()	

def make_patient_fq_dict(in_file):
	patient_dict = {}
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc > 1:
				patient = line[5]
				srr_id = line[2]
				fq1 = srr_id + '.sra_1.fastq.gz'
				fq2 = srr_id + '.sra_2.fastq.gz'
				if patient in patient_dict:
					patient_dict[patient][0].append(fq1)
					patient_dict[patient][1].append(fq2)
				else:
					patient_dict[patient] = [[fq1], [fq2]]
	return(patient_dict)

def combine_fq_file(r1_to_combine, r2_to_combine, r1_fq, r2_fq):
	print r1_fq, r1_to_combine
	print r2_fq, r2_to_combine
	with open(r1_fq, 'w') as r1_fh:
		cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
		cat_files.wait()
	with open(r2_fq, 'w') as r2_fh:
		cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
		cat_files.wait()

##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

def get_vars_around_rs3184504(bam_list, name_prefix):
	final_vcf = name_prefix + '.gatkHC.vcf.gz'
	gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fa_file, '-L', rs3184504_bed, '-I', bam_list, '-dontUseSoftClippedBases', '-stand_call_conf', '20', '--filter_reads_with_N_cigar', '-o', final_vcf])
	gatk_hc.wait()

def genotype_rs3184504_pipeline(info_file):
	patient_fq_dict = make_patient_fq_dict(info_file)
	project = info_file.split('_')[0] + '_0820'
	for patient in patient_fq_dict:
		print(patient, patient_fq_dict[patient])
		fq1s = patient_fq_dict[patient][0]
		fq2s = patient_fq_dict[patient][1]
		fq1 = patient + '_combined.r1.fq.gz'
		fq2 = patient + '_combined.r2.fq.gz'
		##combine fq i.e. 2 per patient
		# combine_fq_file(fq1s, fq2s, fq1, fq2)
		##align with star
		# star_align_paired_end_2_pass(patient, star_index_dir, thread_number, fq1, fq2)
		##process files ? needed for just one snp?
		# proccess_bam_files(patient)
		##just add rg so haplotype caller works
		# add_rg_to_bam(patient)
	##call variants
	make_list_of_bams(patient_fq_dict, '.with_rg.bam', bamlist)
	get_vars_around_rs3184504(bamlist, project)



##run methods
working_dir = '/home/atimms/ngs_data/rnaseq/eric_rnaseq_0720/GSE131411'
os.chdir(working_dir)			
GSE131411_info_file = 'GSE131411_info.txt'

genotype_rs3184504_pipeline(GSE131411_info_file)


