#!/usr/bin/python
import sys
import subprocess
import os
import glob
import dobyns_gemini_pipeline_cybertron_v10

##parameters
delim = '\t'


##file names etc
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.b37.vcf'
indels_1000g = ref_dir + '1000G_phase1.indels.b37.vcf'
rtc_intervals = ref_dir + 'Mills_1000G.intervals'
dbsnp = ref_dir + 'dbsnp_138.b37.vcf'
exome_capture_bed = ref_dir + 'dobyns_exome.in_any_padded_target.1015.bed'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
exome_bed_for_coverage = ref_dir + 'dobyns_exome.in_all_targets.1015.bed'


##programs
samtools = 'samtools'
bedtools = 'bedtools'
picard = 'picard'
bwa = 'bwa'
gatk = '/cm/shared/apps/GATK/3.7/GenomeAnalysisTK.jar'
bcftools = 'bcftools'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
bgzip = 'bgzip'

plink = '/home/atimms/programs/plink'


def relign_target_creater_knowns():
	gatk_rtc = subprocess.Popen(['java', '-Xmx3g', '-jar', gatk, '-T', 'RealignerTargetCreator', '-R', fasta, '-nt', '15', '-o', rtc_intervals, '-known', indels_mills, '-known', indels_1000g])
	gatk_rtc.wait()

##trying to address issue with single end read
'''
def convert_bam_fastq_bedtools(bamfile, r1_fastq, r2_fastq):
	read_sorted_bam = bamfile[:-4] + 'n_sorted.bam'
	st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-u', '-@', '5', '-f', '2', bamfile], stdout=subprocess.PIPE)
	st_n_sort = subprocess.Popen([samtools, 'sort', '-nO', 'bam', '-o', read_sorted_bam, '-@', '10', '-m', '10G', '-T', 'tempy', '-'], stdin=st_sam_bam_pe.stdout)
	st_n_sort.wait()
	bam_fq = subprocess.Popen([bedtools, 'bamtofastq', '-i', read_sorted_bam, '-fq', r1_fastq, '-fq2', r2_fastq])
	bam_fq.wait()
'''

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
	gatk_ir = subprocess.Popen(['java', '-Xmx80g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam, '-L', exome_capture_bed, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
	gatk_ir.wait()
	##bqsr
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-L', exome_capture_bed, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-L', exome_capture_bed, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
	gatk_pr.wait()

def align_with_bwa_one_at_time_single_end(sample, r1_fq):
	# rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
	rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ionxpress"
	post_bwa_bam = sample + '.bwa.bam'
	sort_bam = sample + '.bwa_sort.bam'
	mkdup_bam = sample + '.bwa_mkdup.bam'
	realigned_bam = sample + '.bwa_religned.bam'
	gatk_bam = sample + final_bam_suffix
	# mkdup_bam = sample + '.bwa_mkdup.bam'
	##bwa and convert to bam
	bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq], stdout=subprocess.PIPE)
	st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
	st_sam_bam_pe.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam,  '-@', '10', '-m', '10G', post_bwa_bam])
	st_sort_pe.wait()
	##mark duplicates
	# picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md = subprocess.Popen([picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md.wait()
	##realign around indels
	gatk_ir = subprocess.Popen(['java', '-Xmx80g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam,'-L', exome_capture_bed, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
	gatk_ir.wait()
	##bqsr
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '15', '-R', fasta, '-I', realigned_bam,'-L', exome_capture_bed, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '15', '-R', fasta, '-I', realigned_bam,'-L', exome_capture_bed, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
	gatk_pr.wait()


##delete files supplied in list
def delete_unwated_files(file_extensions):
	files_to_delete = []
	for ext in file_extensions:
		files = glob.glob(ext)
		files_to_delete.extend(files)
	for f in files_to_delete:
		os.remove(f)

##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

##call freebayes in all sample listed in list of bamfiles
def variant_calling_freebayes(bamlist, name_prefix):
	vcf_temp1 = 'temp_fb1.vcf'
	vcf_temp2 = 'temp_fb2.vcf.gz'
	final_vcf = name_prefix + '.freebayes.vcf.gz'
	vcf_fh = open(vcf_temp1, 'w')
	freebayes_run = subprocess.Popen([freebayes,'-f' , fasta, '-L', bamlist, '-t', exome_capture_bed], stdout=vcf_fh)
	freebayes_run.wait()
	vcf_fh.close()
	bgzip_run = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()

##not used!!
def variant_calling_gatk_hc(sample_dict, name_prefix):
	vcf_temp0 = 'temp_gatk0.vcf'
	vcf_raw_snps = 'temp_raw_snps.vcf'
	vcf_filtered_snps = 'temp_filtered_snps.vcf'
	vcf_raw_indels = 'temp_raw_indels.vcf'
	vcf_filtered_indels = 'temp_filtered_indels.vcf'
	vcf_temp1 = 'temp_gatk1.vcf'
	vcf_temp2 = 'temp_gatk2.vcf.gz'
	final_vcf = name_prefix + '.gatkHC.vcf.gz'
	gvcf_files = []
	##get g.vcfs for each sample (add -L)
	for sample in sample_dict:
		gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', sample + final_bam_suffix,'-L', exome_capture_bed, '--emitRefConfidence', 'GVCF', '--variant_index_type', 'LINEAR', '--variant_index_parameter', '128000', '-o', sample + '.g.vcf'])
		gatk_hc.wait()
		gvcf = ['-V', sample + '.g.vcf']
		gvcf_files.extend(gvcf)
	print gvcf_files
	##genotype g.vcfs(add -L)
	command = ['java', '-Xmx20g', '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', fasta, '-nt', '15','-L', exome_capture_bed] + gvcf_files + ['-o', vcf_temp0]
	gatk_gg = subprocess.Popen(command)
	gatk_gg.wait()
	##split data in SNPs and indels and apply manual variant filtering (add -L)
	snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '15', '-V', vcf_temp0, '-o', vcf_raw_snps, '-selectType', 'SNP'])
	snp_cut.wait()
	snp_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_snps, '-o', vcf_filtered_snps, '--filterExpression ', "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "manual_filter"])
	snp_vf.wait()
	indel_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '15', '-V', vcf_temp0, '-o', vcf_raw_indels, '-selectType', 'INDEL'])
	indel_cut.wait()
	indel_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_indels, '-o', vcf_filtered_indels, '--filterExpression ', "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0", '--filterName', "manual_filter"])
	indel_vf.wait()
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '15', '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', vcf_temp1, '-genotypeMergeOptions', 'UNSORTED'])
	combine_var.wait()
	bgzip_cmd = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_cmd.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()


def variant_calling_gatk_hc_no_gvcf(bamlist, name_prefix):
	vcf_temp0 = 'temp_gatk0.vcf'
	vcf_raw_snps = 'temp_raw_snps.vcf'
	vcf_filtered_snps = 'temp_filtered_snps.vcf'
	vcf_raw_indels = 'temp_raw_indels.vcf'
	vcf_filtered_indels = 'temp_filtered_indels.vcf'
	vcf_temp1 = 'temp_gatk1.vcf'
	vcf_temp2 = 'temp_gatk2.vcf.gz'
	final_vcf = name_prefix + '.gatkHC.vcf.gz'
	bam_files = []
	##run haplotype caller
	gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', bamlist,'-L', exome_capture_bed, '-o', vcf_temp0])
	gatk_hc.wait()
	##split data in SNPs and indels and apply manual variant filtering (add -L)
	snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_snps, '-selectType', 'SNP'])
	snp_cut.wait()
	snp_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_snps, '-o', vcf_filtered_snps, '--filterExpression ', "QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "snp_filter"])
	snp_vf.wait()
	indel_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_indels, '-selectType', 'INDEL'])
	indel_cut.wait()
	indel_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_indels, '-o', vcf_filtered_indels, '--filterExpression ', "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0", '--filterName', "indel_filter"])
	indel_vf.wait()
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '5', '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', vcf_temp1, '-genotypeMergeOptions', 'UNSORTED'])
	combine_var.wait()
	bgzip_cmd = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_cmd.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()

def intesect_two_vcf_files(vcf1, vcf2, output_dir):
	bcf_isec = subprocess.Popen([bcftools, 'isec', vcf1, vcf2, '-p', output_dir ])
	bcf_isec.wait()

##get exome coverage for bam files
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

def plink_relatadness_check(vcf, file_prefix):
	##bzgip vcf file
	if os.path.isfile(vcf):
		bgzip_cmd = subprocess.Popen([bgzip, vcf])
		bgzip_cmd.wait()
	##correct filtering??
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(DP)>50", '-V', 'indels', '-o', 'temp_plink.vcf.gz', '-O', 'z', vcf + '.gz'])
	bcftools_filter.wait()
	##generate plink file from vcf
	make_plink = subprocess.Popen([plink, '--vcf', 'temp_plink.vcf.gz', '--vcf-filter', '--const-fid', '--out', 'temp.pass_q50_dp50'])
	make_plink.wait()
	##check sex -results in .sexcheck file
	plink_sex = subprocess.Popen([plink, '--bfile', 'temp.pass_q50_dp50', '--check-sex', '--out', file_prefix + '.pass_q50_dp50'])
	plink_sex.wait()
	##ibd check
	plink_ibd = subprocess.Popen([plink,  '--bfile', 'temp.pass_q50_dp50', '--genome', '--out', file_prefix + '.pass_q50_dp50'])
	plink_ibd.wait()


def run_scalpel(samples, bam_suffix, prefix):
	print samples
	if len(samples) == 3:
		scalpel_trio_disco = subprocess.Popen([scalpel_discovery, '--denovo', '--dad', samples[1] + bam_suffix, '--mom', samples[2] + bam_suffix, '--aff', samples[0] + bam_suffix, '--sib', samples[0] + bam_suffix, '--bed', exome_capture_bed, '--ref', fasta, '--numprocs', '10', '--two-pass'])
		scalpel_trio_disco.wait()
		# scalpel_trio_export = subprocess.Popen([scalpel_discovery, '--denovo', '--bed', exome_capture_bed, '--ref', fasta])
		# scalpel_trio_export.wait()

##make intervals file for indel realinment step
# relign_target_creater_knowns()

##call all methods
##sample list needed if using scalpel (proband, dad, mom, sib(if available)), but don't use it not using that program
# def call_all_exome_methods(working_dir, prefix, bam_or_fastq, sample_dict, sample_list):
def call_all_exome_methods_inc_gemini(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
	os.chdir(working_dir)
	##covert bam to fastq and then process
	# '''
	if bam_or_fastq == 'bam':
		for sample in sample_dict:
			read1_fastq = sample + '.r1.fastq'
			read2_fastq = sample + '.r2.fastq'
			# convert_bam_fastq_bedtools(sample_dict[sample], read1_fastq, read2_fastq)
			convert_bam_fastq_picard(sample_dict[sample], read1_fastq, read2_fastq)
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	elif bam_or_fastq == 'messy_bam':
		for sample in sample_dict:
			read1_fastq = sample + '.r1.fastq'
			read2_fastq = sample + '.r2.fastq'
			convert_bam_fastq_bedtools(sample_dict[sample], read1_fastq, read2_fastq)
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	elif bam_or_fastq == 'fastq':
		for sample in sample_dict:
			read1_fastq = sample_dict[sample][0]
			read2_fastq = sample_dict[sample][1]
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
			##single fq
			# align_with_bwa_one_at_time_single_end(sample,sample_dict[sample])
	else:
		print 'must specify if starting with fastq or bam file'
	# '''
	##variant calling
	make_list_of_bams(sample_dict, final_bam_suffix, bamlist)
	variant_calling_freebayes(bamlist, prefix)
	variant_calling_gatk_hc_no_gvcf(bamlist, prefix)
	intesect_two_vcf_files(prefix + '.gatkHC.vcf.gz', prefix + '.freebayes.vcf.gz', prefix + '.intersected_vcfs')
	##coverage
	calculate_exome_coverage(sample_dict, final_bam_suffix, exome_bed_for_coverage, prefix)
	##plink for sex and relatedness -- work in progress
	plink_relatadness_check(prefix + '.intersected_vcfs/0002.vcf', prefix)
	##delete unwanted files
	delete_unwated_files(['*.r1.fastq', '*.r2.fastq', '*.bwa.bam', '*.bwa_sort.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*', '*bwa_mkdup.bai', '*bwa_mkdup.bam'])
	##run gemini
	dobyns_gemini_pipeline_cybertron_v10.standard_gemini_protocol(working_dir, prefix, ped_type)

def call_just_gemini(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
	os.chdir(working_dir)
	print 'hello'
	##run gemini
	dobyns_gemini_pipeline_cybertron_v10.standard_gemini_protocol(working_dir, prefix, ped_type)

def call_gatk_hc_only(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
	os.chdir(working_dir)
	##variant calling
	make_list_of_bams(sample_dict, final_bam_suffix, bamlist)
	variant_calling_gatk_hc_no_gvcf(bamlist, prefix)

def call_after_mapping(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
	os.chdir(working_dir)
	print 'hello'
	##variant calling
	make_list_of_bams(sample_dict, final_bam_suffix, bamlist)
	variant_calling_freebayes(bamlist, prefix)
	variant_calling_gatk_hc_no_gvcf(bamlist, prefix)
	intesect_two_vcf_files(prefix + '.gatkHC.vcf.gz', prefix + '.freebayes.vcf.gz', prefix + '.intersected_vcfs')
	##coverage
	calculate_exome_coverage(sample_dict, final_bam_suffix, exome_bed_for_coverage, prefix)
	##plink for sex and relatedness -- work in progress
	plink_relatadness_check(prefix + '.intersected_vcfs/0002.vcf', prefix)
	##delete unwanted files
	delete_unwated_files(['*.r1.fastq', '*.r2.fastq', '*.bwa.bam', '*.bwa_sort.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*', '*bwa_mkdup.bai', '*bwa_mkdup.bam'])
	##run gemini
	dobyns_gemini_pipeline_cybertron_v10.standard_gemini_protocol(working_dir, prefix, ped_type)

def call_all_exome_methods_without_gemini(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
	os.chdir(working_dir)
	##covert bam to fastq and then process
	# '''
	if bam_or_fastq == 'bam':
		for sample in sample_dict:
			read1_fastq = sample + '.r1.fastq'
			read2_fastq = sample + '.r2.fastq'
			# convert_bam_fastq_bedtools(sample_dict[sample], read1_fastq, read2_fastq)
			convert_bam_fastq_picard(sample_dict[sample], read1_fastq, read2_fastq)
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	elif bam_or_fastq == 'messy_bam':
		for sample in sample_dict:
			read1_fastq = sample + '.r1.fastq'
			read2_fastq = sample + '.r2.fastq'
			convert_bam_fastq_bedtools(sample_dict[sample], read1_fastq, read2_fastq)
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	elif bam_or_fastq == 'fastq':
		for sample in sample_dict:
			read1_fastq = sample_dict[sample][0]
			read2_fastq = sample_dict[sample][1]
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
			##single fq
			# align_with_bwa_one_at_time_single_end(sample,sample_dict[sample])
	else:
		print 'must specify if starting with fastq or bam file'
	# '''
	##variant calling
	make_list_of_bams(sample_dict, final_bam_suffix, bamlist)
	variant_calling_freebayes(bamlist, prefix)
	variant_calling_gatk_hc_no_gvcf(bamlist, prefix)
	intesect_two_vcf_files(prefix + '.gatkHC.vcf.gz', prefix + '.freebayes.vcf.gz', prefix + '.intersected_vcfs')
	##coverage
	calculate_exome_coverage(sample_dict, final_bam_suffix, exome_bed_for_coverage, prefix)
	##plink for sex and relatedness -- work in progress
	plink_relatadness_check(prefix + '.intersected_vcfs/0002.vcf', prefix)
	##delete unwanted files
	delete_unwated_files(['*.r1.fastq', '*.r2.fastq', '*.bwa.bam', '*.bwa_sort.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*', '*bwa_mkdup.bai', '*bwa_mkdup.bam'])


def intersect_vcfs_and_gemini(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
	os.chdir(working_dir)
	intesect_two_vcf_files(prefix + '.gatkHC.vcf.gz', prefix + '.freebayes.vcf.gz', prefix + '.intersected_vcfs')
	##coverage
	# calculate_exome_coverage(sample_dict, final_bam_suffix, exome_bed_for_coverage, prefix)
	##plink for sex and relatedness -- work in progress
	plink_relatadness_check(prefix + '.intersected_vcfs/0002.vcf', prefix)
	##delete unwanted files
	# delete_unwated_files(['*.r1.fastq', '*.r2.fastq', '*.bwa.bam', '*.bwa_sort.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*', '*bwa_mkdup.bai', '*bwa_mkdup.bam'])
	##run gemini
	dobyns_gemini_pipeline_cybertron_v10.standard_gemini_protocol(working_dir, prefix, ped_type)

def just_plink(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
	os.chdir(working_dir)
	plink_relatadness_check(prefix + '.intersected_vcfs/0002.vcf', prefix)

def call_freebayes_only(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
	os.chdir(working_dir)
	##variant calling
	make_list_of_bams(sample_dict, final_bam_suffix, bamlist)
	variant_calling_freebayes(bamlist, prefix)

def call_just_mapping(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
	os.chdir(working_dir)
	##covert bam to fastq and then process
	# '''
	if bam_or_fastq == 'bam':
		for sample in sample_dict:
			read1_fastq = sample + '.r1.fastq'
			read2_fastq = sample + '.r2.fastq'
			# convert_bam_fastq_bedtools(sample_dict[sample], read1_fastq, read2_fastq)
			convert_bam_fastq_picard(sample_dict[sample], read1_fastq, read2_fastq)
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	elif bam_or_fastq == 'messy_bam':
		for sample in sample_dict:
			read1_fastq = sample + '.r1.fastq'
			read2_fastq = sample + '.r2.fastq'
			convert_bam_fastq_bedtools(sample_dict[sample], read1_fastq, read2_fastq)
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	elif bam_or_fastq == 'fastq':
		for sample in sample_dict:
			read1_fastq = sample_dict[sample][0]
			read2_fastq = sample_dict[sample][1]
			align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
			##single fq
			# align_with_bwa_one_at_time_single_end(sample,sample_dict[sample])
	else:
		print 'must specify if starting with fastq or bam file'


