#!/usr/bin/python
import sys
import subprocess
import os
import glob
# import dobyns_gemini_pipeline_cybertron_v8

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
gatk = '/cm/shared/apps/GATK/3.7/GenomeAnalysisTK.jar'
StrelkaGermlineWorkflow = '/home/atimms/programs/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py'
plink = '/home/atimms/programs/plink'
# samtools = 'samtools'
# bedtools = 'bedtools'
# picard = 'picard'
# bwa = 'bwa'
# gatk = '/cm/shared/apps/GATK/3.7/GenomeAnalysisTK.jar'
# bcftools = 'bcftools'
# freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
# bgzip = 'bgzip'


##methods
def convert_cram_fastq(cramfile, r1_fastq, r2_fastq):
	print(cramfile, r1_fastq, r2_fastq)
	sample_name = cramfile.split('.')[0]
	temp_cram = sample_name + 'temp.cram'
	st_sort_pe = subprocess.Popen(['samtools', 'sort', '-n', '-O', 'cram', '-T', sample_name, '-o', temp_cram, '-@', '10', '-m', '10G', cramfile])
	st_sort_pe.wait()
	# st_index = subprocess.Popen(['samtools', 'index', '-@', '20', temp_cram])
	# st_index.wait()
	picard_sam_fq = subprocess.Popen(['samtools', 'fastq', '-@', '20', '-1', r1_fastq, '-2', r2_fastq, temp_cram])
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
	'''
	bwa_pe = subprocess.Popen(['bwa', 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
	st_sam_bam_pe = subprocess.Popen(['samtools', 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
	st_sam_bam_pe.wait()
	##sort bam
	st_sort_pe = subprocess.Popen(['samtools', 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam, '-@', '10', '-m', '10G', post_bwa_bam])
	st_sort_pe.wait()
	##mark duplicates
	# picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md = subprocess.Popen(['picard', 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md.wait()
	'''
	##realign around indels
	gatk_ir = subprocess.Popen(['java', '-Xmx80g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
	gatk_ir.wait()
	##bqsr
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
	gatk_pr.wait()


##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

def variant_calling_gatk_hc_no_gvcf(bamlist, name_prefix):
	vcf_temp0 = name_prefix + 'temp_gatk0.vcf'
	vcf_raw_snps = name_prefix + 'temp_raw_snps.vcf'
	vcf_filtered_snps = name_prefix + 'temp_filtered_snps.vcf'
	vcf_raw_indels = name_prefix + 'temp_raw_indels.vcf'
	vcf_filtered_indels = name_prefix + 'temp_filtered_indels.vcf'
	vcf_temp1 = name_prefix + 'temp_gatk1.vcf'
	vcf_temp2 = name_prefix + 'temp_gatk2.vcf.gz'
	final_vcf = name_prefix + '.gatkHC.vcf.gz'
	bam_files = []
	##run haplotype caller
	# gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', bamlist, '-o', vcf_temp0])
	##using nct
	# '''
	gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-nct', '10', '-R', fasta, '-I', bamlist, '-o', vcf_temp0])
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
	# '''
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '5', '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', vcf_temp1, '-genotypeMergeOptions', 'UNSORTED'])
	combine_var.wait()
	bgzip_cmd = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip_cmd.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()

def variant_calling_strelka(sample_dict, bam_suffix, ped_name, w_dir):
	# ${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py \
	# --bam NA12878.cram \
	# --bam NA12891.cram \
	# --bam NA12892.cram \
	# --referenceFasta hg19.fa \
	# --runDir ${STRELKA_ANALYSIS_PATH}
	##set up var calling and run
	bam_cmds = []
	for sample in sample_dict:
		bam = sample + bam_suffix
		bam_cmds.extend(['--bam', w_dir + '/' + bam])
	analysis_dir = w_dir + '/' + ped_name + '_strelka'
	'''
	##configure
	st_config = subprocess.Popen([StrelkaGermlineWorkflow] + bam_cmds + ['--referenceFasta', fasta, '--runDir', analysis_dir])
	st_config.wait()	
	##run ?? cores
	#${STRELKA_ANALYSIS_PATH}/runWorkflow.py -m sge -j 36
	# st_run = subprocess.Popen([analysis_dir + '/runWorkflow.py', '-m', 'sge', '-j', '20'])
	st_run = subprocess.Popen([analysis_dir + '/runWorkflow.py', '-m', 'local', '-j', '10'])
	st_run.wait()
	'''
	##get passed variants and normalize
	st_gvcf = analysis_dir + '/results/variants/variants.vcf.gz'
	vcf_temp1 = ped_name + 'temp_st1.vcf.gz'
	vcf_temp2 = ped_name + 'temp_st2.vcf.gz'	
	final_vcf = ped_name + '.strelka.vcf.gz'
	##get filtered
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-o', vcf_temp1, '-O', 'z', st_gvcf])
	bcftools_filter.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()

def intesect_two_vcf_files(vcf1, vcf2, output_dir):
	bcf_isec = subprocess.Popen(['bcftools', 'isec', vcf1, vcf2, '-p', output_dir ])
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

def plink_relatadness_check(vcf, file_prefix):
	##bzgip vcf file
	if os.path.isfile(vcf):
		bgzip_cmd = subprocess.Popen(['bgzip', vcf])
		bgzip_cmd.wait()
	##correct filtering??
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(DP)>50", '-V', 'indels', '-o', 'temp_plink.vcf.gz', '-O', 'z', vcf + '.gz'])
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

##delete files supplied in list
def delete_unwated_files(file_extensions):
	files_to_delete = []
	for ext in file_extensions:
		files = glob.glob(ext)
		files_to_delete.extend(files)
	for f in files_to_delete:
		os.remove(f)

def convert_to_fastq_and_map(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
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
	elif bam_or_fastq == 'cram':
		for sample in sample_dict:
			read1_fastq = sample + '.r1.fastq.gz'
			read2_fastq = sample + '.r2.fastq.gz'
			# convert_cram_fastq(sample_dict[sample], read1_fastq, read2_fastq)
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

def variants_calling_qc_and_gemini(working_dir, prefix, bam_or_fastq, sample_dict, ped_type):
	os.chdir(working_dir)
	'''
	##variant calling
	variant_calling_strelka(sample_dict, final_bam_suffix, prefix,working_dir)
	make_list_of_bams(sample_dict, final_bam_suffix, bamlist)
	variant_calling_gatk_hc_no_gvcf(bamlist, prefix)
	
	intesect_two_vcf_files(prefix + '.gatkHC.vcf.gz', prefix + '.strelka.vcf.gz', prefix + '.intersected_vcfs')
	##coverage
	calculate_exome_coverage(sample_dict, final_bam_suffix, exome_bed_for_coverage, prefix)
	##plink for sex and relatedness -- work in progress
	plink_relatadness_check(prefix + '.intersected_vcfs/0002.vcf', prefix)
	'''
	##delete unwanted files
	delete_unwated_files(['*.r1.fastq', '*.r2.fastq', '*.bwa.bam', '*.bwa_sort.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*', '*bwa_mkdup.bai', '*bwa_mkdup.bam'])
	
	##run gemini
	# dobyns_gemini_pipeline_cybertron_v8.standard_gemini_protocol(working_dir, prefix, ped_type)






##run methods
working_dir = '/home/atimms/ngs_data/genomes/daniela_test_0719'

##for converting to fastq and mappaing split into individual
# convert_to_fastq_and_map(working_dir, 'F003', 'cram', {'F03-00006':'F03-00006.cram'}, 'trio')
# convert_to_fastq_and_map(working_dir, 'F003', 'cram', {'F03-00007':'F03-00007.cram'}, 'trio')
# convert_to_fastq_and_map(working_dir, 'F003', 'cram', {'F03-00008':'F03-00008.cram'}, 'trio')
# convert_to_fastq_and_map(working_dir, 'F005', 'cram', {'F03-00011':'F03-00011.cram'}, 'trio')
# convert_to_fastq_and_map(working_dir, 'F005', 'cram', {'F03-00050':'F03-00050.cram'}, 'trio')
# convert_to_fastq_and_map(working_dir, 'F005', 'cram', {'F03-00012':'F03-00012.cram'}, 'trio')
##variant calling and gemini by ped
variants_calling_qc_and_gemini(working_dir, 'F003', 'cram', {'F03-00006':'F03-00006.cram', 'F03-00007':'F03-00007.cram', 'F03-00008':'F03-00008.cram'}, 'trio')
# variants_calling_qc_and_gemini(working_dir, 'F005', 'cram', {'F03-00011':'F03-00011.cram', 'F03-00050':'F03-00050.cram', 'F03-00012':'F03-00012.cram'}, 'trio')




