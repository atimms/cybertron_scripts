#!/usr/bin/python
import sys
import subprocess
import os
import glob

##info 
'''
combine all daniela's exomes
cp bwa_gatk.bam from:
/home/atimms/ngs_data/exomes/backed_up/daniela_0519
/archive/luquetti_d/microtia_exomes
/home/atimms/ngs_data/exomes/working/daniela_mosaic_exomes_redo_1219


'''


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
gatk = '/cm/shared/apps/GATK/3.7/GenomeAnalysisTK.jar'



##methods
def combine_fq_file(r1_to_combine, r2_to_combine, r1_fq, r2_fq):
	print r1_fq, r1_to_combine
	print r2_fq, r2_to_combine
	with open(r1_fq, 'w') as r1_fh:
		cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
		cat_files.wait()
	with open(r2_fq, 'w') as r2_fh:
		cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
		cat_files.wait()


def convert_bam_fastq_picard(bamfile, r1_fastq, r2_fastq):
	# picard_sam_fq = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
	picard_sam_fq = subprocess.Popen(['picard', 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
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
	##realign around indels
	gatk_ir = subprocess.Popen(['java', '-Xmx80g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam, '-L', exome_capture_bed, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
	gatk_ir.wait()
	##bqsr
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-L', exome_capture_bed, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-L', exome_capture_bed, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
	gatk_pr.wait()

def variant_calling_gatk_hc(samples):
	gvcf_files = []
	##get g.vcfs for each sample (add -L)
	for sample in samples:
		st_index = subprocess.Popen(['samtools', 'index', sample + final_bam_suffix])
		st_index.wait()
		gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', sample + final_bam_suffix, '-L', exome_capture_bed, '--emitRefConfidence', 'GVCF', '--variant_index_type', 'LINEAR', '--variant_index_parameter', '128000', '-o', sample + '.g.vcf'])
		gatk_hc.wait()
		gvcf = ['-V', sample + '.g.vcf']
		gvcf_files.extend(gvcf)
	return gvcf_files


def combine_gvcfs(samples, out_vcf):
	gvcf_files = []
	for sample in samples:
		gvcf = ['-V', sample + '.g.vcf']
		gvcf_files.extend(gvcf)
	print gvcf_files
	gatk_hc = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineGVCFs', '-R', fasta] + gvcf_files + ['-o', out_vcf])
	gatk_hc.wait()

def genotype_gvcfs_and_filter(in_gvcfs, name_prefix):
	vcf_temp0 = name_prefix + 'temp_gatk0.vcf'
	vcf_raw_snps = name_prefix + 'temp_raw_snps.vcf'
	vcf_filtered_snps = name_prefix + 'temp_filtered_snps.vcf'
	vcf_raw_indels = name_prefix + 'temp_raw_indels.vcf'
	vcf_filtered_indels = name_prefix + 'temp_filtered_indels.vcf'
	vcf_temp1 = name_prefix + 'temp_gatk1.vcf'
	vcf_temp2 = name_prefix + 'temp_gatk2.vcf.gz'
	final_vcf = name_prefix + '.gatkHC.vcf.gz'
	gvcf_cmds = []
	for in_gvcf in in_gvcfs:
		gvcf = ['-V', in_gvcf]
		gvcf_cmds.extend(gvcf)
	print gvcf_cmds
	##genotype g.vcfs
	##with nt 
	# command = ['java', '-Xmx100g', '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', fasta, '-nt', '10'] + gvcf_cmds + ['-o', vcf_temp0]
	##straight up
	command = ['java', '-Xmx100g', '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', fasta] + gvcf_cmds + ['-o', vcf_temp0]
	gatk_gg = subprocess.Popen(command)
	gatk_gg.wait()

	# '''
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
	# '''

def renalyze_from_bam(samples):
	for sample in samples:
		bam = glob.glob(sample + '*bam')
		fq1 = sample + '.r1.fastq'
		fq2 = sample + '.r2.fastq'
		print(sample, bam, fq1, fq2)
		convert_bam_fastq_picard(bam[0], fq1, fq2)
		align_with_bwa_one_at_time(sample, fq1, fq2)

def renalyze_from_fastq(samples):
	for sample in samples:
		fq1_files = sorted(glob.glob(sample + '*R1_001.fastq.gz'))
		fq2_files = sorted(glob.glob(sample + '*R2_001.fastq.gz'))
		fq1 = sample + '.r1.fastq.gz'
		fq2 = sample + '.r2.fastq.gz'
		print(fq1_files, fq2_files, fq1, fq2)
		combine_fq_file(fq1_files, fq2_files, fq1, fq2)
		align_with_bwa_one_at_time(sample, fq1, fq2)

##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/daniela_exomes_recall_0720'
os.chdir(working_dir)
project = 'microtia_exomes_recall_0720'

##run in three batches
b1_samples = ['101000201', '101000202', '101000203', '101000301', '101000302', '101000401', '101000402', '101000403', '101000501', '101000502', 
		'101000503', '101000601', '101000602', '101000701', '101000702', '101000703', '101000801', '101000802', '101000803', '101001001', 
		'101001002', '101001003', '101001301', '101001302', '101001303', '101001901', '101001902', '101001903', '101002001', '101002002', 
		'101002003', '101002101', '101002102', '101002103', '101002201', '101002202', '101002203', '101002301', '101002302', '101002303', 
		'101002401', '101002402', '101002403', '101002501', '101002502', '101002503', '101002601', '101002602', '101002603', '101002801', 
		'101002802', '101002803', '101003101', '101003102', '101003103', '101003201', '101003202', '101003203']
gvcf1 = project + '.1.g.vcf'
# variant_calling_gatk_hc(b1_samples)
# combine_gvcfs(b1_samples, gvcf1)

b2_samples = ['103000101', '103000102', '103000103', '103000201', '103000202', '103000203', '103000601', '103000602', '103000603', '103000701', 
		'103000702', '103001101', '103001102', '103001103', '103001301', '103001302', '103001303', '103001401', '103001402', '103001403', 
		'103001501', '103001502', '103001503', '106000401', '106000402', '106000403', '106001201', '106001401', '106001501', '106001502', 
		'106001503', '106001801', '106001805', '106001806a', '106001806b', '106001809', '106002001', '106002002', '106002003', '106002301', 
		'106002302', '106002303', '106002901', '106002902', '106002903', '106002907', '107000101', '107000102', '107000103', '107000201', 
		'107000202', '107000203', '107000301', '107000302', '107000303', '107001101', '107001102', '107001103']
gvcf2 = project + '.2.g.vcf'
# variant_calling_gatk_hc(b2_samples)
# combine_gvcfs(b2_samples, gvcf2)

b3_samples = ['107001201', '107001202', '107001203', '107001301', '107001302', '107001303', '107001501', '107001502', '107001503', '401000201', 
		'401000202', '401000203', 'CFMCLP-01-JV', 'CFMCLP-02-AS', 'CFMCLP-02-AS-MS', 'CFMCLP-03-AA', 'CFMCLP-03-AA-CA', 'CFMCLP-03-AA-JA', 
		'CFMCLP-04-AG', 'CFMCLP-04-AG-HG', 'CFMCLP-04-AG-HGS', 'CFMCLP-04-AG-SCC', 'CFM-MOS-01-01-B', 'CFM-MOS-01-01-T', 'CFM-MOS-01-11-B', 
		'CFM-MOS-01-21-B', 'CFM-MOS-03-01', 'CFM-MOS-03-02', 'CFM-MOS-03-04', 'CFM-MOS-03-11', 'CFM-MOS-03-21', 'CFM-MOS-06-01-B', 
		'CFM-MOS-06-01-T', 'CFM-MOS-06-11-B', 'CFM-MOS-06-21-B', 'CFM-MOS-09-01-B', 'CFM-MOS-09-01-T', 'CFM-MOS-13-01-B', 'CFM-MOS-13-01-T', 
		'CFM-MOS-13-11-B', 'CFM-MOS-13-21-B', 'CFM-MOS-14-01-B', 'CFM-MOS-14-01-T', 'CFM-MOS-14-11-B', 'CFM-MOS-14-21-B', 'CFM-MOS-15-01-B', 
		'CFM-MOS-15-01-T', 'CFM-MOS-15-11-B', 'CFM-MOS-15-21-B', 'CFM-MOS-16-01-B', 'CFM-MOS-16-01-T', 'CFM-MOS-18-01-B', 'CFM-MOS-18-01-T']
gvcf3 = project + '.3.g.vcf'
# variant_calling_gatk_hc(b3_samples)
# combine_gvcfs(b3_samples, gvcf3)

##issue with indexing some of the bams so no gvcfs created, bams corrupted, so need to remake then transfer to make combined vcf and rss
sample_rpt_from_bams = ['101002101', '101002103', '103000201', '103000603', '103000701', '106000401', '107000101', '107000103', '107000303', 
		'107001101', '107001102', '107001103']
sample_rpt_from_fq = ['101003103', '101003202', '103001301', '103001302', '103001303', '107001203', '107001301', '107001501', '107001503']
# working_dir = '/home/atimms/ngs_data/exomes/working/daniela_exomes_recall_0720/rpts'
# os.chdir(working_dir)
##redo commands
# renalyze_from_bam(sample_rpt_from_bams)
# renalyze_from_fastq(sample_rpt_from_fq)
# variant_calling_gatk_hc(sample_rpt_from_bams)
# variant_calling_gatk_hc(sample_rpt_from_fq)


##combine batches
# combine_gvcfs(b1_samples, gvcf1)
# combine_gvcfs(b2_samples, gvcf2)
# combine_gvcfs(b3_samples, gvcf3)

##genotype gvcf
gvcfs_to_genotype = [gvcf1, gvcf2, gvcf3]
genotype_gvcfs_and_filter(gvcfs_to_genotype, project)


