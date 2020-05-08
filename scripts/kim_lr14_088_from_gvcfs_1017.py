#!/usr/bin/python
import sys
import subprocess
import os
import glob
import dobyns_gemini_pipeline_cybertron_v6

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/genomes/kim_lr14-088_1017/working'
os.chdir(working_dir)

##file names etc
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'genome.fa'
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
gatk = '/cm/shared/apps/GATK/GenomeAnalysisTK.jar'
bcftools = 'bcftools'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
bgzip = 'bgzip'

plink = '/home/atimms/programs/plink'



def variant_calling_gatk_hc(samples, name_prefix):
	vcf_temp0 = 'temp_gatk0.vcf'
	vcf_raw_snps = 'temp_raw_snps.vcf'
	vcf_filtered_snps = 'temp_filtered_snps.vcf'
	vcf_raw_indels = 'temp_raw_indels.vcf'
	vcf_filtered_indels = 'temp_filtered_indels.vcf'
	vcf_temp1 = 'temp_gatk1.vcf'
	vcf_temp2 = 'temp_gatk2.vcf.gz'
	final_vcf = name_prefix + '.gatkHC.vcf.gz'
	##get g.vcfs for each sample (add -L)
	gvcf_files = []
	for sample in samples:
		gvcf = ['-V', sample + '.genome.vcf']
		gvcf_files.extend(gvcf)
	##genotype g.vcfs(add -L)
	command = ['java', '-Xmx20g', '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', fasta, '-nt', '15'] + gvcf_files + ['-o', vcf_temp0]
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


##run methods
ped_name = 'LR14-088'
gvcfs = ['PG0006685-BLD', 'PG0006686-BLD', 'PG0006688-BLD']
variant_calling_gatk_hc(gvcfs, ped_name)





