#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys


'''
requires:
##without pisces
qsub -Iq cdbrmq -l mem=50gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
#java 1.8 for gatk
module load java/1.8.0_202 


'''

##set input variables and parameters
delim = '\t'
threads = '5'

##programs
bgzip = '/home/atimms/programs/samtools-1.11/bin/bgzip'
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
slivar = '/home/atimms/programs/slivar_0921/slivar'
# slivar = '/home/atimms/programs/slivar_0121/slivar_dev'
# pslivar = '/home/atimms/programs/slivar_0121/pslivar'
slivar_functions = '/home/atimms/programs/slivar_0921/slivar-0.2.5/js/slivar-functions_at.js'
# slivar_functions = '/home/atimms/programs/slivar_0921/slivar-0.2.5/js/slivar-functions.js'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
freebayes = '/home/atimms/programs/freebayes-1.3.4-linux-static-AMD64'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
tabix = '/home/atimms/programs/samtools-1.11/bin/tabix'
# plink = '/home/atimms/programs/plink'
somalier = '/home/atimms/programs/somalier_0721/somalier'

##files
fasta = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens_assembly38.fasta'
exome_capture_bed = '/home/atimms/ngs_data/references/exome_beds_hg38/hg38_targets_combined_padded_0721.bed'
##sorted specifically for this reference
hg38_refseq_exons = '/home/atimms/ngs_data/references/hg38/hg38_RefSeq_exons.sorted.bed'
ens_gff = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens.GRCh38.103.chr_added.gff3.gz'
gnomad_gnotate = '/home/atimms/ngs_data/references/slivar/gnomad.hg38.genomes.v3.fix.zip'
topmed_gnotate = '/home/atimms/ngs_data/references/slivar/topmed.hg38.dbsnp.151.zip'
#wget -qO - https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | zcat | cut -f 1,21,24 | tail -n+2 | awk '{ printf("%s\tpLI=%.3g;oe_lof=%.5g\n", $1, $2, $3)}' > pli.lookup
pli_lookup = '/home/atimms/ngs_data/references/slivar/pli.lookup'
#grep -v '^#' omim_genemap2_0321.txt | cut -f9,13 | grep -v '^\s' > /home/atimms/ngs_data/references/slivar/omim.lookup
omim_lookup = '/home/atimms/ngs_data/references/slivar/omim.lookup'
#cut -f2,3,4,6,8 dbdb.gene_association.txt | awk 'BEGIN{FS="\t"}{ printf("%s\tinheritance=%s;phenotype=%s;syndrome=%s;loe=%s\n", $1, $2, $3, $4, $5)}' > ~/ngs_data/references/slivar/dbdb.lookup
dbdb_lookup = '/home/atimms/ngs_data/references/slivar/dbdb.lookup'
hg38_selfchain = '/home/atimms/ngs_data/references/slivar/selfchain-LCR.hg38.bed'
hg38_refseq_bed = '/home/atimms/ngs_data/references/slivar/hg38.refseq.bed'
decipher_delevopmental_file = '/home/atimms/ngs_data/references/slivar/DDG2P_11_2_2021.csv'
dbsnp_common_vcf = '/home/atimms/ngs_data/references/hg38_gatk/common_all_20180418.chr_added.vcf.gz'
somalier_sites_vcf = '/home/atimms/ngs_data/references/somalier/sites.hg38.vcf.gz'

##ped types
trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
duo_types = ['duo', 'duo*', 'parent_sibship']
single_types = ['singleton', 'singleton*', 'sibship']
duo_single_types = duo_types + single_types
multiplex_types = ['multiplex']

##file sufffixes
gatk_vcf_suffix = '.gatk38hc.vcf.gz'
fb_vcf_suffix = '.freebayes13.vcf.gz'
int_vcf_dir = '.intersect_vcfs'
int_vcf_suffix = '.intersect_vcfs/0002.vcf'
gatk_bcftools_vcf_suffix = '.gatk38hc.bcftools.GRCh38_103.vcf.gz'
int_bcftools_vcf_suffix = '.intersect.bcftools.GRCh38_103.vcf.gz'

##annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,knownGene,ensGene,dbnsfp41a,clinvar_20210123,gnomad30_genome']
av_operation = ['-operation', 'g,g,g,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput', '-arg', '-splicing 10 ,,,,,' ]


##methods

##make list of all bam files to be analyzed
def make_list_of_bams(samples, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in samples:
			bam = 'Preprocessing/' + sample + '/Recalibrated/' + sample + '.recal.bam'
			outfh.write(bam + '\n')

def variant_calling_gatk_hc_no_gvcf(bamlist, name_prefix):
	vcf_temp0 = name_prefix + 'temp_gatk0.vcf'
	vcf_raw_snps = name_prefix + 'temp_raw_snps.vcf'
	vcf_filtered_snps = name_prefix + 'temp_filtered_snps.vcf'
	vcf_raw_indels = name_prefix + 'temp_raw_indels.vcf'
	vcf_filtered_indels = name_prefix + 'temp_filtered_indels.vcf'
	vcf_temp1 = name_prefix + 'temp_gatk1.vcf'
	final_vcf = name_prefix + gatk_vcf_suffix
	##run haplotype caller
	gatk_hc = subprocess.Popen(['java', '-Xmx50g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', bamlist,'-L', exome_capture_bed, '-o', vcf_temp0])
	gatk_hc.wait()
	##split data in SNPs and indels and apply manual variant filtering (add -L)
	snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-V', vcf_temp0, '-o', vcf_raw_snps, '-selectType', 'SNP'])
	snp_cut.wait()
	snp_vf = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_snps, '-o', vcf_filtered_snps, '--filterExpression ', "QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "snp_filter"])
	snp_vf.wait()
	indel_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-V', vcf_temp0, '-o', vcf_raw_indels, '-selectType', 'INDEL'])
	indel_cut.wait()
	indel_vf = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_indels, '-o', vcf_filtered_indels, '--filterExpression ', "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0", '--filterName', "indel_filter"])
	indel_vf.wait()
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', vcf_temp1, '-genotypeMergeOptions', 'UNSORTED'])
	combine_var.wait()
	bgzip_cmd = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_cmd.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	##decomposing and subsetting vcf - AD looks good with gatk38, so sed not needed
	# zless_vcf = subprocess.Popen(['zless', vcf_temp1 + '.gz'], stdout=subprocess.PIPE)
	# sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
	#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
	# bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', final_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
	bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', final_vcf, '-O', 'z', vcf_temp1 + '.gz'])
	bcftools_norm.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()

def make_sample_dict(input_file):
	ped_dict = {}
	with open(input_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				sample_name = line[0]
				diagnosis = line[1]
				sex = line[2]
				ped_dict[sample_name] = [diagnosis, sex]
	##retrun analysis info
	return(ped_dict)


##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/kim_exomes_0621/'
os.chdir(working_dir)

exome_info_file = 'kim_exome_analysis_0921.txt'
project_name = exome_info_file.split('.')[0]
bams_list = project_name + '.bams.list'
##make dict of samples
sample_dict = make_sample_dict(exome_info_file)
##make list of bams
make_list_of_bams(sample_dict, bams_list)
##gatk hc
variant_calling_gatk_hc_no_gvcf(bams_list, project_name)

