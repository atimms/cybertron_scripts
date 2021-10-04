#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys
import filtering_annotated


'''
requires:
##without pisces
qsub -Iq cdbrmq -l mem=50gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
#java 1.8 for gatk
module load java/1.8.0_202 

##for homozygosity_mapping_cybertron i need python2 so just use:
module load biobuilds/2017.11
module load gcc/6.1.0 ##gcc version used to make bedtools using in homo mapping

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
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'
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
gatk_filtered_vcf_suffix = '.gatk38hc.std_chrs.vcf.gz'
fb_vcf_suffix = '.freebayes13.vcf.gz'
int_vcf_dir = '.intersect_vcfs'
int_vcf_suffix = '.intersect_vcfs/0002.vcf'
gatk_bcftools_vcf_suffix = '.gatk38hc.bcftools.GRCh38_103.vcf.gz'
int_bcftools_vcf_suffix = '.intersect.bcftools.GRCh38_103.vcf.gz'

##annovar parameters... update generic avinputs, order cases and controls
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,knownGene,ensGene,dbnsfp41a,clinvar_20210123,gnomad30_genome,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic', '-genericdbfile', 'B03_12_1.avinput,B05_13_1.avinput,B40_17_1.avinput,N04_15_1.avinput,NBB_13646_1.avinput,NBB_5294_1.avinput,NBB_5419_1.avinput,NBB_5531_1.avinput,NBB_5565_1.avinput,NBB_5841_1.avinput,NBB_5864_1.avinput,NBB_5978_1.avinput,NBB_6033_1.avinput,NBB4999.avinput,NBB5454.avinput,NBB5574.avinput,NBB6041.avinput,B28_17_1.avinput,B36_17_1.avinput,M23_17_1.avinput,NBB_1347_1.avinput,NBB_1714_1.avinput,NBB_4599_1.avinput,NBB_4787_1.avinput,NBB_4916_1.avinput,NBB_5242_1.avinput,NBB_5334_1.avinput,NBB_5936_1.avinput,NBB1862.avinput,NBB5387.avinput,NBB5813.avinput,NBB5893.avinput,U30_17_1.avinput,U38_18_1.avinput' ]
av_operation = ['-operation', 'g,g,g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput', '-arg', '-splicing 10 ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' ]

##filtering
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = ['synonymous SNV', 'unknown']

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


def filter_vcf_non_std_chr(in_vcf, out_vcf):
	vcf_temp1 = 'temp1.vcf.gz'
	vcf_temp2 = 'temp2.vcf.gz'
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-r', 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY', '-o', vcf_temp1, '-O', 'z', in_vcf])
	bcftools_filter.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', out_vcf, vcf_temp2])
	bcf_norm2.wait()

##make avinput files
def convert_to_annovar_move_to_annovar_folder(samples, vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
	con_ann.wait()
	for sample in samples:
		temp_av_file = 'temp.' + sample + '.avinput'
		av_file = sample + '.avinput'
		shutil.move(temp_av_file, av_file)
		shutil.copy(av_file, str(av_ref_dir[0]))

def run_table_annovar(vcf, av_prefix):
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def multianno_to_annotated(av_prefix): 
	extra_header = ['B03_12_1', 'B05_13_1', 'B40_17_1', 'N04_15_1', 'NBB_13646_1', 'NBB_5294_1', 'NBB_5419_1', 'NBB_5531_1', 'NBB_5565_1', 'NBB_5841_1', 'NBB_5864_1', 'NBB_5978_1', 'NBB_6033_1', 'NBB4999', 'NBB5454', 'NBB5574', 'NBB6041', 'B28_17_1', 'B36_17_1', 'M23_17_1', 'NBB_1347_1', 'NBB_1714_1', 'NBB_4599_1', 'NBB_4787_1', 'NBB_4916_1', 'NBB_5242_1', 'NBB_5334_1', 'NBB_5936_1', 'NBB1862', 'NBB5387', 'NBB5813', 'NBB5893', 'U30_17_1', 'U38_18_1', 'chr', 'pos', '', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'B03_12_1_vcf','B05_13_1_vcf','B28_17_1_vcf','B36_17_1_vcf','B40_17_1_vcf','M23_17_1_vcf','N04_15_1_vcf','NBB1862_vcf','NBB4999_vcf','NBB5387_vcf','NBB5454_vcf','NBB5574_vcf','NBB5813_vcf','NBB5893_vcf','NBB6041_vcf','NBB_1347_1_vcf','NBB_13646_1_vcf','NBB_1714_1_vcf','NBB_4599_1_vcf','NBB_4787_1_vcf','NBB_4916_1_vcf','NBB_5242_1_vcf','NBB_5294_1_vcf','NBB_5334_1_vcf','NBB_5419_1_vcf','NBB_5531_1_vcf','NBB_5565_1_vcf','NBB_5841_1_vcf','NBB_5864_1_vcf','NBB_5936_1_vcf','NBB_5978_1_vcf','NBB_6033_1_vcf','U30_17_1_vcf','U38_18_1_vcf']
	multianno = av_prefix + '.hg38_multianno.txt'
	annotated = av_prefix + '.annotated.txt'
	with open(multianno, "r") as multi, open(annotated, "w") as final:
		line_count = 0
		for line in multi:
			line = line.rstrip().split(delim)

			line_count += 1
			if line_count == 1:
				header = line[:89] + extra_header
				final.write(delim.join(header) + '\n')
			else:
				line_out = line[:123] + line[126:]
				final.write(delim.join(line_out) + '\n')


def filter_ann_file_exonic(file_prefix):
	##exonic_variants
	filtering_annotated.filter(working_dir, "or", file_prefix + '.annotated.txt' , file_prefix + "_1.temp", [col_exon, col_exon], ['==','=='], exon_definition)
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", file_prefix + "_1.temp", file_prefix + "_2.temp", [col_function,col_function], ['!=','!='], syn_definition)
	##passed fi;ter
	filtering_annotated.filter(working_dir, "and", file_prefix + "_2.temp", file_prefix + "_3.temp", [130], ['=='], ['PASS'])
	##not in control
	filtering_annotated.filter(working_dir, "and", file_prefix + "_3.temp", file_prefix + "_4.temp", [107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123], ['==','==','==','==','==','==','==','==','==','==','==','==','==','==','==','==','=='], ['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.'])
	##<0.001 in all gnomad
	filtering_annotated.filter(working_dir, "and", file_prefix + "_4.temp", file_prefix + ".exonic.rare.not_in_ctl.xls", [77], ['<='], [0.001])


##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/kim_exomes_0621/'
os.chdir(working_dir)

exome_info_file = 'kim_exome_analysis_0921.txt'
project_name = exome_info_file.split('.')[0]
bams_list = project_name + '.bams.list'
result_vcf = project_name + gatk_vcf_suffix
filtered_vcf= project_name + gatk_filtered_vcf_suffix

##make dict of samples
sample_dict = make_sample_dict(exome_info_file)

##make list of bams
# make_list_of_bams(sample_dict, bams_list)

##gatk hc
# variant_calling_gatk_hc_no_gvcf(bams_list, project_name)

##annotate and format
# filter_vcf_non_std_chr(result_vcf, filtered_vcf)
# convert_to_annovar_move_to_annovar_folder(sample_dict, filtered_vcf)
# run_table_annovar(filtered_vcf, project_name)
multianno_to_annotated(project_name)

##filter rare (absent or <0.001), exonic (not synonympuse), and not in any controls
filter_ann_file_exonic(project_name)




