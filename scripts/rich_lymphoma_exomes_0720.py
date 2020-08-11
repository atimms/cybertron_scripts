#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil
import filtering_annotated

##info
'''
##convert_bam_fastq_picard
module load biobuilds

##xengsort - already made index
qsub -Iq longq -l mem=200gb,ncpus=20
module load local_python/3.7.6
source activate xengsort

##gatk23 pipeline
https://gatk.broadinstitute.org/hc/en-us/articles/360037498992--How-to-Map-reads-to-a-reference-with-alternate-contigs-like-GRCH38
https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently

'''


##parameters
delim = '\t'

##dirs and file
xengsort_index = '/home/atimms/ngs_data/references/xengsort/mouse_human_ref_0720.h5'
##for gatk etc
ref_dir = '/home/atimms/ngs_data/references/hg38_gatk/'
fasta = ref_dir + 'Homo_sapiens_assembly38.fasta'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
indels_1000g = ref_dir + '1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf'
dbsnp = ref_dir + 'Homo_sapiens_assembly38.dbsnp138.vcf'
exome_capture_bed = ref_dir + 'Twist_Exome_Target_hg38_100bp_padding.bed'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
exome_bed_for_coverage = ref_dir + 'Twist_Exome_Target_hg38_sorted.bed' 
mutect_gnomad_vcf = '/home/atimms/ngs_data/references/hg38/af-only-gnomad.hg38.vcf.gz'

##programs
bwa = 'bwa'
samtools = 'samtools'
picard = 'picard'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
bedtools = 'bedtools'
bgzip = 'bgzip'
bcftools = 'bcftools'
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
gatk4 = '/home/atimms/programs/gatk-4.1.3.0/gatk'

##annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,dbnsfp35a,rmsk,genomicSuperDups,gnomad211_genome,gnomad211_exome,gnomad30_genome,clinvar_20200316,cosmic90_noncoding,cosmic90_coding,generic,generic,generic,generic,generic,generic,generic,generic,generic', '-genericdbfile', 'TZ001.avinput,TZ002.avinput,TZ003.avinput,TZ004.avinput,TZ005.avinput,TZ006.avinput,TZ007.avinput,TZ008.avinput,TZ009.avinput']
av_operation = ['-operation', 'g,f,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.','-vcfinput']
av_protocol_mt2 = ['-protocol', 'refGene,dbnsfp35a,rmsk,genomicSuperDups,gnomad211_genome,gnomad211_exome,gnomad30_genome,clinvar_20200316,cosmic90_noncoding,cosmic90_coding']
av_operation_mt2 = ['-operation', 'g,f,r,r,f,f,f,f,f,f']


##filtering parameters
col_exon = 15
exon_definition = ['exonic', 'splicing']
col_function = 18
syn_definition = 'synonymous_SNV'
af_cols = [92,109,126]
freq_req = 0.1
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
	gatk_bam = sample + final_bam_suffix
	'''
	##bwa and convert to bam
	bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
	st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
	st_sam_bam_pe.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam, '-@', '10', '-m', '10G', post_bwa_bam])
	st_sort_pe.wait()
	'''
	##mark duplicates
	# picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md = subprocess.Popen([picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md.wait()
	
	##bqsr
	gatk_br = subprocess.Popen(['java', '-Djava.io.tmpdir=.', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', mkdup_bam, '-L', exome_capture_bed, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Djava.io.tmpdir=.', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', mkdup_bam, '-L', exome_capture_bed, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
	gatk_pr.wait()

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



##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

##delete files supplied in list
def delete_unwated_files(file_extensions):
	files_to_delete = []
	for ext in file_extensions:
		files = glob.glob(ext)
		files_to_delete.extend(files)
	for f in files_to_delete:
		os.remove(f)

def variant_calling_gatk_hc_no_gvcf(bam_list, name_prefix):
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
	gatk_hc = subprocess.Popen(['java', '-Djava.io.tmpdir=.', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', bam_list,'-L', exome_capture_bed, '-o', vcf_temp0])
	gatk_hc.wait()
	##split data in SNPs and indels and apply manual variant filtering (add -L)
	snp_cut = subprocess.Popen(['java', '-Djava.io.tmpdir=.', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_snps, '-selectType', 'SNP'])
	snp_cut.wait()
	snp_vf = subprocess.Popen(['java', '-Djava.io.tmpdir=.', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_snps, '-o', vcf_filtered_snps, '--filterExpression ', "QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "snp_filter"])
	snp_vf.wait()
	indel_cut = subprocess.Popen(['java', '-Djava.io.tmpdir=.', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_indels, '-selectType', 'INDEL'])
	indel_cut.wait()
	indel_vf = subprocess.Popen(['java', '-Djava.io.tmpdir=.', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_indels, '-o', vcf_filtered_indels, '--filterExpression ', "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0", '--filterName', "indel_filter"])
	indel_vf.wait()
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Djava.io.tmpdir=.', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '5', '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', vcf_temp1, '-genotypeMergeOptions', 'UNSORTED'])
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

def split_info_field(info_list):
	indices = [i for i, s in enumerate(info_list) if 'ANNOVAR_DATE' in s]
	# print indices
	i_count = 0
	final_list = []
	for info in info_list:
		# print info
		if i_count > indices[0] and info != 'ALLELE_END':
			info2 = info.split('=')[1]
			#print info2
			final_list.append(info2)
		i_count += 1
	return final_list

def annotate_vars_from_vcf(out_prefix, samples):
	'''
	norm_vcf = out_prefix + '.gatkHC.vcf.gz'
	##make av files for all sample and copy to annovar dir
	# con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', norm_vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', '0620'])
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', norm_vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
	con_ann.wait()
	for sample in samples:
		temp_av_file = 'temp.' + sample + '.avinput'
		av_file = sample + '.avinput'
		shutil.move(temp_av_file, av_file)
		shutil.copy(av_file, str(av_ref_dir[0]))

	##annotate vcf file
	command = [table_annovar] + av_buildver + [norm_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()
	##split into individual samples
	post_annovar_vcf = out_prefix + '.' + av_genome + '_multianno.vcf'
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', post_annovar_vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', out_prefix])
	con_ann.wait()
	'''
	##format files...
	for sample in samples:
		avinput = out_prefix + '.' + sample + '.avinput'
		outfile = out_prefix + '.' + sample + '.annotated.txt'
		head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Zygosity', 'Qual', 'Coverage', 'Filter', 'Pos', 'Ref2', 'Alt2', 'Format', 'Info', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'SIFT_score', 'SIFT_converted_rankscore', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_rankscore', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_converted_rankscore', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_converted_rankscore', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_score_rankscore', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_converted_rankscore', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_converted_rankscore', 'PROVEAN_pred', 'VEST3_score', 'VEST3_rankscore', 'MetaSVM_score', 'MetaSVM_rankscore', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_rankscore', 'MetaLR_pred', 'M-CAP_score', 'M-CAP_rankscore', 'M-CAP_pred', 'REVEL_score', 'REVEL_rankscore', 'MutPred_score', 'MutPred_rankscore', 'CADD_raw', 'CADD_raw_rankscore', 'CADD_phred', 'DANN_score', 'DANN_rankscore', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_rankscore', 'fathmm-MKL_coding_pred', 'Eigen_coding_or_noncoding', 'Eigen-raw', 'Eigen-PC-raw', 'GenoCanyon_score', 'GenoCanyon_score_rankscore', 'integrated_fitCons_score', 'integrated_fitCons_score_rankscore', 'integrated_confidence_value', 'GERP++_RS', 'GERP++_RS_rankscore', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP20way_mammalian', 'phyloP20way_mammalian_rankscore', 'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons20way_mammalian', 'phastCons20way_mammalian_rankscore', 'SiPhy_29way_logOdds', 'SiPhy_29way_logOdds_rankscore', 'Interpro_domain', 'GTEx_V6p_gene', 'GTEx_V6p_tissue', 'rmsk', 'genomicSuperDups', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'AF', 'AF_raw', 'AF_male', 'AF_female', 'AF_afr', 'AF_ami', 'AF_amr', 'AF_asj', 'AF_eas', 'AF_fin', 'AF_nfe', 'AF_oth', 'AF_sas', 'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG', 'cosmic90_noncoding', 'cosmic90_coding', 'TZ001', 'TZ002', 'TZ003', 'TZ004', 'TZ005', 'TZ006', 'TZ007', 'TZ008', 'TZ009']
		head_out = delim.join(head + ['\n'])
		with open(avinput, "r") as av, open(outfile, "w") as final:
			final.write(head_out)
			for line in av:
				line = line.strip('\n').split(delim)
				stuff = line[0:8] + [line[14]] + [line[9]] + line[11:13]
				info = line[15].split(';')
				info_list = split_info_field(info)
				other_stuff = line[16:]
				line_out = delim.join(stuff + other_stuff + info_list +['\n'])
				final.write(line_out)


def filter_ann_file_somatic_only(file_prefix, sample):
	file_prefix = file_prefix + '.' + sample
	##remove if in control, i.e. not in TZ008
	if sample == 'TZ001' or sample == 'TZ002':
		filtering_annotated.filter(working_dir, "and", file_prefix + '.annotated.txt', file_prefix + ".not_in_TZ008.xls", [153], ['=='], ['.'])
	##remove if in control, i.e. not in TZ009
	elif sample == 'TZ003':
		filtering_annotated.filter(working_dir, "and", file_prefix + '.annotated.txt', file_prefix + ".not_in_TZ009.xls", [154], ['=='], ['.'])
	##remove if in control, i.e. not in TZ007
	elif sample == 'TZ004' or sample == 'TZ005' or sample == 'TZ006':
		filtering_annotated.filter(working_dir, "and", file_prefix + '.annotated.txt', file_prefix + ".not_in_TZ007.xls", [152], ['=='], ['.'])
	else:
		print(sample, 'sample name not recognized')




def filter_ann_file_causal(file_prefix, sample):
	file_prefix = file_prefix + '.' + sample
	##if not passed by gatk
	filtering_annotated.filter(working_dir, "and", file_prefix + '.annotated.txt', file_prefix + "_0.temp", [9], ['=='], ['PASS'])

	##remove if in control, i.e. not in TZ008
	if sample == 'TZ001' or sample == 'TZ002':
		filtering_annotated.filter(working_dir, "and", file_prefix + "_0.temp", file_prefix + ".not_in_TZ008.xls", [153], ['=='], ['.'])
		##exonic_variants
		filtering_annotated.filter(working_dir, "or", file_prefix + ".not_in_TZ008.xls" , file_prefix + "_1.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
		##remove synonymous
		filtering_annotated.filter(working_dir, "and", file_prefix + "_1.temp", file_prefix + "_2.temp", [col_function], ['!='], [syn_definition])
		##<10% in all gnomad
		filtering_annotated.filter(working_dir, "and", file_prefix + "_2.temp", file_prefix + ".not_in_TZ008.exonic.rare.xls", af_cols, ['<=','<=','<='], [freq_req,freq_req,freq_req])

	##remove if in control, i.e. not in TZ009
	elif sample == 'TZ003':
		filtering_annotated.filter(working_dir, "and", file_prefix + "_0.temp", file_prefix + ".not_in_TZ009.xls", [154], ['=='], ['.'])
		##exonic_variants
		filtering_annotated.filter(working_dir, "or", file_prefix + ".not_in_TZ009.xls" , file_prefix + "_1.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
		##remove synonymous
		filtering_annotated.filter(working_dir, "and", file_prefix + "_1.temp", file_prefix + "_2.temp", [col_function], ['!='], [syn_definition])
		##<10% in all gnomad
		filtering_annotated.filter(working_dir, "and", file_prefix + "_2.temp", file_prefix + ".not_in_TZ009.exonic.rare.xls", af_cols, ['<=','<=','<='], [freq_req,freq_req,freq_req])

	##remove if in control, i.e. not in TZ007
	elif sample == 'TZ004' or sample == 'TZ005' or sample == 'TZ006':
		filtering_annotated.filter(working_dir, "and", file_prefix + "_0.temp", file_prefix + ".not_in_TZ007.xls", [152], ['=='], ['.'])
		##exonic_variants
		filtering_annotated.filter(working_dir, "or", file_prefix + ".not_in_TZ007.xls" , file_prefix + "_1.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
		##remove synonymous
		filtering_annotated.filter(working_dir, "and", file_prefix + "_1.temp", file_prefix + "_2.temp", [col_function], ['!='], [syn_definition])
		##<10% in all gnomad
		filtering_annotated.filter(working_dir, "and", file_prefix + "_2.temp", file_prefix + ".not_in_TZ007.exonic.rare.xls", af_cols, ['<=','<=','<='], [freq_req,freq_req,freq_req])

	else:
		print(sample, 'sample name not recognized')
		# ##exonic_variants
		# filtering_annotated.filter(working_dir, "or", sample + ".not_in_TZ008.xls" , sample + "_1.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
		# ##remove synonymous
		# filtering_annotated.filter(working_dir, "and", sample + "_1.temp", sample + "_2.temp", [col_function], ['!='], [syn_definition])


def exome_analysis_pipeline(s_dict, project):
	##align reads
	for sample in s_dict:
		read1_fastq = outdir = './' + sample + '_xe/' + sample + '-graft.1.fq'
		read2_fastq = outdir = './' + sample + '_xe/' + sample + '-graft.2.fq'
		# align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	##call using gatk 
	# make_list_of_bams(s_dict, final_bam_suffix, bamlist)
	# variant_calling_gatk_hc_no_gvcf(bamlist, project)
	##annotate using annovar
	# annotate_vars_from_vcf(project, s_dict)
	##filter variants
	for sample in s_dict:
		print(sample)
		# filter_ann_file_somatic_only(project, sample)
		# filter_ann_file_causal(project, sample)
	##get exome coverage
	calculate_exome_coverage(s_dict, final_bam_suffix, exome_bed_for_coverage, project)
	##delete file
	# delete_unwated_files(['*.bwa.bam', '*.bwa_sort.bam', '*.metrics', '*.bwa_religned.bam', '*.bwa_religned.bai', '*.recal_data.table', '*temp*', '*bwa_mkdup.bai', '*bwa_mkdup.bam'])

def call_mutect2_sample(tumor_name, normal_name, bam_suffix, working_dir):
	##issues, no PONs, and germline resources, ? fasta 
	sample_names = [tumor_name, normal_name]
	bams = []
	for s in sample_names:
		bamfile = s + bam_suffix
		ibam = ['-I', bamfile]
		bams.extend(ibam)
	normal = ['-normal', normal_name]

	print(bams, normal)
	##file names
	raw_vcf = tumor_name + '.mt2.temp_raw.vcf.gz'
	filtered_vcf = tumor_name + '.mt2.temp_filtered.vcf.gz'
	vcf_temp1 = tumor_name + '.mt2.temp1.vcf.gz'
	vcf_temp2 = tumor_name + '.mt2.temp2.vcf.gz'
	final_vcf = tumor_name + '.mt2.final.vcf.gz'
	# '''
	##run mt2, use --pcr-indel-model NONE if pcr free library prep and --f1r2-tar-gz for filtering by read bias
	# mt2_cmd = [gatk4, 'Mutect2', '-R', fasta] + bams + normal + ['--germline-resource', mutect_gnomad_vcf, '--tmp-dir', working_dir, '--f1r2-tar-gz', out_prefix + 'f1r2.tar.gz', '-O', raw_vcf]
	mt2_cmd = [gatk4, 'Mutect2', '-R', fasta] + bams + normal + ['--germline-resource', mutect_gnomad_vcf, '--tmp-dir', working_dir, '-O', raw_vcf]
	run_mt2 = subprocess.Popen(mt2_cmd)
	run_mt2.wait()
	##filter calls
	##gatk FilterMutectCalls -R ref.fasta -V unfiltered.vcf -O filtered.vcf
	filter_mt2 = subprocess.Popen([gatk4, 'FilterMutectCalls', '-R', fasta, '-V', raw_vcf, '-O', filtered_vcf])
	filter_mt2.wait()
	##get filtered -- don't use this
	# bcftools_filter = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-o', vcf_temp1, '-O', 'z', filtered_vcf])
	# bcftools_filter.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, filtered_vcf])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()
	# '''

def annotate_mt2_vars_from_vcf(out_prefix, norm_vcf):
	##annotate vcf file
	# command = [table_annovar] + av_buildver + [norm_vcf] + av_ref_dir + av_protocol_mt2 + av_operation_mt2 + av_options + ['-out', out_prefix]
	# annovar = subprocess.Popen(command)
	# annovar.wait()
	##format files...
	multi = out_prefix + '.hg38_multianno.txt'
	outfile = out_prefix + '.mt2.annotated.txt'
	extra_head = ['Filter', 'Info', 'Format', 'tumor', 'normal']
	with open(multi, "r") as in_fh, open(outfile, "w") as out_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc == 1:
				line_out = line[:136] + extra_head
			else:
				line_out = line[:136] + line[145:]
			out_fh.write(delim.join(line_out) + '\n')

def filter_ann_file_mt2(file_prefix):
	##if not passed by mt2
	filtering_annotated.filter(working_dir, "and", file_prefix + '.mt2.annotated.txt', file_prefix + '.mt2.passed.xls', [137], ['=='], ['PASS'])

	##exonic_variants
	filtering_annotated.filter(working_dir, "or", file_prefix + '.mt2.annotated.txt' , file_prefix + "_1.temp", [6, 6], ['==','=='], [exon_definition[0],exon_definition[1]])
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", file_prefix + "_1.temp", file_prefix + "_2.temp", [9], ['!='], ['synonymous SNV'])
	##<10% in all gnomad
	filtering_annotated.filter(working_dir, "and", file_prefix + "_2.temp", file_prefix + ".mt2.exonic.rare.xls", [83,100,117], ['<=','<=','<='], [freq_req,freq_req,freq_req])


def mutect2_pipeline(s_dict, project, bamfile_suffix, work_dir):
	for control in s_dict:
		for sample in s_dict[control]:
			# call_mutect2_sample(sample, control, bamfile_suffix, work_dir)
			##annotate using annovar
			mt2_vcf = sample + '.mt2.final.vcf.gz'
			annotate_mt2_vars_from_vcf(sample, mt2_vcf)
			##filter ann txt file for passed and exonic
			filter_ann_file_mt2(sample)



##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/rich_exomes_0720'
os.chdir(working_dir)

project_name = 'rich_exomes_0720'
mt2_project_name = 'rich_exomes_0720' + '_mt2'
sample_dict = {'TZ001': 'TZ001.368077.bam', 'TZ002': 'TZ002.368078.bam', 'TZ003': 'TZ003.368079.bam', 'TZ004': 'TZ004.368080.bam', 
		'TZ005': 'TZ005.368081.bam', 'TZ006': 'TZ006.368082.bam', 'TZ007': 'TZ007.370762.bam', 'TZ008': 'TZ008.368084.bam', 
		'TZ009': 'TZ009.368085.bam'}
mt2_sample_dict = {'TZ007': ['TZ004', 'TZ005', 'TZ006'], 'TZ008': ['TZ001', 'TZ002'], 'TZ009': ['TZ003'] }
bam_suffix = '.bwa_gatk.bam'

##make fastq files from bams
# convert_bam_fastq_picard(sample_dict)

##run xengsort
# sample_dict = ['TZ001']
# run_xengsort(sample_dict)

##align, call vars and annotate (need to get coverage)
# exome_analysis_pipeline(sample_dict, project_name)

##call and filter mutect2 
mutect2_pipeline(mt2_sample_dict, mt2_project_name, bam_suffix, working_dir)



