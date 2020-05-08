#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import numpy
'''
module load local_python/3.6.4 
'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/genomes/fetal_seq_genome_mh_0518'
os.chdir(working_dir)

##file names etc
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
exome_capture_bed = ref_dir + 'dobyns_exome.in_any_padded_target.1015.bed'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
pisces_fastq_dir = '/home/atimms/ngs_data/references/illumina_human_g1k_v37/'
mh_ref = '/home/atimms/programs/MosaicHunter/resources/'
mh_rpts = mh_ref + 'all_repeats.b37.bed'
mh_common_site = mh_ref + 'WGS.error_prone.b37.bed'
mh_dbsnp = mh_ref + 'dbsnp_137.b37.tsv'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.b37.vcf'
indels_1000g = ref_dir + '1000G_phase1.indels.b37.vcf'
rtc_intervals = ref_dir + 'Mills_1000G.intervals'
dbsnp = ref_dir + 'dbsnp_138.b37.vcf'
# exome_capture_bed = ref_dir + 'dobyns_exome.in_any_padded_target.1015.bed'
exome_bed_for_coverage = ref_dir + 'dobyns_exome.in_all_targets.1015.bed'

acceptable_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
				'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
				'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 
				'chr23', 'chr24', 'chr25', 'chrx', 'chry', 'chrX', 'chrY', 
				'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
				'13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', 
				'24', '25', 'x', 'y', 'X', 'Y']

##programs
mosaic_hunter = '/home/atimms/programs/MosaicHunter/build/mosaichunter.jar'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
ann_var = '/home/atimms/programs/annovar/annotate_variation.pl'
# gatk = '/cm/shared/apps/GATK/GenomeAnalysisTK.jar'
# gatk2 = '/home/atimms/programs/GenomeAnalysisTK.jar'
pisces = '/cm/shared/apps/Pisces/5.1.6.54/Pisces.exe'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
vardictjava = '/home/atimms/programs/VarDictJava/build/install/VarDict/bin/VarDict'
vardict_maf_req = '0.01'
bgzip = '/home/atimms/programs/htslib-1.8/bgzip'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
bcftools = '/home/atimms/programs/bcftools-1.8/bcftools'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147']
av_operation = ['-operation', 'g,r,r,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']



def variant_calling_single_mosaichunter(samples, bam_suffix, sex):
	for sample in samples:
		in_bam = sample + bam_suffix
		clean_bam = sample + '.temp_clean.bam'
		##clean the bam
		##samtools view -h -f 0x2 input.bam | perl -ne 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))' | samtools view -Sb - >cleaner.bam
		'''
		get_paired = subprocess.Popen(['samtools', 'view', '-h', '-f', '0x2', in_bam], stdout=subprocess.PIPE)
		rm_mismatched = subprocess.Popen(['perl', '-ne', 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))'], stdin=get_paired.stdout, stdout=subprocess.PIPE)
		st_view = subprocess.Popen(['samtools', 'view','-b', '-o', clean_bam], stdin=rm_mismatched.stdout)
		st_view.wait()
		st_index_run = subprocess.Popen(['samtools', 'index', clean_bam])
		st_index_run.wait()
		print 'made clean_bam:', clean_bam
		'''
		##run mosaic hunter in single 'mode'
		##default
		outdir = sample + '_mh_10_600_rpt_2'
		run_mh_single = subprocess.Popen(['java', '-jar', mosaic_hunter, 'genome', '-P', 'input_file=' + clean_bam, '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=single',
				 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
				 '-P', 'output_dir=' + outdir, '-P', 'base_number_filter.min_minor_allele_number=2', '-P', 'base_number_filter.min_minor_allele_percentage=1', '-P', 'mosaic_filter.sex', sex, 
				 '-P', 'min_depth=10', '-P', 'max_depth=600', '-P', 'depth_filter.max_depth=600', '-P', 'depth_filter.min_depth=10'])
		run_mh_single.wait()
		'''
		outdir = sample + '_mh_default'
		run_mh_single1 = subprocess.Popen(['java', '-jar', mosaic_hunter, 'genome', '-P', 'input_file=' + clean_bam, '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=single',
				 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
				 '-P', 'output_dir=' + outdir, '-P', 'mosaic_filter.sex', sex])
		run_mh_single1.wait()
		outdir = sample + '_mh_allele_freq_number'
		run_mh_single1b = subprocess.Popen(['java', '-jar', mosaic_hunter, 'genome', '-P', 'input_file=' + clean_bam, '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=single',
				 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
				 '-P', 'output_dir=' + outdir, '-P', 'mosaic_filter.sex', sex, '-P', 'min_minor_allele_number=1', '-P', 'min_minor_allele_percentage=1'])
		run_mh_single1b.wait()

		outdir = sample + '_mh_min_max_depth'
		run_mh_single2 = subprocess.Popen(['java', '-jar', mosaic_hunter, 'genome', '-P', 'input_file=' + clean_bam, '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=single',
				 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
				 '-P', 'output_dir=' + outdir, '-P', 'mosaic_filter.sex', sex, '-P', 'min_depth=10', '-P', 'max_depth=600', '-P', 'depth_filter.max_depth=300', '-P', 'depth_filter.min_depth=10'])
		run_mh_single2.wait()
		outdir = sample + '_mh_depth_mlf_ap'
		run_mh_single3 = subprocess.Popen(['java', '-jar', mosaic_hunter, 'genome', '-P', 'input_file=' + clean_bam, '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=single',
				 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
				 '-P', 'output_dir=' + outdir, '-P', 'mosaic_filter.sex', sex, '-P', 'min_depth=10', '-P', 'max_depth=600', '-P', 'depth_filter.max_depth=300', '-P', 'depth_filter.min_depth=10',
				 '-P', 'mosaic_like_filter.min_minor_allele_percentage=1'])
		run_mh_single3.wait()
		outdir = sample + '_mh_depth_mlf_ap_mrate6'
		run_mh_single4 = subprocess.Popen(['java', '-jar', mosaic_hunter, 'genome', '-P', 'input_file=' + clean_bam, '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=single',
				 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
				 '-P', 'output_dir=' + outdir, '-P', 'mosaic_filter.sex', sex, '-P', 'min_depth=10', '-P', 'max_depth=600', '-P', 'depth_filter.max_depth=300', '-P', 'depth_filter.min_depth=10',
				 '-P', 'mosaic_like_filter.min_minor_allele_percentage=1', '-P', 'mosaic_filter.mosaic_rate = 1e-6'])
		run_mh_single4.wait()
		outdir = sample + '_mh_depth_mlf_ap_mrate8'
		run_mh_single5 = subprocess.Popen(['java', '-jar', mosaic_hunter, 'genome', '-P', 'input_file=' + clean_bam, '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=single',
				 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
				 '-P', 'output_dir=' + outdir, '-P', 'mosaic_filter.sex', sex, '-P', 'min_depth=10', '-P', 'max_depth=600', '-P', 'depth_filter.max_depth=300', '-P', 'depth_filter.min_depth=10',
				 '-P', 'mosaic_like_filter.min_minor_allele_percentage=1', '-P', 'mosaic_filter.mosaic_rate = 1e-8'])
		run_mh_single5.wait()
		'''
		##delete clean bam
		'''
		os.remove(clean_bam)
		os.remove(clean_bam + '.bai')
		print 'deleting', clean_bam, clean_bam + '.bai'
		'''

def make_list_of_bams(bam_files, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')

def genotype_vars_with_freebayes(bedfile, final_vcf):
	vcf_temp1 = 'temp_fb1.vcf'
	vcf_temp2 = 'temp_fb2.vcf'
	with open(vcf_temp1, 'w') as vcf_fh:
		run_freebayes = subprocess.Popen([freebayes, '-f', fasta, '-L', bamlist, '-t', bedfile, '--haplotype-length', '0', '--min-alternate-count', '1', '--min-alternate-fraction', '0', '--pooled-continuous', '--report-monomorphic'], stdout=vcf_fh)
		run_freebayes.wait()
	run_bgzip = subprocess.Popen([bgzip, vcf_temp1])
	run_bgzip.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	os.remove( vcf_temp1 + '.gz')

##convert vcf file to individual avinput file
def convert_to_annovar(vcf, project_prefix):
	av_prefix = project_prefix
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', av_prefix])
	con_ann.wait()

##annotate vcf file
def table_annovar_vcf(vcf, project_prefix):
		out_prefix = project_prefix
		command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()

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

def get_sample_names(vcf):
	with open(vcf, "r") as vcf_fh:
		for line in vcf_fh:
			if line[:2] != '##':
				if line[0] == "#":
					line = line.rstrip().split(delim)
					samples = line[9:]
	return samples

def format_avinput(project_prefix, annovar_vcf):
	multianno = project_prefix + '.hg19_multianno.txt'
	outfile = project_prefix + '.annotated.txt'
	samples = get_sample_names(annovar_vcf)
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred', 'integrated_fitCons_score', 'integrated_confidence_value', 'GERP++_RS', 'phyloP7way_vertebrate', 'phyloP20way_mammalian', 'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'PopFreqMax', '1000G_ALL', '1000G_AFR', '1000G_AMR', '1000G_EAS', '1000G_EUR', '1000G_SAS', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'ESP6500siv2_ALL', 'ESP6500siv2_AA', 'ESP6500siv2_EA', 'CG46', 'avsnp147', 'na', 'na', 'coverage', 'chr', 'pos', 'id', 'ref2', 'alt2', 'qual', 'filter', 'info', 'format'] + samples
	head_out = delim.join(head + ['\n'])
	with open(multianno, "r") as multi_fh, open(outfile, "w") as out_fh:
		out_fh.write(head_out)
		line_count = 0
		for line in multi_fh:
			line_count += 1
			if line_count >1:
				out_fh.write(line)


##calls all annovar methods
def annotate_vcf_file(vcf, project_prefix):
	table_annovar_vcf(vcf, project_prefix)
	post_annovar_vcf = project_prefix + '.' + av_genome + '_multianno.vcf'
	convert_to_annovar(post_annovar_vcf, project_prefix)
	format_avinput(project_prefix, post_annovar_vcf)


def median(l):
	return numpy.median(numpy.array(l))

def combine_mh_results_with_fb_vcf(inbed, inann, outfile, sample_pos):
	##make list from bed file
	bed_var_list = []
	with open(inbed, "r") as inb_fh, open(outfile, "w") as out_fh:
		for line_bed in inb_fh:
			line_bed = line_bed.rstrip().split(delim)
			bed_var = '_'.join([line_bed[0]] + line_bed[2:5])
			# print bed_var
			bed_var_list.append(bed_var)

	with open(outfile, "w") as out_fh:
		with open(inann, "r") as ina_fh:
			line_count = 0
			for line_ann in ina_fh:
				line_ann = line_ann.rstrip().split(delim)
				line_count += 1
				if line_count == 1:
					header = line_ann[:5] + line_ann[10:12] + [line_ann[46], line_ann[65]]
					sample_list = line_ann[sample_pos[0]:sample_pos[1]]
					for sample in sample_list:
						s_cov = sample + ' coverage'
						s_aaf = sample + ' aaf'
						header.extend([s_cov, s_aaf])
					# print header
					out_fh.write(delim.join(header + ['max_aaf', 'median_cov', 'mean_cov' , '\n']))
					##get position of maternal tissue
					# maternal_position = [i for i, s in enumerate(sample_list) if 'maternal' in s]
					# print maternal_position
				else:
					ann_var = '_'.join([line_ann[0]] + line_ann[2:5])
					if ann_var in bed_var_list:
						# print ann_var
						##get all genotype
						genotypes = line_ann[sample_pos[0]:sample_pos[1]]
						# ##maternal genotype
						# maternal_info = genotypes[maternal_position[0]].split(':')
						# print maternal_info
						# mat_alleles = [maternal_info[2], maternal_info[4]]
						##line to write so var + 
						line_out = line_ann[:5] + line_ann[10:12] + [line_ann[46], line_ann[65]]
						##get all genotype and write
						aafs, covs = [], []
						for genotype in genotypes:
							genotype_info = genotype.split(':')
							alleles = genotype_info[2].split(',')
							print(genotype_info, alleles)
							if alleles[0] == '.':
								coverage = 0
							else:
								coverage = int(alleles[0]) + int(alleles[1])
							if coverage == 0:
								allele_freq = 0
							else:
								##get aaf 
								allele_freq = float(alleles[1]) / coverage
							# print coverage, allele_freq
							aafs.append(allele_freq)
							covs.append(coverage)
							line_out.extend([str(coverage), str(allele_freq)])
						cov_median = median(covs)
						cov_mean = sum(covs) / len(covs)
						# print covs, cov_median, cov_mean
						line_out.extend([str(max(aafs)), str(cov_median), str(cov_mean), '\n'])
						out_fh.write(delim.join(line_out))


##run methods
#parameters for this run
project_name = 'fetalseq_genome_0518'
tissue_list = ['Brain', 'Heart', 'Liver', 'Lung', 'Muscle', 'Placenta', 'fg_combined']

gatk_bam_suffix = '.bwa_gatk.bam'
bam_files_for_cov = [t + gatk_bam_suffix for t in tissue_list]

##on all samples
##mosaic hunter -- finally decided just to use combined file-- using '_mh_10_600_rpt_2'
# tissue_list = ['fg_combined']
# variant_calling_single_mosaichunter(tissue_list, gatk_bam_suffix, 'M')

##step 2 take candidate vars and annaotate / check against controls
project_name = 'fetal_variants_0718'
bam_files = glob.glob('*.bam')
bamlist = 'bams.list'
bed_file = project_name + '.bed'
fb_vcf = project_name + '.fb.vcf'
ann_file = project_name + '.annotated.txt'
final_file = project_name + '.added_info.xls'
# make_list_of_bams(bam_files, bamlist)
# genotype_vars_with_freebayes(bed_file, fb_vcf)
# annotate_vcf_file(fb_vcf, project_name)
combine_mh_results_with_fb_vcf(bed_file, ann_file, final_file, [78,88])

