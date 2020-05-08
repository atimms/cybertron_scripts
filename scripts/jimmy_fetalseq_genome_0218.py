#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil
import numpy
import pybedtools as pbt

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/genomes/fetal_seq_genome_0218'
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
mh_common_site = mh_ref + 'WES_Agilent_50M.error_prone.b37.bed'
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
gatk = '/cm/shared/apps/GATK/GenomeAnalysisTK.jar'
gatk2 = '/home/atimms/programs/GenomeAnalysisTK.jar'
pisces = '/cm/shared/apps/Pisces/5.2.0.1/Pisces/Pisces.dll'
stitcher = '/cm/shared/apps/Pisces/5.2.0.1/Stitcher/Stitcher.dll'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
vardictjava = '/home/atimms/programs/VarDictJava/build/install/VarDict/bin/VarDict'
vardict_maf_req = '0.01'
bgzip = 'bgzip'
bcftools = 'bcftools'

#annovar parameters
# av_genome = 'hg19'
# av_buildver = ['-buildver', av_genome]
# av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147']
# av_protocol_pisces = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147,vcf', '-vcfdbfile']
# av_operation = ['-operation', 'g,r,r,f,f,f']
# av_operation_pisces = ['-operation', 'g,r,r,f,f,f,f']
# av_options = ['-otherinfo', '-remove', '-nastring', '.','-arg', '-splicing 10 ,,,,,']
# av_options_vcf = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput','-arg', '-splicing 10 ,,,,,']
# av_options_vcf_pisces = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput','-arg', '-splicing 10 ,,,,,,']
##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,gnomad_genome,kaviar_20150923,avsnp147']
av_operation = ['-operation', 'g,r,r,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']


def run_stitcher_on_bam(samples, bam_suffix):
	bam_files = []
	for sample in samples:
		in_bam = sample + bam_suffix
		out_bam = sample + '.bwa_gatk_stitcher.bam' ##???
		bam_files.append(out_bam)
		##check bai file exists and make if it isn't there
		if os.path.isfile(in_bam + '.bai'):
			print 'bam %s alreaded indexed'%in_bam
		else:
			print 'indexing bam file:', in_bam
			st_index = subprocess.Popen(['samtools', 'index', in_bam])
			st_index.wait()
		##run sticher, then sort and index bam
		run_stitcher = subprocess.Popen(['dotnet', stitcher, '-Bam', in_bam, '-t', '18', '-OutFolder', './'])
		run_stitcher.wait()

def combine_sticher_bams(samples, combined_prefix, bam_suffix):
	outfile = combined_prefix + bam_suffix
	sorted_out_bam = combined_prefix + '.bwa_gatk.stitched.sorted.bam'
	# '''
	bams = []
	for sample in samples:
		in_bam = sample + bam_suffix
		sorted_bam = sample + '.bwa_gatk.stitched.sorted.bam'
		out_bam = sample + '.bwa_gatk.stitched.rg_adjusted.bam'
		bams.extend(['INPUT=', out_bam])
		# '''
		st_sort_pe = subprocess.Popen(['samtools', 'sort', '-O', 'bam', '-T', sample, '-o', sorted_bam, '-@', '10', '-m', '10G', in_bam])
		st_sort_pe.wait()
		st_index = subprocess.Popen(['samtools', 'index', sorted_bam])
		st_index.wait()
		##make new bam file with adjusted RG info
		picard_arrd = subprocess.Popen(['picard', 'AddOrReplaceReadGroups', 'INPUT=' + sorted_bam, 'OUTPUT=' + out_bam, 'RGPL=illumina', 'RGLB=' + sample, 'RGPU=machine1', 'RGID=' + sample, 'RGSM=' + combined_prefix, 'TMP_DIR=' + working_dir, 'CREATE_INDEX=true', 'VALIDATION_STRINGENCY=LENIENT'])
		picard_arrd.wait()
		# '''
	##then merge those sams
	picard_merge = subprocess.Popen(['picard', 'MergeSamFiles'] + bams + ['OUTPUT=' + outfile, 'TMP_DIR=' + working_dir, 'CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'])
	picard_merge.wait()
	# '''
	st_sort_pe = subprocess.Popen(['samtools', 'sort', '-O', 'bam', '-T', combined_prefix, '-o', sorted_out_bam, '-@', '10', '-m', '10G', outfile])
	st_sort_pe.wait()
	st_index = subprocess.Popen(['samtools', 'index', sorted_out_bam])
	st_index.wait()


def variant_calling_single_pisces(ped_name, samples, bam_suffix):
	out_dir = ped_name + '_pisces'
	bam_files = []
	for sample in samples:
		bam = sample + bam_suffix
		bam_files.append(bam)
		##check bai file exists and make if it isn't there
		if os.path.isfile(bam + '.bai'):
			print 'bam %s alreaded indexed'%bam
		else:
			print 'indexing bam file:', bam
			st_index = subprocess.Popen(['samtools', 'index', bam])
			st_index.wait()
	##run picses on all bam 
	run_pisces = subprocess.Popen(['dotnet', pisces, '-G', pisces_fastq_dir, '-BamPaths', ','.join(bam_files), '-MinVF', '0.01', '-gVCF', 'FALSE', '-t', '18', '-OutFolder', out_dir])
	run_pisces.wait()

def merge_vcf_files(vcfs, out_vcf):
	vcfs_to_combine_cmd = []
	for vcf in vcfs:
		vcfs_to_combine_cmd.extend(['--variant', vcf])
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R', fasta, '-nt', '15'] + vcfs_to_combine_cmd + ['-o', out_vcf, '-genotypeMergeOptions', 'UNSORTED', '-filteredAreUncalled'])
	combine_var.wait()


def merge_pisces_results(samples, snp_comb_vcf, indel_comb_vcf):
	snp_vcfs_to_combine, indel_vcfs_to_combine = [], []
	for sample in samples:
		in_vcf = 'fetalseq_genome_0218_pisces/' + sample + '.bwa_gatk.stitched.sorted.vcf'
		snp_vcf = sample + '.pisces_passed_snps.vcf'
		indels_vcf = sample + '.pisces_passed_indels.vcf'
		snp_vcfs_to_combine.append(snp_vcf)
		indel_vcfs_to_combine.append(indels_vcf)
		##remove non passed varianst
		# '''
		snp_cut = subprocess.Popen(['java', '-Xmx200g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '15', '-V', in_vcf, '-o', snp_vcf, '--excludeFiltered','-selectType', 'SNP', '-restrictAllelesTo', 'BIALLELIC'])
		snp_cut.wait()
		indel_cut = subprocess.Popen(['java', '-Xmx200g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '15', '-V', in_vcf, '-o', indels_vcf, '--excludeFiltered', '-selectType', 'INDEL', '-restrictAllelesTo', 'BIALLELIC'])
		indel_cut.wait()
		# '''
	merge_vcf_files(snp_vcfs_to_combine, snp_comb_vcf)
	merge_vcf_files(indel_vcfs_to_combine, indel_comb_vcf)

def make_list_of_bams(samples, bamlist_file):
	bam_files = [s + '.bwa_gatk.bam' for s in samples]
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')

def genotype_vars_with_gatk(bamlist, combined_vcf, gatk_norm_vcf):
	vcf_temp1 = gatk_norm_vcf + 'temp_21.vcf'
	vcf_temp2 = gatk_norm_vcf + 'temp_22.vcf.gz'
	gatk_ug = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'UnifiedGenotyper', '-R', fasta, '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', combined_vcf, '-dcov', '2500', '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', vcf_temp1])
	gatk_ug.wait()
	##split multi-allelic variants calls in separate lines, and left normalize indels
	bgzip_run = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'v', '-o', gatk_norm_vcf, vcf_temp2])
	bcf_norm2.wait()

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
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred', 'integrated_fitCons_score', 'integrated_confidence_value', 'GERP++_RS', 'phyloP7way_vertebrate', 'phyloP20way_mammalian', 'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'gnomAD_genome_ALL', 'gnomAD_genome_AFR', 'gnomAD_genome_AMR', 'gnomAD_genome_ASJ', 'gnomAD_genome_EAS', 'gnomAD_genome_FIN', 'gnomAD_genome_NFE', 'gnomAD_genome_OTH', 'Kaviar_AF', 'Kaviar_AC', 'Kaviar_AN', 'avsnp147', 'na', 'na', 'coverage', 'chr', 'pos', 'id', 'ref2', 'alt2', 'qual', 'filter', 'info', 'format'] + samples
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

def reformat_mutect_gatk_ann_txt(ind, inann_suffix, outfile_suffix, sample_pos):
	inann = ind + inann_suffix
	outfile = ind + outfile_suffix
	##make list from bed file
	with open(outfile, "w") as out_fh, open(inann, "r") as ina_fh:
		line_count = 0
		for line_ann in ina_fh:
			line_ann = line_ann.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				header = line_ann[:12] + [line_ann[46], line_ann[54], line_ann[57]]
				sample_list = line_ann[sample_pos[0]:sample_pos[1]]
				for sample in sample_list:
					s_cov = sample + ' coverage'
					s_aaf = sample + ' aaf'
					header.extend([s_cov, s_aaf])
				# print header
				out_fh.write(delim.join(header + ['\n']))
				##get position of maternal tissue
				maternal_position = [i for i, s in enumerate(sample_list) if 'maternal' in s]
				# print maternal_position
			else:
				##get all genotype
				genotypes = line_ann[sample_pos[0]:sample_pos[1]]
				##line to write so var + 
				line_out = line_ann[:12] + [line_ann[46], line_ann[54], line_ann[57]]
				##get all genotype and write
				aafs, covs = [], []
				for genotype in genotypes:
					genotype_info = genotype.split(':')
					print genotype_info
					if genotype_info == ['./.']:
						alleles = [0,0]
					elif genotype_info[1] == '.' or genotype_info[0] == './.':
						alleles = [0,0]
					else:
						alleles = [genotype_info[1].split(',')[0], genotype_info[1].split(',')[1]]
					print alleles
					# if alleles[0] == '.':
					# 	coverage = 0
					# else:
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
				##remove
				'''
				cov_median = median(covs)
				cov_mean = sum(covs) / len(covs)
				# print covs, cov_median, cov_mean
				line_out.extend([str(max(aafs)), str(cov_median), str(cov_mean), '\n'])
				'''
				out_fh.write(delim.join(line_out +['\n']))

def filter_var_file(ind, infile_suffix, outfile_suffix):
	infile = ind + infile_suffix
	outfile = ind + outfile_suffix
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count == 1:
				out_fh.write(line)
			else:
				line = line.split(delim)
				rmsk = line[10]
				supdup = line[11]
				gnomad = line[12]
				kaviar = line[13]
				chrom = line[0]
				if rmsk == '.' and supdup == '.' and chrom[:2] != 'GL':
					if gnomad == '.'or float(gnomad) <= 0.01:
						if kaviar == '.'or float(kaviar) <= 0.01:
							out_fh.write(delim.join(line))


def filter_var_file_extras(ind, infile_suffix, outfile_suffix, aaf_req, control_reads_req, alt_reads_req, aaf_calc_pos, cov_calc_pos, aaf_cont_pos, cov_cont_pos):
	infile = ind + infile_suffix
	outfile = ind + outfile_suffix
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				header = line + ['average_cov', 'max_aaf', 'sum_alt_reads', 'control_max_alt_reads', '\n']
				out_fh.write(delim.join(header))
			else:
				aafs = [float(line[i]) for i in tuple(aaf_calc_pos)]
				covs = [int(line[i]) for i in tuple(cov_calc_pos)]
				alt_reads = numpy.multiply(aafs,covs)
				max_aaf = max(aafs)
				sum_alt_reads = sum(alt_reads)
				average_cov = sum(covs) / len(covs)
				cont_aafs = [float(line[i]) for i in tuple(aaf_cont_pos)]
				cont_covs = [int(line[i]) for i in tuple(cov_cont_pos)]
				cont_alt_reads = numpy.multiply(cont_aafs,cont_covs)
				cont_max_alt_reads = max(cont_alt_reads)
				chrom = line[0]
				# print aafs, covs, alt_reads
				# print max_alt_reads, max_aaf
				if max_aaf < aaf_req and sum_alt_reads >= alt_reads_req and cont_max_alt_reads <= control_reads_req:
					# print 'kept'
					lineout = line + [str(average_cov), str(max_aaf), str(sum_alt_reads), str(cont_max_alt_reads), '\n']
					out_fh.write(delim.join(lineout))

def make_file_for_graphing_in_r(ind, infile_suffix, outfile_suffix, tissue_dict):
	infile = ind + infile_suffix
	var_info_dict = {}
	for tissue in tissue_dict:
		outfile = tissue + outfile_suffix
		with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
			out_fh.write(delim.join(['var', 'coverage', 'aaf', '\n']))
			lc = 0
			for line in in_fh:
				lc += 1
				if lc >1:
					line = line.split(delim)
					var = '_'.join(line[:5])
					cov = line[tissue_dict[tissue][0]]
					aaf = line[tissue_dict[tissue][1]]
					out_fh.write(delim.join([var, cov, aaf, '\n']))

def vars_specific_to_one_sample(ind, infile_suffix, outfile_suffix, tissue_dict, reads_req, aaf_req):
	infile = ind + infile_suffix
	for tissue in tissue_dict:
		outfile = ind + '.' + tissue + outfile_suffix
		with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
			line_count = 0
			for line in in_fh:
				line_count += 1
				line = line.rstrip().split(delim)
				if line_count == 1:
					header = line + ['\n']
					out_fh.write(delim.join(header))
				else:
					t_cov = int(line[tissue_dict[tissue][0]])
					t_aaf = float(line[tissue_dict[tissue][1]])
					t_reads = round(t_cov * t_aaf)
					other_cov = [int(line[tissue_dict[t][0]]) for t in tissue_dict if t is not tissue]
					other_aaf = [float(line[tissue_dict[t][1]]) for t in tissue_dict if t is not tissue]
					other_max_aaf = max(other_aaf)
					other_max_alt_reads = round(max(numpy.multiply(other_aaf,other_cov)))
					print tissue
					print(t_cov, t_aaf, t_reads)
					# print other_aaf, other_max_aaf, other_max_alt_reads
					if t_reads >= reads_req and t_aaf >= aaf_req and other_max_alt_reads < reads_req and other_max_aaf < aaf_req:
						lineout = line + ['\n']
						out_fh.write(delim.join(lineout))

def vars_in_all_samples(ind, infile_suffix, outfile_suffix, reads_req, aaf_req, aaf_calc_pos, cov_calc_pos):
	infile = ind + infile_suffix
	outfile = ind + outfile_suffix
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				header = line + ['\n']
				out_fh.write(delim.join(header))
			else:
				aafs = [float(line[i]) for i in tuple(aaf_calc_pos)]
				covs = [int(line[i]) for i in tuple(cov_calc_pos)]
				alt_reads = numpy.multiply(aafs,covs)
				# print aafs, alt_reads
				if all(i >= aaf_req for i in aafs) and all(i >= reads_req for i in alt_reads):
					# print 'passed'
					lineout = line + ['\n']
					out_fh.write(delim.join(lineout))


##run methods
#parameters for this run
project_name = 'fetalseq_genome_0218'
samples_for_sticher = ['Brain', 'Heart', 'Liver', 'Lung', 'Muscle', 'Placenta']
samples_for_variant_calling = ['Brain', 'Heart', 'Liver', 'Lung', 'Muscle', 'Placenta', 'fg_combined']
samples_with_controls = ['Brain', 'Heart', 'Liver', 'Lung', 'Muscle', 'Placenta', 'fg_combined', 'LR16-302', 'LR16-302f', 'LR16-302m']
gatk_bam_suffix = '.bwa_gatk.bam'
stitched_suffix = '.bwa_gatk.stitched.bam'
stiched_sorted_suffix = '.bwa_gatk.stitched.sorted.bam'
combined_name = 'fg_combined'
all_sample_snp_vcf = project_name + '.all.pisces_snps.vcf'
all_sample_indels_vcf = project_name + '.all.pisces_indels.vcf'
gatk_snp_genotyes_vcf = project_name + '.gatk_snp_genotyes.vcf'
gatk_snp_prefix = gatk_snp_genotyes_vcf.split('.')[0]
gatk_indel_genotyes_vcf =  project_name + '.gatk_indel_genotyes.vcf'

##run sticher on individual then combine bams
# run_stitcher_on_bam(samples_for_sticher, gatk_bam_suffix)
##repeat just brain sample
# run_stitcher_on_bam(['Brain'], gatk_bam_suffix)
# combine_sticher_bams(samples_for_sticher, combined_name, stitched_suffix)
# variant_calling_single_pisces(project_name, samples_for_variant_calling, stiched_sorted_suffix)
##combine
# merge_pisces_results(samples_for_variant_calling, all_sample_snp_vcf, all_sample_indels_vcf)
##genotype snps and indels
# make_list_of_bams(samples_with_controls, bamlist)
# genotype_vars_with_gatk(bamlist, all_sample_snp_vcf, gatk_snp_genotyes_vcf)
##not done
# genotype_vars_with_gatk(bamlist, all_sample_indels_vcf, gatk_indel_genotyes_vcf)


##annotate and filter snps
# annotate_vcf_file(gatk_snp_genotyes_vcf, gatk_snp_prefix)
# reformat_mutect_gatk_ann_txt(gatk_snp_prefix, '.annotated.txt', '.pisces.all_vars.xls', [70,80])
# filter_var_file(gatk_snp_prefix, '.pisces.all_vars.xls', '.pisces.filtered_vars.xls')
##positions for getting tissue and controls aafs/converage
pos_for_aaf_calcs = [16,18,26,28,30,32]
pos_for_cov_calcs = [15,17,25,27,29,31]
pos_for_aaf_controls = [20,22,24]
pos_for_cov_controls = [19,21,23]
##number or greater reads requires when combining all sample
alt_reads_reqs = [3,5,10,20]
##if this many reads in controls remove
control_reads_to_filters = [0,1]
##maximum aaf in any tissue --remove if greater
aaf_reqs = [0.35,0.75]
##list of positions of cov/aaf of ind tissues
tissue_dict = {'Brain':[15,16], 'Heart':[17,18], 'Liver':[25,26], 'Lung':[27,28], 'Muscle':[29,30], 'Placenta':[31,32]}


'''
##for many tests
for alt_reads_req in alt_reads_reqs:
	for control_reads_to_filter in control_reads_to_filters:
		for aaf_req in aaf_reqs:
			filtered_prefix = '.pisces.filtered_vars_0218.' + str(aaf_req) + '_aaf.' + str(control_reads_to_filter) + '_cont_reads.' + str(alt_reads_req) + '_alt_reads.xls'
			graphing_prefix = '.pisces.filtered_vars_0218.' + str(aaf_req) + '_aaf.' + str(control_reads_to_filter) + '_cont_reads.' + str(alt_reads_req) + '_alt_reads.for_r.txt'
			##filter vars
			filter_var_file_extras(gatk_snp_prefix, '.pisces.filtered_vars.xls', filtered_prefix, aaf_req, control_reads_to_filter, alt_reads_req, pos_for_aaf_calcs, pos_for_cov_calcs, pos_for_aaf_controls, pos_for_cov_controls)
			##take var file and split for graphing
			if aaf_req == 0.75:
				make_file_for_graphing_in_r(gatk_snp_prefix, filtered_prefix, graphing_prefix, tissue_dict)
'''




##filter in different ways ie in all samples
##vars just in one sample
# vars_specific_to_one_sample(gatk_snp_prefix, '.pisces.filtered_vars_0218.0.35_aaf.1_cont_reads.5_alt_reads.xls', '_specific.xls', tissue_dict, 2, 0.01)
# vars_specific_to_one_sample(gatk_snp_prefix, '.test.xls', '_specific.xls', tissue_dict, 2, 0.01)
vars_in_all_samples(gatk_snp_prefix, '.pisces.filtered_vars_0218.0.35_aaf.1_cont_reads.5_alt_reads.xls', '_in_all_samples.xls', 2, 0.01, pos_for_aaf_calcs, pos_for_cov_calcs)



