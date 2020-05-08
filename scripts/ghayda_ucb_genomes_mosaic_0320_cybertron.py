#!/usr/bin/env python
import os
import subprocess
import glob
import shutil


'''
load these:
module load mono/5.10.1.47
module load biobuilds
'''

##set input variables and parameters
delim = '\t'

##program
gatk4 = '/home/atimms/programs/gatk-4.1.3.0/gatk'
pisces = '/cm/shared/apps/Pisces/5.1.6.54/Pisces.exe'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
ann_var = '/home/atimms/programs/annovar_1019/annotate_variation.pl'
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
gatk2 = '/home/atimms/programs/GenomeAnalysisTK.jar'

##files etc
pisces_fastq_dir = '/home/atimms/ngs_data/references/illumina_hg38/'
fasta = '/home/atimms/ngs_data/references/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta'
mutect_gnomad_vcf = '/home/atimms/ngs_data/references/hg38/af-only-gnomad.hg38.vcf.gz'



##annovar parameters --- add cosmic
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp35a,dbnsfp31a_interpro,avsnp150,gnomad211_genome,gnomad211_exome,gnomad30_genome,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']

##methods
def make_analysis_dict(input_file):
	analysis_dict = {}
	with open(input_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				ped_name = line[0]
				sample_name = line[1]
				id_name = line[6]
				anal_type = line[7]
				affection = line[5]
			##add data to analysis_dict
				if ped_name in analysis_dict:
					if id_name in analysis_dict[ped_name][0]:
						print('issue with sample', id_name, sample_name)
					else:
						analysis_dict[ped_name][0][id_name] = sample_name
						if affection == '2':
							analysis_dict[ped_name][1].append(id_name)
						elif affection == '1':
							analysis_dict[ped_name][2].append(id_name)
				else:
					##if affected sample
					if affection == '2':
						analysis_dict[ped_name] = [{id_name:sample_name}, [id_name], [], anal_type]
					##or parent
					elif affection == '1':
						analysis_dict[ped_name] = [{id_name:sample_name}, [], [id_name], anal_type]
	return analysis_dict


def variant_calling_single_pisces(ped_name, samples, bam_suffix):
	out_dir = ped_name + '_pisces'
	# out_dir_gvcf = ped_name + '_pisces_gvcf'
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
	##run picses on all bam - ? used vcf or gvcf
	# run_pisces_gvcf = subprocess.Popen(['mono', pisces, '-G', pisces_fastq_dir, '-BamPaths', ','.join(bam_files), '-MinVF', '0.01', '-gVCF', 'FALSE', '-t', '18', '-OutFolder', out_dir_gvcf])
	# run_pisces_gvcf.wait()
	run_pisces = subprocess.Popen(['mono', pisces, '-G', pisces_fastq_dir, '-BamPaths', ','.join(bam_files), '-MinVF', '0.01', '-t', '18', '-OutFolder', out_dir])
	run_pisces.wait()

def make_list_of_bams(bam_files, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')

def filter_norm_pisces_vcf(pisces_vcf,pisces_filtered_vcf, sample_id):
	vcf_temp1 = sample_id + 'temp_01.vcf.gz'
	vcf_temp2 = sample_id + 'temp_02.vcf.gz'
	##compress and index vcf
	bgzip_run = subprocess.Popen(['bgzip', pisces_vcf])
	bgzip_run.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', pisces_vcf + '.gz'])
	bcf_index.wait()
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-r', 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM', '-o', vcf_temp1, '-O', 'z', pisces_vcf + '.gz'])
	bcftools_filter.wait()
	##split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_index2 = subprocess.Popen(['bcftools', 'index', vcf_temp1])
	bcf_index2.wait()
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'v', '-o', pisces_filtered_vcf, vcf_temp2])
	bcf_norm2.wait()

def genotype_vars_with_gatk(bamlist, combined_vcf, gatk_norm_vcf, sample_id):
	vcf_temp1 = sample_id + 'temp_21.vcf'
	vcf_temp2 = sample_id + 'temp_22.vcf.gz'
	# gatk_ug = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk2, '-T', 'UnifiedGenotyper', '-R', fasta, '-nct', '8', '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', combined_vcf, '-dcov', '2500', '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', vcf_temp1])
	# gatk_ug = subprocess.Popen(['java', '-Djava.io.tmpdir=.', '-Xmx120g', '-jar', gatk2, '-T', 'UnifiedGenotyper', '-R', fasta, '-nct', '8', '-nt', '4', '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', combined_vcf, '-dcov', '2500', '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', vcf_temp1])
	##this works most of the time
	# gatk_ug = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk2, '-T', 'UnifiedGenotyper', '-R', fasta, '-nct', '3', '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', combined_vcf, '-dcov', '2500', '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', vcf_temp1])
	##temp without nct
	gatk_ug = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk2, '-T', 'UnifiedGenotyper', '-R', fasta, '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', combined_vcf, '-dcov', '2500', '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', vcf_temp1])
	gatk_ug.wait()
	##split multi-allelic variants calls in separate lines, and left normalize indels
	bgzip_run = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'v', '-o', gatk_norm_vcf, vcf_temp2])
	bcf_norm2.wait()

def annotate_filtered_pisces_vcf(in_vcf, sample_id):
	##run_table_annovar
	command = [table_annovar] + av_buildver + [in_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', sample_id]
	annovar = subprocess.Popen(command)
	annovar.wait()

def add_allele_counts_to_annotation(genotyped_var_vcf, in_file, outfile):
	##make dict with new genotype info
	genotyped_dict = {}
	samples = []
	with open(genotyped_var_vcf, "r") as gvv_fh:
		for line in gvv_fh:
			if line[1] != '#':
				if line[0] == '#':
					line = line.rstrip().split(delim)
					samples = line[9:]
				else:
					line = line.rstrip().split(delim)
					var = '_'.join(line[:5])
					vcf_format = line[8]
					genotypes = line[9:]
					# print genotypes
					all_alleles = []
					for genotype in genotypes:
						if genotype == './.':
							alleles = '0,0'
						else:
							alleles = genotype.split(":")[1]
						all_alleles.append(alleles)
					# alleles = [i.split(":")[1] for i in genotypes]
					genotyped_dict[var] =  all_alleles
	# for i in genotyped_dict:
	# 	print i, genotyped_dict[i]
	# print samples
	##add genotyped info to ann.txt and parse
	# '''
	with open(in_file, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				out_fh.write(delim.join(line + samples) + '\n')
	
			else:
				v = '_'.join(line[141:146])
				if v in genotyped_dict:
					out_fh.write(delim.join(line + genotyped_dict[v] + ['\n']))
				else:
					print 'variant %s not genotyped'%v
					nas = ['na'] * len(samples)
					out_fh.write(delim.join(line + nas) + '\n')
	# '''


def format_pisces_results_per_proband(ped, sample, parents, id_dict):
	sample_parents_bamfiles = [sample + '.bam'] + [p + '.bam' for p in parents]
	bamlist = sample + 'bams.list'
	pisces_vcf = ped + '_pisces/' + sample + '.vcf'
	pisces_filtered_vcf = sample + '.pisces_filtered.vcf'
	genotyped_var_vcf = sample + '.pisces.gatk_genotyes.vcf'
	ann_file = sample + '.pisces.hg38_multianno.txt'
	ann_plus_geno_file = sample + '.pisces.hg38_annotated.xls'
	filtered_file = sample + '.pisces.hg38_filtered.xls'
	##genotype all vars on parents after filtering for passed
	make_list_of_bams(sample_parents_bamfiles, bamlist)
	# filter_norm_pisces_vcf(pisces_vcf,pisces_filtered_vcf, sample)
	genotype_vars_with_gatk(bamlist, pisces_filtered_vcf, genotyped_var_vcf, sample)
	##annotate
	# annotate_filtered_pisces_vcf(pisces_filtered_vcf, sample + '.pisces')
	##add genotype data to annotation --- missing genotypes? and 
	# add_allele_counts_to_annotation(genotyped_var_vcf, ann_file, ann_plus_geno_file)
	##filter vars
	# filter_split_pisces_vars(ann_plus_geno_file, filtered_file_suffix)

def variant_calling_pisces(ped_dict, bam_suffix):
	for ped in ped_dict:
		samples = ped_dict[ped][0].keys()
		print(samples)
		##call variants individually
		# variant_calling_single_pisces(ped, samples, bam_suffix)
		##combine and change ids and split affected
		affected_samples = ped_dict[ped][1]
		parent_ids = ped_dict[ped][2]
		new_id_dict = ped_dict[ped][0]
		for affected_sample in affected_samples:
			format_pisces_results_per_proband(ped, affected_sample, parent_ids, new_id_dict)

def call_mutect2_sample_parents(sample, parents, bam_suffix, working_dir):
	##issues, no PONs, and germline resources, ? fasta 
	sample_names = [sample] + parents
	bams, normals = [], []
	for s in sample_names:
		bamfile = s + bam_suffix
		ibam = ['-I', bamfile]
		bams.extend(ibam)
	for p in parents:
		normal = ['-normal', p]
		normals.extend(normal)
	# print(bams, normals)
	##file names
	raw_vcf = sample + '.mt2.temp_raw.vcf.gz'
	filtered_vcf = sample + '.mt2.temp_filtered.vcf.gz'
	vcf_temp1 = sample + '.mt2.temp1.vcf.gz'
	vcf_temp2 = sample + '.mt2.temp2.vcf.gz'
	final_vcf = sample + '.mt2.final.vcf.gz'
	# '''
	##run mt2, use --pcr-indel-model NONE if pcr free library prep and --f1r2-tar-gz for filtering by read bias
	# mt2_cmd = [gatk4, 'Mutect2', '-R', fasta] + bams + normals + ['--germline-resource', mutect_gnomad_vcf, '--tmp-dir', working_dir, '--f1r2-tar-gz', out_prefix + 'f1r2.tar.gz', '-O', raw_vcf]
	mt2_cmd = [gatk4, 'Mutect2', '-R', fasta] + bams + normals + ['--germline-resource', mutect_gnomad_vcf, '--tmp-dir', working_dir, '-O', raw_vcf]
	run_mt2 = subprocess.Popen(mt2_cmd)
	run_mt2.wait()
	##filter calls
	##gatk FilterMutectCalls -R ref.fasta -V unfiltered.vcf -O filtered.vcf
	filter_mt2 = subprocess.Popen([gatk4, 'FilterMutectCalls', '-R', fasta, '-V', raw_vcf, '-O', filtered_vcf])
	filter_mt2.wait()
	##get filtered
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-o', vcf_temp1, '-O', 'z', filtered_vcf])
	bcftools_filter.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()
	# '''
	


def variant_calling_mutect2(ped_dict, bam_suffix, work_dir):
	for ped in ped_dict:
		##call variants individually
		affected_samples = ped_dict[ped][1]
		parent_ids = ped_dict[ped][2]
		new_id_dict = ped_dict[ped][0]
		for affected_sample in affected_samples:
			##get filtered mutect calls
			call_mutect2_sample_parents(affected_sample, parent_ids, bam_suffix, work_dir)
			##annotate mutect results

##master method
def run_analysis(working_dir, analysis_file):
	os.chdir(working_dir)
	project = analysis_file.split('.')[0]
	##read in analysis file and make 2 dicts
	analysis_dict = make_analysis_dict(analysis_file)
	##check...
	for p in analysis_dict:
		print(p, analysis_dict[p])
	##call pisces and filter etc
	variant_calling_pisces(analysis_dict, '.bam')
	##call mutect and filter etc
	# variant_calling_mutect2(analysis_dict, '.bam', working_dir)






##run methods
work_dir = '/home/atimms/ngs_data/genomes/ghayda_ucb_0320'
##test on two peds
sample_file = 'ghayda_ucb_0320_1.txt'
##all samples
# sample_file = 'ghayda_ucb_0320_all.txt'
##run master method
# run_analysis(work_dir, sample_file)


##run mutect and pisces steps by ped
##040720 mutect
# sample_file = 'LR10-246.txt' 
# sample_file = 'LR12-101.txt' #273096 050420
# sample_file = 'LR12-123.txt'
# sample_file = 'LR12-249.txt'
# sample_file = 'LR12-257.txt'
# sample_file = 'LR12-259.txt'
# sample_file = 'LR12-265.txt'
# sample_file = 'LR12-269.txt'
# sample_file = 'LR13-282.txt'
# sample_file = 'LR14-155.txt'
# sample_file = 'LR16-214.txt'
## to do (swtich to mutect and change mem)
# sample_file = 'LR16-432.txt'
# sample_file = 'LR17-337.txt'
# sample_file = 'LR17-408.txt'
##run master method
# run_analysis(work_dir, sample_file)

##for pisces formatting
##273267
# sample_file = 'LR10-246.txt'
# run_analysis(work_dir, sample_file)
# sample_file = 'LR12-101.txt'
# run_analysis(work_dir, sample_file)
# sample_file = 'LR12-123.txt'
# run_analysis(work_dir, sample_file)
##273268
# sample_file = 'LR12-249.txt'
# run_analysis(work_dir, sample_file)
# sample_file = 'LR12-257.txt'
# run_analysis(work_dir, sample_file)
# sample_file = 'LR12-259.txt'
# run_analysis(work_dir, sample_file)
##273269
# sample_file = 'LR12-265.txt'
# run_analysis(work_dir, sample_file)
# sample_file = 'LR12-269.txt'
# run_analysis(work_dir, sample_file)
##273270
# sample_file = 'LR13-282.txt'
# run_analysis(work_dir, sample_file)
# sample_file = 'LR14-155.txt'
# run_analysis(work_dir, sample_file)
# sample_file = 'LR16-214.txt'
# run_analysis(work_dir, sample_file)
##273271
# sample_file = 'LR16-432.txt'
# run_analysis(work_dir, sample_file)
# sample_file = 'LR17-337.txt'
# run_analysis(work_dir, sample_file)
# sample_file = 'LR17-408.txt'
# run_analysis(work_dir, sample_file)

##rpts with no nct
##273428
# sample_file = 'LR12-249.txt'
# run_analysis(work_dir, sample_file)
##273489
# sample_file = 'LR12-269.txt'
# run_analysis(work_dir, sample_file)
##
sample_file = 'LR12-101.txt'
run_analysis(work_dir, sample_file)

