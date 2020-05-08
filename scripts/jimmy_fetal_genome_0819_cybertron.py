#!/usr/bin/python
import os
import subprocess
import glob
import shutil

##set input variables and parameters
delim = '\t'

'''
load these:
module load java/1.8.0_121 
module load biobuilds/2016.11
module load vt/0.5772-60f436c3 
module load mono/5.10.1.47
module load Pisces/5.1.6.54

'''

##programs
gatk = '/cm/shared/apps/GATK/3.7/GenomeAnalysisTK.jar'
StrelkaGermlineWorkflow = '/home/atimms/programs/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py'
plink = '/home/atimms/programs/plink'
mosaic_hunter = '/home/atimms/programs/MosaicHunter/build/mosaichunter.jar'
# gatk4 = '/home/atimms/programs/gatk-4.1.3.0/gatk'
gatk4 = '/home/atimms/programs/gatk-4.1.4.1/gatk'
pisces = '/cm/shared/apps/Pisces/5.1.6.54/Pisces.exe'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
ann_var = '/home/atimms/programs/annovar_1019/annotate_variation.pl'
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'

##files
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
bedtools_genome_file = ref_dir + 'human_g1k_v37.fasta.genome'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.b37.vcf'
indels_1000g = ref_dir + '1000G_phase1.indels.b37.vcf'
rtc_intervals = ref_dir + 'Mills_1000G.intervals'
dbsnp = ref_dir + 'dbsnp_138.b37.vcf'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
mh_ref = '/home/atimms/programs/MosaicHunter/resources/'
mh_rpts = mh_ref + 'all_repeats.b37.bed'
mh_common_site = mh_ref + 'WGS.error_prone.b37.bed'
mh_dbsnp = mh_ref + 'dbsnp_137.b37.tsv'
# mutect_gnomad_vcf = ref_dir + 'af-only-gnomad.raw.sites.b37.vcf.gz'
mutect_gnomad_vcf = ref_dir + 'small_exac_common_3_b37.vcf.gz'
pisces_fastq_dir = '/home/atimms/ngs_data/references/illumina_human_g1k_v37/'


##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,gnomad211_genome,avsnp147']
av_operation = ['-operation', 'g,r,r,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']

##methods

def combine_fq_files(work_dir, r1_to_combine, r2_to_combine, r1_fq, r2_fq):
	os.chdir(work_dir)
	print r1_fq, r1_to_combine
	print r2_fq, r2_to_combine
	with open(r1_fq, 'w') as r1_fh:
		cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
		cat_files.wait()
	with open(r2_fq, 'w') as r2_fh:
		cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
		cat_files.wait()


def align_with_bwa_one_at_time(sample_dict, working_dir):
	for sample in sample_dict:
		print(sample, sample_dict[sample])
		r1_fq = sample_dict[sample][0]
		r2_fq = sample_dict[sample][1]
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
		post_bwa_bam = sample + '.bwa.bam'
		sort_bam = sample + '.bwa_sort.bam'
		mkdup_bam = sample + '.bwa_mkdup.bam'
		realigned_bam = sample + '.bwa_religned.bam'
		gatk_bam = sample + final_bam_suffix
		# mkdup_bam = sample + '.bwa_mkdup.bam'
		##bwa and convert to bam
		# '''
		bwa_pe = subprocess.Popen(['bwa', 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen(['samtools', 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()
		##sort bam
		st_sort_pe = subprocess.Popen(['samtools', 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam, '-@', '10', '-m', '10G', post_bwa_bam])
		st_sort_pe.wait()
		# '''
		##mark duplicates
		# picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
		picard_md = subprocess.Popen(['picard', 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics', 'TMP_DIR=' + working_dir])
		picard_md.wait()
		# '''
		##realign around indels
		gatk_ir = subprocess.Popen(['java', '-Xmx80g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
		gatk_ir.wait()
		##bqsr
		gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
		gatk_br.wait()
		gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
		gatk_pr.wait()


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
	with open(multianno, "r") as multi_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in multi_fh:
			line_count += 1
			line = line.split(delim)
			if line_count ==1:
				head = line[:30] + ['chr', 'pos', 'id', 'ref2', 'alt2', 'qual', 'filter', 'info', 'format'] + samples
				out_fh.write(delim.join(head) + '\n')
			else:
				out_fh.write(delim.join(line[:30] + line[33:]))


##calls all annovar methods
def annotate_vcf_file(vcf, project_prefix):
	table_annovar_vcf(vcf, project_prefix)
	annovar_vcf = project_prefix + '.hg19_multianno.vcf'
	format_avinput(project_prefix, annovar_vcf)



def calculate_genome_coverage(samples, bam_suffix, genome_file, prefix, work_dir):
	os.chdir(work_dir)
	##parameters and header
	coverage_required = [1,5,10,20,30,40,50,100,200]
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
			##just the file
			# bedtools_cov = subprocess.Popen([bedtools, 'genomecov', '-ibam', bam, '-g', genome_file], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen(['bedtools', 'genomecov', '-ibam', bam, '-g', genome_file], stdout=subprocess.PIPE)
			grep_cmd = subprocess.Popen(['grep', '^genome'], stdin=bedtools_cov.stdout, stdout=hist_fh)
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



def variant_calling_single_mosaichunter(samples, bam_suffix):
	for sample in samples:
		##sex for input into MH, works as fetus is female
		if sample == 'Dad':
			sex = 'M'
		else:
			sex = 'F'
		in_bam = sample + bam_suffix
		clean_bam = sample + '.temp_clean.bam'
		##clean the bam
		##samtools view -h -f 0x2 input.bam | perl -ne 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))' | samtools view -Sb - >cleaner.bam
		# '''
		get_paired = subprocess.Popen(['samtools', 'view', '-h', '-f', '0x2', in_bam], stdout=subprocess.PIPE)
		rm_mismatched = subprocess.Popen(['perl', '-ne', 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))'], stdin=get_paired.stdout, stdout=subprocess.PIPE)
		st_view = subprocess.Popen(['samtools', 'view','-b', '-o', clean_bam], stdin=rm_mismatched.stdout)
		st_view.wait()
		st_index_run = subprocess.Popen(['samtools', 'index', clean_bam])
		st_index_run.wait()
		print 'made clean_bam:', clean_bam
		# '''
		##run mosaic hunter in single 'mode'
		##default
		outdir = sample + '_mh_10_600'
		run_mh_single = subprocess.Popen(['java', '-jar', mosaic_hunter, 'genome', '-P', 'input_file=' + clean_bam, '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=single',
				 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
				 '-P', 'output_dir=' + outdir, '-P', 'base_number_filter.min_minor_allele_number=2', '-P', 'base_number_filter.min_minor_allele_percentage=1', '-P', 'mosaic_filter.sex', sex, 
				 '-P', 'min_depth=10', '-P', 'max_depth=600', '-P', 'depth_filter.max_depth=600', '-P', 'depth_filter.min_depth=10'])
		run_mh_single.wait()
		##delete clean bam
		'''
		os.remove(clean_bam)
		os.remove(clean_bam + '.bai')
		print 'deleting', clean_bam, clean_bam + '.bai'
		'''

def variant_calling_trio_mosaichunter(samples, pro_sex):
	##presumes clean bams are already made
	for sample in samples:
		if sample != 'Mom' and sample != 'Dad':
			clean_bams = [sample + '.temp_clean.bam', 'Dad.temp_clean.bam', 'Mom.temp_clean.bam']
			##default
			outdir = sample + '_mh_trio_results_default'
			run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'genome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
					clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
					 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
					 '-P', 'output_dir=' + outdir, '-P', 'mosaic_filter.sex', pro_sex, '-P', 'min_depth=10', '-P', 'max_depth=600', '-P', 'depth_filter.max_depth=600', '-P', 'depth_filter.min_depth=10'])
			run_mh_trio.wait()
			##all 3
			outdir = sample + '_mh_trio_results_all3'
			run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'genome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
					clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
					 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
					 '-P', 'output_dir=' + outdir, '-P', 'depth_filter.min_depth=10', '-P', 'min_minor_allele_number=2', 
					 '-P', 'min_minor_allele_percentage=1', '-P', 'min_p_value=0', '-P', 'base_number_filter.min_minor_allele_number=2', '-P', 
					 'base_number_filter.min_minor_allele_percentage=1', '-P', 'strand_bias_filter.p_value_cutoff=0', '-P', 'within_read_position_filter.p_value_cutoff=0', '-P', 'mosaic_filter.sex', pro_sex, 
					 '-P', 'min_depth=10', '-P', 'max_depth=600', '-P', 'depth_filter.max_depth=600', '-P', 'depth_filter.min_depth=10'])
			run_mh_trio.wait()



def variant_calling_mutect2(samples, bam_suffix, working_dir, out_prefix):
	##issues, no PONs
	## from https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php#--native-pair-hmm-threads and 
	## https://gatkforums.broadinstitute.org/gatk/discussion/24057/how-to-call-somatic-mutations-using-gatk4-mutect2#latest
	bams, tbams = [], []
	for sample in samples:
		ibam = ['-I', sample + bam_suffix]
		bams.extend(ibam)
		if sample != 'Mom' and sample != 'Dad':
			tbams.extend(ibam)
	##get vars
	all_chrs_to_analyze = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
	chr_to_analyze = '21'
	##vcf names
	# raw_vcf = out_prefix + '.' + chr_to_analyze + '.mt2.raw.vcf.gz'
	all_filtered_vcf = out_prefix + '.mt2.filtered.vcf.gz'
	out_vcf = out_prefix + '.mt2.vcf.gz'
	vcf_temp = out_prefix + '.mt2.temp.vcf.gz'
	final_vcf = out_prefix + '.mt2.norm.vcf.gz'
	'''
	##run manually per chromosome and then combine and move to next steps
	##run mt2, using --pcr-indel-model NONE as pcr free library prep and --f1r2-tar-gz for filtering by read bias
	mt2_cmd = [gatk4, 'Mutect2', '-R', fasta] + bams + ['-normal', 'Mom', '-normal', 'Dad', '--germline-resource', mutect_gnomad_vcf, '--pcr-indel-model', 'NONE',
			'--tmp-dir', working_dir, '--f1r2-tar-gz', chr_to_analyze + 'f1r2.tar.gz', '-L', chr_to_analyze, '-O', raw_vcf]
	run_mt2 = subprocess.Popen(mt2_cmd)
	run_mt2.wait()
	##merge raw vcfs
	in_vcfs = glob.glob('*' + '.mt2.raw.vcf.gz')
	
	bcf_cc = subprocess.Popen(['bcftools', 'concat'] + in_vcfs + ['-O', 'z', '-o', out_vcf])
	bcf_cc.wait()
	
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp, out_vcf])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()
	'''
	##not run as major issues !!!
	# ##get list of commands to add files for all chromosomes
	# all_f1r2_files = []
	# for chr_to_anal in all_chrs_to_analyze:
	# 	f1r2 = ['-I', chr_to_anal + 'f1r2.tar.gz']
	# 	all_f1r2_files.extend(f1r2)
	# ##LearnReadOrientationModel using all the f1r2.tar.gz, from https://software.broadinstitute.org/gatk/documentation/article?id=24057
	# #gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
	# run_lrom = subprocess.Popen([gatk4, 'LearnReadOrientationModel'] + all_f1r2_files + ['-O', 'read-orientation-model.tar.gz'])
	# run_lrom.wait()
	# ##these need to be done on the genome not per chromosome, so just do once
	# ##Run GetPileupSummaries to summarize read support for a set number of known variant sites, for the read oriantation
	# #gatk GetPileupSummaries -I tumor.bam -V chr17_small_exac_common_3_grch38.vcf.gz -L chr17_small_exac_common_3_grch38.vcf.gz -O getpileupsummaries.table
	# run_gps = subprocess.Popen([gatk4, 'GetPileupSummaries'] + tbams + ['-V', mutect_gnomad_vcf, '-L', mutect_gnomad_vcf, '-O', 'getpileupsummaries.table'])
	# run_gps.wait()
	# ##Estimate contamination with CalculateContamination
	# #gatk CalculateContamination -I getpileupsummaries.table -tumor-segmentation segments.table -O calculatecontamination.table
	# run_cc = subprocess.Popen([gatk4, 'CalculateContamination', '-I', 'getpileupsummaries.table', '-tumor-segmentation', 'segments.table', '-O', 'contamination.table'])
	# run_cc.wait()
	# ##Finally, pass the learned read orientation model to FilterMutectCallswith the -ob-priors argument:
	# #gatk FilterMutectCalls -V unfiltered.vcf [--tumor-segmentation segments.table] [--contamination-table contamination.table] --ob-priors read-orientation-model.tar.gz -O filtered.vcf
	# ##get list of all raw vcfs
	# raw_vcf_files = []
	# for chr_to_anal in all_chrs_to_analyze:
	# 	raw_vcf = ['-V', out_prefix + '.' + chr_to_anal + '.mt2.raw.vcf.gz']
	# 	raw_vcf_files.extend(raw_vcf)
	# run_fmc = subprocess.Popen([gatk4, 'FilterMutectCalls'] + raw_vcf_files + ['--tumor-segmentation', 'segments.table', '--contamination-table', 'contamination.table', 
	# 		'--ob-priors', 'read-orientation-model.tar.gz', '-O', filtered_vcf])
	# run_fmc.wait()
	'''
	##appempt2 -- not using LearnReadOrientationModel
	##Run GetPileupSummaries to summarize read support for a set number of known variant sites, for the read oriantation
	#gatk GetPileupSummaries -I tumor.bam -V chr17_small_exac_common_3_grch38.vcf.gz -L chr17_small_exac_common_3_grch38.vcf.gz -O getpileupsummaries.table
	run_gps = subprocess.Popen([gatk4, 'GetPileupSummaries', '--java-options', "-Xmx200g"] + tbams + ['-V', mutect_gnomad_vcf, '-L', mutect_gnomad_vcf, '-O', 'getpileupsummaries.table', '--tmp-dir', '.'])
	run_gps.wait()
	#gatk CalculateContamination -I getpileupsummaries.table -tumor-segmentation segments.table -O calculatecontamination.table
	run_cc = subprocess.Popen([gatk4, 'CalculateContamination', '-I', 'getpileupsummaries.table', '-O', 'contamination.table', '--tmp-dir', '.'])
	run_cc.wait()

	##filter using contaimination table per chrom
	in_vcfs = []
	for chrom in all_chrs_to_analyze:
		in_vcf = out_prefix + '.' + chrom + '.mt2.raw.vcf.gz'
		filtered_vcf = out_prefix + '.' + chrom + '.mt2.filtered.vcf.gz'
		in_vcfs.append(filtered_vcf)
		run_fmc = subprocess.Popen([gatk4, 'FilterMutectCalls', '-R', fasta, '-V', in_vcf, '--contamination-table', 'contamination.table', '-O', filtered_vcf])
		run_fmc.wait()
	bcf_cc = subprocess.Popen(['bcftools', 'concat'] + in_vcfs + ['-O', 'z', '-o', out_vcf])
	bcf_cc.wait()
	
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp, out_vcf])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', all_filtered_vcf, vcf_temp])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', all_filtered_vcf])
	bcf_index.wait()
	'''
	##for testing
	# out_prefix = 'test'
	##annotate vcf file
	annotate_vcf_file(all_filtered_vcf, out_prefix + '_mt2')
	##reformat (add allele counts)and filter multiple ways
	ann_txt = out_prefix + '_mt2.annotated.txt'
	reformated_txt = out_prefix + '_mt2.all_vars.xls'
	reformat_ann_txt(ann_txt, reformated_txt, [39,47], 'mutect2')
	filter_var_file(reformated_txt, out_prefix + '_mt2', 'mutect2')




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
	run_pisces = subprocess.Popen(['mono', pisces, '-G', pisces_fastq_dir, '-BamPaths', ','.join(bam_files), '-MinVF', '0.01', '-gVCF', 'FALSE', '-t', '18', '-OutFolder', out_dir])
	run_pisces.wait()


def merge_vcf_files(vcfs, out_vcf):
	vcfs_to_combine_cmd = []
	for vcf in vcfs:
		vcfs_to_combine_cmd.extend(['--variant', vcf])
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R', fasta, '-nt', '15'] + vcfs_to_combine_cmd + ['-o', out_vcf, '-genotypeMergeOptions', 'UNSORTED', '-filteredAreUncalled'])
	combine_var.wait()

def merge_pisces_results(samples, snp_comb_vcf):
	snp_vcfs_to_combine = []
	for sample in samples:
		in_vcf = 'fetal_seq_genome_0819_pisces/' + sample + '.bwa_gatk.vcf'
		snp_vcf = sample + '.pisces_passed_snps.vcf'
		snp_vcfs_to_combine.append(snp_vcf)
		##remove non passed varianst
		# '''
		snp_cut = subprocess.Popen(['java', '-Xmx200g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '15', '-V', in_vcf, '-o', snp_vcf, '--excludeFiltered','-selectType', 'SNP', '-restrictAllelesTo', 'BIALLELIC'])
		snp_cut.wait()

		# '''
	merge_vcf_files(snp_vcfs_to_combine, snp_comb_vcf)


def make_list_of_bams(samples, bam_suffix, bamlist_file):
	bam_files = [s + bam_suffix for s in samples]
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')

def genotype_vars_with_gatk(bamlist, combined_vcf, gatk_norm_vcf):
	vcf_temp1 = gatk_norm_vcf + 'temp_21.vcf'
	vcf_temp2 = gatk_norm_vcf + 'temp_22.vcf.gz'
	gatk_ug = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'UnifiedGenotyper', '-R', fasta, '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', combined_vcf, '-dcov', '2500', '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', vcf_temp1])
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

def variant_calling_pisces(samples, bam_suffix, out_prefix):
	##call variants individually
	# variant_calling_single_pisces(out_prefix, samples, bam_suffix)
	##combine
	all_sample_snp_vcf = out_prefix + '.pisces.combined.snp.vcf'
	gatk_snp_genotyes_vcf = out_prefix + '.gatk_genotyes.snp.vcf'
	# merge_pisces_results(samples, all_sample_snp_vcf)
	##genotype snps and indels
	# make_list_of_bams(samples, bam_suffix, bamlist)
	# genotype_vars_with_gatk(bamlist, all_sample_snp_vcf, gatk_snp_genotyes_vcf)
	##annotate vcf file
	# annotate_vcf_file(gatk_snp_genotyes_vcf, out_prefix + '_pisces')
	##reformat (add allele counts)and filtermultiple ways
	ann_txt = out_prefix + '_pisces.annotated.txt'
	reformated_txt = out_prefix + '_pisces.all_vars.xls'
	# reformat_ann_txt(ann_txt, reformated_txt, [39,47], 'pisces')
	filter_var_file(reformated_txt, out_prefix + '_pisces', 'pisces')


def reformat_ann_txt(inann, outfile, sample_pos, var_caller):
	##make list from bed file
	with open(outfile, "w") as out_fh, open(inann, "r") as ina_fh:
		line_count = 0
		for line_ann in ina_fh:
			line_ann = line_ann.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				header = line_ann
				sample_list = line_ann[sample_pos[0]:sample_pos[1]]
				for sample in sample_list:
					s_cov = sample + ' coverage'
					s_aaf = sample + ' aaf'
					s_alt = sample + ' alt reads'
					header.extend([s_cov, s_aaf, s_alt])
				# print header
				out_fh.write(delim.join(header + ['\n']))
				##get position of maternal tissue
				# maternal_position = [i for i, s in enumerate(sample_list) if 'Mom' in s]
				# print maternal_position
			else:
				keep_var = 'yes'
				##get all genotype
				genotypes = line_ann[sample_pos[0]:sample_pos[1]]
				format = line_ann[38].split(':')
				##line to write so var + 
				line_out = line_ann
				##get all genotype and write
				for genotype in genotypes:
					genotype_info = genotype.split(':')
					if var_caller == 'strelka':
						# print genotype_info
						if format[4] == 'AD':
							all_alleles = genotype_info[4].split(',')
						elif format[5] == 'AD':
							all_alleles = genotype_info[5].split(',')
						alleles = [all_alleles[0], all_alleles[1]]
						# print alleles
						coverage = int(alleles[0]) + int(alleles[1])
						if coverage == 0:
							allele_freq = 0
						else:
							##get aaf 
							allele_freq = float(alleles[1]) / coverage
						# print coverage, allele_freq
						line_out.extend([str(coverage), str(allele_freq), str(alleles[1])])

					elif var_caller == 'mutect2':
						all_alleles = genotype_info[1].split(',')
						alleles = [all_alleles[0], all_alleles[1]]
						##remove multiallelic
						# print(all_alleles, len(all_alleles))
						if len(all_alleles) !=2:
							keep_var == 'no'
							print(all_alleles)
						coverage = int(alleles[0]) + int(alleles[1])
						if coverage == 0:
							allele_freq = 0
						else:
							##get aaf 
							allele_freq = float(alleles[1]) / coverage
						# print coverage, allele_freq
						line_out.extend([str(coverage), str(allele_freq), str(alleles[1])])
					elif var_caller == 'pisces':
						if len(genotype_info) != 5:
							line_out.extend(['0','0','0'])
							# print(genotype_info)
						else:
							if genotype_info[1] == '.':
								line_out.extend(['0','0','0'])
							else:
								# print(genotype_info)
								all_alleles = genotype_info[1].split(',')
								alleles = [all_alleles[0], all_alleles[1]]
								##remove multiallelic
								# print(all_alleles, len(all_alleles))
								if len(all_alleles) !=2:
									keep_var == 'no'
									print(all_alleles, '>2 allles')
								coverage = int(alleles[0]) + int(alleles[1])
								if coverage == 0:
									allele_freq = 0
								else:
									##get aaf 
									allele_freq = float(alleles[1]) / coverage
								# print coverage, allele_freq
								line_out.extend([str(coverage), str(allele_freq), str(alleles[1])])
				if keep_var == 'yes':
					out_fh.write(delim.join(line_out +['\n']))

def filter_var_file(infile, out_prefix, var_caller):
	##make list from bed file
	not_in_parents = out_prefix + '.not_called_in_parents.xls'
	not_in_parents_aaf = out_prefix + '.not_called_in_parents.rare.xls'
	not_in_parents_aaf_3reads = out_prefix + '.not_called_in_parents.rare.3reads.xls'
	not_in_parents_3reads = out_prefix + '.not_called_in_parents.3reads.xls'
	lt3_in_parents = out_prefix + '.lt3_in_parents.xls'
	lt3_in_parents_aaf = out_prefix + '.lt3_in_parents.rare.xls'
	lt3_in_parents_aaf_3reads = out_prefix + '.lt3_in_parents.rare.3reads.xls'
	lt3_in_parents_3reads = out_prefix + '.lt3_in_parents.3reads.xls'
	with open(infile, "r") as in_fh, open(not_in_parents, "w") as nip_fh, open(not_in_parents_aaf, "w") as nipa_fh, \
		open(not_in_parents_aaf_3reads, "w") as nipa3_fh, open(not_in_parents_3reads, "w") as nip3_fh, \
		open(lt3_in_parents, "w") as lip_fh, open(lt3_in_parents_aaf, "w") as lipa_fh, open(lt3_in_parents_aaf_3reads, "w") as lipa3_fh, \
		open(lt3_in_parents_3reads, "w") as lip3_fh:
		line_count, passed_count, nip_count,nipa_count,nipa3_count, nip3_count, lip_count,lipa_count,lipa3_count, lip3_count = 0, 0, 0, 0, 0, 0,0,0,0,0
		for line in in_fh:
			line_count += 1
			if line_count == 1:
				nip_fh.write(line)
				nipa_fh.write(line)
				nipa3_fh.write(line)
				nip3_fh.write(line)
				lip_fh.write(line)
				lipa_fh.write(line)
				lipa3_fh.write(line)
				lip3_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				chrom = line[0]
				rmsk = line[10]
				segdup = line[11]
				ref = line[3]
				alt = line[4]
				if var_caller == 'strelka':
					all_child_coverage = int(line[59])
					dad_genotype = line[41].split(':')[0]
					mom_genotype = line[44].split(':')[0]
					dad_alt_read_count = int(line[55])
					mom_alt_read_count = int(line[64])
					sample_alt_pos = [49,52,58,67,70]
				elif var_caller == 'mutect2':
					all_child_coverage = int(line[53])
					dad_genotype = line[42].split(':')[0]
					mom_genotype = line[46].split(':')[0]
					dad_alt_read_count = int(line[58])
					mom_alt_read_count = int(line[70])
					sample_alt_pos = [49,52,61,64,67]
				elif var_caller == 'mosaichunter':
					all_child_coverage = int(line[62])
					dad_genotype = line[46].split(':')[0]
					mom_genotype = line[39].split(':')[0]
					dad_alt_read_count = int(line[70])
					mom_alt_read_count = int(line[49])
					sample_alt_pos = [52,55,58,61,67]
				elif var_caller == 'pisces':
					# print(line)
					all_child_coverage = int(line[53])
					dad_genotype = line[42].split(':')[0]
					mom_genotype = line[46].split(':')[0]
					dad_alt_read_count = int(line[58])
					mom_alt_read_count = int(line[70])
					sample_alt_pos = [49,52,61,64,67]
				gnomad_af = line[12]
				if gnomad_af == '.':
					gnomad_af = 0
				sample_read_counts = [line[i] for i in sample_alt_pos]
				sample_read_counts_int = [int(i) for i in sample_read_counts]
				##remove if in repeat region, or high/low coverage in combined child (remove indels)
				# if not chrom.startswith('GL') and rmsk == '.' and segdup == '.' and all_child_coverage >= 100 <= 1000:
				# print(rmsk, segdup, all_child_coverage, ref, alt)
				if not chrom.startswith('GL') and rmsk == '.' and segdup == '.' and all_child_coverage >= 100 and all_child_coverage <= 1000 and ref != '-' and alt != '-' and len(ref) == 1 and len(alt) == 1:
					passed_count += 1
					if mom_genotype == '0/0' and dad_genotype == '0/0':
						nip_fh.write(delim.join(line) + '\n')
						nip_count += 1
						if gnomad_af <= 0.001:
							nipa_count += 1
							nipa_fh.write(delim.join(line) + '\n')
							if max(sample_read_counts_int) >= 3:
								nipa3_count += 1
								nipa3_fh.write(delim.join(line) + '\n')
						if max(sample_read_counts_int) >= 3:
							nip3_fh.write(delim.join(line) + '\n')
							nip3_count += 1
					if (dad_alt_read_count + mom_alt_read_count) < 3:
						lip_fh.write(delim.join(line) + '\n')
						lip_count += 1
						if gnomad_af <= 0.001:
							lipa_count += 1
							lipa_fh.write(delim.join(line) + '\n')
							if max(sample_read_counts_int) >= 3:
								lipa3_count += 1
								lipa3_fh.write(delim.join(line) + '\n')
						if max(sample_read_counts_int) >= 3:
							lip3_fh.write(delim.join(line) + '\n')
							lip3_count += 1

	print(line_count, passed_count, nip_count,nipa_count,nipa3_count, nip3_count,lip_count,lipa_count,lipa3_count, lip3_count)




def variant_calling_strelka(sample_dict, bam_suffix, w_dir ,ped_name):
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
	final_vcf = ped_name + '.strelka.vcf.gz'
	# '''
	##configure
	st_config = subprocess.Popen([StrelkaGermlineWorkflow] + bam_cmds + ['--referenceFasta', fasta, '--runDir', analysis_dir])
	st_config.wait()	
	##run ?? cores
	#${STRELKA_ANALYSIS_PATH}/runWorkflow.py -m sge -j 36
	# st_run = subprocess.Popen([analysis_dir + '/runWorkflow.py', '-m', 'sge', '-j', '20'])
	st_run = subprocess.Popen([analysis_dir + '/runWorkflow.py', '-m', 'local', '-j', '10'])
	st_run.wait()
	
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
	# '''
	##for testing
	# ped_name = 'test'
	# final_vcf = ped_name + '.strelka.vcf.gz'
	##annotate vcf file
	# annotate_vcf_file(final_vcf, ped_name + '_strelka')
	##reformat (add allele counts)and filtermultiple ways
	ann_txt = ped_name + '_strelka.annotated.txt'
	reformated_txt = ped_name + '_strelka.all_vars.xls'
	# reformat_ann_txt(ann_txt, reformated_txt, [39,47], 'strelka')
	filter_var_file(reformated_txt, ped_name + '_strelka', 'strelka')


def genotype_vars_with_freebayes(bedfile, final_vcf):
	vcf_temp1 = 'temp_fb1.vcf'
	vcf_temp2 = 'temp_fb2.vcf'
	with open(vcf_temp1, 'w') as vcf_fh:
		run_freebayes = subprocess.Popen([freebayes, '-f', fasta, '-L', bamlist, '-t', bedfile, '--haplotype-length', '0', '--min-alternate-count', '1', '--min-alternate-fraction', '0', '--pooled-continuous', '--report-monomorphic'], stdout=vcf_fh)
		run_freebayes.wait()
	run_bgzip = subprocess.Popen(['bgzip', vcf_temp1])
	run_bgzip.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	os.remove( vcf_temp1 + '.gz')


def combine_and_format_mh(mh_vars, inann, outfile, sample_pos):
	mhvar_dict = {}
	with open(mh_vars, "r") as mhv_fh:
		for line in mhv_fh:
			line = line.rstrip().split(delim)
			chrom = line[2]
			pos = line[3]
			ref = line[4]
			if line[4] == line[8]:
				alt = line[10]
			else:
				alt = line[8]
			var = '_'.join([chrom,pos,ref,alt])
			anal = '_'.join(line[:2])
			if var in mhvar_dict:
				mhvar_dict[var].append(anal)
			else:
				mhvar_dict[var] = [anal]
	for v in mhvar_dict:
		print(v, mhvar_dict[v])


	##make list from bed file
	v_list = []
	with open(outfile, "w") as out_fh, open(inann, "r") as ina_fh:
		line_count = 0
		for line_ann in ina_fh:
			line_ann = line_ann.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				header = line_ann
				sample_list = line_ann[sample_pos[0]:sample_pos[1]]
				for sample in sample_list:
					s_cov = sample + ' coverage'
					s_aaf = sample + ' aaf'
					s_alt = sample + ' alt reads'
					header.extend([s_cov, s_aaf, s_alt])
				# print header
				out_fh.write(delim.join(header + ['analysis', '\n']))
				##get position of maternal tissue
				# maternal_position = [i for i, s in enumerate(sample_list) if 'Mom' in s]
				# print maternal_position
			else:
				variant = '_'.join([line_ann[0]] + line_ann[2:5])
				# print(variant)
				if variant in mhvar_dict and variant not in v_list:
					v_list.append(variant)
					analysis = ','.join(mhvar_dict[variant])
					##get all genotype
					genotypes = line_ann[sample_pos[0]:sample_pos[1]]
					format = line_ann[38].split(':')
					##line to write so var + 
					line_out = line_ann
					##get all genotype and write
					for genotype in genotypes:
						genotype_info = genotype.split(':')
						all_alleles = genotype_info[2].split(',')
						if len(all_alleles) == 2:
							alleles = [all_alleles[0], all_alleles[1]]
							coverage = int(alleles[0]) + int(alleles[1])
							##get aaf 
							allele_freq = float(alleles[1]) / coverage
							# print(all_alleles)
							line_out.extend([str(coverage), str(allele_freq), str(alleles[1])])
						else:
							print(line)
					out_fh.write(delim.join(line_out +[analysis, '\n']))


def analyze_mosaichunter_results(sample_dict, file_prefix):
	##filnames
	combined_mh_file = file_prefix + '.all_mh.xls'
	combined_mh_bed = file_prefix + '.all_mh.bed'
	freebayes_vcf = file_prefix + '.freebayes_geno.vcf'
	freebayes_annotated = file_prefix + '_mh.annotated.txt'
	reformated_txt = file_prefix + '_mosaichunter.all_vars.xls'
	##combine in a single file, and make bed
	with open(combined_mh_file, "w") as out_fh, open(combined_mh_bed, "w") as outbed_fh:
		for sample in sample_dict:
			single_mh_final = sample + '_mh_10_600/final.passed.tsv'
			trio_mh_final = sample + '_mh_trio_results_default/final.passed.tsv'
			trio_a3_mh_final = sample + '_mh_trio_results_all3/final.passed.tsv'
			if sample == 'Mom' or sample == 'Dad':
				with open(single_mh_final, "r") as smf_fh:
					for line in smf_fh:
						out_fh.write(sample + '\tsingle\t' + line)
						linesplit = line.split(delim)
						if linesplit[2] == linesplit[6]:
							alt = linesplit[8]
						else:
							alt = linesplit[6]
						bed_out = [linesplit[0], str(int(linesplit[1]) - 1), linesplit[1], linesplit[2], alt, sample, 'single']
						outbed_fh.write(delim.join(bed_out) + '\n')
			else:
				with open(single_mh_final, "r") as smf_fh, open(trio_mh_final, "r") as tmf_fh, open(trio_a3_mh_final, "r") as tamf_fh:
					for line in smf_fh:
						out_fh.write(sample + '\tsingle\t' + line)
						linesplit = line.split(delim)
						if linesplit[2] == linesplit[6]:
							alt = linesplit[8]
						else:
							alt = linesplit[6]
						bed_out = [linesplit[0], str(int(linesplit[1]) - 1), linesplit[1], linesplit[2], alt, sample, 'single']
						outbed_fh.write(delim.join(bed_out) + '\n')
					for line in tmf_fh:
						out_fh.write(sample + '\ttrio\t' + line)
						linesplit = line.split(delim)
						if linesplit[2] == linesplit[6]:
							alt = linesplit[8]
						else:
							alt = linesplit[6]
						bed_out = [linesplit[0], str(int(linesplit[1]) - 1), linesplit[1], linesplit[2], alt, sample, 'trio']
						outbed_fh.write(delim.join(bed_out) + '\n')
					for line in tamf_fh:
						out_fh.write(sample + '\ttrio_all3\t' + line)
						linesplit = line.split(delim)
						if linesplit[2] == linesplit[6]:
							alt = linesplit[8]
						else:
							alt = linesplit[6]
						bed_out = [linesplit[0], str(int(linesplit[1]) - 1), linesplit[1], linesplit[2], alt, sample, 'trio_all3']
						outbed_fh.write(delim.join(bed_out) + '\n')
	##genotype, combine all data together and filter
	# genotype_vars_with_freebayes(combined_mh_bed, freebayes_vcf)
	##annotate vars
	# annotate_vcf_file(freebayes_vcf, file_prefix + '_mh')
	##combine data and filter
	# combine_and_format_mh(combined_mh_file, freebayes_annotated, reformated_txt, [39,47])
	##filter
	filter_var_file(reformated_txt, file_prefix + '_mh', 'mosaichunter')

def intersect_var_callers(var_callers, infile_suffix, out_prefix):
	big_dict = {}
	for var_caller in var_callers:
		infile = 'fetal_seq_genome_0819_' + var_caller + infile_suffix
		with open(infile, "r") as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.rstrip().split(delim)
				if lc == 1:
					header = line
				else:
					var = '_'.join(line[:5])
					if var in big_dict:
						# big_dict[var][1].append(var_caller)
						big_dict[var][var_caller] = [header, line]
					else:
						# big_dict[var] = [line, [var_caller], header]
						big_dict[var] = {var_caller: [header, line]}
	##get counts/vars
	in_all_4, pis_mt2, pis_mt2_st, pis_mt2_mh = 0, 0, 0, 0
	in_all_4_file = out_prefix + '.in_all_4' + infile_suffix
	pis_mt2_file = out_prefix + '.pisces_mt2' + infile_suffix
	pis_mt2_st_file = out_prefix + '.pisces_mt2_strekla' + infile_suffix
	pis_mt2_mh_file = out_prefix + '.pisces_mt2_mh' + infile_suffix
	with open(in_all_4_file, "w") as ia4_fh, open(pis_mt2_file, "w") as pmt_fh, open(pis_mt2_st_file, "w") as pmts_fh, open(pis_mt2_mh_file, "w") as pmth_fh:

		for v in big_dict:
			var_caller_names = big_dict[v].keys()
			print(v, var_caller_names)

			if len(var_caller_names) == 4:
				in_all_4 += 1
				if in_all_4 == 1:
					ia4_fh.write(delim.join(big_dict[v]['mh'][0]) + '\n')
				ia4_fh.write(delim.join(big_dict[v]['mh'][1]) + '\n')
			if 'pisces' in var_caller_names and 'mt2' in var_caller_names:
				pis_mt2 += 1
				if pis_mt2 == 1:
					pmt_fh.write(delim.join(big_dict[v]['mt2'][0]) + '\n')
				pmt_fh.write(delim.join(big_dict[v]['mt2'][1]) + '\n')
			if 'pisces' in var_caller_names and 'mt2' in var_caller_names and 'strelka' in var_caller_names:
				pis_mt2_st += 1
				if pis_mt2_st == 1:
					pmts_fh.write(delim.join(big_dict[v]['mt2'][0]) + '\n')
				pmts_fh.write(delim.join(big_dict[v]['mt2'][1]) + '\n')
			if 'pisces' in var_caller_names and 'mt2' in var_caller_names and 'mh' in var_caller_names:
				pis_mt2_mh += 1
				if pis_mt2_mh == 1:
					pmth_fh.write(delim.join(big_dict[v]['mh'][0]) + '\n')
				pmth_fh.write(delim.join(big_dict[v]['mh'][1]) + '\n')
	print(in_all_4, pis_mt2, pis_mt2_st, pis_mt2_mh)


def run_all_methods(work_dir, info_dict, project):
	os.chdir(work_dir)
	##map
	# align_with_bwa_one_at_time(info_dict, work_dir)
	##mosaic hunter - single mode
	# variant_calling_single_mosaichunter(info_dict, final_bam_suffix)
	##mosaic hunter - trio mode, have clean bams from single mode
	# variant_calling_trio_mosaichunter(info_dict, 'F')
	##mosaic hunter - combine data and genotype 
	# analyze_mosaichunter_results(info_dict, project)
	##mutect2
	# variant_calling_mutect2(info_dict, final_bam_suffix, work_dir, project)
	##pisces
	# variant_calling_pisces(info_dict, final_bam_suffix, project)
	##strelka
	# variant_calling_strelka(info_dict, final_bam_suffix, work_dir, project)
	##compare variant callers, get insections etc
	var_callers = ['mh', 'mt2', 'pisces', 'strelka']
	file_suffix = '.lt3_in_parents.xls'
	intersect_var_callers(var_callers, file_suffix, project)



##run methods
working_dir = '/home/atimms/ngs_data/genomes/fetal_seq_genome_0819'
project_name = 'fetal_seq_genome_0819'
child_r1_fqs_to_combine = ['S1_S0_L001_R1_001.fastq.gz', 'S2_S0_L001_R1_001.fastq.gz', 'S3_S0_L001_R1_001.fastq.gz', 'S4_S0_L001_R1_001.fastq.gz', 'S5_S0_L001_R1_001.fastq.gz']
child_r2_fqs_to_combine = ['S1_S0_L001_R2_001.fastq.gz', 'S2_S0_L001_R2_001.fastq.gz', 'S3_S0_L001_R2_001.fastq.gz', 'S4_S0_L001_R2_001.fastq.gz', 'S5_S0_L001_R2_001.fastq.gz']
combined_name = 'Child_combined'
combined_r1_fq = combined_name + '_R1.fastq.gz'
combined_r2_fq = combined_name + '_R2.fastq.gz'

fastq_dict = {'Brain': ['S1_S0_L001_R1_001.fastq.gz', 'S1_S0_L001_R2_001.fastq.gz'], 'Cerebellum': ['S2_S0_L001_R1_001.fastq.gz', 'S2_S0_L001_R2_001.fastq.gz'], 
		'Heart': ['S3_S0_L001_R1_001.fastq.gz', 'S3_S0_L001_R2_001.fastq.gz'], 'Kidney': ['S4_S0_L001_R1_001.fastq.gz', 'S4_S0_L001_R2_001.fastq.gz'], 
		'Lung': ['S5_S0_L001_R1_001.fastq.gz', 'S5_S0_L001_R2_001.fastq.gz'], 'Mom': ['S6_S0_L001_R1_001.fastq.gz', 'S6_S0_L001_R2_001.fastq.gz'], 
		'Dad': ['S7_S0_L001_R1_001.fastq.gz', 'S7_S0_L001_R2_001.fastq.gz'], combined_name:[combined_r1_fq, combined_r2_fq]}

##combine all tissue fqs to make child combined fq
# combine_fq_files(working_dir, child_r1_fqs_to_combine, child_r2_fqs_to_combine, combined_r1_fq, combined_r2_fq)

##temp catch up
# fastq_dict = {combined_name:[combined_r1_fq, combined_r2_fq]}
##map and variant calling
run_all_methods(working_dir, fastq_dict, project_name)

##coverage
# calculate_genome_coverage(fastq_dict, '.bwa_gatk.bam', bedtools_genome_file, project_name, working_dir)


