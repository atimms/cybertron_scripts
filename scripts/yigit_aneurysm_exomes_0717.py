#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'
# working_dir = '/home/atimms/ngs_data/exomes/aneurysm_exomes_0717'
working_dir = '/home/atimms/ngs_data/exomes/aneurysm_exomes_1017'
os.chdir(working_dir)
##file names etc
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.b37.vcf'
indels_1000g = ref_dir + '1000G_phase1.indels.b37.vcf'
rtc_intervals = ref_dir + 'Mills_1000G.intervals'
dbsnp = ref_dir + 'dbsnp_138.b37.vcf'
exome_bed_for_coverage = ref_dir + 'dobyns_exome.in_all_targets.1015.bed'
exome_capture_bed = ref_dir + 'dobyns_exome.in_any_padded_target.1015.bed'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
pisces_fastq_dir = '/home/atimms/ngs_data/references/illumina_human_g1k_v37/'
mh_ref = '/home/atimms/programs/MosaicHunter/resources/'
mh_rpts = mh_ref + 'all_repeats.b37.bed'
mh_common_site = mh_ref + 'WES_Agilent_50M.error_prone.b37.bed'
mh_dbsnp = mh_ref + 'dbsnp_137.b37.tsv'

##programs
samtools = 'samtools'
bedtools = 'bedtools'
picard = 'picard'
bwa = 'bwa'
bcftools = 'bcftools'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
bgzip = 'bgzip'
plink = '/home/atimms/programs/plink'
mosaic_hunter = '/home/atimms/programs/MosaicHunter/build/mosaichunter.jar'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
# gatk = '/cm/shared/apps/GATK/GenomeAnalysisTK.jar'
gatk = '/home/atimms/programs/GenomeAnalysisTK.jar'
pisces = '/cm/shared/apps/Pisces/5.1.6.54/Pisces.exe'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
vardictjava = '/home/atimms/programs/VarDictJava/build/install/VarDict/bin/VarDict'
vardict_maf_req = '0.01'

#annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147']
# av_protocol_pisces = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147,vcf', '-vcfdbfile']
av_operation = ['-operation', 'g,r,r,f,f,f']
# av_operation_pisces = ['-operation', 'g,r,r,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.','-arg', '-splicing 10 ,,,,,']
av_options_vcf = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput','-arg', '-splicing 10 ,,,,,']
# av_options_vcf_pisces = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput','-arg', '-splicing 10 ,,,,,,']

def convert_bam_fastq_picard(bamfile, r1_fastq, r2_fastq):
	# picard_sam_fq = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
	picard_sam_fq = subprocess.Popen([picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
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
	bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
	st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
	st_sam_bam_pe.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam, '-@', '10', '-m', '10G', post_bwa_bam])
	st_sort_pe.wait()
	##mark duplicates
	# picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md = subprocess.Popen([picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
	picard_md.wait()
	##realign around indels
	gatk_ir = subprocess.Popen(['java', '-Xmx80g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam, '-L', exome_capture_bed, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
	gatk_ir.wait()
	##bqsr
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-L', exome_capture_bed, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-L', exome_capture_bed, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
	gatk_pr.wait()

def variant_calling_paired_mosaichunter(pro_name, bams, delete_clean_bams):
	outfiles = []
	##clean bam files
	clean_bams = []
	for bam in bams:
		sample_name = bam.split('.')[0]
		clean_bam = sample_name + '.temp_clean.bam'
		clean_bams.append(clean_bam)
		##clean the bam
		# '''
		##samtools view -h -f 0x2 input.bam | perl -ne 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))' | samtools view -Sb - >cleaner.bam
		get_paired = subprocess.Popen(['samtools', 'view', '-h', '-f', '0x2', bam], stdout=subprocess.PIPE)
		rm_mismatched = subprocess.Popen(['perl', '-ne', 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))'], stdin=get_paired.stdout, stdout=subprocess.PIPE)
		st_view = subprocess.Popen(['samtools', 'view','-b', '-o', clean_bam], stdin=rm_mismatched.stdout)
		st_view.wait()
		bgzip_run = subprocess.Popen(['samtools', 'index', clean_bam])
		bgzip_run.wait()
		# '''
	print 'clean_bams:', clean_bams
	##run mosaic hunter in trio 'mode'
	##_mh_paired_naive
	outdir = pro_name + '_mh_paired_naive'
	run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[1], '-P', 'mosaic_filter.control_bam_file=' + 
			clean_bams[0], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=paired_naive',
			 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
			 '-P', 'output_dir=' + outdir, '-P', 'depth_filter.min_depth=10', '-P', 'min_minor_allele_number=2', 
			 '-P', 'min_minor_allele_percentage=1', '-P', 'min_p_value=0', '-P', 'base_number_filter.min_minor_allele_number=2', '-P', 
			 'base_number_filter.min_minor_allele_percentage=1', '-P', 'strand_bias_filter.p_value_cutoff=0', '-P', 'within_read_position_filter.p_value_cutoff=0'])
	run_mh_trio.wait()
	outfiles.append(outdir + '/final.passed.tsv')
	##_mh_paired_naive
	outdir = pro_name + '_mh_paired_fisher'
	run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[1], '-P', 'mosaic_filter.control_bam_file=' + 
			clean_bams[0], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=paired_fisher',
			 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
			 '-P', 'output_dir=' + outdir, '-P', 'depth_filter.min_depth=10', '-P', 'min_minor_allele_number=2', 
			 '-P', 'min_minor_allele_percentage=1', '-P', 'min_p_value=0', '-P', 'base_number_filter.min_minor_allele_number=2', '-P', 
			 'base_number_filter.min_minor_allele_percentage=1', '-P', 'strand_bias_filter.p_value_cutoff=0', '-P', 'within_read_position_filter.p_value_cutoff=0'])
	run_mh_trio.wait()
	outfiles.append(outdir + '/final.passed.tsv')
	##delete clean bams
	if delete_clean_bams == 'yes':
		for b in clean_bams:
			os.remove(b)
			os.remove(b + '.bai')
			print 'deleting', b, b + '.bai'
		return outfiles


def variant_calling_vardict_paired(ped_name, bam_files):
	##<path_to_vardict_folder>/build/install/VarDict/bin/VarDict -G /path/to/hg19.fa -f $AF_THR -N tumor_sample_name -b "/path/to/tumor.bam|/path/to/normal.bam" -z -F -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | VarDict/testsomatic.R | VarDict/var2vcf_paired.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR
	for bam in bam_files:
		##check bai file exists and make if it isn't there
		if os.path.isfile(bam + '.bai'):
			print 'bam %s alreaded indexed'%bam
		else:
			print 'indexing bam file:', bam
			st_index = subprocess.Popen(['samtools', 'index', bam])
			st_index.wait()
	##make list of vcfs
	vcf_file = ped_name + '.vardict.vcf'
	out_name = ped_name + '.vd'
	with open(vcf_file, 'w') as out_fh:
		# run_vardict = subprocess.Popen([vardictjava, '-G', fasta, '-f', vardict_maf_req, '-N', out_name, '-b', bam_files[1] + '|' + bam_files[0], '-c', '1', '-S', '2', '-E', '3', exome_capture_bed], stdout=subprocess.PIPE)
		# run_teststrandbias = subprocess.Popen(['teststrandbias.R'], stdin=run_vardict.stdout, stdout=subprocess.PIPE)
		# run_var2vcf = subprocess.Popen(['var2vcf_valid.pl', '-N', out_name, '-E', '-f', vardict_maf_req], stdin=run_teststrandbias.stdout, stdout=out_fh)
		# run_var2vcf.wait()
		run_vardict = subprocess.Popen([vardictjava, '-G', fasta, '-f', vardict_maf_req, '-N', out_name, '-b', bam_files[1] + '|' + bam_files[0], '-c', '1', '-S', '2', '-E', '3', exome_capture_bed], stdout=subprocess.PIPE)
		run_teststrandbias = subprocess.Popen(['testsomatic.R'], stdin=run_vardict.stdout, stdout=subprocess.PIPE)
		run_var2vcf = subprocess.Popen(['var2vcf_paired.pl', '-N', out_name, '-f', vardict_maf_req], stdin=run_teststrandbias.stdout, stdout=out_fh)
		run_var2vcf.wait()
	'''
	##filter 
	filter_vardict_vcf(temp_vcf, vcf)
	##copy parents vcf to annovar ref dir
	cp_vcf = subprocess.Popen(['cp', parents_vcf, str(av_ref_dir[0])])
	cp_vcf.wait()
	##annoate results
	annotate_vardict_output(vcfs[:-2], parents_vcf, bam_files[-2:])
	'''

def combine_vcf_files_gatk(vcfs, out_vcf):
	vcfs_with_v = []
	for vcf in vcfs:
		vcf_with_v = ['-V', vcf]
		vcfs_with_v.extend(vcf_with_v)
	print vcfs_with_v
	# combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fa_file, '-nt', '15', '--variant', vcfs[0],'--variant', vcfs[1],'-o', out_vcf, '-genotypeMergeOptions', 'UNIQUIFY'])
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '15'] + vcfs_with_v + ['-o', out_vcf, '-genotypeMergeOptions', 'REQUIRE_UNIQUE'])
	combine_var.wait()

def make_vcf_with_vars_just_in_tumor(in_vcf, out_vcf):
	with open(in_vcf, "r") as in_fh, open(out_vcf, "w") as out_fh:
		line_count, passed_count = 0,0
		for line in in_fh:
			if line[0] == '#':
				out_fh.write(line)
			else:
				line_count += 1
				line = line.rstrip().split(delim)
				filter_col = line[6]
				gl_genotye = line[9]
				an_genotype = line[10]
				if filter_col == 'PASS' and gl_genotye == './.' and an_genotype != './.':
					passed_count += 1
					out_fh.write(delim.join(line) + '\n')
		print in_vcf, line_count, passed_count

def make_list_of_bams(bam_files, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')

def genotype_vars_with_gatk(bamlist, in_vcf, out_vcf, ped):
	vcf_temp1 = ped + 'temp_21.vcf'
	vcf_temp2 = ped + 'temp_22.vcf.gz'
	gatk_ug = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'UnifiedGenotyper', '-R', fasta, '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', in_vcf, '-dcov', '2500', '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', vcf_temp1])
	gatk_ug.wait()
	##split multi-allelic variants calls in separate lines, and left normalize indels
	bgzip_run = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'v', '-o', out_vcf, vcf_temp2])
	bcf_norm2.wait()

def add_gatk_genotypes_to_ann_txt(in_file, out_file, filter_file, gatk_vcf):
	##make dict with new genotype info
	ann_temp = out_file + '.temp.xls'
	genotyped_dict = {}
	with open(gatk_vcf, "r") as gvv_fh:
		for line in gvv_fh:
			if line[0] != '#':
				line = line.rstrip().split(delim)
				var = '_'.join(line[:5])
				genotypes = line[9:]
				# print genotypes

				##germline_reads
				gl_genotye = genotypes[0]

				if gl_genotye == './.':
					gl_alt_reads = 'not_covered'
					gl_aaf = 'na'
				elif len(gl_genotye.split(':')[1].split(',')) == 1:
					gl_alt_reads = 'not_covered'
					gl_aaf = 'not_covered'
				elif len(gl_genotye.split(':')[1].split(',')) >2:
					gl_alt_reads = 'multi-allelic'
					gl_aaf = 'multi-allelic'
				else:
					gl_ref_reads = float(gl_genotye.split(':')[1].split(',')[0])
					gl_alt_reads = int(gl_genotye.split(':')[1].split(',')[1])
					gl_aaf = gl_alt_reads / (gl_ref_reads + gl_alt_reads)
				##tumor reads
				an_genotye = genotypes[1]
				if an_genotye == './.':
					an_alt_reads = 'not_covered'
					an_aaf = 'na'
				elif len(an_genotye.split(':')[1].split(',')) == 1:
					an_alt_reads = 'not_covered'
					an_aaf = 'not_covered'
				elif len(an_genotye.split(':')[1].split(',')) > 2:
					an_alt_reads = 'multi-allelic'
					an_aaf = 'multi-allelic'
				else:
					an_ref_reads = float(an_genotye.split(':')[1].split(',')[0])
					an_alt_reads = int(an_genotye.split(':')[1].split(',')[1])
					an_aaf = an_alt_reads / (an_ref_reads + an_alt_reads)
				genotyped_dict[var] = genotypes + [str(gl_aaf), str(gl_alt_reads), str(an_aaf), str(an_alt_reads)]
	# for i in genotyped_dict:
	# 	print i, genotyped_dict[i]
	##add genotyped info to ann.txt and parse
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				out_fh.write(delim.join(line[:66] + ['germline_pisces', 'aneurysm_pisces', 'germline_genotyped', 'aneurysm_genotypes', 'germline_aaf', 'germline_alt_reads', 'aneurysm_aaf', 'aneurysm_alt_reads', '\n']))
			else:
				v = '_'.join(line[69:74])
				if v in genotyped_dict:
					out_fh.write(delim.join(line[:66] + line[78:80] + genotyped_dict[v] + ['\n']))
				else:
					print 'varaint %s not genotyped'%v
					out_fh.write(delim.join(line[:66] + line[78:80] + ['not_genotyped','not_genotyped', 'not_genotyped', 'not_genotyped', 'not_genotyped', 'not_genotyped', '\n']))
	with open(out_file, "r") as in_fh, open(filter_file, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				out_fh.write(delim.join(line + ['\n']))
			else:
				an_alt_reads = line[73]
				gl_alt_reads = line[71]
				##rmove ungenotype and no coverage
				if an_alt_reads != 'not_covered' and an_alt_reads != 'not_genotyped' and gl_alt_reads != 'not_covered' and gl_alt_reads != 'not_genotyped':
					##if any reads in gl and must be 2 or reads in anuresym
					if int(an_alt_reads) >= 2 and int(gl_alt_reads) <1:
						##not in repeat region
						if line[10] == '.' and line[11] == '.':
							out_fh.write(delim.join(line + ['\n']))


def annotate_pisces_output(proband_vcf, gatk_genos):
	##add parent vcf to annovar parameters
	proband = proband_vcf.split('.')[0]
	print proband_vcf, proband
	# '''
	##annotate
	command = [table_annovar] + av_buildver + [proband_vcf] + av_ref_dir + av_protocol + av_operation + av_options_vcf + ['-out', proband]
	annovar = subprocess.Popen(command)
	annovar.wait()
	# '''
	multianno = proband + '.hg19_multianno.txt'
	multianno_vcf = proband + '.hg19_multianno.vcf'
	annotated = proband + '.pisces.xls'
	annotated_filtered = proband + '.pisces_filtered.xls'
	print multianno, annotated
	add_gatk_genotypes_to_ann_txt(multianno, annotated, annotated_filtered, gatk_genos)


def variant_calling_pisces(ped_name, bam_files):
	out_dir = ped_name + '_pisces'
	vcfs = []
	combined_vcf = ped_name + '.combined_pisces.vcf'
	aneur_var_vcf = ped_name + '.vars_in_an.pisces.vcf'
	gatk_geno_vcf = ped_name + '.gatk_genotyed.pisces.vcf'
	# '''
	for bam in bam_files:
		##check bai file exists and make if it isn't there
		if os.path.isfile(bam + '.bai'):
			print 'bam %s alreaded indexed'%bam
		else:
			print 'indexing bam file:', bam
			st_index = subprocess.Popen(['samtools', 'index', bam])
			st_index.wait()
		##make list of vcfs
		vcf = out_dir + '/' + bam.rsplit('.', 1)[0] + '.vcf'
		vcfs.append(vcf)
	print vcfs
	##run picses on all bam

	run_pisces = subprocess.Popen(['mono', pisces, '-G', pisces_fastq_dir, '-BamPaths', ','.join(bam_files), '-MinVF', '0.01', '-i', exome_capture_bed, '-t', '16', '-OutFolder', out_dir])
	run_pisces.wait()
	##combine vcf files
	combine_vcf_files_gatk(vcfs, combined_vcf)
	##filter vcf to get only variants that appear in the anurism
	filer_vcf = make_vcf_with_vars_just_in_tumor(combined_vcf, aneur_var_vcf)
	make_list_of_bams(bam_files, 'bams.list')
	genotype_vars_with_gatk('bams.list', aneur_var_vcf, gatk_geno_vcf, ped_name)
	# '''
	##annaotate tumor and the combine with 
	annotate_pisces_output(aneur_var_vcf, gatk_geno_vcf)


##0717 batch
'''
##bam_dict so bam files for each test so normal and tumor
bam_dict = {'468':['Ferreira.Undefined.100760.Undefined.bina.alignment.recal.bam', 'Ferreira.Undefined.100768.Undefined.bina.alignment.recal.bam'], 
		'498a':['Ferreira.Undefined.100762.Undefined.bina.alignment.recal.bam', 'Ferreira.Undefined.100770.Undefined.bina.alignment.recal.bam'],
		'498b':['Ferreira.Undefined.100762.Undefined.bina.alignment.recal.bam', 'Ferreira.Undefined.1499.Undefined.bina.alignment.recal.bam']}
# bam_dict = {'468':['Ferreira.Undefined.100760.Undefined.bina.alignment.recal.bam', 'Ferreira.Undefined.100768.Undefined.bina.alignment.recal.bam']}


##run methods
for ped_name in bam_dict:
	bam_files = []
	for bam in bam_dict[ped_name]:
		sample = bam.split('.')[2]
		read1_fastq = sample + '.r1.fastq'
		read2_fastq = sample + '.r2.fastq'
		final_bam = sample + final_bam_suffix
		bam_files.append(final_bam)
		# convert_bam_fastq_picard(bam, read1_fastq, read2_fastq)
		# align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	# variant_calling_paired_mosaichunter(ped_name, bam_files, 'yes')
	# variant_calling_vardict_paired(ped_name, bam_files)
	variant_calling_pisces(ped_name, bam_files)
'''

##1017 batch --- 21094
# '''
##bam_dict so bam files for each test so normal and tumor
bam_dict = {'485':['Ferreira.Undefined.100761.Undefined.bina.alignment.recal.bam', 'Ferreira.Undefined.100769.Undefined.bina.alignment.recal.bam'], 
		'507':['Ferreira.Undefined.100763.Undefined.bina.alignment.recal.bam', 'Ferreira.Undefined.100772.Undefined.bina.alignment.recal.bam'],
		'509':['Ferreira.Undefined.100764.Undefined.bina.alignment.recal.bam', 'Ferreira.Undefined.100773.Undefined.bina.alignment.recal.bam'],
		'VH912':['17-90076s2.bwa.gatk.bam', '17-90076s1.bwa.gatk.bam']}
# bam_dict = {'507':['Ferreira.Undefined.100763.Undefined.bina.alignment.recal.bam', 'Ferreira.Undefined.100772.Undefined.bina.alignment.recal.bam']}


##run methods 
for ped_name in bam_dict:
	bam_files = []
	for bam in bam_dict[ped_name]:
		if bam.split('.')[0] == 'Ferreira':
			sample = bam.split('.')[2]
		else:
			sample = bam.split('.')[0]
		read1_fastq = sample + '.r1.fastq'
		read2_fastq = sample + '.r2.fastq'
		final_bam = sample + final_bam_suffix
		bam_files.append(final_bam)
		convert_bam_fastq_picard(bam, read1_fastq, read2_fastq)
		align_with_bwa_one_at_time(sample, read1_fastq, read2_fastq)
	variant_calling_paired_mosaichunter(ped_name, bam_files, 'yes')
	variant_calling_vardict_paired(ped_name, bam_files)
	variant_calling_pisces(ped_name, bam_files)
# '''


