#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'


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


##programs
mosaic_hunter = '/home/atimms/programs/MosaicHunter/build/mosaichunter.jar'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
gatk = '/cm/shared/apps/GATK/GenomeAnalysisTK.jar'
pisces = '/cm/shared/apps/Pisces/5.1.6.54/Pisces.exe'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
vardictjava = '/home/atimms/programs/VarDictJava/build/install/VarDict/bin/VarDict'
vardict_maf_req = '0.01'


#annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147']
av_protocol_pisces = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp30a,popfreq_all_20150413,avsnp147,vcf', '-vcfdbfile']
av_operation = ['-operation', 'g,r,r,f,f,f']
av_operation_pisces = ['-operation', 'g,r,r,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.','-arg', '-splicing 10 ,,,,,']
av_options_vcf = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput','-arg', '-splicing 10 ,,,,,']
av_options_vcf_pisces = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput','-arg', '-splicing 10 ,,,,,,']

def variant_calling_paired_mutect2(bams, name_prefix):
	tumor_bam = bams[0]
	normal_bam = bams[1]
	vcf_temp1 = name_prefix + '.temp_1.vcf'
	vcf_temp2 = name_prefix + '.temp_2.vcf.gz'
	final_vcf = name_prefix + '.mutect.vcf.gz'
	##run mutect caller
	call_mutect = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'MuTect2', '-R', fasta, '-I:tumor', tumor_bam, '-I:normal', normal_bam, '-L', exome_capture_bed, '-o', vcf_temp1])
	call_mutect.wait()
	bgzip_run = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()

##annotate vcf file
def table_annovar_vcf(vcf, project_prefix):
		out_prefix = project_prefix
		command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()


def variant_calling_trio_mosaichunter(pro_name, bams, pro_sex, pro_cov):
	outfiles = []
	##clean bam files
	clean_bams = []
	for bam in bams:
		sample_name = bam.split('.')[0]
		clean_bam = sample_name + '.temp_clean.bam'
		clean_bams.append(clean_bam)
		##clean the bam
		##samtools view -h -f 0x2 input.bam | perl -ne 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))' | samtools view -Sb - >cleaner.bam
		get_paired = subprocess.Popen(['samtools', 'view', '-h', '-f', '0x2', bam], stdout=subprocess.PIPE)
		rm_mismatched = subprocess.Popen(['perl', '-ne', 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))'], stdin=get_paired.stdout, stdout=subprocess.PIPE)
		st_view = subprocess.Popen(['samtools', 'view','-b', '-o', clean_bam], stdin=rm_mismatched.stdout)
		st_view.wait()
		bgzip_run = subprocess.Popen(['samtools', 'index', clean_bam])
		bgzip_run.wait()
	print 'clean_bams:', clean_bams
	# clean_bams = ['/data/atimms/dobyns_mosaic_test_0317/LR16-173_avm_lesion.temp_clean.bam', '/data/atimms/dobyns_mosaic_test_0317/LR16-173f.temp_clean.bam', '/data/atimms/dobyns_mosaic_test_0317/LR16-173m.temp_clean.bam']
	##run mosaic hunter in trio 'mode'
	##default
	outdir = pro_name + '_mh_trio_results_default'
	run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
			clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
			 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
			 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov])
	run_mh_trio.wait()
	outfiles.append(outdir + '/final.passed.tsv')
	# ##lower cov
	# outdir = pro_name + '_mh_trio_results_lowcov'
	# run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
	# 		clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
	# 		 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
	# 		 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov, '-P', 'depth_filter.min_depth=10'])
	# run_mh_trio.wait()
	# outfiles.append(outdir + '/final.passed.tsv')
	# ##lower aaf
	# outdir = pro_name + '_mh_trio_results_lowaaf'
	# run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
	# 		clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
	# 		 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
	# 		 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov, '-P', 'min_minor_allele_number=2', '-P', 'min_minor_allele_percentage=1', 
	# 		 '-P', 'base_number_filter.min_minor_allele_number=2', '-P', 'base_number_filter.min_minor_allele_percentage=1'])
	# run_mh_trio.wait()
	# outfiles.append(outdir + '/final.passed.tsv')
	# ##no filter
	# outdir = pro_name + '_mh_trio_results_nofilter'
	# run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
	# 		clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
	# 		 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
	# 		 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov, '-P', 'min_p_value=0', '-P', ' strand_bias_filter.p_value_cutoff=0', 
	# 		 '-P', 'within_read_position_filter.p_value_cutoff=0'])
	# run_mh_trio.wait()
	# outfiles.append(outdir + '/final.passed.tsv')
	##all 3
	outdir = pro_name + '_mh_trio_results_all3'
	run_mh_trio = subprocess.Popen(['java', '-jar', mosaic_hunter, 'exome', '-P', 'input_file=' + clean_bams[0], '-P', 'mosaic_filter.father_bam_file=' + 
			clean_bams[1], '-P', 'mosaic_filter.mother_bam_file=' + clean_bams[2], '-P', 'reference_file=' + fasta, '-P', 'mosaic_filter.mode=trio',
			 '-P', 'mosaic_filter.dbsnp_file=' + mh_dbsnp, '-P', 'repetitive_region_filter.bed_file=' + mh_rpts, '-P', 'common_site_filter.bed_file=' + mh_common_site,
			 '-P', 'output_dir=' + outdir, '-P', 'syscall_filter.depth=' + pro_cov, '-P', 'depth_filter.min_depth=10', '-P', 'min_minor_allele_number=2', 
			 '-P', 'min_minor_allele_percentage=1', '-P', 'min_p_value=0', '-P', 'base_number_filter.min_minor_allele_number=2', '-P', 
			 'base_number_filter.min_minor_allele_percentage=1', '-P', 'strand_bias_filter.p_value_cutoff=0', '-P', 'within_read_position_filter.p_value_cutoff=0'])
	run_mh_trio.wait()
	outfiles.append(outdir + '/final.passed.tsv')
	##delete clean bams
	for b in clean_bams:
		os.remove(b)
		os.remove(b + '.bai')
		print 'deleting', b, b + '.bai'
	return outfiles

def reformat_annovar(multianno, annotated):
	with open(multianno, "r") as infh, open(annotated, "w") as outfh:
		line_count = 0
		for line in infh:
			line_count += 1
			line = line.strip('\n').split(delim)
			if line_count == 1:
				header = line[:66] + ['ref/alt', 'aaf', 'prior probability', '\n']
				outfh.write(delim.join(header))
			else:
				location = line[5]
				if location == 'exonic' or location == 'splicing':
					max_maf = line[46]
					if max_maf == '.':
						max_maf = 0
					if float(max_maf) <= 0.01:
						alleles = line[66].split(',')
						aaf = float(alleles[1])/ (int(alleles[1]) + int(alleles[0]))
						print location, alleles, aaf
						line_out = line[:67] + [str(aaf), line[67], '\n']
						outfh.write(delim.join(line_out))

def annotate_mh_output(s_name, files_to_annotate):
	for mh_var_file in files_to_annotate:
		print mh_var_file
		analysis = mh_var_file.split('/')[0].split('_')[-1]
		print analysis
		avinput = s_name + '.avinput'
		##convert mh output to an avinput file
		with open(mh_var_file, "r") as infh, open(avinput, "w") as outfh:
			for line in infh:
				line = line.strip('\n').split(delim)
				chrom = line[0]
				pos = line[1]
				ref = line[2]
				if ref == line[6]:
					alt = line[8]
					alleles = line[7] + ',' + line[9]
				else:
					alt = line[6]
					alleles = line[9] + ',' + line[7]
					# print line
					# print chrom, pos, ref, alt, alleles
				post_prob = line[29]
				# print line
				print chrom, pos, ref, alt, alleles, post_prob
				outfh.write(delim.join([chrom, pos, pos, ref, alt, alleles, post_prob, '\n']))
		##run annovar
		command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', s_name]
		annovar = subprocess.Popen(command)
		annovar.wait()
		#filter annovar output
		multianno = s_name + '.hg19_multianno.txt'
		annotated = s_name + '.mosaichunter_trio.' + analysis +  '.xls'
		reformat_annovar(multianno, annotated)


def reformat_annovar_mutect(multianno, annotated):
	with open(multianno, "r") as infh, open(annotated, "w") as outfh:
		line_count = 0
		for line in infh:
			line_count += 1
			line = line.strip('\n').split(delim)
			if line_count == 1:
				header = line[:66] + ['filter', 'affected_cov', 'affected_aaf', 'unaffected_cov', 'unaffected_aaf', '\n']
				outfh.write(delim.join(header))
			else:
				location = line[5]
				if location == 'exonic' or location == 'splicing':
					max_maf = line[46]
					if max_maf == '.':
						max_maf = 0
					if float(max_maf) <= 0.01:
						a_alleles = line[78].split(':')[1].split(',')
						print a_alleles
						a_cov = int(a_alleles[1]) + int(a_alleles[0])
						if a_cov == 0:
							a_aaf = 0
						else:
							a_aaf = float(a_alleles[1])/ (int(a_alleles[1]) + int(a_alleles[0]))
						# print 'a', location, a_alleles, a_cov, a_aaf
						b_alleles = line[79].split(':')[1].split(',')
						# print b_alleles
						b_cov = int(b_alleles[1]) + int(b_alleles[0])
						if b_cov == 0:
							b_aaf = 'no coverage'
						else:
							b_aaf = float(b_alleles[1])/ (int(b_alleles[1]) + int(b_alleles[0]))
						# print 'b', location, b_alleles, b_cov, b_aaf
						if line[75] == 'PASS':
							line_out = line[:66] + [line[75]] +  [str(a_cov), str(a_aaf),str(b_cov), str(b_aaf), '\n']
							outfh.write(delim.join(line_out))


def reformat_annovar_pisces(multianno, annotated):
	with open(multianno, "r") as infh, open(annotated, "w") as outfh:
		line_count = 0
		for line in infh:
			line_count += 1
			line = line.strip('\n').split(delim)
			if line_count == 1:
				header = line[:66] + ['q', 'coverage', 'filter', 'format', 'info', '\n']
				outfh.write(delim.join(header))
			else:
				location = line[5]
				if location == 'exonic' or location == 'splicing':
					max_maf = line[46]
					if max_maf == '.':
						max_maf = 0
					##if <1% maf passed by pisces and not in parents, print line
					if float(max_maf) <= 0.01 and line[76] == 'PASS' and line[66] == '.':
						if line[76] == 'PASS':
							line_out = line[:66] + line[68:70] + [line[76]] + line[78:80] + ['\n']
							outfh.write(delim.join(line_out))

def annotate_mutect_output(s_name):
	vcf = s_name + '.mutect.vcf.gz'
	##run annovar
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options_vcf + ['-out', s_name]
	annovar = subprocess.Popen(command)
	annovar.wait()
	#filter annovar output
	multianno = s_name + '.hg19_multianno.txt'
	annotated = s_name + '.mutect.xls'
	reformat_annovar_mutect(multianno, annotated)


def annotate_pisces_output(proband_vcfs, parents_vcf):
	##add parent vcf to annovar parameters
	av_protocol_pisces_plus = av_protocol_pisces + [parents_vcf]
	# av_protocol_pisces.append(parents_vcf)
	# print av_protocol_pisces
	for proband_vcf in proband_vcfs:
		proband = proband_vcf.split('.')[0]
		print proband_vcf, proband
		##annotate
		command = [table_annovar] + av_buildver + [proband_vcf] + av_ref_dir + av_protocol_pisces_plus + av_operation_pisces + av_options_vcf_pisces + ['-out', proband]
		annovar = subprocess.Popen(command)
		annovar.wait()
		multianno = proband + '.hg19_multianno.txt'
		annotated = proband.split('/')[1] + '.pisces.xls'
		print multianno, annotated
		reformat_annovar_pisces(multianno, annotated)



def combine_vcf_files_gatk(vcfs, out_vcf):
	vcfs_with_v = []
	for vcf in vcfs:
		vcf_with_v = ['-V', vcf]
		vcfs_with_v.extend(vcf_with_v)
	print vcfs_with_v
	# combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fa_file, '-nt', '15', '--variant', vcfs[0],'--variant', vcfs[1],'-o', out_vcf, '-genotypeMergeOptions', 'UNIQUIFY'])
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '15'] + vcfs_with_v + ['-o', out_vcf, '-genotypeMergeOptions', 'REQUIRE_UNIQUE'])
	combine_var.wait()


def variant_calling_pisces(ped_name, bam_files, ped_type):
	out_dir = ped_name + '_pisces'
	vcfs = []
	combined_vcf = ped_name + '.combined_pisces.vcf'
	parents_vcf = ped_name + '.parents_pisces.vcf'
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
	##loaded pisces module
	# run_pisces = subprocess.Popen(['pisces', '-G', pisces_fastq_dir, '-BamPaths', ','.join(bam_files), '-MinVF', '0.01', '-i', exome_capture_bed, '-t', '16', '-OutFolder', out_dir])
	run_pisces.wait()
	##combine vcf files
	combine_vcf_files_gatk(vcfs, combined_vcf)
	combine_vcf_files_gatk(vcfs[-2:], parents_vcf)
	##copy parents vcf to annovar ref dir
	# shutil.copy(parents_vcf, str(av_ref_dir[0]))
	cp_vcf = subprocess.Popen(['cp', parents_vcf, str(av_ref_dir[0])])
	cp_vcf.wait()
	##annoatate proband vcf giving parents info
	annotate_pisces_output(vcfs[:-2], parents_vcf)


def _safe_to_float(x):
	if x is None:
		return None
	else:
		try:
			return float(x)
		except ValueError:
			return None

def filter_vardict_vcf(in_vcf, out_vcf):
	with open(in_vcf, "r") as in_fh, open(out_vcf, "w") as out_fh:
		line_count, lad_count, lfq_count, passed_count = 0,0,0,0
		for line in in_fh:
			if line[0] == '#':
				out_fh.write(line)
			else:
				passed = 'pass'
				line_count += 1
				parts = line.rstrip().split(delim)
				# print parts, parts[8], parts[9]
				sample_ft = {a: v for (a, v) in zip(parts[8].split(":"), parts[9].split(":"))}
				ft2_names = [b.split('=')[0] for b in parts[7].split(';')]
				ft2_values = [b.split('=')[1] for b in parts[7].split(';')]
				sample_ft2 = {a: v for (a, v) in zip(ft2_names, ft2_values)}
				qual = _safe_to_float(parts[5])
				dp = _safe_to_float(sample_ft.get("DP"))
				af = _safe_to_float(sample_ft.get("AF"))
				nm = _safe_to_float(sample_ft2.get("NM"))
				mq = _safe_to_float(sample_ft2.get("MQ"))
				# print dp, af, nm, mq, qual
				if dp is not None and af is not None:
					if dp * af < 6:
						if nm is not None and mq is not None:
							if (mq < 55.0 and nm > 1.0) or (mq < 60.0 and nm > 2.0):
								lad_count += 1
								passed = 'LowAlleleDepth'
					if dp < 10:
						lad_count += 1
						passed = 'LowAlleleDepth'
					if qual is not None and qual < 45:
						lad_count += 1
						passed = 'LowAlleleDepth'
				if af is not None and qual is not None:
					if af < 0.2 and qual < 55 :
						lfq_count += 1
						passed = 'LowFreqQuality'
				if passed == 'pass':
					passed_count += 1
					out_fh.write(line)
	print line_count, lad_count, lfq_count, passed_count

def annotate_vardict_output(proband_vcfs, parents_vcf):
	##add parent vcf to annovar parameters
	av_protocol_pisces_plus = av_protocol_pisces + [parents_vcf]
	# av_protocol_pisces.append(parents_vcf)
	# print av_protocol_pisces
	for proband_vcf in proband_vcfs:
		proband = proband_vcf.split('.')[0]
		print proband_vcf, proband
		##annotate
		command = [table_annovar] + av_buildver + [proband_vcf] + av_ref_dir + av_protocol_pisces_plus + av_operation_pisces + av_options_vcf_pisces + ['-out', proband]
		annovar = subprocess.Popen(command)
		annovar.wait()
		multianno = proband + '.hg19_multianno.txt'
		annotated = proband + '.vardict.xls'
		print multianno, annotated
		reformat_annovar_pisces(multianno, annotated)

def variant_calling_vardict_trio(ped_name, bam_files, ped_type):
	vcfs = []
	combined_vcf = ped_name + '.combined_vardict.vcf'
	parents_vcf = ped_name + '.parents_vardict.vcf'
	for bam in bam_files:
		##check bai file exists and make if it isn't there
		if os.path.isfile(bam + '.bai'):
			print 'bam %s alreaded indexed'%bam
		else:
			print 'indexing bam file:', bam
			st_index = subprocess.Popen(['samtools', 'index', bam])
			st_index.wait()
		##make list of vcfs
		temp_vcf = bam.rsplit('.', 2)[0] + '.vd_temp.vcf'
		vcf = bam.rsplit('.', 2)[0] + '.vd.vcf'
		vcfs.append(vcf)
		out_name = bam.rsplit('.',2)[0] + '.vd'
		with open(temp_vcf, 'w') as out_fh:
			run_vardict = subprocess.Popen([vardictjava, '-G', fasta, '-f', vardict_maf_req, '-N', out_name, '-b', bam, '-c', '1', '-S', '2', '-E', '3', exome_capture_bed], stdout=subprocess.PIPE)
			run_teststrandbias = subprocess.Popen(['teststrandbias.R'], stdin=run_vardict.stdout, stdout=subprocess.PIPE)
			run_var2vcf = subprocess.Popen(['var2vcf_valid.pl', '-N', out_name, '-E', '-f', vardict_maf_req], stdin=run_teststrandbias.stdout, stdout=out_fh)
			run_var2vcf.wait()
		filter_vardict_vcf(temp_vcf, vcf)
	print vcfs
	##combine vcf files
	combine_vcf_files_gatk(vcfs, combined_vcf)
	combine_vcf_files_gatk(vcfs[-2:], parents_vcf)
	##copy parents vcf to annovar ref dir
	cp_vcf = subprocess.Popen(['cp', parents_vcf, str(av_ref_dir[0])])
	cp_vcf.wait()
	##annoate results
	annotate_vardict_output(vcfs[:-2], parents_vcf)

def run_mosaic_variant_calling(work_dir, ped_name, bam_files, proband_sex, proband_coverage, analysis_type):
	os.chdir(work_dir)
	if analysis_type == 'paired':
		variant_calling_paired_mutect2(bam_files, ped_name)
		annotate_mutect_output(ped_name)
	elif analysis_type == 'trio':
		sample_name = bam_files[0].split('.')[0]
		files_to_annotate = variant_calling_trio_mosaichunter(sample_name, bam_files, proband_sex, proband_coverage)
		annotate_mh_output(sample_name, files_to_annotate)
	elif analysis_type == 'pisces_trio':
		variant_calling_pisces(ped_name, bam_files, 'trio')
	elif analysis_type == 'vardict_trio':
		variant_calling_vardict_trio(ped_name, bam_files, 'trio')


##run methods


##paired analysis with mutect
##what the fields are
##run_mosaic_variant_calling(working_dir, 'name', ['affected_bam', 'unaffected_bam'], 'sex', 'pro_cov', 'analysis')
##empty one
##run_mosaic_variant_calling(working_dir, '', ['', ''], '', '', 'paired')
##ghayda peds 0617
# working_dir = '/data/atimms/novartis_fcd_exomes_0317'
# run_mosaic_variant_calling(working_dir, 'LR11-347_paired', ['LR11-347_1209066.bwa_gatk.bam', 'LR11-347_1112100.bwa_gatk.bam'], 'F', '', 'paired')

##trio analysis with mosaic humnter
##what the fields are
##run_mosaic_variant_calling(working_dir, 'name', ['pro_bam', 'dad_bam', 'mom_bam'], 'sex', 'pro_cov', 'analysis')
##empty one
##run_mosaic_variant_calling(working_dir, '', ['', '', ''], '', '', 'trio')

##jimmys files
# working_dir = '/data/atimms/dobyns_mosaic_test_0317'
# run_mosaic_variant_calling(working_dir, 'LR16-173', ['LR16-173_avm_lesion.bwa_gatk.bam', 'LR16-173f.bwa_gatk.bam', 'LR16-173m.bwa_gatk.bam'], 'F', '320', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR13-356', ['LR13-356_avm_lesion.bwa_gatk.bam', 'LR13-356f.bwa_gatk.bam', 'LR13-356m.bwa_gatk.bam'], 'M', '385', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR16-030', ['LR16-030.bwa_gatk.bam', 'LR16-030f.bwa_gatk.bam', 'LR16-030m.bwa_gatk.bam'], 'M', '360', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR16-031', ['LR16-031.bwa_gatk.bam', 'LR16-031f.bwa_gatk.bam', 'LR16-031m.bwa_gatk.bam'], 'M', '980', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR16-172', ['LR16-172.bwa_gatk.bam', 'LR16-172f.bwa_gatk.bam', 'LR16-172m.bwa_gatk.bam'], 'M', '355', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR14-252', ['LR14-252_frozen_scalp.bwa_gatk.bam', 'LR14-252f.bwa_gatk.bam', 'LR14-252m.bwa_gatk.bam'], 'F', '475', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR14-252', ['LR14-252_skin.bwa_gatk.bam', 'LR14-252f.bwa_gatk.bam', 'LR14-252m.bwa_gatk.bam'], 'F', '190', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR16-065', ['LR16-065_eye.bwa_gatk.bam', 'LR16-065f.bwa_gatk.bam', 'LR16-065m.bwa_gatk.bam'], 'M', '490', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR16-065', ['LR16-065_saliva.bwa_gatk.bam', 'LR16-065f.bwa_gatk.bam', 'LR16-065m.bwa_gatk.bam'], 'M', '160', 'trio')
##kims peds
# working_dir = '/data/atimms/dobyns_mosaic_test_0317'
# run_mosaic_variant_calling(working_dir, 'LR04-350', ['LR04-350.bwa_gatk.bam', 'LR04-350f.bwa_gatk.bam', 'LR04-350m.bwa_gatk.bam'], 'M', '72', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR08-108', ['LR08-108.bwa_gatk.bam', 'LR08-108f.bwa_gatk.bam', 'LR08-108m.bwa_gatk.bam'], 'M', '109', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR08-337', ['LR08-337.bwa_gatk.bam', 'LR08-337f.bwa_gatk.bam', 'LR08-337m.bwa_gatk.bam'], 'M', '81', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR09-203', ['LR09-203.bwa_gatk.bam', 'LR09-203f.bwa_gatk.bam', 'LR09-203m.bwa_gatk.bam'], 'M', '69', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR09-387', ['LR09-387.bwa_gatk.bam', 'LR09-387f.bwa_gatk.bam', 'LR09-387m.bwa_gatk.bam'], 'M', '191', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR14-048', ['LR14-048.bwa_gatk.bam', 'LR14-048f.bwa_gatk.bam', 'LR14-048m.bwa_gatk.bam'], 'M', '132', 'trio')
##ghayda's data --- need average cov
# working_dir = '/data/atimms/novartis_fcd_exomes_0317'
# run_mosaic_variant_calling(working_dir, 'LR12-243', ['LR12-243_brain.bwa_gatk.bam', 'LR12-243f.bwa_gatk.bam', 'LR12-243m.bwa_gatk.bam'], 'M', '297', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-245', ['LR12-245_brain.bwa_gatk.bam', 'LR12-245f.bwa_gatk.bam', 'LR12-245m.bwa_gatk.bam'], 'F', '301', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-250', ['LR12-250_brain.bwa_gatk.bam', 'LR12-250f.bwa_gatk.bam', 'LR12-250m.bwa_gatk.bam'], 'F', '291', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-255', ['LR12-255_brain.bwa_gatk.bam', 'LR12-255f.bwa_gatk.bam', 'LR12-255m.bwa_gatk.bam'], 'M', '233', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-259', ['LR12-259_brain.bwa_gatk.bam', 'LR12-259f.bwa_gatk.bam', 'LR12-259m.bwa_gatk.bam'], 'F', '272', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-260', ['LR12-260_brain.bwa_gatk.bam', 'LR12-260f.bwa_gatk.bam', 'LR12-260m.bwa_gatk.bam'], 'F', '344', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR12-269', ['LR12-269_brain.bwa_gatk.bam', 'LR12-269f.bwa_gatk.bam', 'LR12-269m.bwa_gatk.bam'], 'M', '203', 'trio')
##mom not realated so don't analyze this way.
# run_mosaic_variant_calling(working_dir, 'LR13-354', ['LR13-354_brain.bwa_gatk.bam', 'LR13-354f.bwa_gatk.bam', 'LR13-354m.bwa_gatk.bam'], 'F', '', 'trio')
##kims new peds
# working_dir = '/data/atimms/dobyns_mosaic_test_0317'
# run_mosaic_variant_calling(working_dir, 'LR01-079', ['LR01-079.bwa_gatk.bam', 'LR01-079f.bwa_gatk.bam', 'LR01-079m.bwa_gatk.bam'], 'M', '49', 'trio')
# run_mosaic_variant_calling(working_dir, 'LR04-376', ['LR04-376.bwa_gatk.bam', 'LR04-376f.bwa_gatk.bam', 'LR04-376m.bwa_gatk.bam'], 'M', '54', 'trio')

##analysis with Pisces
# working_dir = '/data/atimms/dobyns_mosaic_test_0317'
# run_mosaic_variant_calling(working_dir, 'LR16-173', ['LR16-173_avm_lesion.bwa_gatk.bam', 'LR16-173_saliva.bwa_gatk.bam', 'LR16-173f.bwa_gatk.bam', 'LR16-173m.bwa_gatk.bam'], 'F', '320', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR13-356', ['LR13-356_avm_lesion.bwa_gatk.bam', 'LR13-356_saliva.bwa_gatk.bam', 'LR13-356f.bwa_gatk.bam', 'LR13-356m.bwa_gatk.bam'], 'M', '385', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR16-030', ['LR16-030.bwa_gatk.bam', 'LR16-030f.bwa_gatk.bam', 'LR16-030m.bwa_gatk.bam'], 'M', '360', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR16-031', ['LR16-031.bwa_gatk.bam', 'LR16-031f.bwa_gatk.bam', 'LR16-031m.bwa_gatk.bam'], 'M', '980', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR16-172', ['LR16-172.bwa_gatk.bam', 'LR16-172f.bwa_gatk.bam', 'LR16-172m.bwa_gatk.bam'], 'M', '355', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR14-252', ['LR14-252_frozen_scalp.bwa_gatk.bam', 'LR14-252_skin.bwa_gatk.bam', 'LR14-252f.bwa_gatk.bam', 'LR14-252m.bwa_gatk.bam'], 'F', '475', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR16-065', ['LR16-065_eye.bwa_gatk.bam', 'LR16-065_saliva.bwa_gatk.bam', 'LR16-065f.bwa_gatk.bam', 'LR16-065m.bwa_gatk.bam'], 'M', '490', 'pisces_trio')
# #kims peds
# working_dir = '/data/atimms/dobyns_mosaic_test_0317'
# run_mosaic_variant_calling(working_dir, 'LR04-350', ['LR04-350.bwa_gatk.bam', 'LR04-350f.bwa_gatk.bam', 'LR04-350m.bwa_gatk.bam'], 'M', '72', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR08-108', ['LR08-108.bwa_gatk.bam', 'LR08-108f.bwa_gatk.bam', 'LR08-108m.bwa_gatk.bam'], 'M', '109', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR08-337', ['LR08-337.bwa_gatk.bam', 'LR08-337f.bwa_gatk.bam', 'LR08-337m.bwa_gatk.bam'], 'M', '81', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR09-203', ['LR09-203.bwa_gatk.bam', 'LR09-203f.bwa_gatk.bam', 'LR09-203m.bwa_gatk.bam'], 'M', '69', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR09-387', ['LR09-387.bwa_gatk.bam', 'LR09-387f.bwa_gatk.bam', 'LR09-387m.bwa_gatk.bam'], 'M', '191', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR14-048', ['LR14-048.bwa_gatk.bam', 'LR14-048f.bwa_gatk.bam', 'LR14-048m.bwa_gatk.bam'], 'M', '132', 'pisces_trio')
# #ghayda's data --- need average cov
# working_dir = '/data/atimms/novartis_fcd_exomes_0317'
# run_mosaic_variant_calling(working_dir, 'LR12-243', ['LR12-243_brain.bwa_gatk.bam', 'LR12-243f.bwa_gatk.bam', 'LR12-243m.bwa_gatk.bam'], 'M', '297', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR12-245', ['LR12-245_brain.bwa_gatk.bam', 'LR12-245f.bwa_gatk.bam', 'LR12-245m.bwa_gatk.bam'], 'F', '301', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR12-250', ['LR12-250_brain.bwa_gatk.bam', 'LR12-250f.bwa_gatk.bam', 'LR12-250m.bwa_gatk.bam'], 'F', '291', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR12-255', ['LR12-255_brain.bwa_gatk.bam', 'LR12-255f.bwa_gatk.bam', 'LR12-255m.bwa_gatk.bam'], 'M', '233', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR12-259', ['LR12-259_brain.bwa_gatk.bam', 'LR12-259f.bwa_gatk.bam', 'LR12-259m.bwa_gatk.bam'], 'F', '272', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR12-260', ['LR12-260_brain.bwa_gatk.bam', 'LR12-260f.bwa_gatk.bam', 'LR12-260m.bwa_gatk.bam'], 'F', '344', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR12-269', ['LR12-269_brain.bwa_gatk.bam', 'LR12-269f.bwa_gatk.bam', 'LR12-269m.bwa_gatk.bam'], 'M', '203', 'pisces_trio')
# ##kims new peds
# working_dir = '/data/atimms/dobyns_mosaic_test_0317'
# run_mosaic_variant_calling(working_dir, 'LR01-079', ['LR01-079.bwa_gatk.bam', 'LR01-079f.bwa_gatk.bam', 'LR01-079m.bwa_gatk.bam'], '', '', 'pisces_trio')
# run_mosaic_variant_calling(working_dir, 'LR04-376', ['LR04-376.bwa_gatk.bam', 'LR04-376f.bwa_gatk.bam', 'LR04-376m.bwa_gatk.bam'], 'M', '54', 'pisces_trio')


