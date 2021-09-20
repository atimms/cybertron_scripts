#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'

'''
requires:

##for pisces
qsub -Iq cdbrmq -l mem=200gb,ncpus=20 -P 19833a08-f6fb-4bea-8526-8a79069da878
##for other
qsub -Iq cdbrmq -l mem=30gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
#java 1.8 for gatk
module load java/1.8.0_202 


'''


##file names etc
fasta = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens_assembly38.fasta'
exome_capture_bed = '/home/atimms/ngs_data/references/exome_beds_hg38/hg38_targets_combined_padded_0721_noM.bed'
bamlist = 'bams.list'
pisces_fastq_dir = '/home/atimms/ngs_data/references/illumina_hg38/'

##programs
pisces = '/home/atimms/programs/pisces_all/Pisces'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
bgzip = '/home/atimms/programs/samtools-1.11/bin/bgzip'
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
samtools = '/home/atimms/programs/samtools-1.11/bin/samtools'

##annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol_pisces_no_parents = ['-protocol', 'refGene,knownGene,ensGene,rmsk,genomicSuperDups,dbnsfp41a,clinvar_20210123,gnomad30_genome']
av_operation_pisces_no_parents = ['-operation', 'g,g,g,r,r,f,f,f']
av_options_vcf_pisces_no_parents = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput', '-arg', '-splicing 10 ,,,,,,,' ]
av_protocol_pisces = ['-protocol', 'refGene,knownGene,ensGene,rmsk,genomicSuperDups,dbnsfp41a,clinvar_20210123,gnomad30_genome,vcf', '-vcfdbfile']
av_operation_pisces = ['-operation', 'g,g,g,r,r,f,f,f,f']
av_options_vcf_pisces = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput', '-arg', '-splicing 10 ,,,,,,,,' ]



##methods
def genotype_vars_with_gatk(bamlist, combined_vcf, gatk_norm_vcf):
	vcf_temp1 = combined_vcf.split('.')[0] + 'temp_21.vcf'
	vcf_temp2 = combined_vcf.split('.')[0] + 'temp_22.vcf.gz'
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

def make_list_of_bams(bam_files, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')

def reformat_annovar_pisces_with_parents(multianno, annotated, parents_bamfiles, multianno_vcf, var_caller):
	ann_temp = annotated + '.temp.xls'
	ann_temp2 = annotated + '.temp2.xls'
	bamlist = 'bams.list'
	multianno_vcf_temp = multianno_vcf.rsplit('.',1)[0] + 'temp.vcf'
	genotyped_var_vcf = multianno.split('.')[0] + '.gatk_genotyped.temp.vcf'
	kept_vars = []
	##filter varaints for locations and more
	with open(multianno, "r") as infh, open(ann_temp, "w") as outfh:
		line_count = 0
		for line in infh:
			line_count += 1
			line = line.strip('\n').split(delim)
			if line_count == 1:
				header = line[:91] + ['q', 'coverage', 'filter', 'format', 'info', 'var_combined', 'sample_alt_reads', '\n']
				outfh.write(delim.join(header))
			else:
				location = line[5]
				exonic_func = line[8]
				# rmsk = line[10]
				# supdup = line[11]
				##if exonic (not synon)
				if (location == 'exonic' or location == 'splicing') and exonic_func != 'synonymous SNV':
					# if rmsk == '.' and supdup == '.':
					max_maf = line[78]
					if max_maf == '.':
						max_maf = 0
					##if <1% maf,  passed by pisces and not in parents, print line
					if float(max_maf) <= 0.01 and line[101] == 'PASS' and line[91] == '.':
						variant = '_'.join(line[95:100])
						if var_caller == 'vardict':
							pro_alt_reads = line[104].split(':')[2]
						if var_caller == 'pisces':
							# if len(line[79].split(':')[2].split(',')[1]) != 2:
							# 	print line[79].split(':')[2]
							pro_alt_reads = line[104].split(':')[2].split(',')[1]
							# pro_alt_reads = line[79].split(':')[2]
						# line_out = line[:66] + line[68:70] + [line[76]] + line[78:80] + [variant, pro_alt_reads] + ['\n']
						line_out = line[:91] + line[93:95] + [line[101]] + line[103:] + [variant, pro_alt_reads] + ['\n']
						kept_vars.append(variant)
						# print variant
						outfh.write(delim.join(line_out))
	##parse multianno_vcf for those varaints we kept (for use in genotyping)
	with open(multianno_vcf, "r") as ma_fh, open(multianno_vcf_temp, "w") as out_fh:
		for line in ma_fh:
			if line[0] == '#':
				out_fh.write(line)
			else:
				line = line.split(delim)
				var = '_'.join(line[:5])
				if var in kept_vars:
					out_fh.write(delim.join(line))
	##genotype all vars on parents
	make_list_of_bams(parents_bamfiles, bamlist)
	# '''
	genotype_vars_with_gatk(bamlist, multianno_vcf_temp, genotyped_var_vcf)
	# '''
	##make dict with new genotype info
	genotyped_dict = {}
	with open(genotyped_var_vcf, "r") as gvv_fh:
		for line in gvv_fh:
			if line[0] != '#':
				line = line.rstrip().split(delim)
				var = '_'.join(line[:5])
				genotypes = line[9:]
				print(genotypes)
				all_alt_reads = []
				for genotype in genotypes:
					if genotype == './.':
						alt_reads = 'not_covered'
					elif len(genotype.split(':')[1].split(',')) == 1:
						alt_reads = 0
					else:
						alt_reads = int(genotype.split(':')[1].split(',')[1])
					all_alt_reads.append(alt_reads)
				genotyped_dict[var] = genotypes + [str(max(all_alt_reads))]
	for i in genotyped_dict:
		print(i, genotyped_dict[i])
	##add genotyped info to ann.txt and parse
	with open(ann_temp, "r") as in_fh, open(ann_temp2, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				if len(parents_bamfiles) == 2:
					out_fh.write(delim.join(line + ['parent1', 'parent2', 'max_parents_alt_reads', '\n']))
				elif len(parents_bamfiles) == 1:
					out_fh.write(delim.join(line + ['parent', 'max_parents_alt_reads', '\n']))
			else:
				v = line[96]
				if v in genotyped_dict:
					out_fh.write(delim.join(line + genotyped_dict[v] + ['\n']))
				else:
					print('varaint %s not genotyped'%v)
					out_fh.write(delim.join(line + ['na', 'na', 'na', '\n']))
	with open(ann_temp2, "r") as in_fh, open(annotated, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				out_fh.write(delim.join(line + ['\n']))
			else:
				sample_alt_reads = int(line[97])
				if len(parents_bamfiles) == 2:
					parent_alt_reads = line[100]
				elif len(parents_bamfiles) == 1:
					parent_alt_reads = line[99]
				##remove if not covered in parents
				if sample_alt_reads >= 2 and parent_alt_reads != 'not_covered' and parent_alt_reads != 'na':
					if int(parent_alt_reads) <= 2:
						out_fh.write(delim.join(line + ['\n']))

def reformat_annovar_pisces_single(multianno, annotated, var_caller):
	##filter varaints for locations and more
	with open(multianno, "r") as infh, open(annotated, "w") as outfh:
		line_count = 0
		for line in infh:
			line_count += 1
			line = line.strip('\n').split(delim)
			if line_count == 1:
				header = line[:91] + ['q', 'coverage', 'filter', 'format', 'info', 'var_combined', 'sample_alt_reads', '\n']
				outfh.write(delim.join(header))
			else:
				location = line[5]
				exonic_func = line[8]
				# rmsk = line[10]
				# supdup = line[11]
				##if exonic (not synon)
				if (location == 'exonic' or location == 'splicing') and exonic_func != 'synonymous SNV':
					# if rmsk == '.' and supdup == '.':
					max_maf = line[78]
					if max_maf == '.':
						max_maf = 0
					##if <1% maf,  passed by pisces and not in parents, print line
					if float(max_maf) <= 0.01 and line[100] == 'PASS':
						variant = '_'.join(line[94:99])
						if var_caller == 'vardict':
							pro_alt_reads = line[103].split(':')[2]
						if var_caller == 'pisces':
							# if len(line[79].split(':')[2].split(',')[1]) != 2:
							# 	print line[79].split(':')[2]
							pro_alt_reads = line[103].split(':')[2].split(',')[1]
							# pro_alt_reads = line[79].split(':')[2]
						# line_out = line[:66] + line[68:70] + [line[76]] + line[78:80] + [variant, pro_alt_reads] + ['\n']
						line_out = line[:91] + line[92:94] + [line[100]] + line[102:] + [variant, pro_alt_reads] + ['\n']
						# kept_vars.append(variant)
						# print variant
						##if >=2 reads in proband keep
						if int(pro_alt_reads)>= 2:
							outfh.write(delim.join(line_out))

def annotate_pisces_output(proband_vcfs, parents_vcf, parents_bams, parents_in_ped):
	##add parent vcf to annovar parameters
	av_protocol_pisces_plus = av_protocol_pisces + [parents_vcf]
	# av_protocol_pisces.append(parents_vcf)
	# print av_protocol_pisces
	for proband_vcf in proband_vcfs:
		proband = proband_vcf.split('.')[0]
		multianno = proband + '.hg38_multianno.txt'
		multianno_vcf = proband + '.hg38_multianno.vcf'
		annotated = proband.split('/')[1] + '.pisces.xls'
		# print multianno, annotated
		print(proband_vcf, proband)
		print(parents_vcf, parents_bams)
		print(parents_in_ped)
		##annotate including parents if we have them
		if parents_in_ped == 'no':
			# '''
			command = [table_annovar] + av_buildver + [proband_vcf] + av_ref_dir + av_protocol_pisces_no_parents + av_operation_pisces_no_parents + av_options_vcf_pisces_no_parents + ['-out', proband]
			annovar = subprocess.Popen(command)
			annovar.wait()
			# '''
			reformat_annovar_pisces_single(multianno, annotated, 'pisces')
		else:
			# '''
			command = [table_annovar] + av_buildver + [proband_vcf] + av_ref_dir + av_protocol_pisces_plus + av_operation_pisces + av_options_vcf_pisces + ['-out', proband]
			annovar = subprocess.Popen(command)
			annovar.wait()
			# '''
			reformat_annovar_pisces_with_parents(multianno, annotated, parents_bams, multianno_vcf, 'pisces')

def combine_vcf_files_gatk(vcfs, out_vcf):
	vcfs_with_v = []
	for vcf in vcfs:
		vcf_with_v = ['-V', vcf]
		vcfs_with_v.extend(vcf_with_v)
	print(vcfs_with_v)
	# combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fa_file, '-nt', '15', '--variant', vcfs[0],'--variant', vcfs[1],'-o', out_vcf, '-genotypeMergeOptions', 'UNIQUIFY'])
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '15'] + vcfs_with_v + ['-o', out_vcf, '-genotypeMergeOptions', 'REQUIRE_UNIQUE'])
	combine_var.wait()

def get_vars_in_all_files(in_files, out_file):
	var_dict = {}
	with open(out_file, "w") as outfh:
		##put all vars in a dict
		file_count = 0
		for infile in in_files:
			file_count += 1
			with open(infile, "r") as infh:
				line_count = 0
				for line in infh:
					line_count += 1
					line = line.strip('\n').split(delim)
					if line_count == 1:
						if file_count == 1:
							outfh.write(delim.join(line[:73] + ['\n']))
					else:
						variant = '_'.join(line[:5])
						if variant in var_dict:
							var_dict[variant][0] += 1
						else:
							var_dict[variant] = [1, line[:73]]
		for v in var_dict:
			##if var in all files write line
			# print var_dict[v]
			if var_dict[v][0] == len(in_files):
				outfh.write(delim.join(var_dict[v][1] + ['\n']))




def variant_calling_pisces(ped_name, pro_bam_files, parent_bam_files):
	out_dir = ped_name + '_pisces'
	all_bams = pro_bam_files + parent_bam_files
	pro_vcfs = [out_dir + '/' + bam.rsplit('/', 1)[1].rsplit('.', 1)[0] + '.vcf' for bam in pro_bam_files]
	parent_vcfs = [out_dir + '/' + bam.rsplit('/', 1)[1].rsplit('.', 1)[0] + '.vcf' for bam in parent_bam_files]
	parent_comined_vcf = ped_name + '.parents_pisces.vcf'
	does_ped_have_parents = 'no'
	##check bai file exists and make if it isn't there
	for bam in all_bams:
		if os.path.isfile(bam + '.bai'):
			print('bam %s alreaded indexed'%bam)
		else:
			print('indexing bam file:', bam)
			st_index = subprocess.Popen([samtools, 'index', bam])
			st_index.wait()
	print(pro_bam_files, parent_bam_files)
	print(pro_vcfs, parent_vcfs)
	##run pisces on all samples
	# '''
	# run_pisces = subprocess.Popen(['mono', pisces, '-G', pisces_fastq_dir, '-BamPaths', ','.join(all_bams), '-MinVF', '0.01', '-i', exome_capture_bed, '-t', '16', '-OutFolder', out_dir])
	run_pisces = subprocess.Popen([pisces, '-g', pisces_fastq_dir, '-BamPaths', ','.join(all_bams), '-MinVF', '0.01', '-i', exome_capture_bed, '-t', '18', '-gVCF', 'false', '-ThreadByChr', 'True', '-FilterDuplicates', 'true', '-OutFolder', out_dir])
	run_pisces.wait()
	# '''
	##if have parents vcf combine
	if len(parent_vcfs) == 1:
		cp_vcf1 = subprocess.Popen(['cp', parent_vcfs[0], parent_comined_vcf])
		cp_vcf1.wait()
		cp_vcf2 = subprocess.Popen(['cp', parent_comined_vcf, str(av_ref_dir[0])])
		cp_vcf2.wait()
		does_ped_have_parents = 'yes'
	elif len(parent_vcfs) >= 2:
		combine_vcf_files_gatk(parent_vcfs, parent_comined_vcf)
		cp_vcf2 = subprocess.Popen(['cp', parent_comined_vcf, str(av_ref_dir[0])])
		cp_vcf2.wait()
		does_ped_have_parents = 'yes'
	##annotate and filter the variants
	annotate_pisces_output(pro_vcfs, parent_comined_vcf, parent_bam_files, does_ped_have_parents)
	##if 2 or more probands / samples get intersecion
	if len(pro_bam_files) > 1:
		pro_files = [bam.rsplit('.', 2)[0] + '.pisces.xls' for bam in pro_bam_files]
		in_all_file = ped_name + '.in_all_affected.pisces.xls'
		get_vars_in_all_files(pro_files, in_all_file)




def run_mosaic_variant_calling(work_dir, ped_dict):
	os.chdir(work_dir)
	for ped in ped_dict:
		##if mosaic analysis wanted i.e. not na for all samples
		cc_count = ped_dict[ped][2].count('case') +  ped_dict[ped][2].count('control')
		# print(ped, cc_count)
		##get case/control bams
		if cc_count > 0:
			case_bams, control_bams = [], []
			mos_index = 0
			for mos_type in ped_dict[ped][2]:
				sample = ped_dict[ped][0][mos_index]
				if mos_type == 'case':
					case_bams.append('Preprocessing/' + sample + '/Recalibrated/' + sample + '.recal.bam')
				elif mos_type == 'control':
					control_bams.append('Preprocessing/' + sample + '/Recalibrated/' + sample + '.recal.bam')
				else:
					print(sample, ' is defined as ', mos_type, " which means it isn't anayzed")
				mos_index += 1
			print(ped, case_bams, control_bams)
			variant_calling_pisces(ped, case_bams, control_bams)

##test
'''
ped_dict = {'NBB4999': [['NBB4999'], 'singleton', ['na']], 'LR10-064': [['LR10-064', 'LR10-064f', 'LR10-064m'], 'trio', ['case', 'control', 'control']], 'LR20-210': [['LR20-210'], 'singleton', ['case']], 'LR21-027': [['LR21-027'], 'singleton', ['case']]}
# ped_dict = {'NBB4999': [['NBB4999'], 'singleton', ['na']], 'LR10-064': [['LR10-064', 'LR10-064f', 'LR10-064m'], 'trio', ['case', 'control', 'control']]}

working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_genedx_0821/'
run_mosaic_variant_calling(working_dir, ped_dict)
'''


