#!/usr/bin/env python
import os
import subprocess
import glob
import shutil


'''
load these:
module load java/1.8.0_121 
module load mono/5.10.1.47
module load biobuilds

'''

##set input variables and parameters
delim = '\t'

##programs
gatk2 = '/home/atimms/programs/GenomeAnalysisTK.jar'
gatk4 = '/home/atimms/programs/gatk-4.1.3.0/gatk'
pisces = '/cm/shared/apps/Pisces/5.1.6.54/Pisces.exe'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
manta_config = '/home/atimms/programs/manta-1.6.0.centos6_x86_64/bin/configManta.py'

##files
fasta = '/home/atimms/ngs_data/references/mm10/mm10_alt.fa'
pisces_fastq_dir = '/home/atimms/ngs_data/references/illumina_mm10_alt/'

##sample dict 
sample_dict = {'378561':'KRT_RBX_88s', '378562':'KRT_RBX_138s', '378563':'KRT_RBX_233s', '378564':'KRT_R7B_255s', 
		'378565':'KRT_RKZ_66s', '378566':'KRT_RKZ_116s', '378567':'KRT_RKZ_156s', '378568':'KRT_R79_99s', 
		'378569':'KRT_RBX_233t', '378570':'KRT_R7B_255t', '378571':'KRT_RKZ_66t', '378572':'KRT_R79_99t'}


##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic', '-genericdbfile', 'mgp.v5.merged.snps_all.dbSNP142.avinput']
av_operation = ['-operation', 'g,r,r,f,f']
av_options_vcf = ['-otherinfo', '-remove', '-vcfinput']
av_options = ['-otherinfo', '-remove']


##methods
def make_analysis_dict(input_file):
	analysis_dict = {}
	with open(input_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				gt_name = line[3]
				sample_name = line[1]
				bam = line[2]
				tissue = line[6]
			##add data to analysis_dict
				if gt_name in analysis_dict:
					if tissue == 'sorted_tumor':
						analysis_dict[gt_name][1].append(sample_name)
						analysis_dict[gt_name][3].append(bam)
					elif tissue == 'tail_gDNA':
						analysis_dict[gt_name][0].append(sample_name)
						analysis_dict[gt_name][2].append(bam)
					else:
						print(tissue + ' not recognized')
				else:
					if tissue == 'sorted_tumor':
						analysis_dict[gt_name] = [[], [sample_name], [], [bam]]
					elif tissue == 'tail_gDNA':
						analysis_dict[gt_name] = [[sample_name], [], [bam], []]
					else:
						print(tissue + ' not recognized')

	return analysis_dict

##make avinput files
def convert_to_annovar_move_to_annovar_folder(samples, vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'mt2'])
	con_ann.wait()
	avfiles = []
	for sample in samples:
		av_file = 'mt2.' + sample + '.avinput'
		avfiles.append(av_file)
		# os.rename(temp_av_file, av_file)	
		shutil.copy(av_file, str(av_ref_dir[0]))
	return(avfiles)
	
def run_table_annovar(avinputs, av_prefix):
	print(avinputs, av_prefix)
	for avinput in avinputs:
		av_prefix_ext = av_prefix + '.' + avinput.split('.')[1]
		command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix_ext]
		runannovar = subprocess.Popen(command)
		runannovar.wait()

def run_table_annovar_vcf(vcf):
	av_prefix_ext = vcf.rsplit('.',3)[0]
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options_vcf + ['-out', av_prefix_ext]
	runannovar = subprocess.Popen(command)
	runannovar.wait()

def multianno_to_annotated(avinputs, av_prefix):
	for avinput in avinputs:
		file_prefix = av_prefix + '.' + avinput.split('.')[1]
		head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142', 'mgp.v5', 'zygosity', 'info', 'format', 'vcf_info']
		multianno = file_prefix + '.mm10_multianno.txt'
		annotated = file_prefix + '.annotated.txt'
		with open(multianno, "r") as multi, open(annotated, "w") as final:
			final.write(delim.join(head) + '\n')
			line_count = 0
			for line in multi:
				line_count += 1
				if line_count > 1:
					line = line.split(delim)
					line_out = line[:15] + line[24:]
					final.write(delim.join(line_out))

def multianno_to_annotated_vcf(av_prefix):
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142', 'mgp.v5', 'zygosity', 'info', 'format', 'vcf_info']
	multianno = av_prefix + '.mm10_multianno.txt'
	annotated = av_prefix + '.annotated.txt'
	with open(multianno, "r") as multi, open(annotated, "w") as final:
		final.write(delim.join(head) + '\n')
		line_count = 0
		for line in multi:
			line_count += 1
			if line_count > 1:
				line = line.split(delim)
				line_out = line[:15] + line[24:]
				final.write(delim.join(line_out))

def variant_calling_mutect2(out_prefix, n_bams, t_bams, working_dir):
	##issues, no PONs, and germline resources, ? fasta 
	all_bams = n_bams + t_bams
	bams, normals, sample_names = [], [], []
	for bamfile in all_bams:
		ibam = ['-I', bamfile]
		bams.extend(ibam)

		if bamfile in n_bams:
			sample = ['-normal', bamfile.split('.')[1]]
			normals.extend(sample)
		##get ids for annotation
		else:
			sample_id = bamfile.split('.')[1]
			sample_names.append(sample_id)

	# print(bams, normals)
	##file names
	raw_vcf = out_prefix + '.mt2.temp_raw.vcf.gz'
	filtered_vcf = out_prefix + '.mt2.temp_filtered.vcf.gz'
	vcf_temp1 = out_prefix + '.mt2.temp1.vcf.gz'
	vcf_temp2 = out_prefix + '.mt2.temp2.vcf.gz'
	final_vcf = out_prefix + '.mt2.final.vcf.gz'
	'''
	##run mt2, use --pcr-indel-model NONE if pcr free library prep and --f1r2-tar-gz for filtering by read bias
	# mt2_cmd = [gatk4, 'Mutect2', '-R', fasta] + bams + normals + ['--germline-resource', mutect_gnomad_vcf, '--tmp-dir', working_dir, '--f1r2-tar-gz', out_prefix + 'f1r2.tar.gz', '-O', raw_vcf]
	mt2_cmd = [gatk4, 'Mutect2', '-R', fasta] + bams + normals + ['--tmp-dir', working_dir, '-O', raw_vcf]
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
	'''
	##annotate mutect results
	print(out_prefix, sample_names, final_vcf)
	avs = convert_to_annovar_move_to_annovar_folder(sample_names, final_vcf)
	run_table_annovar(avs, out_prefix + '.mt2')
	multianno_to_annotated(avs, out_prefix + '.mt2')

def filter_mt2_vars(infile, outfile):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		# out_fh.write(delim.join(head) + '\n')
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count == 1:
				out_fh.write(line)
			else:
				line = line.split(delim)
				func = line[5]
				rmsk = line[10]
				segdup = line[11]
				alleles = [i.split(':')[1].split(',') for i in line[17:19]]
				t_alleles_sum = sum([int(i) for i in alleles[0]])
				n_alleles_sum = sum([int(i) for i in alleles[1]])
				# print(alleles, t_alleles_sum, n_alleles_sum)
				if t_alleles_sum >= 10 and n_alleles_sum >= 10:
					if func == 'intronic' or func == 'intergenic':
						if rmsk == '.' and segdup == '.':
							out_fh.write(delim.join(line))
					else:
						out_fh.write(delim.join(line))




def variant_calling_mutect2_ind(n_bams, t_bams, working_dir):
	##one tumor bam at a time
	all_tumor_bam = []
	for bam in t_bams:
		##get bams to analyze
		bams = []
		bams.extend(['-I', bam])
		bams.extend(['-I', n_bams[0]])
		##info on samples etc
		sample_id = bam.split('.')[1]
		out_prefix = bam.rsplit('.', 1)[0]
		normal = ['-normal', n_bams[0].split('.')[1]]
		print(bams, sample_id, out_prefix, normal)
		##file names
		raw_vcf = out_prefix + '.mt2.temp_raw.vcf.gz'
		filtered_vcf = out_prefix + '.mt2.temp_filtered.vcf.gz'
		vcf_temp1 = out_prefix + '.mt2.temp1.vcf.gz'
		vcf_temp2 = out_prefix + '.mt2.temp2.vcf.gz'
		final_vcf = out_prefix + '.mt2.final.vcf.gz'
		'''
		##run mt2, use --pcr-indel-model NONE if pcr free library prep and --f1r2-tar-gz for filtering by read bias
		# mt2_cmd = [gatk4, 'Mutect2', '-R', fasta] + bams + normals + ['--germline-resource', mutect_gnomad_vcf, '--tmp-dir', working_dir, '--f1r2-tar-gz', out_prefix + 'f1r2.tar.gz', '-O', raw_vcf]
		mt2_cmd = [gatk4, 'Mutect2', '-R', fasta] + bams + normal + ['--tmp-dir', working_dir, '-O', raw_vcf]
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
		##annotate mutect results
		print(out_prefix, sample_id, final_vcf)
		# avs = convert_to_annovar_move_to_annovar_folder([sample_id], final_vcf)
		# run_table_annovar_vcf(final_vcf)
		# multianno_to_annotated_vcf(out_prefix + '.mt2')
		##filter vars
		annotated = out_prefix + '.mt2.annotated.txt'
		filtered = out_prefix + '.mt2.filtered.xls'
		filter_mt2_vars(annotated, filtered)

def merge_filter_pisces_vcfs(in_vcfs, out_vcf, group_id):
	vcf_temp1 = group_id + 'temp_01.vcf.gz'
	vcf_temp2 = group_id + 'temp_02.vcf.gz'
	##get vcfs ready
	in_vcfs_gz = []
	for pisces_vcf in in_vcfs:
		in_vcfs_gz.append(pisces_vcf + '.gz')
		##compress and index vcf
		bgzip_run = subprocess.Popen(['bgzip', pisces_vcf])
		bgzip_run.wait()
		bcf_index = subprocess.Popen(['bcftools', 'index', pisces_vcf + '.gz'])
		bcf_index.wait()		

	print(in_vcfs_gz)
	##combine vcfs
	if len(in_vcfs_gz) == 1:
		bcftools_filter = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-r', '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y', '-o', vcf_temp1, '-O', 'z'] + in_vcfs_gz)
		bcftools_filter.wait()
	else:
		bcf_merge = subprocess.Popen(['bcftools', 'merge', '-f', 'PASS', '-r', '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y', '-O', 'z', '-o', vcf_temp1] + in_vcfs_gz)
		bcf_merge.wait()

	##split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_index2 = subprocess.Popen(['bcftools', 'index', vcf_temp1])
	bcf_index2.wait()
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'v', '-o', out_vcf, vcf_temp2])
	bcf_norm2.wait()

def make_list_of_bams(bam_files, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')

def genotype_vars_with_gatk(bamlist, combined_vcf, gatk_norm_vcf, sample_id):
	vcf_temp1 = sample_id + 'temp_21.vcf'
	vcf_temp2 = sample_id + 'temp_22.vcf.gz'
	# gatk_ug = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk2, '-T', 'UnifiedGenotyper', '-R', fasta, '-nct', '8', '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', combined_vcf, '-dcov', '2500', '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', vcf_temp1])
	gatk_ug = subprocess.Popen(['java', '-Xmx120g', '-jar', gatk2, '-T', 'UnifiedGenotyper', '-R', fasta, '-gt_mode', 'GENOTYPE_GIVEN_ALLELES', '-alleles', combined_vcf, '-dcov', '2500', '-out_mode', 'EMIT_ALL_SITES', '-I', bamlist, '-o', vcf_temp1])
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
	# '''
	command = [table_annovar] + av_buildver + [in_vcf] + av_ref_dir + av_protocol + av_operation + av_options_vcf + ['-out', sample_id]
	annovar = subprocess.Popen(command)
	annovar.wait()
	# '''
	##reformat annotated
	##get sample names
	sample_file = sample_id + 'samples_temp.txt'
	with open(sample_file, "w") as sf_fh:
		bgzip_run = subprocess.Popen(['bcftools', 'query', '-l', in_vcf], stdout=sf_fh)
		bgzip_run.wait()
	sample_names = []
	with open(sample_file, "r") as sf_fh:
		for line in sf_fh:
			sample_name = line.rsplit('.',1)[0]
			sample_names.append(sample_name)
	##change mutlianno to annotate
	multianno = sample_id + '.mm10_multianno.txt'
	annotated =  sample_id + '.pisces.annotated.xls'
	with open(multianno, "r") as multi_fh, open(annotated, "w") as ann_fh:
		line_count = 0
		for line in multi_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				head = line[:13] + ['mgp_v5', 'vcf_var', 'pisces_vcf_format'] + sample_names
				ann_fh.write(delim.join(head) + '\n')
			else:
				var = '_'.join(line[17:22])
				line_out = line[:14] + [var] + line[25:]
				ann_fh.write(delim.join(line_out) + '\n')			
				



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
				v = line[14]
				if v in genotyped_dict:
					out_fh.write(delim.join(line + genotyped_dict[v] + ['\n']))
				else:
					print 'varaint %s not genotyped'%v
					nas = ['na'] * len(samples)
					out_fh.write(delim.join(line + nas) + '\n')
	# '''

def filter_split_pisces_vars(in_file, out_suffix):
	##first time just get samples
	with open(in_file, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				vcf_info = line[16:]
				if len(vcf_info) == 3:
					samples = [line[16]]
				else:
					samples = line[16:19]
	print(samples)
	##second time split and filter
	sc = 0
	for sample in samples:
		sc += 1
		out_file = sample +  out_suffix
		print(in_file, out_file)
		with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
			line_count = 0
			for line in in_fh:
				line_count += 1
				line = line.rstrip().split(delim)
				vcf_info = line[16:]
				pisces_sample_info = line[15 + sc]
				normal_gatk_info = line[-1]
				if len(vcf_info) == 3:
					gatk_sample_info = line[16 + sc]
				else:
					gatk_sample_info = line[18 + sc]


				##get header
				if line_count == 1:
					header = line[:16] + [pisces_sample_info, gatk_sample_info, normal_gatk_info]
					out_fh.write(delim.join(header) + '\n')
				else:
					# print(line)
					##if var not genotyped it mean it is multiallelic, or no info for either sample don't use
					# if gatk_sample_info != 'na' and normal_gatk_info != '.' and gatk_sample_info != '.':
					if len(gatk_sample_info.split(',')) == 2  and len(normal_gatk_info.split(',')) == 2:
						t_alleles_sum = sum([int(i) for i in gatk_sample_info.split(',')])
						n_alleles_sum = sum([int(i) for i in normal_gatk_info.split(',')])
						# print(normal_gatk_info)
						n_alt = int(normal_gatk_info.split(',')[1])
						t_alt = int(gatk_sample_info.split(',')[1])
						func = line[5]
						rmsk = line[10]
						segdup = line[11]
						# print(t_alleles_sum, n_alleles_sum, n_alt, func, rmsk, segdup)
						##filter by coverage (10x in tumor and normal)
						if t_alleles_sum >= 10 and n_alleles_sum >= 10:
							##filter by having 0 or 1 reads in the normal sample and 1 or more alt reads in the tumor
							if n_alt <= 1 and t_alt >=1:
								line_out = line[:16] + [pisces_sample_info, gatk_sample_info, normal_gatk_info]
								if func == 'intronic' or func == 'intergenic':
									if rmsk == '.' and segdup == '.':
										out_fh.write(delim.join(line_out) + '\n')
								else:
									out_fh.write(delim.join(line_out) + '\n')



def variant_calling_pisces(ped_name, bam_files, tumor_bams, normal_bams):
	print(ped_name, bam_files)
	out_dir = ped_name + '_pisces'
	for bam in bam_files:
		##check bai file exists and make if it isn't there
		if os.path.isfile(bam + '.bai'):
			print 'bam %s alreaded indexed'%bam
		else:
			print 'indexing bam file:', bam
			st_index = subprocess.Popen(['samtools', 'index', bam])
			st_index.wait()
	##run picses on all bam
	'''
	run_pisces = subprocess.Popen(['mono', pisces, '-G', pisces_fastq_dir, '-BamPaths', ','.join(bam_files), '-MinVF', '0.01', '-t', '8', '-OutFolder', out_dir])
	run_pisces.wait()
	'''
	##merge tumor vcfs, then genotype on all samples and annotate
	tumor_vcfs = [ped_name + '_pisces/' + t.rsplit('.', 1)[0] + '.vcf' for t in  tumor_bams]
	pisces_filtered_vcf = ped_name + '.pisces.combined_filtered.vcf'
	genotyped_var_vcf = ped_name + '.pisces.gatk_geno.vcf'
	ann_file = ped_name + '.pisces.annotated.xls'
	ann_plus_geno_file = ped_name + '.pisces.annotated_genos.xls'
	filtered_file_suffix = '.pisces.filtered.xls'
	bamlist = 'bams.list'
	##merge vcfs
	# merge_filter_pisces_vcfs(tumor_vcfs, pisces_filtered_vcf, ped_name)
	##genotype all vars on parents after filtering for passed
	# make_list_of_bams(bam_files, bamlist)
	# genotype_vars_with_gatk(bamlist, pisces_filtered_vcf, genotyped_var_vcf, ped_name)
	##annotate
	# annotate_filtered_pisces_vcf(pisces_filtered_vcf, ped_name)
	##add genotype data to annotation
	# add_allele_counts_to_annotation(genotyped_var_vcf, ann_file, ann_plus_geno_file)
	##filter vars
	filter_split_pisces_vars(ann_plus_geno_file, filtered_file_suffix)

def call_manta_somatic(ped_name, tumor_bams, normal_bams):
	print(ped_name, tumor_bams, normal_bams)
	##get all normal bams
	nbams = []
	for normal_bam in normal_bams:
		nbam = ['--normalBam', normal_bam]
		nbams.extend(nbam)
	##run one analysis per tumor bam
	for tumor_bam in tumor_bams:
		tbam = ['--tumorBam', tumor_bam]
		outdir = tumor_bam.split('.')[0] + '_manta'
		##call svs
		make_manta_config = subprocess.Popen([manta_config] + nbams + tbam + ['--referenceFasta', fasta, '--exome', '--generateEvidenceBam', '--runDir', outdir])
		make_manta_config.wait()
		manta_run = subprocess.Popen([outdir + '/runWorkflow.py', '-j', '8'])
		manta_run.wait()

def run_manta_on_all(all_bams):
	##get all normal bams
	bams = []
	for bam_file in all_bams:
		bam = ['--bam', bam_file]
		bams.extend(bam)
	outdir = 'all_samples_manta'
	##call svs
	make_manta_config = subprocess.Popen([manta_config] + bams + ['--referenceFasta', fasta, '--exome', '--generateEvidenceBam', '--runDir', outdir])
	make_manta_config.wait()
	manta_run = subprocess.Popen([outdir + '/runWorkflow.py', '-j', '8'])
	manta_run.wait()

def filter_tail_pisces_vcf(in_vcf, out_vcf, group_id):
	vcf_temp1 = group_id + 'temp_11.vcf.gz'
	vcf_temp2 = group_id + 'temp_12.vcf.gz'
	##compress and index vcf
	bgzip_run = subprocess.Popen(['bgzip', in_vcf])
	bgzip_run.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', in_vcf + '.gz'])
	bcf_index.wait()		
	##filter vcf i.e. passed, het, in main chr
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-g', 'het', '-r', '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y', '-o', vcf_temp1, '-O', 'z', in_vcf + '.gz'])
	bcftools_filter.wait()


	##split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_index2 = subprocess.Popen(['bcftools', 'index', vcf_temp1])
	bcf_index2.wait()
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'v', '-o', out_vcf, vcf_temp2])
	bcf_norm2.wait()


def filter_split_pisces_cnv_vars(in_file, out_suffix):
	##first time just get samples
	with open(in_file, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				vcf_info = line[16:]
				if len(vcf_info) == 3:
					samples = [line[17]]
				else:
					samples = line[17:20]
	print(samples)
	##second time split and filter
	sc = 0
	for sample in samples:
		sc += 1
		out_file = sample_dict[sample] +  out_suffix
		print(in_file, out_file)
		with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
			line_count = 0
			for line in in_fh:
				line_count += 1
				line = line.rstrip().split(delim)
				vcf_info = line[16:]
				pisces_sample_info = line[16]
				normal_gatk_info = line[-1]
				gatk_sample_info = line[16 + sc]

				##get header
				if line_count == 1:
					header = line[:16] + [pisces_sample_info, gatk_sample_info, normal_gatk_info]
					out_fh.write(delim.join(header) + '\n')
				else:
					# print(line)
					##if var not genotyped it mean it is multiallelic, or no info for either sample don't use
					# if gatk_sample_info != 'na' and normal_gatk_info != '.' and gatk_sample_info != '.':
					if len(gatk_sample_info.split(',')) == 2  and len(normal_gatk_info.split(',')) == 2:
						t_alleles_sum = sum([float(i) for i in gatk_sample_info.split(',')])
						n_alleles_sum = sum([float(i) for i in normal_gatk_info.split(',')])
						# print(normal_gatk_info)
						n_alt = int(normal_gatk_info.split(',')[1])
						t_alt = int(gatk_sample_info.split(',')[1])
						if t_alt == 0:
							t_af = 0
						else:
							t_af = t_alt/t_alleles_sum
						if n_alt == 0:
							n_af = 0
						else:
							n_af = n_alt/n_alleles_sum
						func = line[5]
						rmsk = line[10]
						segdup = line[11]
						# print(t_alleles_sum, n_alleles_sum, n_alt, func, rmsk, segdup)
						##filter by coverage (10x in tumor and normal)
						if t_alleles_sum >= 20 and n_alleles_sum >= 20:
							##filter by having hom in tumor and het in tail (normal)
							if t_af <= 0.1 or t_af >=0.9:
								if n_af >= 0.3 and n_af <= 0.7:
									if rmsk == '.' and segdup == '.':
										line_out = line[:16] + [pisces_sample_info, gatk_sample_info, normal_gatk_info, str(t_af), str(n_af)]
										out_fh.write(delim.join(line_out) + '\n')
									# if func == 'intronic' or func == 'intergenic':
									# 	if rmsk == '.' and segdup == '.':
									# 		out_fh.write(delim.join(line_out) + '\n')
									# else:
									# 	out_fh.write(delim.join(line_out) + '\n')


def pisces_cnv_hunt(ped_name, bam_files, tumor_bams, normal_bams):
	print(ped_name, bam_files)
	##merge tumor vcfs, then genotype on all samples and annotate
	tail_vcf = ped_name + '_pisces/' + normal_bams[0].rsplit('.', 1)[0] + '.vcf'
	pisces_filtered_vcf = ped_name + '.pisces_tail.filtered_het.vcf'
	genotyped_var_vcf = ped_name + '.pisces_tail.gatk_geno.vcf'
	ann_file = ped_name + '_tail.pisces.annotated.xls'
	ann_plus_geno_file = ped_name + '.pisces_tail.annotated_genos.xls'
	filtered_file_suffix =  '.pisces.het_tail_not_het_tumor.xls'
	bamlist = 'bams.list'
	'''
	##filter tail vcf, get passed het vars
	filter_tail_pisces_vcf(tail_vcf, pisces_filtered_vcf, ped_name)
	##genotype all vars on parents after filtering for passed
	make_list_of_bams(bam_files, bamlist)
	genotype_vars_with_gatk(bamlist, pisces_filtered_vcf, genotyped_var_vcf, ped_name)
	##annotate
	annotate_filtered_pisces_vcf(pisces_filtered_vcf, ped_name + '_tail')
	##add genotype data to annotation
	add_allele_counts_to_annotation(genotyped_var_vcf, ann_file, ann_plus_geno_file)
	'''
	##filter and split i.e. get all vars that are hom in the tumor
	filter_split_pisces_cnv_vars(ann_plus_geno_file, filtered_file_suffix)


##master method
def run_analysis(working_dir, analysis_file):
	os.chdir(working_dir)
	project = analysis_file.split('.')[0]
	##read in analysis file and make 2 dicts
	analysis_dict = make_analysis_dict(analysis_file)
	##check...
	every_bam = []
	for group in analysis_dict:
		# print(group, analysis_dict[group])
		normal_bams = analysis_dict[group][2]
		tumor_bams = analysis_dict[group][3]
		all_bams = normal_bams + tumor_bams
		every_bam.extend(all_bams)
		all_samples = analysis_dict[group][0] + analysis_dict[group][1]
		##run mutect
		##combined -gave weird results so use individual
		# variant_calling_mutect2(group, normal_bams, tumor_bams, working_dir)
		##ran individually
		# variant_calling_mutect2_ind(normal_bams, tumor_bams, working_dir)
		##call pisces
		# variant_calling_pisces(group, all_bams, tumor_bams, normal_bams)
		##use pisces results to look for CNV
		pisces_cnv_hunt(group, all_bams, tumor_bams, normal_bams)
		##call manta on tumor normal pairs
		# call_manta_somatic(group, tumor_bams, normal_bams)
	##run manta on all the bams in gl
	# run_manta_on_all(every_bam)





##run methods
work_dir = '/home/atimms/ngs_data/exomes/working/eric_mouse_exomes_0320'
sample_file = 'exomes_0320_all.txt'
# sample_file = 'exomes_0320_1.txt'
run_analysis(work_dir, sample_file)



