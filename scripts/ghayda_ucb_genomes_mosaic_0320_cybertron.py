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

##sample dict
sample_dict = {'UC1402089': 'LR13-282_UC1402089', 'UC1309030': 'LR13-282_UC1309030', 'UC1403115': 'LR13-282f_UC1403115', 
		'UC1309021': 'LR13-282m_UC1309021', 'UC1405044': 'LR14-155_UC1405044', 'UC1406081': 'LR14-155_UC1406081', 
		'UC1406080': 'LR14-155f_UC1406080', 'UC1406082': 'LR14-155m_UC1406082', 'UC1711041': 'LR17-408_UC1711041', 
		'UC1711046': 'LR17-408_UC1711046', 'UC1711032': 'LR17-408f_UC1711032', 'UC1711033': 'LR17-408m_UC1711033', 
		'UC1210009': 'LR12-123_UC1210009', 'UC1401004': 'LR12-123_UC1401004', 'UC1410166': 'LR12-123f_UC1410166', 
		'UC1205034': 'LR12-123m_UC1205034', 'UC1401017': 'LR12-269_UC1401017', 'UC1109008': 'LR12-269_UC1109008', 
		'UC1305096': 'LR12-269f_UC1305096', 'UC1305097': 'LR12-269m_UC1305097', 'UC1401005': 'LR12-257_UC1401005', 
		'UC1906025': 'LR12-257_UC1906025', 'UC1906053': 'LR12-257f_UC1906053', 'UC1906054': 'LR12-257m_UC1906054', 
		'UC1812014': 'LR16-432a1_UC1812014', 'UC1611017': 'LR16-432a1_UC1611017', 'UC1908020': 'LR16-432f_UC1908020', 
		'UC1908021': 'LR16-432m_UC1908021', 'UC1401007': 'LR12-259_UC1401007', 'UC1305007': 'LR12-259_UC1305007', 'UC1305008': 'LR12-259f_UC1305008', 'UC1305009': 'LR12-259m_UC1305009', 'UC1401013_2': 'LR12-265_UC1401013', 'UC1311084': 'LR12-265_UC1311084', 'UC1302036': 'LR12-265f_UC1302036', 'UC1302037': 'LR12-265m_UC1302037', 'UC1701038-2': 'LR12-101_UC1701038', 'UC1206006': 'LR12-101_UC1206006', 'UC1206140': 'LR12-101f_UC1206140', 'UC1206141': 'LR12-101m_UC1206141', 'UC1401039': 'LR12-249_UC1401039', 'UC1908029': 'LR12-249_UC1908029', 'UC1908030': 'LR12-249f_UC1908030', 'UC1908031': 'LR12-249m_UC1908031', 'UC1708005': 'LR17-337_UC1708005', 'UC1906037': 'LR17-337_UC1906037', 'UC1906038': 'LR17-337f_UC1906038', 'UC1906039': 'LR17-337m_UC1906039', 'UC1308074': 'LR10-246_UC1308074', 'UC1504075': 'LR10-246f_UC1504075', 'UC1308075': 'LR10-246m_UC1308075', 'UC1604151': 'LR16-214_UC1604151', 'UC1605024': 'LR16-214f_UC1605024', 'UC1604152': 'LR16-214m_UC1604152'}


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

def filter_pisces_vars(infile, outfile, sample_name):
	##filter the txt file and add new header
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		# out_fh.write(delim.join(head) + '\n')
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				s_names = line[151:]
				header = line[:12] + [line[15], line[18], line[53], line[67]] + line[83:138] + ['format', 'pisces_vcf_info'] + s_names
				out_fh.write(delim.join(header) + '\n')

			else:
				func = line[5]
				ex_func = line[8]
				rmsk = line[10]
				segdup = line[11]
				gnomad_afs = [line[84], line[101], line[118]]
				gnomad_afs_format = []
				##change . to zero and make into floats
				for af in gnomad_afs:
					if af == '.':
						new_af = 0
					else:
						new_af = float(af)
					gnomad_afs_format.append(new_af)
				# print(gnomad_afs, gnomad_afs_format)
				##get all alleles
				# all_alleles = line[151:]
				all_alleles = [i.split(',') for i in line[151:]]
				allele_counts = len(sum(all_alleles, []))
				# print(all_alleles, allele_counts)
				##only keep if have 2 alleles per sample
				if allele_counts == 6:
					##turn into ints
					alleles_int = [[int(i) for i in j] for j in all_alleles]
					##get sums
					alleles_sum = [sum(i) for i in alleles_int]
					# print(alleles_int, alleles_sum)
					##get pro and parents alleles
					pro_alleles = alleles_int.pop(s_names.index(sample_name))
					pro_alt = pro_alleles[1]
					parent_alleles = alleles_int
					parent_max_alt = max(parent_alleles[0][1], parent_alleles[1][1])
					# print(pro_alleles,parent_alleles,parent_max_alt)
					##coverage >=10 in all samples and <=1% in gnomad (3 versions)
					if min(alleles_sum) >= 10 and max(gnomad_afs_format) <= 0.01:
						if parent_max_alt <= 1 and pro_alt >= 2:
							##get exonic/splicing but not silent
							if func == 'exonic' or func == 'splicing':
								if ex_func != 'synonymous SNV':
									line_out = line[:12] + [line[15], line[18], line[53], line[67]] + line[83:138] + line[149:]
									out_fh.write(delim.join(line_out) + '\n')

def format_pisces_results_per_proband(ped, sample, parents, id_dict):
	sample_parents_bamfiles = [sample + '.bam'] + [p + '.bam' for p in parents]
	bamlist = sample + 'bams.list'
	pisces_vcf = ped + '_pisces/' + sample + '.vcf'
	pisces_filtered_vcf = sample + '.pisces_filtered.vcf'
	genotyped_var_vcf = sample + '.pisces.gatk_genotyes.vcf'
	ann_file = sample + '.pisces.hg38_multianno.txt'
	ann_plus_geno_file = sample + '.pisces.hg38_annotated.xls'
	filtered_file = sample_dict[sample] + '.pisces.hg38.exonic_filtered.xls'
	##genotype all vars on parents after filtering for passed
	# make_list_of_bams(sample_parents_bamfiles, bamlist)
	# filter_norm_pisces_vcf(pisces_vcf,pisces_filtered_vcf, sample)
	# genotype_vars_with_gatk(bamlist, pisces_filtered_vcf, genotyped_var_vcf, sample)
	##annotate
	# annotate_filtered_pisces_vcf(pisces_filtered_vcf, sample + '.pisces')
	##add genotype data to annotation --- missing genotypes? and 
	# add_allele_counts_to_annotation(genotyped_var_vcf, ann_file, ann_plus_geno_file)
	##filter vars
	filter_pisces_vars(ann_plus_geno_file, filtered_file, sample)

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
	
def filter_mt2_vars(infile, in_vcf, outfile):
	##get vcf samples from vcf file
	with open(in_vcf, "r") as vcf_fh:
		for line in vcf_fh:
			if line.startswith('#CHROM'):
				samples = line.rstrip().split(delim)[9:]
	print(infile, in_vcf, outfile, samples)
	##filter the txt file and add new header
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		# out_fh.write(delim.join(head) + '\n')
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				header = line[:12] + [line[15], line[18], line[53], line[67]] + line[83:138] + ['format'] + samples
				out_fh.write(delim.join(header) + '\n')
			else:
				func = line[5]
				ex_func = line[8]
				rmsk = line[10]
				segdup = line[11]
				gnomad_afs = [line[84], line[101], line[118]]
				gnomad_afs_format = []
				##change . to zero and make into floats
				for af in gnomad_afs:
					if af == '.':
						new_af = 0
					else:
						new_af = float(af)
					gnomad_afs_format.append(new_af)
				# print(gnomad_afs, gnomad_afs_format)
				##get all alleles
				alleles = [i.split(':')[1].split(',') for i in line[150:]]
				##turn into ints
				alleles_int = [[int(i) for i in j] for j in alleles]
				##get sums
				alleles_sum = [sum(i) for i in alleles_int]
				# print(alleles, alleles_sum)
				##coverage >=10 in all samples and <=1% in gnomad (3 versions)
				if min(alleles_sum) >= 10 and max(gnomad_afs_format) <= 0.01:
					##get exonic/splicing but not silent
					if func == 'exonic' or func == 'splicing':
						if ex_func != 'synonymous SNV':
							line_out = line[:12] + [line[15], line[18], line[53], line[67]] + line[83:138] + line[149:]
							out_fh.write(delim.join(line_out) + '\n')


def variant_calling_mutect2(ped_dict, bam_suffix, work_dir):
	for ped in ped_dict:
		##call variants individually
		affected_samples = ped_dict[ped][1]
		parent_ids = ped_dict[ped][2]
		new_id_dict = ped_dict[ped][0]
		for affected_sample in affected_samples:
			mt2_vcf = affected_sample + '.mt2.final.vcf.gz'
			multianno = affected_sample +  '.mt2.hg38_multianno.txt'
			multianno_vcf = affected_sample +  '.mt2.hg38_multianno.vcf'
			filtered_file = sample_dict[affected_sample] + '.mt2.hg38.exonic_filtered.xls'
			##get filtered mutect calls
			# call_mutect2_sample_parents(affected_sample, parent_ids, bam_suffix, work_dir)
			# print(affected_sample, parent_ids, bam_suffix, work_dir)
			##annotate mutect results
			# annotate_filtered_pisces_vcf(mt2_vcf, affected_sample + '.mt2')
			##filter mutect results
			filter_mt2_vars(multianno, multianno_vcf, filtered_file)

			
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
# sample_file = 'ghayda_ucb_0320_1.txt'
##all samples
sample_file = 'ghayda_ucb_0320_all.txt'
##run master method
run_analysis(work_dir, sample_file)


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
##run UC1206006 indivdually
# os.chdir(work_dir)
# call_mutect2_sample_parents('UC1206006', ['UC1206140', 'UC1206141'], '.bam', work_dir)


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
##273490
# sample_file = 'LR12-101.txt'
# run_analysis(work_dir, sample_file)
##273509
# sample_file = 'LR12-123.txt'
# run_analysis(work_dir, sample_file)
##
# sample_file = 'LR17-408.txt'
# run_analysis(work_dir, sample_file)








