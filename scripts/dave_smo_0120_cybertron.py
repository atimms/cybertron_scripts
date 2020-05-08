#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import filtering_annotated
import homozygosity_mapping_cybertron

##parameters
delim = '\t'
thread_number = '20'

##modules 
'''
module load biobuilds/2017.11
module load gcc/6.1.0 ##for homozygosity_mapping_cybertron
module load singularity ##for SVE (SV calling)
module load local_python/2.7.14 ##for SVE (SV calling)
'''


##working dir
# working_dir = '/data/atimms/timon_0317'
working_dir = '/home/atimms/ngs_data/enu_mapping/dave_smo_0120'
os.chdir(working_dir)

##programs and files
gatk = '/cm/shared/apps/GATK/GenomeAnalysisTK.jar'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
fasta = '/home/atimms/ngs_data/references/mm10/mm10.fa'
##first 2 columns of fasta.fai file
bedtools_genome_file  = '/home/atimms/ngs_data/references/mm10/mm10.fa.genome'
bwa = 'bwa'
samtools = 'samtools'
bcftools = 'bcftools'
picard = 'picard'
bedtools = 'bedtools'
# bedtools = '/home/atimms/programs/bedtools2-master/bin/bedtools'
bgzip = 'bgzip'
exome_bed_file = '/home/atimms/ngs_data/references/mm10/mm10_refGene_coding_exons.bed'
sve_singularity = '/home/atimms/singularity_containers/sve.simg'
fusorsv = '/home/atimms/programs/SVE/scripts/FusorSV/FusorSV.py'
fusorsv_pickle = '/home/atimms/programs/SVE/scripts/FusorSV/data/models/default.pickle'
manta_config = '/home/atimms/programs/manta-1.6.0.centos6_x86_64/bin/configManta.py'


##files
bamslist_file = 'bams.list'
# st_avinputs = ['K50000022.avinput', 'K50000024.avinput', 'K50000095.avinput', 'combined.avinput']
# 

##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic', 
	'-genericdbfile', 'Smo90.avinput,Smo92.avinput,Smo93.avinput,Smo99.avinput,SmoCre32.avinput,SmoCre33.avinput,SmoCre58.avinput,SmoCre74.avinput,SmoCre79.avinput,SmoCre80.avinput,SmoCr427.avinput,SmoCr470.avinput,SmoCr477.avinput,SmoCr487.avinput,mgp.v5.merged.snps_all.dbSNP142.avinput,129P2_OlaHsd.mgp.v5.snps.dbSNP142.avinput,129S1_SvImJ.mgp.v5.snps.dbSNP142.avinput,129S5SvEvBrd.mgp.v5.snps.dbSNP142.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,,,,,,,,,,,,,']
#av_options = ['-otherinfo', '-remove', '-vcfinput']

##annovar parameters - chr15 region
av_protocol_c15 = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic', 
	'-genericdbfile', 'chr15.Smo90.avinput,chr15.Smo92.avinput,chr15.Smo93.avinput,chr15.Smo99.avinput,chr15.SmoCre32.avinput,chr15.SmoCre33.avinput,chr15.SmoCre58.avinput,chr15.SmoCre74.avinput,chr15.SmoCre79.avinput,chr15.SmoCre80.avinput,chr15.SmoCr427.avinput,chr15.SmoCr470.avinput,chr15.SmoCr477.avinput,chr15.SmoCr487.avinput,mgp.v5.merged.snps_all.dbSNP142.avinput,129P2_OlaHsd.mgp.v5.snps.dbSNP142.avinput,129S1_SvImJ.mgp.v5.snps.dbSNP142.avinput,129S5SvEvBrd.mgp.v5.snps.dbSNP142.avinput']
av_operation_c15 = ['-operation', 'g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f']
av_options_c15 = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,,,,,,,,,,,,,', '-vcfinput']
av_protocol_redo = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic', 
	'-genericdbfile', 'redo.Smo90.avinput,redo.Smo92.avinput,redo.Smo93.avinput,redo.Smo99.avinput,redo.SmoCre32.avinput,redo.SmoCre33.avinput,redo.SmoCre58.avinput,redo.SmoCre74.avinput,redo.SmoCre79.avinput,redo.SmoCre80.avinput,redo.SmoCr427.avinput,redo.SmoCr470.avinput,redo.SmoCr477.avinput,redo.SmoCr487.avinput,mgp.v5.merged.snps_all.dbSNP142.avinput,129P2_OlaHsd.mgp.v5.snps.dbSNP142.avinput,129S1_SvImJ.mgp.v5.snps.dbSNP142.avinput,129S5SvEvBrd.mgp.v5.snps.dbSNP142.avinput']


##align with bwa
def align_with_bwa(sample_dict):
	for sample in sample_dict:
		r1_fq = sample_dict[sample][0]
		r2_fq = sample_dict[sample][1]
		print sample, r1_fq, r2_fq
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
		pe_bam = sample + post_bwa_bam
		sort_bam = sample + sorted_bam
		pic_dup_bam = sample + mkdup_bam
		# '''
		bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-q', '20', '-b', '-@', '5', '-o', pe_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()
		st_sort_pe = subprocess.Popen([samtools, 'sort','-O', 'bam', '-o', sort_bam, '-T', sample,'-@', '10', '-m', '10G', pe_bam])
		st_sort_pe.wait()
		# '''
		##mark duplicates
		picard_md = subprocess.Popen([picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + pic_dup_bam, 'METRICS_FILE=' + sample + '.metrics', 'TMP_DIR=' + working_dir])
		picard_md.wait()

##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

##call samtools on bamfiles
def variant_calling_samtools(samples, bam_suffix, final_vcf_suffix):
	for sample in samples:
		bamfile = sample + bam_suffix
		vcf_temp1 = sample + '.temp_st.vcf.gz'
		vcf_temp2 = sample + '.temp_st2.vcf.gz'
		final_vcf = sample + final_vcf_suffix
		# '''
		stmp = subprocess.Popen([samtools,'mpileup', '-ug','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, bamfile], stdout=subprocess.PIPE)
		bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '19', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
		bcft.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1])
		bcf_index.wait()
		#split multi-allelic variants calls in separate lines
		bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1])
		bcf_norm1.wait()
		# '''
		bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf +'.gz', vcf_temp2])
		bcf_norm2.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
		bcf_index.wait()

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
	output_file = prefix + '.exome.coverage_data.txt'
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


def calculate_genome_coverage(samples, bam_suffix, genome_file, prefix):
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
			##just the file
			# bedtools_cov = subprocess.Popen([bedtools, 'genomecov', '-ibam', bam, '-g', genome_file], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen([bedtools, 'genomecov', '-ibam', bam, '-g', genome_file], stdout=subprocess.PIPE)
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



##make avinput files
def convert_to_annovar(samples, vcf_suffix):
	for sample in samples:
		vcf = sample + vcf_suffix
		con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
		con_ann.wait()
	temp_files = glob.glob('temp*.avinput')
	for temp_file in temp_files:
		real_file = temp_file[5:]
		os.rename(temp_file, real_file)
		shutil.copy(real_file, str(av_ref_dir[0]))

def run_table_annovar(samples):
	for sample in samples:
		avinput = sample + '.avinput'
		av_prefix = avinput.rsplit('.',1)[0]
		command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()

def multianno_to_annotated(samples): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 
	'genomicSuperDups', 'snp142', 'Smo90.avinput','Smo92.avinput','Smo93.avinput','Smo99.avinput','SmoCre32.avinput','SmoCre33.avinput','SmoCre58.avinput',
	'SmoCre74.avinput','SmoCre79.avinput','SmoCre80.avinput', 'SmoCr427.avinput','SmoCr470.avinput','SmoCr477.avinput','SmoCr487.avinput', 'mgp.v5.snps', 
	'129P2_OlaHsd.mgp.v5', '129S1_SvImJ.mgp.v5', '129S5SvEvBrd.mgp.v5' , 'Zygosity', 'Qual', 'Coverage', 'Chr', 'Pos', 'Filter', 'Ref2', 'Alt2','Qual', 
	'Filter', 'GT_info', 'Format', 'Info']
	head_out = delim.join(head + ['\n'])
	for sample in samples:
		avinput = sample + '.avinput'
		av_prefix = avinput.rsplit('.',1)[0]
		multianno = av_prefix + '.mm10_multianno.txt'
		annotated = av_prefix + '.annotated.txt'
		with open(multianno, "r") as multi, open(annotated, "w") as final:
			final.write(head_out)
			line_count = 0
			for line in multi:
				line_count += 1
				if line_count > 1:
					final.write(line)

def make_bed_from_ann_txt(infile, outfile):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				out_fh.write(delim.join(line[:3]) + '\n')

def run_manta_svs(samples):
	bams = []
	for sample in samples:
		in_bam = sample + '.bwa_mkdup.bam'
		bams.append(in_bam)
	##call svs
	manta_config = subprocess.Popen([manta_config, '--bam', bams[0], '--bam', bams[1], '--referenceFasta', fasta, '--runDir', '.'])
	manta_config.wait()
	manta_run = subprocess.Popen(['./runWorkflow.py'])
	manta_run.wait()

def make_table_of_genotypes(samples, file_suffix, all_snp_suffix, out_file):
	snp_dict = {}
	##get all vars and put in a dict to get counts
	for sample in samples:
		in_file = sample + file_suffix
		with open(in_file, "r") as in_fh:
			line_count = 0
			for line in in_fh:
				line = line.split(delim)
				line_count += 1
				if line_count > 1:
					snp = '_'.join(line[:5])
					snp_dict[snp] = []
	print(len(snp_dict))
	##add genotypes, so make dict per sample then add if seen atall
	for sample in samples:
		sample_dict = {}
		in_file = sample + all_snp_suffix
		with open(in_file, "r") as in_fh:
			line_count = 0
			for line in in_fh:
				line = line.split(delim)
				line_count += 1
				if line_count > 1:
					snp2 = '_'.join(line[:5])
					sample_dict[snp2] = line[27]
		for s in snp_dict:
			if s in sample_dict:
				snp_dict[s].append(sample_dict[s])
			else:
				snp_dict[s].append('wt')

	##make final file
	with open(out_file, "w") as outfh:
		header = ['Chr','Start','End','Ref','Alt'] + samples
		outfh.write(delim.join(header) + '\n')
		for snppy in snp_dict:
			line_out = snppy.split('_') + snp_dict[snppy]
			outfh.write(delim.join(line_out) + '\n')


def make_list_of_bams(bam_files, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')


##call samtools on bamfiles
def variant_calling_samtools_region(project, final_vcf_suffix, bamlist_file, region):
	vcf_temp1 = project + '.temp_st.vcf.gz'
	vcf_temp2 = project + '.temp_st2.vcf.gz'
	final_vcf = project + final_vcf_suffix
	# '''
	stmp = subprocess.Popen([samtools,'mpileup', '-ug','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist_file, '-r', region], stdout=subprocess.PIPE)
	bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '19', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
	bcft.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	# '''
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf +'.gz', vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()

##make avinput files
def convert_to_annovar_move_to_annovar_folder(samples, vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'chr15'])
	con_ann.wait()
	for sample in samples:
		av_file = 'chr15.' + sample + '.avinput'
		shutil.copy(av_file, str(av_ref_dir[0]))

def convert_to_annovar_move_to_annovar_folder_redo(samples, vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'redo'])
	con_ann.wait()
	for sample in samples:
		av_file = 'redo.' + sample + '.avinput'
		shutil.copy(av_file, str(av_ref_dir[0]))


def run_table_annovar_chr15(vcf, av_prefix):
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol_c15 + av_operation_c15 + av_options_c15 + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def run_table_annovar_genome(vcf, av_prefix):
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol_redo + av_operation_c15 + av_options_c15 + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()


def multianno_to_annotated_chr15(av_prefix): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 
	'genomicSuperDups', 'snp142', 'Smo90','Smo92','Smo93','Smo99','SmoCre32','SmoCre33','SmoCre58', 'SmoCre74','SmoCre79','SmoCre80', 'SmoCr427','SmoCr470',
	'SmoCr477','SmoCr487', 'mgp.v5.snps', '129P2_OlaHsd.mgp.v5', '129S1_SvImJ.mgp.v5', '129S5SvEvBrd.mgp.v5' , 'Zygosity', 'Qual', 'Coverage', 'Chr', 'Pos', 
	'Filter', 'Ref2', 'Alt2','Qual', 'Filter', 'GT_info', 'Format', 'Smo90', 'Smo93', 'Smo99', 'SmoCre58', 'SmoCre74', 'SmoCre79', 'SmoCre80', 'SmoCr427', 
	'SmoCr470', 'SmoCr477', 'SmoCr487', 'Smo92', 'SmoCre32', 'SmoCre33']
	head_out = delim.join(head + ['\n'])
	multianno = av_prefix + '.mm10_multianno.txt'
	annotated = av_prefix + '.annotated.txt'
	with open(multianno, "r") as multi, open(annotated, "w") as final:
		final.write(head_out)
		line_count = 0
		for line in multi:
			line_count += 1
			if line_count > 1:
				final.write(line)


def filter_rpt(file_prefix):
	##remove if in rmsk, segdup
	filtering_annotated.filter(working_dir, "and", file_prefix + '.annotated.txt', file_prefix + ".c15_vars.xls", [11,12], ['==','=='], ['.','.'])

def filter_for_coverage(file_prefix):
	infile = file_prefix + ".c15_vars.xls"
	outfile = file_prefix + ".rpts_coverage.xls"
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count == 1:
				out_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				genotypes = line[43:]
				# print(genotypes, len(genotypes))
				depths = [int(i.split(":")[2]) for i in genotypes]
				print(depths, min(depths), sum(depths))
				if min(depths) > 5 and sum(depths) < 281:
					out_fh.write(delim.join(line) + '\n')

def genotype_snps_in_chr15_region(samples, prefix):
	bams = [i + '.bwa_mkdup.bam' for i in samples]
	list_bam_file = 'bams.list'
	st_vcf_suffix = '.st.vcf.gz'
	roi = 'chr15:0-50000000'
	# make_list_of_bams(bams, list_bam_file)
	# variant_calling_samtools_region(prefix, st_vcf_suffix, list_bam_file, roi)
	# convert_to_annovar_move_to_annovar_folder(samples, prefix + st_vcf_suffix)
	# run_table_annovar_chr15(prefix + st_vcf_suffix, prefix)
	# multianno_to_annotated_chr15(prefix)
	# filter_rpt(prefix)
	filter_for_coverage(prefix)

def variant_calling_samtools_combined(project, final_vcf_suffix, bamlist_file):
	vcf_temp1 = project + '.temp_st.vcf.gz'
	vcf_temp2 = project + '.temp_st2.vcf.gz'
	final_vcf = project + final_vcf_suffix
	'''
	stmp = subprocess.Popen([samtools,'mpileup', '-ug','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist_file], stdout=subprocess.PIPE)
	bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '19', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
	bcft.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	'''
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()


def regenotype_snps_in_genome(samples, prefix):
	bams = [i + '.bwa_mkdup.bam' for i in samples]
	list_bam_file = 'bams.list'
	st_vcf_suffix = '.st.vcf.gz'
	roi = 'chr15:0-50000000'
	# make_list_of_bams(bams, list_bam_file)
	# variant_calling_samtools_combined(prefix, st_vcf_suffix, list_bam_file)
	# convert_to_annovar_move_to_annovar_folder_redo(samples, prefix + st_vcf_suffix)
	run_table_annovar_genome(prefix + st_vcf_suffix, prefix)
	multianno_to_annotated_chr15(prefix)
	filter_rpt(prefix)
	filter_for_coverage(prefix)


##call methods
##parameters
project_name = 'dave_smo_0120'
# samples = ['K50000022']
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
st_vcf_suffix = '.st.vcf'
fq_dict = {'Smo90': ['Smo90_USD16094073L_HL2NMDSXX_L1_1.fq.gz', 'Smo90_USD16094073L_HL2NMDSXX_L1_2.fq.gz'],
		'Smo92': ['Smo92_USD16094074L_HL2NMDSXX_L1_1.fq.gz', 'Smo92_USD16094074L_HL2NMDSXX_L1_2.fq.gz'],
		'Smo93': ['Smo93_USD16094075L_HL2NMDSXX_L1_1.fq.gz', 'Smo93_USD16094075L_HL2NMDSXX_L1_2.fq.gz'],
		'Smo99': ['Smo99_USD16094076L_HL2NMDSXX_L1_1.fq.gz', 'Smo99_USD16094076L_HL2NMDSXX_L1_2.fq.gz'],
		'SmoCre32': ['SmoCre32_USD16094067L_HL2NMDSXX_L1_1.fq.gz', 'SmoCre32_USD16094067L_HL2NMDSXX_L1_2.fq.gz'],
		'SmoCre33': ['SmoCre33_USD16094068L_HL2NMDSXX_L1_1.fq.gz', 'SmoCre33_USD16094068L_HL2NMDSXX_L1_2.fq.gz'],
		'SmoCre58': ['SmoCre58_USD16094069L_HL2NMDSXX_L1_1.fq.gz', 'SmoCre58_USD16094069L_HL2NMDSXX_L1_2.fq.gz'],
		'SmoCre74': ['SmoCre74_USD16094070L_HL2NMDSXX_L1_1.fq.gz', 'SmoCre74_USD16094070L_HL2NMDSXX_L1_2.fq.gz'],
		'SmoCre79': ['SmoCre79_USD16094071L_HL2NMDSXX_L1_1.fq.gz', 'SmoCre79_USD16094071L_HL2NMDSXX_L1_2.fq.gz'],
		'SmoCre80': ['SmoCre80_USD16094072L_HL2NMDSXX_L1_1.fq.gz', 'SmoCre80_USD16094072L_HL2NMDSXX_L1_2.fq.gz']}

new_fq_dict = {'SmoCr427': ['SmoCr427_CSFP190050239-1a_HTF77DSXX_L1_1.fq.gz', 'SmoCr427_CSFP190050239-1a_HTF77DSXX_L1_2.fq.gz'],
		'SmoCr470': ['SmoCr470_CSFP190050236-1a_HTF77DSXX_L1_1.fq.gz', 'SmoCr470_CSFP190050236-1a_HTF77DSXX_L1_2.fq.gz'],
		'SmoCr477': ['SmoCr477_CSFP190050237-1a_HTF77DSXX_L1_1.fq.gz', 'SmoCr477_CSFP190050237-1a_HTF77DSXX_L1_2.fq.gz'],
		'SmoCr487': ['SmoCr487_CSFP190050238-1a_HTF77DSXX_L1_1.fq.gz', 'SmoCr487_CSFP190050238-1a_HTF77DSXX_L1_2.fq.gz']}

##combine dicts
all_fq_dict = fq_dict.copy()
all_fq_dict.update(new_fq_dict)


affected = ['Smo90', 'Smo93', 'Smo99', 'SmoCre58', 'SmoCre74', 'SmoCre79', 'SmoCre80', 'SmoCr427', 'SmoCr470', 'SmoCr477', 'SmoCr487']
unaffected = ['Smo92', 'SmoCre32', 'SmoCre33']
one_affected = [affected[0]]
other_samples = affected[1:] + unaffected
samples = fq_dict.keys()

##map with bwa and process with samtools etc
# '''
# align_with_bwa(new_fq_dict)
# variant_calling_samtools(new_fq_dict, mkdup_bam, st_vcf_suffix)
# convert_to_annovar(all_fq_dict, st_vcf_suffix + '.gz')
# ##just need one sample to find variants all affected
# run_table_annovar(one_affected)
# multianno_to_annotated(one_affected)
# ##then repeat for completness
# run_table_annovar(other_samples)
# multianno_to_annotated(other_samples)
# '''

##coverage
'''
calculate_genome_coverage(samples, mkdup_bam, bedtools_genome_file, project_name)
calculate_exome_coverage(samples, mkdup_bam, exome_bed_file, project_name)
'''

##align combined fq files for SV analysis
# comb_fq_dict = {'smo_unaffected': ['smo_unaffected.r1.fq.gz', 'smo_unaffected.r2.fq.gz'],
# 		'smo_affected': ['smo_affected.r1.fq.gz', 'smo_affected.r2.fq.gz']}
# align_with_bwa(comb_fq_dict)
# run_sve_from_singularity(comb_fq_dict)
# run_manta_svs(comb_fq_dict)





##filter and graphing data
##parameters
##filtering -- column nubers start at 1
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = 'synonymous SNV'
zygosity_col = 32
cov_col = 34
cov_definition = [8,50]
qual_col = 40
qual_definition = 30
##het snp mapping 
genome_fai = '/home/atimms/ngs_data/references/mm10/mm10.fa.fai'
window_size = [1000000,5000000,2000000]
# window_size = [10000000]
step_size = 1000000
info_col = 40
naf_values = [0.8,0.9,0.95]

##filter vars
'''
for sample in one_affected:
	##remove if not in all 11 affected
	filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "1.temp", [14,16,17,20,21,22,23,24,25,26,27], ['!=','!=','!=','!=','!=','!=','!=','!=','!=','!=','!='], ['','','','','','','','','','',''])
	##filter variants by coverage and quality 
	filtering_annotated.filter(working_dir, "and", sample + "1.temp", sample + "2.temp", [cov_col,cov_col,qual_col], ['>=','<=','>='], [cov_definition[0], cov_definition[1],qual_definition])
	##remove if in  rpt region
	filtering_annotated.filter(working_dir, "and", sample + '2.temp', sample + '.in_all_aff.all_vars.xls', [11,12], ['==','=='], ['',''])
	##filter variants by all affected being hom
	filtering_annotated.filter(working_dir, "and", sample + '.in_all_aff.all_vars.xls', sample + '.in_all_aff.hom_vars.xls', [14,16,17,20,21,22,23,24,25,26,27], ['==','==','==','==','==','==','==','==','==','==','=='], ['hom','hom','hom','hom','hom','hom','hom','hom','hom','hom','hom'])
	##get 129 variants for all and hom vars
	filtering_annotated.filter(working_dir, "or", sample + '.in_all_aff.all_vars.xls', sample +'.in_all_aff.all_vars_129.xls', [29,30,31], ['!=','!=','!='], ['','',''])
	filtering_annotated.filter(working_dir, "or", sample + '.in_all_aff.hom_vars.xls', sample +'.in_all_aff.hom_vars_129.xls', [29,30,31], ['!=','!=','!='], ['','',''])
	##all vars not in any unaffected
	filtering_annotated.filter(working_dir, "and", sample + '.in_all_aff.all_vars.xls', sample +'.in_all_aff_not_unaffected.all_vars.xls', [15,18,19], ['==','==','=='], ['','',''])
	filtering_annotated.filter(working_dir, "and", sample + '.in_all_aff.all_vars_129.xls', sample +'.in_all_aff_not_unaffected.all_vars_129.xls', [15,18,19], ['==','==','=='], ['','',''])
	#hom vars not hom in unaffected
	filtering_annotated.filter(working_dir, "and", sample + '.in_all_aff.hom_vars.xls', sample +'.in_all_aff_not_hom_unaffected.hom_vars.xls', [15,18,19], ['!=','!=','!='], ['hom','hom','hom'])
	filtering_annotated.filter(working_dir, "and", sample + '.in_all_aff.hom_vars_129.xls', sample +'.in_all_aff_not_hom_unaffected.hom_vars_129.xls', [15,18,19], ['!=','!=','!='], ['hom','hom','hom'])
'''

##make files to graph the data
files_to_graph = ['Smo90.in_all_aff.all_vars.xls', 'Smo90.in_all_aff.hom_vars.xls', 'Smo90.in_all_aff.all_vars_129.xls', 'Smo90.in_all_aff.hom_vars_129.xls', 'Smo90.in_all_aff_not_unaffected.all_vars.xls', 
		'Smo90.in_all_aff_not_unaffected.all_vars_129.xls', 'Smo90.in_all_aff_not_hom_unaffected.hom_vars.xls', 'Smo90.in_all_aff_not_hom_unaffected.hom_vars_129.xls']
# files_to_graph = ['Smo90.in_all_aff_not_unaffected.all_vars.xls', 'Smo90.in_all_aff_not_hom_unaffected.hom_vars.xls']

##make bed from enu vars
'''
beds_to_graph = []
for file_to_graph in files_to_graph:
	bed_to_graph = file_to_graph.rsplit('.', 1)[0] + '.bed'
	beds_to_graph.append(bed_to_graph)
	make_bed_from_ann_txt(file_to_graph, bed_to_graph)

for bed_to_graph in beds_to_graph:
	for ws in window_size:
		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping_cybertron.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]
		print genome_and_window
		window_bed = genome_and_window + '.bed'
		##all shared snps
		out_bed = bed_to_graph.rsplit('.', 1)[0] + '.' + genome_and_window + '.bed'
		##bedtools intersect 
		with open('temp.bed', "w") as naf_fh: 
			hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', window_bed, '-b', bed_to_graph, '-c'], stdout=naf_fh)
			hom_bt_intersect.wait()
		##add header
		with open(out_bed, "w") as out_fh, open('temp.bed', "r") as in_fh:
			out_fh.write(delim.join(['chr', 'start', 'end', 'snp_number', '\n']))
			for line in in_fh:
				out_fh.write(line)
'''

c15_project = project_name + '_chr15'
all_samples = affected + unaffected

##look at chr15 snps
# genotype_snps_in_chr15_region(all_samples, c15_project)

##re genotype all snps
redo_project = project_name + '_redo'
regenotype_snps_in_genome(all_samples, redo_project)


