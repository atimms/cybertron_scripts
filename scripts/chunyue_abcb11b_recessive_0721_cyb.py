#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import filtering_annotated
import homozygosity_mapping_cybertron
import math

##parameters
delim = '\t'
thread_number = '20'

##modules 
'''
qsub -Iq cdbrmq -l mem=20gb,ncpus=1 -P 2cf4ba3b-f9ef-4cfc-a2c1-be85f5db6730
module load biobuilds/2017.11
module load gcc/6.1.0 ##for homozygosity_mapping_cybertron

graphing done in enu_hom_mapping_cyb.R
'''


##files
fasta = '/home/atimms/ngs_data/references/danRer11/danRer11.fa'
fasta_fai = '/home/atimms/ngs_data/references/danRer11/danRer11.fa.fai'
##first 2 columns of fasta.fai file
bedtools_genome_file  = '/home/atimms/ngs_data/references/danRer11/danRer11.fa.genome'

##programs
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'

##annovar parameters
av_genome = 'danRer11'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,ensGene,rmsk,generic,generic,generic,generic', '-genericdbfile', 'TP1FRes.avinput,TP1FUnr.avinput,ResP1.avinput,ResP2.avinput']
av_operation = ['-operation', 'g,g,r,f,f,f,f']
av_options_vcf = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,', '-vcfinput']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,', '--nopolish']


##methods

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
	print(output_file)
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


##align with bwa
def align_with_bwa(sample_dict, working_dir):
	for sample in sample_dict:
		r1_fq = sample_dict[sample][0]
		r2_fq = sample_dict[sample][1]
		print(sample, r1_fq, r2_fq)
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
		pe_bam = sample + post_bwa_bam
		sort_bam = sample + sorted_bam
		pic_dup_bam = sample + mkdup_bam
		# '''
		bwa_pe = subprocess.Popen(['bwa', 'mem', '-M', '-t', '5', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen(['samtools', 'view', '-q', '20', '-b', '-@', '5', '-o', pe_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()
		st_sort_pe = subprocess.Popen(['samtools', 'sort','-O', 'bam', '-o', sort_bam, '-T', sample,'-@', '10', '-m', '10G', pe_bam])
		st_sort_pe.wait()
		# '''
		##mark duplicates
		picard_md = subprocess.Popen(['picard', 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + pic_dup_bam, 'METRICS_FILE=' + sample + '.metrics', 'TMP_DIR=' + working_dir])
		picard_md.wait()

##make list of all bam files to be analyzed
def make_list_of_bams(sample_names, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample_name in sample_names:
			outfh.write(sample_name + bam_suffix + '\n')


##call samtools on bamfiles
def variant_calling_samtools(project, final_vcf_suffix, bamlist_file):
	vcf_temp1 = project + '.temp_st.vcf.gz'
	vcf_temp2 = project + '.temp_st2.vcf.gz'
	final_vcf = project + final_vcf_suffix
	# '''
	# stmp = subprocess.Popen([samtools,'mpileup', '-ug','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist_file], stdout=subprocess.PIPE)
	# bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '19', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
	stmp = subprocess.Popen(['samtools','mpileup', '-Ou','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist_file], stdout=subprocess.PIPE)
	bcft = subprocess.Popen(['bcftools','call', '-vmO', 'z', '--threads', '9', '-o', vcf_temp1], stdin=stmp.stdout)
	bcft.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	# '''
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()

##make avinput files
def convert_to_annovar_move_to_annovar_folder(samples, vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
	con_ann.wait()
	for sample in samples:
		temp_av_file = 'temp.' + sample + '.avinput'
		av_file = sample + '.avinput'
		shutil.move(temp_av_file, av_file)
		shutil.copy(av_file, str(av_ref_dir[0]))


def convert_to_annovar_sample_combined(vcf, out_prefix):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withfreq', '-allsample', '-outfile', out_prefix + '.avinput'])
	con_ann.wait()

def run_table_annovar_combined(av_prefix):
	avinput = av_prefix + '.avinput'
	command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()


def run_table_annovar(samples):
	for sample in samples:
		avinput = sample + '.avinput'
		av_prefix = avinput.rsplit('.',1)[0]
		command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
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

def multianno_to_annotated(samples):
	for sample in samples:
		multianno = sample + '.danRer11_multianno.txt'
		outfile = sample + '.annotated.xls'
		with open(multianno, "r") as av, open(outfile, "w") as final:
			lc = 0
			for line in av:
				lc += 1
				line = line.strip('\n').split(delim)
				if lc == 1:
					header = line[:16] + ['TP1FRes','TP1FUnr', 'ResP1', 'ResP2', 'zygosity', 'qual', 'cov', 'info', 'format', 'vcf']
					final.write(delim.join(header) + '\n')
				else:
					line_out = line[:23] + line[30:] 
					final.write(delim.join(line_out) + '\n')


def filter_vcf_non_std_chr(in_vcf, out_vcf):
	vcf_temp1 = 'temp1.vcf.gz'
	vcf_temp2 = 'temp2.vcf.gz'
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-r', 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25', '-o', vcf_temp1, '-O', 'z', in_vcf])
	bcftools_filter.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', out_vcf, vcf_temp2])
	bcf_norm2.wait()

def make_bed_from_ann_txt(infile, outfile):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				out_fh.write(delim.join(line[:3]) + '\n')


##run methods
work_dir = '/home/atimms/ngs_data/enu_mapping/chunyue_abcb11b_recessive_0721'
os.chdir(work_dir)


##params
project_name = 'chunyue_abcb11b_recessive_0721'
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
bamlist = 'bams.list'
st_vcf_suffix = '.st.vcf.gz'
result_vcf = project_name + st_vcf_suffix
filtered_vcf = project_name + '.st_filtered_chr.vcf.gz'
fq_dict = {'TP1FRes': ['TP1FRes_USD16091700L_HL2H3DSXX_L2_1.fq.gz', 'TP1FRes_USD16091700L_HL2H3DSXX_L2_2.fq.gz'],
	'TP1FUnr':['TP1FUnr_USD16091701L_HL2H3DSXX_L2_1.fq.gz', 'TP1FUnr_USD16091701L_HL2H3DSXX_L2_2.fq.gz'],
	'ResP1':['ResP1.r1.fastq.gz', 'ResP1.r2.fastq.gz'],'ResP2':['ResP2.r1.fastq.gz', 'ResP2.r2.fastq.gz']}
parent_fq_dict = {'ResP1':['ResP1.r1.fastq.gz', 'ResP1.r2.fastq.gz'],'ResP2':['ResP2.r1.fastq.gz', 'ResP2.r2.fastq.gz']}


##map with bwa and process with samtools etc
# align_with_bwa(parent_fq_dict, work_dir)
# make_list_of_bams(fq_dict, mkdup_bam, bamlist)
# variant_calling_samtools(project_name, st_vcf_suffix, bamlist)

##coverage
# calculate_genome_coverage(fq_dict, mkdup_bam, bedtools_genome_file, project_name)

##annotate with annovar
# filter_vcf_non_std_chr(result_vcf, filtered_vcf)
# convert_to_annovar_move_to_annovar_folder(fq_dict, filtered_vcf)
# run_table_annovar(fq_dict)
# multianno_to_annotated(fq_dict)

##map counts of snps hom in rescued and het/wt in unrescued and het in both parents
hom_snp_suffix = '.resc_hom.xls'
hom_bed_suffix = '.resc_hom.bed'
samples = ['TP1FRes']
window_size = [1000000, 500000, 100000]
step_size = 100000
# window_size = [100000]
# step_size = 10000
working_dir = work_dir
genome_fai = fasta_fai

##filtering
col_exons = [6, 11]
exon_definition = ['exonic', 'splicing']
col_functions = [9, 14]
# syn_definition = 'synonymous SNV'
zygosity_col = 21
cov_col = 23
cov_definition = 10
qual_col = 22
qual_definition = 30
##make more stringent
# cov_definition = 20
# qual_definition = 50

'''
for sample in samples:
	##remove if in repeat region or indel and cov/qual
	filtering_annotated.filter(working_dir, "and", sample + '.annotated.xls', sample + "21.temp", [16], ['=='], [''])
	filtering_annotated.filter(working_dir, "and", sample + "21.temp", sample + "22.temp", [4,5], ['!=','!='], ['-','-'])
	filtering_annotated.filter(working_dir, "and", sample + "22.temp", sample + "23.temp", [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
	##keep if hom in rescued
	filtering_annotated.filter(working_dir, "and", sample + "23.temp", sample + "24.temp", [zygosity_col], ['=='], ['hom'])
	##keep if not hom in unrescued and het in both parents
	filtering_annotated.filter(working_dir, "and", sample + "24.temp", sample + "25.temp", [19,20], ['==','=='], ['het','het'])
	filtering_annotated.filter(working_dir, "and", sample + "25.temp", sample + hom_snp_suffix, [18], ['!='], ['hom'])

##make bed from enu vars
for sample in samples:
	make_bed_from_ann_txt(sample + hom_snp_suffix, sample + hom_bed_suffix)
##make file for graphing
for sample in samples:
	for ws in window_size:
		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping_cybertron.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]
		print genome_and_window
		window_bed = genome_and_window + '.bed'
		##all shared snps
		out_bed = sample + '.resc_hom.' + genome_and_window + '.bed'
		##bedtools intersect 
		with open('temp.bed', "w") as naf_fh: 
			hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', window_bed, '-b', sample + hom_bed_suffix, '-c'], stdout=naf_fh)
			hom_bt_intersect.wait()
		##add header
		with open(out_bed, "w") as out_fh, open('temp.bed', "r") as in_fh:
			out_fh.write(delim.join(['chr', 'start', 'end', 'snp_number']) + '\n')
			for line in in_fh:
				##removing the chr from start of the line
				out_fh.write(line[3:])
'''


##map counts of snps hom wt in rescued and het/wt in unrescued and het in both parents
hom_snp_suffix = '.resc_wt.xls'
hom_bed_suffix = '.resc_wt.bed'
samples = ['TP1FUnr']
window_size = [1000000, 500000, 100000]
step_size = 100000
# window_size = [100000]
# step_size = 10000
working_dir = work_dir
genome_fai = fasta_fai

##filtering
col_exons = [6, 11]
exon_definition = ['exonic', 'splicing']
col_functions = [9, 14]
# syn_definition = 'synonymous SNV'
zygosity_col = 21
cov_col = 23
cov_definition = 10
qual_col = 22
qual_definition = 30
##make more stringent
# cov_definition = 20
# qual_definition = 50

# '''
for sample in samples:
	##remove if in repeat region or indel and cov/qual
	filtering_annotated.filter(working_dir, "and", sample + '.annotated.xls', sample + "21.temp", [16], ['=='], [''])
	filtering_annotated.filter(working_dir, "and", sample + "21.temp", sample + "22.temp", [4,5], ['!=','!='], ['-','-'])
	filtering_annotated.filter(working_dir, "and", sample + "22.temp", sample + "23.temp", [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
	##keep if wt i.e. no call in rescued
	filtering_annotated.filter(working_dir, "and", sample + "23.temp", sample + "24.temp", [17], ['=='], [''])
	##keep if hom or wt in unrescued and het in both parents
	filtering_annotated.filter(working_dir, "and", sample + "24.temp", sample + "25.temp", [19,20], ['==','=='], ['het','het'])
	filtering_annotated.filter(working_dir, "and", sample + "25.temp", sample + hom_snp_suffix, [18], ['!='], [''])

##make bed from enu vars
for sample in samples:
	make_bed_from_ann_txt(sample + hom_snp_suffix, sample + hom_bed_suffix)
##make file for graphing
for sample in samples:
	for ws in window_size:
		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping_cybertron.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]
		print genome_and_window
		window_bed = genome_and_window + '.bed'
		##all shared snps
		out_bed = sample + '.resc_wt.' + genome_and_window + '.bed'
		##bedtools intersect 
		with open('temp.bed', "w") as naf_fh: 
			hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', window_bed, '-b', sample + hom_bed_suffix, '-c'], stdout=naf_fh)
			hom_bt_intersect.wait()
		##add header
		with open(out_bed, "w") as out_fh, open('temp.bed', "r") as in_fh:
			out_fh.write(delim.join(['chr', 'start', 'end', 'snp_number']) + '\n')
			for line in in_fh:
				##removing the chr from start of the line
				out_fh.write(line[3:])
# '''
