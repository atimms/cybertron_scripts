#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import filtering_annotated
import homozygosity_mapping_cybertron
from collections import OrderedDict

##parameters
delim = '\t'
thread_number = '20'

##modules 
'''
use dave_analysis.pbs or
qsub -Iq cdbrmq -l mem=200gb,ncpus=20 -P 2cf4ba3b-f9ef-4cfc-a2c1-be85f5db6730
qsub -Iq cdbrmq -l mem=40gb,ncpus=2 -P 2cf4ba3b-f9ef-4cfc-a2c1-be85f5db6730

module load biobuilds/2017.11
module load gcc/6.1.0 ##gcc version used to make bedtools using in homo mapping
'''

##programs
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
bwa = 'bwa'
samtools = 'samtools'
bcftools = 'bcftools'
picard = 'picard'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
bgzip = 'bgzip'



##files
fasta = '/home/atimms/ngs_data/references/mm10/mm10.fa'
bamlist = 'bams.list'
# st_avinputs = ['K50000022.avinput', 'K50000024.avinput', 'K50000095.avinput', 'combined.avinput']
bedtools_genome_file  = '/home/atimms/ngs_data/references/mm10/mm10.fa.genome'

##file suffix
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
st_vcf_suffix = '.st.vcf.gz'

##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic', '-genericdbfile', 'mgp.v5.merged.snps_all.dbSNP142.avinput,xiaoer19.avinput,xiaoer26.avinput,xiaoer27.avinput,xiaoer28.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,', '-vcfinput']
#av_options = ['-otherinfo', '-remove', '-vcfinput']

##filtering and hom mapping parameters
##filtering
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = ['synonymous SNV', 'unknown']
cov_definition = 10
qual_definition = 30
##het snp mapping 
genome_fai = '/home/atimms/ngs_data/references/mm10/mm10.fa.fai'
# window_size = [1000000,5000000,2000000]
window_size = [10000000]
step_size = 1000000

##methods

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
	stmp = subprocess.Popen([samtools,'mpileup', '-Ou','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist_file], stdout=subprocess.PIPE)
	# bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '10', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
	bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '10', '-o', vcf_temp1], stdin=stmp.stdout)
	bcft.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	# '''
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()

def filter_vcf_non_std_chr(in_vcf, out_vcf):
	vcf_temp1 = 'temp1.vcf.gz'
	vcf_temp2 = 'temp2.vcf.gz'
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-r', 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19', '-o', vcf_temp1, '-O', 'z', in_vcf])
	bcftools_filter.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', out_vcf, vcf_temp2])
	bcf_norm2.wait()

##make avinput files
def convert_to_annovar_move_to_annovar_folder(samples, vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
	con_ann.wait()
	for sample in samples:
		temp_av_file = 'temp.' + sample + '.avinput'
		av_file = sample + '.avinput'
		shutil.move(temp_av_file, av_file)
		shutil.copy(av_file, str(av_ref_dir[0]))

def run_table_annovar(vcf, av_prefix):
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def multianno_to_annotated(av_prefix): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142', 'mgp.v5.snps', 'xiaoer19', 'xiaoer26', 'xiaoer27', 'xiaoer28',  'QUAL', 'vcf_info', 'vcf_format', 'vcf_26', 'vcf_27', 'vcf_19', 'vcf_28', 'dp_26', 'dp_27', 'dp_19', 'dp_28']
	head_out = delim.join(head + ['\n'])
	multianno = av_prefix + '.mm10_multianno.txt'
	annotated = av_prefix + '.annotated.txt'
	with open(multianno, "r") as multi, open(annotated, "w") as final:
		final.write(head_out)
		line_count = 0
		for line in multi:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count > 1:
				# print(line[29:])
				dp1 = line[30].split(":")[2]
				dp2 = line[31].split(":")[2]
				dp3 = line[32].split(":")[2]
				dp4 = line[33].split(":")[2]
				line_out = line[:18] + [line[19]] + line[28:] + [dp1, dp2, dp3, dp4]
				final.write(delim.join(line_out) + '\n')

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

def filter_ann_file_exonic(file_prefix):
	##exonic_variants
	filtering_annotated.filter(working_dir, "or", file_prefix + '.annotated.txt' , file_prefix + "_1.temp", [col_exon, col_exon], ['==','=='], exon_definition)
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", file_prefix + "_1.temp", file_prefix + "_2.temp", [col_function,col_function], ['!=','!='], syn_definition)
	##in all 4 samples
	filtering_annotated.filter(working_dir, "and", file_prefix + "_2.temp", file_prefix + "_3.temp", [15,16,17,18], ['!=','!=','!=','!='], ['.','.','.','.'])
	filtering_annotated.filter(working_dir, "and", file_prefix + "_2.temp", file_prefix + "_31.temp", [15,16,17,18], ['==','==','==','=='], ['hom','hom','hom','hom'])

	##not in dbSNP or mgp 
	filtering_annotated.filter(working_dir, "and", file_prefix + "_3.temp", file_prefix + ".exonic_rare.xls", [13,14], ['==','=='], ['.','.'])
	filtering_annotated.filter(working_dir, "and", file_prefix + "_31.temp", file_prefix + ".hom_exonic_rare.xls", [13,14], ['==','=='], ['.','.'])

def hom_mapping(file_prefix, sample_names):
	##remove if in rmsk or segdups
	filtering_annotated.filter(working_dir, "and", file_prefix + '.annotated.txt', file_prefix + "11.temp", [11,12], ['==','=='], ['.','.'])
	##remove indels
	filtering_annotated.filter(working_dir, "and", file_prefix + "11.temp", file_prefix + "12.temp", [4,5], ['!=','!='], ['-','-'])
	##filter variants by coverage and quality 
	filtering_annotated.filter(working_dir, "and", file_prefix + "12.temp", file_prefix + '13.temp', [19], ['>='], [qual_definition])
	filtering_annotated.filter(working_dir, "and", file_prefix + "13.temp", file_prefix + '.hom_temp.txt', [26,27,28,29], ['>=','>=','>=','>='], [cov_definition,cov_definition,cov_definition,cov_definition])
	for sample in sample_names:
		shutil.copy(file_prefix + '.hom_temp.txt', sample + '.hom_temp.txt')
		for ws in window_size:
			#make bed file with windows and returns genome name and window size variable
			genome_and_window = homozygosity_mapping_cybertron.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]
			if sample == 'xiaoer19':
				zygosity_col = 15
				info_col = 24
			if sample == 'xiaoer26':
				zygosity_col = 16
				info_col = 22
			if sample == 'xiaoer27':
				zygosity_col = 17
				info_col = 23
			if sample == 'xiaoer28':
				zygosity_col = 18
				info_col = 25
			##make bed file from variants
			homozygosity_mapping_cybertron.make_bed_from_ann(working_dir, 'samtools', sample + '.hom_temp.txt', zygosity_col, info_col)
			##hom and het count and hom percentage
			homozygosity_mapping_cybertron.count_and_percentage(working_dir, genome_and_window, sample + '.bed')
			##naf
			homozygosity_mapping_cybertron.naf_in_window(working_dir, genome_and_window, sample + '.bed')
			##total snp number
			homozygosity_mapping_cybertron.total_snp_in_window(working_dir, genome_and_window, sample + '.bed')
			##combine bedgraphs for graphing in r
			homozygosity_mapping_cybertron.combine_bedgraphs_for_r(working_dir, sample, genome_and_window)

##run_methods
working_dir = '/home/atimms/ngs_data/enu_mapping/dave_xiaoer_1221'
os.chdir(working_dir)

project_name = 'dave_xiaoer_1221'
fq_dict = {'xiaoer19':['Xiaoer_19_CSFP210019943-1a_H23WLDSX3_L2_1.fq.gz', 'Xiaoer_19_CSFP210019943-1a_H23WLDSX3_L2_2.fq.gz'],
		'xiaoer26':['Xiaoer_26_CSFP210019944-1a_H23WLDSX3_L2_1.fq.gz', 'Xiaoer_26_CSFP210019944-1a_H23WLDSX3_L2_2.fq.gz'],
		'xiaoer27':['Xiaoer_27_CSFP210019946-1a_H23WLDSX3_L2_1.fq.gz', 'Xiaoer_27_CSFP210019946-1a_H23WLDSX3_L2_2.fq.gz'],
		'xiaoer28':['Xiaoer_28_CSFP210019945-1a_H23WLDSX3_L2_1.fq.gz', 'Xiaoer_28_CSFP210019945-1a_H23WLDSX3_L2_2.fq.gz']}
samples = fq_dict.keys()
result_vcf = project_name + st_vcf_suffix
filtered_vcf = project_name + '.st.autosomes.vcf.gz'

##analysis wanted 
## exonic variants that were shared by all the samples
## common region of homozygosity in all samples


##align all samples
# align_with_bwa(fq_dict)

##call vars
# make_list_of_bams(samples, mkdup_bam, bamlist)
# variant_calling_samtools(project_name, st_vcf_suffix, bamlist)

##coverage
# calculate_genome_coverage(samples, mkdup_bam, bedtools_genome_file, project_name)

##annotate and format
# filter_vcf_non_std_chr(result_vcf, filtered_vcf)
# convert_to_annovar_move_to_annovar_folder(samples, filtered_vcf)
# run_table_annovar(filtered_vcf, project_name)
# multianno_to_annotated(project_name)

##get exonic variants in all 4 samples
# filter_ann_file_exonic(project_name)

##homozygosity mapping for the 4 samples
# hom_mapping(project_name, samples)


def get_shared_hom_regions(file_prefix, sample_names, infile_suffix):
	tmp_beds = []
	for sample in sample_names:
		in_bed = sample + infile_suffix
		filtered_bed = sample + '_naf90.temp.bed'
		tmp_beds.append(filtered_bed)
		
		with open(in_bed, "r") as in_fh, open(filtered_bed, "w") as out_fh:
			for line in in_fh:
				line = line.rstrip().split(delim)
				naf = float(line[3])
				count = int(line[4])
				if naf > 0.9 and count >= 10:
					out_fh.write(delim.join(line) + '\n')
	bt_int_bed = file_prefix + '_naf90.bt_int.temp.bed'
	with open(bt_int_bed, "w") as naf_fh: 
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', tmp_beds[0], '-b'] + tmp_beds[1:] + ['-C', '-filenames'], stdout=naf_fh)
		hom_bt_intersect.wait()
	bed_dict = {}
	with open(bt_int_bed, "r") as in_fh:
		 for line in in_fh:
		 	line = line.rstrip().split(delim)
		 	region = '_'.join(line[:3])
		 	naf = line[3]
		 	count = line[4]
		 	intersect = line[6]
		 	if region in bed_dict:
		 		if intersect != '0':
		 			bed_dict[region][2] += 1
		 	else:
		 		if intersect == '0':
		 			bed_dict[region] = [naf, count, 0]
		 		else:
		 			bed_dict[region] = [naf, count, 1]
	shared_bed = file_prefix + '_naf90.bt_int.in_all.bed'
	with open(shared_bed, "w") as out_fh:
		out_fh.write(delim.join(['chr', 'start', 'finish', 'naf', 'snp_count']) + '\n')
		for r in bed_dict:
			if bed_dict[r][2] == 3:
				out_fh.write(delim.join(r.split('_') + bed_dict[r][:2]) + '\n')





##get shared regions i.e. naf >90 in all 4 samples
bed_suffix = '_mm10_10000kb_1000kb_naf.bedgraph'
get_shared_hom_regions(project_name, samples, bed_suffix)


