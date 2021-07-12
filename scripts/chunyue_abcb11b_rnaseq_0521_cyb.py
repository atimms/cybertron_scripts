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
av_protocol = ['-protocol', 'refGene,ensGene,rmsk,generic,generic,generic,generic,generic,generic,generic,generic,generic,bed', '-genericdbfile', 'RM1.avinput,RM2.avinput,RM3.avinput,UM1.avinput,UM2.avinput,UM3.avinput,WT1.avinput,WT2.avinput,WT3.avinput', '-bedfile', 'rm_snps_0521.bed']
av_operation = ['-operation', 'g,g,r,f,f,f,f,f,f,f,f,f,r']
# av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,', '-vcfinput']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,,,,', '--nopolish']



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
	stmp = subprocess.Popen(['samtools','mpileup', '-ug','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist_file], stdout=subprocess.PIPE)
	bcft = subprocess.Popen(['bcftools','call', '-vmO', 'z', '--threads', '9', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
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

##make avinput files
def convert_to_annovar_move_to_annovar_folder(samples, vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
	con_ann.wait()
	for sample in samples:
		temp_av_file = 'temp.' + sample + '.Aligned.sortedByCoord.out.bam.avinput'
		av_file = sample + '.avinput'
		shutil.move(temp_av_file, av_file)
		shutil.copy(av_file, str(av_ref_dir[0]))

def run_table_annovar(samples):
	for sample in samples:
		avinput = sample + '.avinput'
		av_prefix = avinput.rsplit('.',1)[0]
		command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()

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
					header = line[:16] + ['RM1', 'RM2', 'RM3', 'UM1', 'UM2', 'UM3', 'WT1', 'WT2', 'WT3', 'WGS_SNPs', 'zygosity', 'qual', 'cov', 'info', 'format', 'vcf']
					final.write(delim.join(header) + '\n')
				else:
					line_out = line[:29] + line[36:] 
					final.write(delim.join(line_out) + '\n')


def bcftools_int_wgs_snps(in_vcf, bed, out_vcf):
	##bedtools intersect 
	with open(out_vcf, "w") as out_fh: 
		hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', in_vcf, '-b', bed], stdout=out_fh)
		hom_bt_intersect.wait()

def make_bed_from_ann_txt(infile, outfile):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				out_fh.write(delim.join(line[:3]) + '\n')


##run methods
work_dir = '/home/atimms/ngs_data/rnaseq/chunyue_abcb11b_rnaseq_0521'
os.chdir(work_dir)


##params
project_name = 'chunyue_abcb11b_rnaseq_0521'
star_bam = '.Aligned.sortedByCoord.out.bam'
bamlist = 'bams.list'
st_vcf_suffix = '.st.vcf.gz'
result_vcf = project_name + st_vcf_suffix
filtered_vcf = project_name + '.st_filtered_chr.vcf.gz'
samples = ['RM1', 'RM2', 'RM3', 'UM1', 'UM2', 'UM3', 'WT1', 'WT2', 'WT3']
samples_to_annotate = ['RM1', 'RM2', 'RM3']
wgs_snp_bed = '/home/atimms/ngs_data/references/annovar/danRer11/rm_snps_0521.bed'
wgs_snp_vcf = project_name + '.st.wgs_snps.vcf'

##map with bwa and process with samtools etc
# make_list_of_bams(samples, star_bam, bamlist)
# variant_calling_samtools(project_name, st_vcf_suffix, bamlist)

##annotate with annovar, use bed with snps from wgs
# filter_vcf_non_std_chr(result_vcf, filtered_vcf)
# convert_to_annovar_move_to_annovar_folder(samples, filtered_vcf)
# run_table_annovar(samples_to_annotate)
# multianno_to_annotated(samples_to_annotate)
##filter vcf for wgs snps
# bcftools_int_wgs_snps(filtered_vcf, wgs_snp_bed, wgs_snp_vcf)



##hom mapping for the three RMs
hom_snp_suffix = '.resc_hom.xls'
hom_bed_suffix = '.resc_hom.bed'
# window_size = [1000000, 500000, 100000]
# step_size = 100000
window_size = [100000]
step_size = 10000
zygosity_col = 27
cov_col = 29
cov_definition = 20
qual_col = 28
qual_definition = 30
working_dir = work_dir
genome_fai = fasta_fai

for sample in samples_to_annotate:
	##remove if in repeat region or indel and cov/qual
	filtering_annotated.filter(working_dir, "and", sample + '.annotated.xls', sample + "21.temp", [16], ['=='], [''])
	filtering_annotated.filter(working_dir, "and", sample + "21.temp", sample + "22.temp", [4,5], ['!=','!='], ['-','-'])
	filtering_annotated.filter(working_dir, "and", sample + "22.temp", sample + "23.temp", [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
	##keep if hom not hom in unrescued
	filtering_annotated.filter(working_dir, "and", sample + "23.temp", sample + "24.temp", [zygosity_col], ['=='], ['hom'])
	##keep if not hom in unrescued
	filtering_annotated.filter(working_dir, "and", sample + "24.temp", sample + hom_snp_suffix, [20,21,22], ['!=','!=','!='], ['hom','hom','hom'])

##make bed from enu vars
for sample in samples_to_annotate:
	make_bed_from_ann_txt(sample + hom_snp_suffix, sample + hom_bed_suffix)
##make file for graphing
for sample in samples_to_annotate:
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
