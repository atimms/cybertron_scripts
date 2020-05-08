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
module load gcc/6.1.0 ##gcc version used to make bedtools using in homo mapping
'''


##working dir
working_dir = '/home/atimms/ngs_data/enu_mapping/sergei_enu_0120'
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
bgzip = 'bgzip'
delly = '/home/atimms/programs/delly/src/delly'


##files
bamslist_file = 'bams.list'
# st_avinputs = ['K50000022.avinput', 'K50000024.avinput', 'K50000095.avinput', 'combined.avinput']
exome_bed_file = '/home/atimms/ngs_data/references/mm10/mm10_refGene_coding_exons.bed'
delly_exclude_regions = '/home/atimms/programs/delly/excludeTemplates/mouse.mm10.excl.tsv'

##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic', '-genericdbfile', 'mgp.v5.merged.snps_all.dbSNP142.avinput,C3H_HeH.mgp.v5.snps.dbSNP142.avinput,C3H_HeJ.mgp.v5.snps.dbSNP142.avinput,0120.M392.avinput,0120.M398.avinput,0120.M409.avinput,0120.M558.avinput,0120.1493.avinput,0120.1524.avinput,0120.1525.avinput,0120.1526.avinput,0120.C3H.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,,,,,,,', '-vcfinput']
#av_options = ['-otherinfo', '-remove', '-vcfinput']


##filtering and hom mapping parameters
##filtering
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = 'synonymous SNV'
zygosity_col = 28
cov_col = 30
cov_definition = 10
qual_col = 29
qual_definition = 30
##het snp mapping 
genome_fai = '/home/atimms/ngs_data/references/mm10/mm10.fa.fai'
window_size = [1000000,5000000,2000000]
# window_size = [10000000]
step_size = 1000000
info_col = 40

##methods
##make list of all bam files to be analyzed
def make_list_of_bams(bam_files, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')


##call samtools on bamfiles
def variant_calling_samtools(project, final_vcf_suffix, bamlist_file, region):
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
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', '0120'])
	con_ann.wait()
	for sample in samples:
		av_file = '0120.' + sample + '.avinput'
		shutil.copy(av_file, str(av_ref_dir[0]))

def run_table_annovar(vcf, av_prefix):
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def multianno_to_annotated(av_prefix): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142','mgp.v5.snps', 'C3H_HeH', 'C3H_HeJ', 'M392', 'M398', 'M409', 'M558', '1493', '1524', '1525', '1526', 'C3H', 'Zygosity', 'Qual', 'Coverage', 'Chr', 'Pos', 'Filter', 'Ref2', 'Alt2','Qual', 'Filter', 'GT_info', 'Format','1493', '1524', '1525', '1526', 'C3H', 'M392', 'M398', 'M409', 'M558']
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


def filter_rpt_c3h(file_prefix):
	##remove if in rmsk, segdup
	filtering_annotated.filter(working_dir, "and", file_prefix + '.annotated.txt', file_prefix + "11.temp", [11,12], ['==','=='], ['.','.'])
	##in c3hr mouse
	filtering_annotated.filter(working_dir, "and", file_prefix + '11.temp', file_prefix + ".c3h.xls", [25], ['!='], ['.'])

def filter_for_enu(file_prefix):
	##remove if in rmsk, segdup, dbsnp, mgp, c3h
	filtering_annotated.filter(working_dir, "and", file_prefix + '.annotated.txt', file_prefix + ".enu_rpts.xls", [11,12,13,14,25], ['==','==','==','==','=='], ['.','.','.','.','.'])
	##remove if in  dbsnp, mgp, c3h
	filtering_annotated.filter(working_dir, "and", file_prefix + '.annotated.txt', file_prefix + ".enu.xls", [13,14,25], ['==','==','=='], ['.','.','.'])

##run_methods
project_name = 'sergei_0120_c3'
bams = ['1493.bwa_mkdup.bam', '1524.bwa_mkdup.bam', '1525.bwa_mkdup.bam', '1526.bwa_mkdup.bam', 'C3H.bwa_mkdup.bam', 
		'M392.bwa_mkdup.bam', 'M398.bwa_mkdup.bam', 'M409.bwa_mkdup.bam', 'M558.bwa_mkdup.bam']
sample_names = [i.split('.')[0] for i in bams]
bamlist = 'bams.list'
roi = 'chr3:100000000-160039680'
result_vcf_suffix = '.vcf.gz'
result_vcf = project_name + result_vcf_suffix
mutants = ['1524', '1525', '1526']
##call vars in chr3 region
# make_list_of_bams(bams, bamlist)
# variant_calling_samtools(project_name, result_vcf_suffix, bamlist, roi)

##annotate and filter
# convert_to_annovar_move_to_annovar_folder(sample_names, result_vcf)
# run_table_annovar(result_vcf, project_name)
multianno_to_annotated(project_name)
filter_rpt_c3h(project_name)
filter_for_enu(project_name)
