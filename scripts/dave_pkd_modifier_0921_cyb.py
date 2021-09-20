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
bedtools = 'bedtools'
bgzip = 'bgzip'



##files
fasta = '/home/atimms/ngs_data/references/mm10/mm10.fa'
bamlist = 'bams.list'
# st_avinputs = ['K50000022.avinput', 'K50000024.avinput', 'K50000095.avinput', 'combined.avinput']

##file suffix
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
st_vcf_suffix = '.st.vcf.gz'

##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic,generic', '-genericdbfile', 'mgp.v5.merged.snps_all.dbSNP142.avinput,0620.kc1_wt.avinput,0620.kc2_ko.avinput,129P2_OlaHsd.mgp.v5.snps.dbSNP142.avinput,129S1_SvImJ.mgp.v5.snps.dbSNP142.avinput,129S5SvEvBrd.mgp.v5.snps.dbSNP142.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,', '-vcfinput']
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
	bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '10', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
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


##make avinput files
def convert_to_annovar_move_to_annovar_folder(samples, vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', '0620'])
	con_ann.wait()
	for sample in samples:
		av_file = '0620.' + sample + '.avinput'
		shutil.copy(av_file, str(av_ref_dir[0]))

def run_table_annovar(vcf, av_prefix):
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def multianno_to_annotated(av_prefix): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142', 'mgp.v5.snps', 'kc1_wt', 'kc2_ko', '129P2_OlaHsd', '129S1_SvImJ', '129S5SvEvBrd', 'vcf_info', 'vcf_format', 'vcf_kc1_wt', 'vcf_kc2_ko']
	head_out = delim.join(head + ['\n'])

	multianno = av_prefix + '.mm10_multianno.txt'
	annotated = av_prefix + '.annotated.txt'
	with open(multianno, "r") as multi, open(annotated, "w") as final:
		final.write(head_out)
		line_count = 0
		for line in multi:
			line = line.split(delim)
			line_count += 1
			if line_count > 1:
				line_out = line[:18] + line[28:]
				final.write(delim.join(line_out))


##run_methods
working_dir = '/home/atimms/ngs_data/enu_mapping/dave_pkd_modifier_0921'
os.chdir(working_dir)

project_name = 'dave_pkd1rc_modifier_0921'
fq_dict = {'EP0014-0259-0008':['EP0014-0259-0008_R1_001.fastq.gz', 'EP0014-0259-0008_R2_001.fastq.gz'],
		'EP0014-0259-0009':['EP0014-0259-0009_R1_001.fastq.gz', 'EP0014-0259-0009_R2_001.fastq.gz'],
		'EP0014-0259-0010':['EP0014-0259-0010_R1_001.fastq.gz', 'EP0014-0259-0010_R2_001.fastq.gz']}
samples = fq_dict.keys()
result_vcf = project_name + st_vcf_suffix


##analysis wanted 
## exonic variants that were shared by all the sample
## common region of homozygosity in all samples


##align all samples
align_with_bwa(fq_dict)


##call vars
make_list_of_bams(samples, mkdup_bam, bamlist)
variant_calling_samtools(project_name, st_vcf_suffix, bamlist)

##annotate and filter
# convert_to_annovar_move_to_annovar_folder(samples, result_vcf)
# run_table_annovar(result_vcf, project_name)
# multianno_to_annotated(project_name)
# filter_ann_file(project_name)