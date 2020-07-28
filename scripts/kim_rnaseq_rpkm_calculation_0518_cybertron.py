#!/usr/bin/env python
import sys
import subprocess
import os
import glob


"""
to run:
module load java/1.8.0_121
module load biobuilds/2016.11
"""

##parameters
delim = '\t'
threads = '16'

##setup working directory where results will be
working_dir = '/home/atimms/ngs_data/misc/kim_rpkm_rnaseq_0518'
os.chdir(working_dir)

##programs
rnaseqc = '/home/atimms/programs/RNA-SeQC_v1.1.8.jar'
##ref files etc
genome_name = 'hg19'
# genome_name = 'GRCh38'
fa_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name +  '/genes.gtf'
rrna_gtf = 'hg19_rRNA.gtf'



def add_rg_sort_index_bam(sample_names, in_bam_suffix, out_bam_suffix):
	for sample_name in sample_names:
		in_bam = sample_name + in_bam_suffix
		arg_bam = sample_name + '.with_rg.bam'
		out_bam = sample_name + out_bam_suffix
		picard_arg = subprocess.Popen(['picard', 'AddOrReplaceReadGroups', 'INPUT=' + in_bam, 'OUTPUT=' + arg_bam, 'SO=coordinate', 'RGID=' + sample_name, 'RGLB=' + sample_name, 'RGPL=Illumina', 'RGPU=machine1', 'RGSM=' + sample_name ])
		picard_arg.wait()	
		# picard_rs = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'ReorderSam', 'INPUT=' + in_bam, 'OUTPUT=' + out_bam, 'CREATE_INDEX=true', 'REFERENCE=' + fa_file])
		picard_rs = subprocess.Popen(['picard', 'SortSam', 'INPUT=' + arg_bam, 'OUTPUT=' + out_bam, 'CREATE_INDEX=true', 'SORT_ORDER=coordinate'])
		picard_rs.wait()

def run_rnaseqc(sample_file, out_dir):
	# module_rm = subprocess.Popen(['module', 'rm', 'java/1.8.0_121'])
	# module_rm.wait()
	rnaseqc_cmd = subprocess.Popen(['java', '-jar', rnaseqc, '-s', sample_file, '-t', gtf_file, '-r', fa_file, '-o', out_dir, '-rRNA', rrna_gtf])
	rnaseqc_cmd.wait()

def run_rnaseqc_no_rrna_gtf(sample_file, out_dir):
	# module_rm = subprocess.Popen(['module', 'rm', 'java/1.8.0_121'])
	# module_rm.wait()
	rnaseqc_cmd = subprocess.Popen(['java', '-jar', rnaseqc, '-s', sample_file, '-t', gtf_file, '-r', fa_file, '-o', out_dir])
	rnaseqc_cmd.wait()

def combine_rqc_rpkm_files(sample_names, outfile, results_dir):
	gene_super_dict = {}
	for sample in sample_names:
		ind_gene_dict = {}
		infile = results_dir + sample + '/' + sample + '.metrics.tmp.txt.rpkm.gct'
		line_count = 0
		with open(infile, "r") as in_fh:
			for line in in_fh:
				line_count += 1
				if line_count > 3:
					line = line.rstrip().split(delim)
					gene = line[1]
					rpkm = float(line[2])
					if gene in ind_gene_dict:
						ind_gene_dict[gene].append(rpkm)
					else:
						ind_gene_dict[gene] = [rpkm]
		for g in ind_gene_dict:
			# print g, ind_gene_dict[g]
			ave_rpkm = sum(ind_gene_dict[g]) / len(ind_gene_dict[g])
			if g in gene_super_dict:
				gene_super_dict[g].append(str(ave_rpkm))
			else:
				gene_super_dict[g] = [str(ave_rpkm)]

	print 'dict contains %s genes'%len(gene_super_dict)


	with open(outfile, "w") as out_fh:
		out_fh.write(delim.join(['gene'] + sample_names + ['\n']))
		for g in gene_super_dict:
			if g != "":
				print g, gene_super_dict[g]
				out_fh.write(delim.join([g] + gene_super_dict[g] + ['\n']))


def add_rg_index_bam(sample_names, in_bam_suffix, out_bam_suffix):
	for sample_name in sample_names:
		in_bam = sample_name + in_bam_suffix
		out_bam = sample_name + out_bam_suffix
		picard_arg = subprocess.Popen(['picard', 'AddOrReplaceReadGroups', 'INPUT=' + in_bam, 'OUTPUT=' + out_bam, 'SO=coordinate', 'RGID=' + sample_name, 'RGLB=' + sample_name, 'RGPL=Illumina', 'RGPU=machine1', 'RGSM=' + sample_name ])
		picard_arg.wait()
		samtools_index = subprocess.Popen(['samtools', 'index', out_bam])
		samtools_index.wait()

##run methods
##parameters

##0518
'''
samples = ['H24551_15_3pcw_EGL', 'H24616_21_7pcw_EGL', 'H25174_EGL1', 'H25174_EGL2', 'H26122_EGL1', 'H26122_EGL2', 'H26360_Gr2', 'H26360_GrA_2', 
		'H26360_GrB_1', 'H26362_EGL1', 'H26362_GrB_1', 'H26374_Gr1', 'H26374_GrB_1', 'H26499_PKEGL1', 'H26566_EGL1', 'H26566_EGL2', 'H26857_EGL1', 
		'H26938_EGL1rep', 'H24551_15_3pcw_PK', 'H24616_21_7pcw_PK', 'H25174_PK1_1', 'H25174_PK2_1', 'H26122_PK1_1', 'H26122_PK2_2', 'H26360_Pk2', 
		'H26360_PkA_2', 'H26360_PkB_1', 'H26362_PK1_1', 'H26362_PkA_1', 'H26374_Pk1', 'H26499_PK1_1', 'H26499_PK1_1rep', 'H26566_PK1_1', 'H26566_PK1_1rep', 
		'H26589_PK1_1', 'H26857_PK1_1', 'H26938_PK1_1', 'H26938_PK1_1rep', 'H24551_15_3pcw_WS', 'H24616_21_7pcw_WS', 'H25174_Wver', 'H26122_Wver', 
		'H26360_W2', 'H26360_Whm', 'H26374_W1', 'H26566_Wver', 'H26589_Wver', 'H26712_A', 'H26712_B', 'H26857_Whm', 'H26857_Wver', 'H26938_Wver', 
		'H24551_RL_combined', 'H24616_RL_combined', 'H26362_RL_combined', 'H26566_RL_combined', 'H26857_RL_combined', 'K13521_RL_combined', 'K13532_RL_combined', 
		'K13855_RL_combined', 'K13875_RL_combined']
# samples = ['H24551_15_3pcw_EGL']

sample_info = 'rnaseqc_sample_info.txt'
outdir = 'rnaseqc_results'
star_bam_suffix = '.Aligned.out.bam'
sorted_bam_suffix = '.rg_sorted.bam'
rnaseqc_rpkm_file = 'kim_rnaseq_0518.rpkm.txt'

##sort and add read group -- must have java 8
# add_rg_sort_index_bam(samples, star_bam_suffix, sorted_bam_suffix)

##manually make rna-seqc sample file, 3 column sample name, bam, notes

##run rna-seqc on bam to get rpkm (can't use java 8)
# run_rnaseqc(sample_info, outdir)
# combine_rqc_rpkm_files(samples, rnaseqc_rpkm_file)
'''


##0618 hipscs_rnaseq
'''
samples = ['control78_13', 'control78_14', 'control78_15', 'GSM2071619', 'GSM2071620', 'null77_7', 'null77_8', 'null77_9', 'patient77_10', 'patient77_11', 'patient77_12', 'SAH01px_1', 'SAH01px_2', 'SAH01px_3', 'SAH02_control_4', 'SAH02_control_5', 'SAH02_control_6']
# samples = ['H24551_15_3pcw_EGL']

sample_info = 'rnaseqc_sample_info_0618.txt'
outdir = 'rnaseqc_results_0618/'
star_bam_suffix = '.Aligned.out.bam'
sorted_bam_suffix = '.rg_sorted.bam'
rnaseqc_rpkm_file = 'kim_hipscs_rnaseq_0618.rpkm.txt'

##sort and add read group -- must have java 8
# add_rg_sort_index_bam(samples, star_bam_suffix, sorted_bam_suffix)

##manually make rna-seqc sample file, 3 column sample name, bam, notes

##run rna-seqc on bam to get rpkm (can't use java 8)
# run_rnaseqc(sample_info, outdir)
# combine_rqc_rpkm_files(samples, rnaseqc_rpkm_file, outdir)
'''

##1018 kim one sample
# '''

samples = ['H26362-Wver-tube10']

sample_info = 'rnaseqc_sample_info_1018.txt'
outdir = 'rnaseqc_results_1018/'
star_bam_suffix = '.Aligned.out.bam'
sorted_bam_suffix = '.rg_sorted.bam'
rnaseqc_rpkm_file = 'kim_rnaseq_1018.rpkm.txt'

##sort and add read group -- must have java 8
# add_rg_sort_index_bam(samples, star_bam_suffix, sorted_bam_suffix)

##manually make rna-seqc sample file, 3 column sample name, bam, notes

##run rna-seqc on bam to get rpkm (can't use java 8)
# run_rnaseqc(sample_info, outdir)
# combine_rqc_rpkm_files(samples, rnaseqc_rpkm_file, outdir)


##1018 kim new set of rnaseq sample
'''
##setup working directory where results will be
working_dir = '/home/atimms/ngs_data/rnaseq/kim_rnaseq_hb_1018'
os.chdir(working_dir)

samples = ['ASD_1347_LMD_PK_S2', 'ASD_1347_Section_S1', 'ASD_1714_LMD_PK_S4', 'ASD_1714_Section_S3', 
		'ASD_5419_LMD_PK_S6', 'ASD_5419_Section_S5', 'ASD_5565_LMD_PK_S8', 'ASD_5565_Section_S7']

sample_info = 'kim_1018.rnaseqc_sample_info.txt'
outdir = 'rnaseqc_results_kim_hb_1018/'
star_bam_suffix = '.Aligned.out.bam'
sorted_bam_suffix = '.rg_sorted.bam'
rnaseqc_rpkm_file = 'kim_rnaseq_hb_1018.rpkm.txt'

##sort and add read group -- must have java 8
# add_rg_sort_index_bam(samples, star_bam_suffix, sorted_bam_suffix)

##manually make rna-seqc sample file, 3 column sample name, bam, notes

##run rna-seqc on bam to get rpkm (can't use java 8)
run_rnaseqc(sample_info, outdir)
combine_rqc_rpkm_files(samples, rnaseqc_rpkm_file, outdir)
'''

##0720 kim new set of rnaseq samples
##setup working directory where results will be
working_dir = '/home/atimms/ngs_data/rnaseq/kim_rnaseq_0720/hg19_alignment'
os.chdir(working_dir)

samples = ['Ctrl_13086_Sec', 'Ctrl_13086_Gr', 'Ctrl_13086_PK', 'Ctrl_13086_RL', 'DW22840EE_Sec', 'DW22840EE_Gr', 
		'DW22840EE_RL', 'DW05950_Sec', 'DW05950_CP', 'DW05950_Gr', 'DW05950_Mes', 'DW05950_PK', 'DW05950_RL', 
		'Ctrl_H26074_Sec', 'Ctrl_H26074_RL']

sample_info = 'kim_0720.rnaseqc_sample_info.txt'
outdir = 'rnaseqc_results_kim_0720/'
star_bam_suffix = '.Aligned.sortedByCoord.out.bam'
sorted_bam_suffix = '.rg_sorted.bam'
rnaseqc_rpkm_file = 'kim_rnaseq_0720_b6.rpkm.txt'

##sort and add read group -- must have java 8
# add_rg_sort_index_bam(samples, star_bam_suffix, sorted_bam_suffix)

##manually make rna-seqc sample file, 3 column sample name, bam, notes

##run rna-seqc on bam to get rpkm (can't use java 8)
# run_rnaseqc_no_rrna_gtf(sample_info, outdir)
combine_rqc_rpkm_files(samples, rnaseqc_rpkm_file, outdir)



