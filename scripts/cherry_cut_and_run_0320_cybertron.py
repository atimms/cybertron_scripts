#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'
threads = '20'

'''
loaded biobuids as need access to:
bedtools
R
bowtie2
samtools

'''


##programs
seacr = '/home/atimms/programs/SEACR-master/SEACR_1.1.sh'
bg_to_bw = '/home/atimms/programs/bedGraphToBigWig'
picard = '/home/atimms/programs/picard.jar'
##files
bt_fa = '/home/atimms/ngs_data/references/hg38/hg38.fa'
bt_idx = '/home/atimms/ngs_data/references/hg38/hg38'
fa_info = '/home/atimms/ngs_data/references/hg38/hg38.chrom.sizes'
##methods


def map_using_bowtie2(outbam, fq_list, bowtie_idx):
	#options: --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700.
	paired_end_bowtie2 = subprocess.Popen(['bowtie2', '-x', bowtie_idx, '-p', '10', '-1', fq_list[0], '-2', fq_list[1], '--local', '--very-sensitive-local', 
		'--no-unal', '--no-mixed', '--no-discordant', '--phred33', '-I', '10', '-X', '700'], stdout=subprocess.PIPE)
	st_sam_bam = subprocess.Popen(['samtools', 'view', '-b', '-@', '5', '-o', outbam, '-'], stdin=paired_end_bowtie2.stdout)
	st_sam_bam.wait()

def bam_to_bedgraph(bam, out_bed, fasta_info):
	##combined
	with open(out_bed, 'w') as out_fh:
		bedtools_btb = subprocess.Popen(['bedtools', 'bamtobed', '-bedpe', '-i', bam], stdout=subprocess.PIPE)
		bedtools_sort = subprocess.Popen(['bedtools', 'sort', '-i', 'stdin'], stdin=bedtools_btb.stdout, stdout=subprocess.PIPE)
		bedtools_gc = subprocess.Popen(['bedtools', 'genomecov', '-bg', '-g', fasta_info,  '-i', 'stdin'], stdin=bedtools_sort.stdout, stdout=out_fh)
		bedtools_gc.wait()
	##individual steps
	# temp1_bed = out_bed.rsplit('.', 1)[0] + '.temp1.bed'
	# temp2_bed = out_bed.rsplit('.', 1)[0] + '.temp2.bed'
	# with open(temp1_bed, 'w') as t1_fh:
	# 	bedtools_btb = subprocess.Popen(['bedtools', 'bamtobed', '-bedpe', '-i', bam], stdout=t1_fh)
	# 	bedtools_btb.wait()
	# with open(temp2_bed, 'w') as t2_fh:
	# 	bedtools_sort = subprocess.Popen(['bedtools', 'sort', '-i',temp1_bed ], stdout=t2_fh)
	# 	bedtools_sort.wait()
	# with open(out_bed, 'w') as out_fh:
	# 	bedtools_gc = subprocess.Popen(['bedtools', 'genomecov', '-bg', '-g', fasta_info, '-i', temp2_bed], stdout=out_fh)
	# 	bedtools_gc.wait()		


def make_bigwig_file(bed, chrom_sizes, out_bw):
	# out_bw = bam.split('.')[0] + '.bigwig'
	with open('temp.bed', 'w') as t1_fh:
		bedtools_btb = subprocess.Popen(['cut', '-f1-4', bed], stdout=t1_fh)
		bedtools_btb.wait()	
	bgbw = subprocess.Popen([bg_to_bw, 'temp.bed', chrom_sizes, out_bw])
	bgbw.wait()

def run_seacr_ind(bedgraph, out_prefix, chrom_sizes):
	#bash SEACR_1.1.sh target.bedgraph 0.01 non stringent output
	run_seacr = subprocess.Popen(['bash', seacr, bedgraph, '0.01', 'non', 'stringent', out_prefix])
	run_seacr.wait()
	out_bed = out_prefix + '.stringent.bed'
	out_bw = out_prefix + '.stringent.bigwig'
	make_bigwig_file(out_bed, chrom_sizes ,out_bw)

def run_seacr_igg_control(bedgraph, igg_bg, out_prefix, chrom_sizes):
	#bash SEACR_1.1.sh target.bedgraph 0.01 non stringent output
	run_seacr = subprocess.Popen(['bash', seacr, bedgraph, igg_bg, 'norm', 'stringent', out_prefix])
	run_seacr.wait()
	out_bed = out_prefix + '.stringent.bed'
	out_bw = out_prefix + '.stringent.bigwig'
	make_bigwig_file(out_bed, chrom_sizes ,out_bw)
	run_seacr = subprocess.Popen(['bash', seacr, bedgraph, igg_bg, 'norm', 'relaxed', out_prefix])
	run_seacr.wait()
	out_bed = out_prefix + '.relaxed.bed'
	out_bw = out_prefix + '.relaxed.bigwig'
	make_bigwig_file(out_bed, chrom_sizes ,out_bw)

def run_cut_n_run(infile, work_dir):
	os.chdir(work_dir)
	with open(infile, 'r') as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc >1:
				line = line.rstrip().split(delim)
				sample_name = line[0]
				fq1 = line[1]
				fq2 = line[2]
				fq_files = [fq1, fq2]
				control = False
				if 'IgG' in sample_name:
					control = True
				print(sample_name, fq1, fq2, control)
				'''
				##map reads using bt2
				bt2_bam = sample_name + '.bt2.bam'
				map_using_bowtie2(bt2_bam, fq_files, bt_idx)
				##convert to bed then bedgraph
				sample_bedgraph = sample_name + '.bedGraph'
				bam_to_bedgraph(bt2_bam, sample_bedgraph, fa_info)
				##run seacr in different ways
				##all samples selecting the top 1% of regions by AUC
				seacr_out_auc1 = sample_name + '.auc1'
				seacr_out_igg = sample_name + '.IgG'
				run_seacr_ind(sample_bedgraph, seacr_out_auc1, fa_info)
				##using igg as control
				if 'IgG' not in sample_name:
					print(sample_bedgraph, 'IgG.bedGraph', seacr_out_igg)
					run_seacr_igg_control(sample_bedgraph, 'IgG.bedGraph', seacr_out_igg, fa_info)
				'''
				bt2_bam = sample_name + '.bt2.bam'
				sorted_bam = sample_name + '.bt2_sorted.bam'
				sort_index_bam(bt2_bam, sorted_bam)

				
def build_bowtie2_indexes(fa_file, name):
	bowtie2_index = subprocess.Popen(['bowtie2-build', fa_file, name])
	bowtie2_index.wait()

def sort_index_bam(in_bam, out_bam):
	# picard_rs = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'ReorderSam', 'INPUT=' + in_bam, 'OUTPUT=' + out_bam, 'CREATE_INDEX=true', 'REFERENCE=' + fa_file])
	picard_rs = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'SortSam', 'INPUT=' + in_bam, 'OUTPUT=' + out_bam, 'CREATE_INDEX=true', 'SORT_ORDER=coordinate'])
	picard_rs.wait()

##run methods

working_dir = '/home/atimms/ngs_data/misc/cherry_hu_ret_cutandrun_0320'
sample_file = 'cherry_hu_ret_cutandrun_0320.txt'
##biuld index for hg38
# build_bowtie2_indexes(bt_fa, bt_idx)
##run analysis
run_cut_n_run(sample_file, working_dir)



