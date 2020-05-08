#!/usr/bin/env python
import sys
import subprocess
import os
import glob

'''
##not used now
module load biobuilds/2017.11
loaded biobuids as need access to:
bedtools

if using strand aware i.e. 2 files


0420:
module load local_python/2.7.14
conda create --name make_bigwig
source activate make_bigwig
and make sure you have deeptools i.e. 
conda install -c bioconda deeptools

from now on:
module load local_python/2.7.14
source activate make_bigwig


'''


##parameters
delim = '\t'
threads = '20'

##programs
bg_to_bw = '/home/atimms/programs/bedGraphToBigWig'

##refs
##had to make this myself from fa file headers (weird chr names used)
hg38_chrom_sizes = '/home/atimms/ngs_data/references/igenomes/GRCh38/genome.chrom_sizes'
hg38_fa = '/home/atimms/ngs_data/references/igenomes/GRCh38/genome.fa'

##this is not strand aware i.e. one file per 
def convert_bams_bigwig(work_dir, bams, chrom_sizes, fasta):
	os.chdir(work_dir)
	for bam in bams:
		out_bed = bam.split('.')[0] + '.bedGraph'
		out_bw = bam.split('.')[0] + '.bigwig'
		with open(out_bed, 'w') as out_fh:
			bedtools_gc = subprocess.Popen(['bedtools', 'genomecov', '-g', chrom_sizes,  '-bg', '-ibam', bam], stdout=subprocess.PIPE)
			bedtools_sort = subprocess.Popen(['bedtools', 'sort', '-i', 'stdin'], stdin=bedtools_gc.stdout, stdout=out_fh)
			bedtools_sort.wait()
		bgbw = subprocess.Popen([bg_to_bw, out_bed, chrom_sizes, out_bw])
		bgbw.wait()

def convert_bams_bigwig_stranded(work_dir, bams, single_paired_end):
	# Forward strand > bamCoverage -b a.bam -o a.fwd.bw --samFlagExclude 16
	# Reverse strand > bamCoverage -b a.bam -o a.rev.bw --samFlagInclude 16
	os.chdir(work_dir)
	for bam in bams:
		samtools_index = subprocess.Popen(['samtools', 'index', bam])
		samtools_index.wait()
		forward_bw = bam.split('.')[0] + '.fwd.bigwig'
		reverse_bw = bam.split('.')[0] + '.rev.bigwig'
		if single_paired_end == 'single':
			bc_f = subprocess.Popen(['bamCoverage', '-b', bam, '-o', forward_bw, '--samFlagExclude', '16'])
			bc_f.wait()
			bc_r = subprocess.Popen(['bamCoverage', '-b', bam, '-o', reverse_bw, '--samFlagInclude', '16'])
			bc_r.wait()
		else:
			print('ahhh paired end not set up yet')


def bamcoverage_on_region(work_dir, bams, region, out_suffix):
	os.chdir(work_dir)
	for bam in bams:
		samtools_index = subprocess.Popen(['samtools', 'index', bam])
		samtools_index.wait()
		out_bigwig = bam.split('.')[0] + out_suffix
		bc_f = subprocess.Popen(['bamCoverage', '-b', bam, '-o', out_bigwig, '--region', region])
		bc_f.wait()		


##run methods
working_dir = '/home/atimms/ngs_data/rnaseq/kevin_rnaseq_0619/kevin_organoid_tc_0619'
bams = ['I4_1.Aligned.sortedByCoord.out.bam', 'I4_2.Aligned.sortedByCoord.out.bam', 'I4_3.Aligned.sortedByCoord.out.bam', 'I5_1.Aligned.sortedByCoord.out.bam', 
	'I5_2.Aligned.sortedByCoord.out.bam', 'I5_3.Aligned.sortedByCoord.out.bam', 'RC20_1.Aligned.sortedByCoord.out.bam', 'RC20_2.Aligned.sortedByCoord.out.bam', 
	'RC20_3.Aligned.sortedByCoord.out.bam', 'RC28_1.Aligned.sortedByCoord.out.bam', 'RC28_2.Aligned.sortedByCoord.out.bam', 'RC28_3.Aligned.sortedByCoord.out.bam', 
	'RC5_1.Aligned.sortedByCoord.out.bam', 'RC5_2.Aligned.sortedByCoord.out.bam', 'RC5_3.Aligned.sortedByCoord.out.bam']
# bams = ['I4_1.Aligned.sortedByCoord.out.bam']

##one file i.e. both strands
# convert_bams_bigwig(working_dir, bams, hg38_chrom_sizes, hg38_fa)

##two files i.e. forward and reverse
# convert_bams_bigwig_stranded(working_dir, bams, 'single')

working_dir = '/home/atimms/ngs_data/rnaseq/kim_rett_rnaseq_0119'
bam_files = ['Kim12.sorted.bam', 'Kim14.sorted.bam', 'Kim15.sorted.bam', 'Kim16.sorted.bam', 'Kim18.sorted.bam', 'Kim22.sorted.bam', 'Kim23.sorted.bam', 'Kim25.sorted.bam', 'Kim27.sorted.bam', 'Kim28.sorted.bam', 'Kim30.sorted.bam', 'Kim3.sorted.bam', 'Kim4.sorted.bam', 'Kim5.sorted.bam', 'Kim7.sorted.bam', 'Kim9.sorted.bam']
mecp2_region = 'chrX:153282179:153376247'
bigwig_suffix = '.MECP2.bigwig'
bamcoverage_on_region(working_dir, bam_files, mecp2_region, bigwig_suffix)

