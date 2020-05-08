#!/usr/bin/env python
import sys
import subprocess
import os




##methods
def convert_bam_fastq_bedtools(bamfile, r1_fastq, r2_fastq, work_dir):
	os.chdir(work_dir)
	read_sorted_bam = bamfile[:-4] + 'n_sorted.bam'
	st_n_sort = subprocess.Popen(['samtools', 'sort', '-nO', 'bam', '-o', read_sorted_bam, '-@', '10', '-m', '10G', '-T', 'tempy', bamfile])
	st_n_sort.wait()
	bam_fq = subprocess.Popen(['bedtools', 'bamtofastq', '-i', read_sorted_bam, '-fq', r1_fastq, '-fq2', r2_fastq])
	bam_fq.wait()




##run methods
##working dir
working_dir = '/gpfs/research/rit/atimms/genomes/LR16-302'

print 'poopy eyeballs'
convert_bam_fastq_bedtools('UDN222967-HNVH7CCXX_s8_SL154673-recal.bam', 'UDN222967_s8_r1.fq', 'UDN222967_s8_r2.fq', working_dir)



