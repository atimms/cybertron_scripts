#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'

##methods

def bedtools_sort_merge(in_bed, out_bed):
	with open(out_bed, "w") as out_fh:
		sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n', in_bed], stdout=subprocess.PIPE)
		bt_merge = subprocess.Popen(['bedtools', 'merge', '-i', '-'], stdin=sort_file.stdout, stdout=out_fh)
		bt_merge.wait()


def bedtools_multicov(bams, bed, out_file):
	with open(out_file, "w") as out_fh:
		bt_mc = subprocess.Popen(['bedtools', 'multicov', '-bams'] + bams + ['-bed', bed], stdout=out_fh)
		bt_mc.wait() 





##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/cunn_exome_1120'
os.chdir(working_dir)
##info files etx
bed_file = 'hg19.refseq_exons.chr18_15-130mb.bed'
sorted_bed = 'hg19.refseq_exons.chr18_15-130mb.sorted.bed'
coverage_results = 'AK-01-111.bt_cov.chr18_15-130mb.txt'
bam_files = ['AK-01-111-01a.bwa_gatk.bam', 'AK-01-111-01u.bwa_gatk.bam', 'AK-01-111-02u.bwa_gatk.bam', 'AK-01-111-03u.bwa_gatk.bam']

##1. sort and merge bed file
bedtools_sort_merge(bed_file, sorted_bed)

##2.multicov on all bams
bedtools_multicov(bam_files, sorted_bed, coverage_results)
