#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import homozygosity_mapping_cybertron


'''
module load biobuilds/2017.11
module load gcc/6.1.0 ##for homozygosity_mapping_cybertron
'''

##parameters
delim = '\t'


def get_files_to_graph(ten_counts, window_size, step_size, fai):
	for ws in window_size:
		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping_cybertron.make_windows(working_dir, fai, ws, step_size).split('.')[0]
		print(genome_and_window)
		window_bed = genome_and_window + '.bed'
		##count tens per window
		out_bed = ten_counts.rsplit('.', 1)[0] + '.' + genome_and_window + '.bed'
		##bedtools intersect 
		with open('temp.bed', "w") as naf_fh: 
			hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', window_bed, '-b', ten_counts, '-c'], stdout=naf_fh)
			hom_bt_intersect.wait()
		##add header
		with open(out_bed, "w") as out_fh, open('temp.bed', "r") as in_fh:
			out_fh.write(delim.join(['chr', 'start', 'end', 'ten_count', '\n']))
			for line in in_fh:
				out_fh.write(line)












##working dir
# working_dir = '/data/atimms/timon_0317'
working_dir = '/home/atimms/ngs_data/misc/dave_muga_0420'
os.chdir(working_dir)

genome_fai = '/home/atimms/ngs_data/references/mm10/mm10.fa.fai'
het_count_bed = 'smo_het_count.bed'
ten_bed = 'smo_tens.bed'
fifteen_bed = 'smo_15s.bed'
window_sizes = [1000000,5000000,2000000,10000000]
step_sizes = 1000000

##make graphs
# get_files_to_graph(ten_bed, window_sizes, step_sizes, genome_fai)

get_files_to_graph(fifteen_bed, window_sizes, step_sizes, genome_fai)
##copy whole dir to /archive/beier_d/dave_data so can graph in R








