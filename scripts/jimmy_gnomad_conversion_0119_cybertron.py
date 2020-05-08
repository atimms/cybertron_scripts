#!/usr/bin/env python
import subprocess
import re
import os
import numpy as np
##note
'''
module load local_python/3.6.4
module load biobuilds/2017.11
'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/references/gnomad_data'
os.chdir(working_dir)


histogram_bins = [0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975]




def convert_vcf_annovar_with_pab_max(in_vcf, out_file):
	##bcftools query -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%pab_max\n' gnomad.exomes.r2.1.sites.chr22.vcf.gz
	bt_query = subprocess.Popen(['bcftools','query', '-f', '%CHROM\t%POS\t%END\t%REF\t%ALT\t%pab_max\t%ab_hist_alt_bin_freq\n', '-o', out_file, in_vcf])
	bt_query.wait()


def get_median_from_hist(counts):
	values = []
	counts = counts.split('|')
	counts = [int(i) for i in counts]
	# print(counts)
	count_count = 0
	for count in counts:
		if count != 0:
			# print(count, count_count, histogram_bins[count_count])
			to_add = [histogram_bins[count_count]] * count
			values.extend(to_add)
			# print(to_add, values)
		count_count += 1
	# print(values)
	if len(values) >=1:
		##get median and mean
		h_mean = np.mean(values)
		h_median = np.median(values)
		h_min = min(values)
		h_max = max(values)
		h_count = len(values)
		# print(h_median, h_min, h_max, h_count)
		return([str(h_median), str(h_min), str(h_max), str(h_count)])
	else:
		return(['na','na','na','na' ])




def query_pab_max(avinput, infile, outfile):
	pab_dict = {}
	with open(avinput, "r") as av_ph:
		for line in av_ph:
			line = line.strip('\n').split(delim)
			var = '_'.join(line[:5])
			# print(var)
			pab_dict[var] = line
	with open(outfile, "w") as outph, open(infile, "r") as inph:
		for line in inph:
			line = line.strip('\n').split(delim)
			var2 = '_'.join(line)
			# print(var2)
			if var2 in pab_dict:
				print(pab_dict[var2])
				ab_median_mean = get_median_from_hist(pab_dict[var2][6])
				outph.write(delim.join(pab_dict[var2] + ab_median_mean) + '\n')
			else:
				print('var not found:', var2)

def add_info_to_avinput_get_subset(infile, outfile):
	with open(infile, "r") as in_ph, open(outfile, "w") as outph:
		for line in in_ph:
			line = line.strip('\n').split(delim)
			hist = line[6]
			ab_median_mean = get_median_from_hist(hist)
			c_median = ab_median_mean[0]
			c_count = ab_median_mean[3]
			if c_median != 'na':
				if float(c_median) < 0.25 and int(c_count) > 1000:
					outph.write(delim.join(line + ab_median_mean) + '\n')

##run methods
# gnomad_vcf = 'gnomad.exomes.r2.1.sites.chr22.vcf.gz'
gnomad_vcf = 'gnomad.genomes.r2.1.sites.vcf.gz'
gnomad_avinput = 'gnomad.genomes.r2.1.pab_max.avinput'
gnomad_info = 'gnomad.genomes.r2.1.aaf_candidates.txt'
# gnomad_avinput = 'test.avinput'
test_txt = 'dana_test_0119.txt'
pab_txt = 'dana_test_0119_pab.txt'
##convert the vcf to tab delimted format with pab_max value
# convert_vcf_annovar_with_pab_max(gnomad_vcf, gnomad_avinput)


##get pab_max values for set of queries
# query_pab_max(gnomad_avinput, test_txt, pab_txt)

##just get all vars
add_info_to_avinput_get_subset(gnomad_avinput, gnomad_info)



