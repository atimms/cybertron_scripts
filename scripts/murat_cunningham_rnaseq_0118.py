#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil
import pybedtools as pbt

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/murat_cunningham_0118'
os.chdir(working_dir)

##file names etc
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'

def run_bedtools_genomecov(bam_files):
	#bedtools genomecov -ibam $bam'_'$map_prog'_proccessed.bam' -bg  > 'temp.bed'
	for bam_file in bam_files:
		sample = bam_file.split('.')[0]
		outfile = sample + '.bt_gc.temp.bed'
		with open(outfile, "w") as outfh:
			bt_gc = subprocess.Popen(['bedtools', 'genomecov', '-ibam', bam_file, '-bg', '-split'], stdout=outfh)
			bt_gc.wait()
def filter_by_coverage(bam_files, covs_req):
	for bam_file in bam_files:
		sample = bam_file.split('.')[0]
		bed_file = sample + '.bt_gc.temp.bed'
		for cov_req in covs_req:
			outfile = sample + '.' + str(cov_req) + 'cov.bed'
			with open(outfile, "w") as outfh, open(bed_file, "r") as infh:
				for line in infh:
					line = line.split(delim)
					coverage = int(line[3])
					chrom = line[0]
					if coverage >= cov_req and chrom[0] != 'G':
						outfh.write(delim.join(line))
def calculate_coverage_in_all(bam_files, covs_req):
	for cov_req in covs_req:
		bed_files = []
		union_bed = 'unionbedg.' + str(cov_req) + 'cov.bed'
		ub_sorted = 'unionbedg.' + str(cov_req) + 'cov.merged.bed'
		with open(union_bed, "w") as outfh:
			for bam_file in bam_files:
				sample = bam_file.split('.')[0]
				bed_file = sample + '.' + str(cov_req) + 'cov.bed'
				bed_files.append(bed_file)
			bt_ub = subprocess.Popen(['bedtools', 'unionbedg', '-i'] + bed_files + ['-filler', 'none'], stdout=subprocess.PIPE)
			grep_cmd = subprocess.Popen(['grep', '-v', 'none'], stdin=bt_ub.stdout, stdout=outfh)
			grep_cmd.wait()
		c = pbt.BedTool(union_bed)
		d = c.sort().merge().moveto(ub_sorted)

def add_coverage(covs_req):
	for cov_req in covs_req:
		total_cov = 0
		ub_sorted = 'unionbedg.' + str(cov_req) + 'cov.merged.bed'
		union_bed = 'unionbedg.' + str(cov_req) + 'cov.bed'
		# with open(ub_sorted, "r") as infh:
		with open(union_bed, "r") as infh:
			for line in infh:
				line = line.rstrip().split(delim)
				start = int(line[1])
				end = int(line[2])
				cov = end -start
				total_cov += cov
		print "for coverage of %s or greater total of %s bps are covered in all samples"%(cov_req, total_cov)

	##run ub_sorted
bams = glob.glob('*.accepted_hits.merged.nodups.realigned.recal.bam')
# bams = ['95704.accepted_hits.merged.nodups.realigned.recal.bam', '95703.accepted_hits.merged.nodups.realigned.recal.bam']
# bams = ['95445.accepted_hits.merged.nodups.realigned.recal.bam', '95446.accepted_hits.merged.nodups.realigned.recal.bam', '95447.accepted_hits.merged.nodups.realigned.recal.bam', '95448.accepted_hits.merged.nodups.realigned.recal.bam', '95449.accepted_hits.merged.nodups.realigned.recal.bam', '95450.accepted_hits.merged.nodups.realigned.recal.bam', '95451.accepted_hits.merged.nodups.realigned.recal.bam', '95452.accepted_hits.merged.nodups.realigned.recal.bam', '95453.accepted_hits.merged.nodups.realigned.recal.bam', '95454.accepted_hits.merged.nodups.realigned.recal.bam']
# coverage_wanted = [5,10,20,50]
coverage_wanted = [10,20,50]
# run_bedtools_genomecov(bams)
# filter_by_coverage(bams, coverage_wanted)
# calculate_coverage_in_all(bams, coverage_wanted)
add_coverage(coverage_wanted)


