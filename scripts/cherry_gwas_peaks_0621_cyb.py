#!/usr/bin/env python
import subprocess
import os
from collections import Counter
import glob

'''
info:
qsub -Iq cdbrmq -l mem=20gb,ncpus=1 -P 19833a08-f6fb-4bea-8526-8a79069da878
'''

##parameters
delim = '\t'
##programs
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
##files etc
bt_genome = '/home/atimms/ngs_data/references/10x/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.chrom_sizes'


##methods
def make_combined_snp_bed(in_files, out_file):
	with open(out_file, "w") as out_fh:
		for in_file in in_files:
			with open(in_file, "r") as in_fh:
				study = in_file.split('.')[0]
				for line in in_fh:
					line = line.rstrip().split(delim)
					line_out = line[:3] + [study]
					out_fh.write(delim.join(line_out) + '\n')


def bt_intersect_snps_bed(snp_bed, peak_beds, out_bed):
	##bedtools intersect 
	with open(out_bed, "w") as out_fh: 
		# hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', region_bed, '-b'] + summit_peak_beds + ['-wa', '-wb', '-filenames'], stdout=out_fh)
		# hom_bt_intersect.wait()
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', snp_bed, '-b'] + peak_beds + ['-C', '-filenames'], stdout=out_fh)
		hom_bt_intersect.wait()

def compare_snps_vs_cc_peaks(in_file, out_file):
	snp_count_dict = {}
	peak_names = []
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			# print(line)
			line = line.rstrip().split(delim)
			snp = '_'.join(line[:3])
			study = line[3]
			peak = line[4].split('.')[0] + '.' + line[4].split('.')[1]
			if peak not in peak_names:
				peak_names.append(peak)
			count = int(line[5])
			
			if study in snp_count_dict:
				##if snp seen before
				if snp in snp_count_dict[study]:
					if count == 0:
						snp_count_dict[study][snp][1].append(peak)
					else:
						snp_count_dict[study][snp][0] += 1
						snp_count_dict[study][snp][2].append(peak)
				else:
					if count == 0:
						snp_count_dict[study][snp] = [0, [peak], []]
					else:
						snp_count_dict[study][snp] = [1, [], [peak]]	
			##if study not seen before 
			else:
				if count == 0:
					snp_count_dict[study] = {snp:[0, [peak], []]}
				else:
					snp_count_dict[study] = {snp:[1, [], [peak]]}
	# print(snp_count_dict)
	cc_study_count_dict = {}
	header = ['cell class']
	for st in snp_count_dict:
		header.append(st)
		print('total snps in', st, "=", len(snp_count_dict[st]))
		total_snps = len(snp_count_dict[st])
		all_study_peaks_with_snps = []
		for s in snp_count_dict[st]:
			# print(st, s, snp_count_dict[st][s][0], len(snp_count_dict[st][s][1]), len(snp_count_dict[st][s][2]))
			all_study_peaks_with_snps.extend(snp_count_dict[st][s][2])
		# print(st, len(Counter(all_study_peaks_with_snps)), Counter(all_study_peaks_with_snps))
		##get counts of all snp in cc peaks
		snp_peak_counter = Counter(all_study_peaks_with_snps)
		# print(snp_peak_counter)
		for peak in peak_names:
			# print(peak, snp_peak_counter[peak])
			intersect_count = snp_peak_counter[peak]
			intersect_fraction = (intersect_count/total_snps) * 100
			# print(peak, intersect_count, intersect_fraction, total_snps, st)
			if peak in cc_study_count_dict:
				cc_study_count_dict[peak].append(intersect_fraction)
			else:
				cc_study_count_dict[peak] = [intersect_fraction]
	# print(header)
	# print(cc_study_count_dict)
	with open(out_file, "w") as out_fh:
		out_fh.write(delim.join(header) + '\n')
		for s in cc_study_count_dict:
			print(s, cc_study_count_dict[s])
			line_out = [s] + [str(i) for i in cc_study_count_dict[s]]
			out_fh.write(delim.join(line_out) + '\n')


def make_region_bed(in_file, out_bed):
	region_dict = {}
	lc = 0
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			lc +=1
			if lc >1:
				# print(line)
				line = line.rstrip().split(delim)
				chrom = 'chr' + line[4]	
				start = int(line[5]) -1
				end = int(line[5])
				region = line[0]
				if region in region_dict:
					if chrom == region_dict[region][0]:
						region_dict[region][3] += 1
						if start < region_dict[region][1]:
							region_dict[region][1] = start
						if end > region_dict[region][2]:
							region_dict[region][2] = end
					else:
						print('issue with:', line)						

				else:
					region_dict[region] = [chrom, start, end, 1]
	with open(out_bed, "w") as out_fh:
		for r in region_dict:
			print(region_dict[r] + [r])
			size = region_dict[r][2] - region_dict[r][1]
			if size <2000:
				region_dict[r][1] -= 1000
				region_dict[r][2] += 1000
			line_out = [str(i) for i in region_dict[r] + [r]]
			out_fh.write(delim.join(line_out) + '\n')


def compare_regions_vs_cc_peaks(in_file, out_file_prefix):
	region_dict = {}
	peak_names = []
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			# print(line)
			line = line.rstrip().split(delim)
			region = line[4]
			peak = line[5].split('.')[0] + '.' + line[5].split('.')[1]
			snp_count = int(line[3])
			region_size = int(line[2]) - int(line[1])
			count = int(line[6])
			##get normalized counts
			if count == 0:
				count_snp, count_size = 0, 0
			else:
				count_snp = count/snp_count
				count_size = count/region_size
			##get list/order of peaks 
			if peak not in peak_names:
				peak_names.append(peak)
			##populate dict
			if region in region_dict:
				region_dict[region][0].append(count)
				region_dict[region][1].append(count_snp)
				region_dict[region][2].append(count_size)
			else:
				region_dict[region] = [[count], [count_snp], [count_size]]
	c_out_file = out_file_prefix + '.count.txt'
	c_snp_out_file = out_file_prefix + '.count_snp.txt'
	c_size_out_file = out_file_prefix + '.count_size.txt'
	with open(c_out_file, "w") as out1_fh, open(c_snp_out_file, "w") as out2_fh, open(c_size_out_file, "w") as out3_fh:
		out1_fh.write(delim.join(['region'] + peak_names) + '\n')
		out2_fh.write(delim.join(['region'] + peak_names) + '\n')
		out3_fh.write(delim.join(['region'] + peak_names) + '\n')
		for r in region_dict:
			counts = [str(x) for x in region_dict[r][0]]
			counts_snps = [str(x) for x in region_dict[r][1]]
			counts_size = [str(x) for x in region_dict[r][2]]
			out1_fh.write(delim.join([r] + counts) + '\n')
			out2_fh.write(delim.join([r] + counts_snps) + '\n')
			out3_fh.write(delim.join([r] + counts_size) + '\n')

def compare_lead_snp_vs_cc_peaks(in_file, out_file):
	region_dict = {}
	peak_names = []
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			# print(line)
			line = line.rstrip().split(delim)
			region = line[3]
			peak = line[4].split('.')[0] + '.' + line[4].split('.')[1]
			count = int(line[5])
			##get list/order of peaks 
			if peak not in peak_names:
				peak_names.append(peak)
			##populate dict
			if region in region_dict:
				region_dict[region].append(count)

			else:
				region_dict[region] = [count]
	with open(out_file, "w") as out_fh:
		out_fh.write(delim.join(['region'] + peak_names) + '\n')
		for r in region_dict:
			counts = [str(x) for x in region_dict[r]]
			out_fh.write(delim.join([r] + counts) + '\n')


def make_dbsnp_dict(dbsnp_bed):
	snp_dict = {}
	with open(dbsnp_bed, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			rs = line[3]
			chr_pos = line[0] + '_' + line[2]
			snp_dict[rs] = chr_pos
	return(snp_dict)


def get_lead_snps(infiles, amd2016_data, snp_locus, snp_set):
	amd2016_dict = {}
	##get 'real' pvalues
	# '''
	with open(amd2016_data, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				rs = line[0]
				p_value = line[7]
				amd2016_dict[rs] = p_value
				# print(rs, p_value)
	# '''
	##get SNP per locus/set
	set_dict, locus_dict = {}, {}
	for infile in infiles:
		locus = infile.rsplit('/', 1)[1].split('.')[0]
		with open(infile, "r") as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				if lc > 1:
					line = line.rstrip().replace('"','').split(',')
					locus_set = locus + '_' + line[0]
					# print(locus, locus_set)
					##these files had one extra column
					if locus.startswith('MacTel2020'):
						if len(line) == 10:
							chrom = 'chr' + line[4]
							pos = int(line[5])
							p_value = float(line[9])
						elif len(line) == 11:
							chrom = 'chr' + line[5]
							pos = int(line[6])
							p_value = float(line[10])							
					##need to get real p value
					elif locus.startswith('AMD2016'):
						rs = line[2]
						chrom = 'chr' + line[3]
						pos = int(line[4])
						p_value = amd2016_dict[rs]
					else:
						chrom = 'chr' + line[3]
						pos = int(line[4])
						p_value = float(line[8])
					##add snp to locus dict if the most significant
					if locus in locus_dict:
						if p_value == locus_dict[locus][2]:
							print(p_value, locus_set)
							locus_dict[locus][1].append(pos)
						elif p_value < locus_dict[locus][2]:
							locus_dict[locus][1] = [pos]
							locus_dict[locus][2] = p_value
					else:
						locus_dict[locus] = [chrom,[pos], p_value]
					##then the same for the sets
					if locus_set in set_dict:
						if p_value == set_dict[locus_set][2]:
							# print(p_value, locus_set)
							set_dict[locus_set][1].append(pos)
						elif p_value < set_dict[locus_set][2]:
							set_dict[locus_set][1] = [pos]
							set_dict[locus_set][2] = p_value
					else:
						set_dict[locus_set] = [chrom,[pos], p_value]
					# print(locus, chrom, pos, p_value, locus_dict[locus], set_dict[locus_set])
	##print files
	with open(snp_locus, "w") as out_fh:
		for l in locus_dict:
			chrom = locus_dict[l][0]
			##get mid positon between snps
			if len(locus_dict[l][1]) == 1:
				pos = locus_dict[l][1][0]
			else:
				pos = int(sum(locus_dict[l][1]) / len(locus_dict[l][1]))
			start = str(pos - 1001)
			end = str(pos + 1000)
			##use all lead snps and extend
			# if len(locus_dict[l][1]) == 1:
			# 	pos = locus_dict[l][1][0]
			# 	start = str(pos - 1001)
			# 	end = str(pos + 1000)
			# else:
			# 	start = str(min(locus_dict[l][1]) - 1001)
			# 	end = str(max(locus_dict[l][1]) + 1000)
			out_fh.write(delim.join([chrom, start, end, l]) + '\n')
	with open(snp_set, "w") as out_fh:
		for s in set_dict:
			chrom = set_dict[s][0]
			##get mid positon between snps
			if len(set_dict[s][1]) == 1:
				pos = set_dict[s][1][0]
			else:
				pos = int(sum(set_dict[s][1]) / len(set_dict[s][1]))
			start = str(pos - 1001)
			end = str(pos + 1000)
			##use all lead snps and extend
			# if len(set_dict[s][1]) == 1:
			# 	pos = set_dict[s][1][0]
			# 	start = str(pos - 1001)
			# 	end = str(pos + 1000)
			# else:
			# 	start = str(min(set_dict[s][1]) - 1001)
			# 	end = str(max(set_dict[s][1]) + 1000)
			out_fh.write(delim.join([chrom, start, end, s]) + '\n')




##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_gwas_atac_peaks_0621'
os.chdir(working_dir)

##params
##peaks copied from /active/cherry_t/OrgManuscript_SingleCell_Data/peak_calls/human_cell_class_macs2_0621
cell_class_narrowpeaks = ['human.AC_HC_GC_Precursors.macs2_q0.000001_peaks.narrowPeak', 'human.Developing_Amacrines.macs2_q0.000001_peaks.narrowPeak', 
		'human.Developing_Bipolars.macs2_q0.000001_peaks.narrowPeak', 'human.Developing_Cones.macs2_q0.000001_peaks.narrowPeak', 
		'human.Developing_Ganglions.macs2_q0.000001_peaks.narrowPeak', 'human.Developing_Horizontals.macs2_q0.000001_peaks.narrowPeak', 
		'human.Developing_Rods.macs2_q0.000001_peaks.narrowPeak', 'human.Early_Progenitors.macs2_q0.000001_peaks.narrowPeak', 
		'human.Ganglion_Precursors.macs2_q0.000001_peaks.narrowPeak', 'human.Late_Progenitors.macs2_q0.000001_peaks.narrowPeak', 
		'human.Mature_Amacrines.macs2_q0.000001_peaks.narrowPeak', 'human.Mature_Bipolars.macs2_q0.000001_peaks.narrowPeak', 
		'human.Mature_Cones.macs2_q0.000001_peaks.narrowPeak', 'human.Mature_Horizontals.macs2_q0.000001_peaks.narrowPeak', 
		'human.Mature_Mullers.macs2_q0.000001_peaks.narrowPeak', 'human.Mature_Rods.macs2_q0.000001_peaks.narrowPeak', 
		'human.Photoreceptor_Bipolar_Precursors.macs2_q0.000001_peaks.narrowPeak']
##lifted over narrow peaks to hg19
cell_class_hg19_beds = ['human.AC_HC_GC_Precursors.macs2_q0.000001_peaks.hg19.bed', 'human.Developing_Amacrines.macs2_q0.000001_peaks.hg19.bed', 
		'human.Developing_Bipolars.macs2_q0.000001_peaks.hg19.bed', 'human.Developing_Cones.macs2_q0.000001_peaks.hg19.bed', 
		'human.Developing_Ganglions.macs2_q0.000001_peaks.hg19.bed', 'human.Developing_Horizontals.macs2_q0.000001_peaks.hg19.bed', 
		'human.Developing_Rods.macs2_q0.000001_peaks.hg19.bed', 'human.Early_Progenitors.macs2_q0.000001_peaks.hg19.bed', 
		'human.Ganglion_Precursors.macs2_q0.000001_peaks.hg19.bed', 'human.Late_Progenitors.macs2_q0.000001_peaks.hg19.bed', 
		'human.Mature_Amacrines.macs2_q0.000001_peaks.hg19.bed', 'human.Mature_Bipolars.macs2_q0.000001_peaks.hg19.bed', 
		'human.Mature_Cones.macs2_q0.000001_peaks.hg19.bed', 'human.Mature_Horizontals.macs2_q0.000001_peaks.hg19.bed', 
		'human.Mature_Mullers.macs2_q0.000001_peaks.hg19.bed', 'human.Mature_Rods.macs2_q0.000001_peaks.hg19.bed', 
		'human.Photoreceptor_Bipolar_Precursors.macs2_q0.000001_peaks.hg19.bed']

hg38_dbsnp_bed = 'hg38_dbsnp.bed'
gwas_beds = ['amd2013.snps.chr.bed', 'amd2016.snps.chr.bed', 'early_amd2020.snps.chr.bed', 'mactel.snps.chr.bed', 
		'progression_amd2018.chr.bed', 'progression_baseline_amd2018.chr.bed', 'progression_cnv2018.chr.bed', 
		'retinitis_pigmentosa.chr.bed']
combined_gwas_snp_bed  = 'combined_gwas_snps.bed'
snps_peaks_bt_bed = 'combined_snps.cc_peaks.bt_int.bed'
snps_peaks_sum = 'combined_snps.cc_peaks.counts_for_heatmap.txt'


##cell class peaks per SNP

## 1. combine all snp files
# make_combined_snp_bed(gwas_beds, combined_gwas_snp_bed)

## 2. intersect/compare summit peaks with gwas SNPs
bt_intersect_snps_bed(combined_gwas_snp_bed, cell_class_hg19_beds, snps_peaks_bt_bed)

## 3. compare studies
compare_snps_vs_cc_peaks(snps_peaks_bt_bed, snps_peaks_sum)


##cell class peaks per regio
combined_gwas_snp_region_bed = 'combined_gwas_snp_regions.txt'
gwas_regions_bed = 'combined_gwas_regions.bed'
regions_peaks_bt_bed = 'combined_gwas_regions.cc_peaks.bt_int.bed'
regions_peaks_heatmap_prefix = 'combined_gwas_regions.cc_peaks.for_heatmap'

## 1. make bed per region
# make_region_bed(combined_gwas_snp_region_bed, gwas_regions_bed)

## 2. intersect/compare summit peaks with gwas region
bt_intersect_snps_bed(gwas_regions_bed, cell_class_hg19_beds, regions_peaks_bt_bed)

## 3. compare studies
compare_regions_vs_cc_peaks(regions_peaks_bt_bed, regions_peaks_heatmap_prefix)


##cell class peaks per lead SNP i.e. one per region
gwas_peaks = glob.glob(working_dir + '/gwas_peaks/*csv')
# gwas_peaks = [working_dir + '/gwas_peaks/MacTel2020_locus3.csv']
extra_amd2016_data = working_dir + '/gwas_peaks/Fritsche-26691988.txt'
gwas_lead_snp_region = 'gwas_lead_snp_per_region_0721.bed'
gwas_lead_snp_set = 'gwas_lead_snp_per_set_0721.bed'
gwas_lead_snp_region_bt_bed = 'gwas_lead_snp_per_region_0721.cc_peaks.bt_int.bed'
gwas_lead_snp_region_sum = 'gwas_lead_snp_per_region_0721.cc_peaks.counts_for_heatmap.txt'
gwas_lead_snp_set_bt_bed = 'gwas_lead_snp_per_set_0721.cc_peaks.bt_int.bed'
gwas_lead_snp_set_sum = 'gwas_lead_snp_per_set_0721.cc_peaks.counts_for_heatmap.txt'


##not needed
# dbsnp_dict = make_dbsnp_dict(hg38_dbsnp_bed)


## 1. get lead SNP per region and by snp
##need to get pvalue for amd2016 and get aggregate position for amd2013
# get_lead_snps(gwas_peaks, extra_amd2016_data, gwas_lead_snp_region, gwas_lead_snp_set)

## 2. intersect/compare summit peaks with gwas SNPs
# bt_intersect_snps_bed(gwas_lead_snp_region, cell_class_hg19_beds, gwas_lead_snp_region_bt_bed)
# bt_intersect_snps_bed(gwas_lead_snp_set, cell_class_hg19_beds, gwas_lead_snp_set_bt_bed)

## 3. compare studies
# compare_lead_snp_vs_cc_peaks(gwas_lead_snp_region_bt_bed, gwas_lead_snp_region_sum)
# compare_lead_snp_vs_cc_peaks(gwas_lead_snp_set_bt_bed, gwas_lead_snp_set_sum)

