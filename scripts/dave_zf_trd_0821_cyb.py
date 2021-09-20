#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import filtering_annotated
import homozygosity_mapping_cybertron
import math

##parameters
delim = '\t'
thread_number = '20'

##modules 
'''
qsub -Iq cdbrmq -l mem=30gb,ncpus=2 -P 2cf4ba3b-f9ef-4cfc-a2c1-be85f5db6730
module load biobuilds/2017.11
module load gcc/6.1.0 ##for homozygosity_mapping_cybertron

'''


##files
fasta = '/home/atimms/ngs_data/references/danRer11/danRer11.fa'
fasta_fai = '/home/atimms/ngs_data/references/danRer11/danRer11.fa.fai'
##first 2 columns of fasta.fai file
bedtools_genome_file  = '/home/atimms/ngs_data/references/danRer11/danRer11.fa.genome'

##programs
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'

##annovar parameters
av_genome = 'danRer11'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,ensGene,rmsk,generic,generic,generic,generic', '-genericdbfile', 'SH_1_A5.avinput,SH_2_B3.avinput,SH_3_m.avinput,SH_4_f.avinput']
av_operation = ['-operation', 'g,g,r,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,', '-vcfinput']
#av_options = ['-otherinfo', '-remove', '-vcfinput']


##methods

def filter_parent_snps(in_file, female_snps, male_snps):
	female_snp_outfile = female_snps + '.q50_cov10.txt'
	male_snp_outfile = male_snps + '.q50_cov10.txt'
	female_snp_noalt_outfile = female_snps + '.q50_cov10.no_alt.txt'
	male_snp_noalt_outfile = male_snps + '.q50_cov10.no_alt.txt'
	with open(in_file, 'r') as in_fh, open(female_snp_outfile, 'w') as fout_fh, open(male_snp_outfile, 'w') as mout_fh, open(female_snp_noalt_outfile, 'w') as f2out_fh, open(male_snp_noalt_outfile, 'w') as m2out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc == 1:
				mout_fh.write(delim.join(line) + '\n')
				fout_fh.write(delim.join(line) + '\n')
				m2out_fh.write(delim.join(line) + '\n')
				f2out_fh.write(delim.join(line) + '\n')
			else:
				qual = float(line[5])
				rmsk = line[21]
				male_parent_vcf_info = line[7]
				female_parent_vcf_info = line[8]
				male_pool_vcf_info = line[9]
				female_pool_vcf_info = line[10]
				male_parent_genotype = line[22]
				female_parent_genotype = line[23]
				male_pool_genotype = line[24]
				female_pool_genotype = line[25]
				vcf_infos = line[7:11]
				# coverages = [int(i.split(':')[2]) for i in vcf_infos]
				coverage_info = [i.split(':')[4].split(',') for i in vcf_infos]
				coverages =[]
				for c in coverage_info:
					cov = [int(i) for i in c]
					sum_cov = sum(cov)
					coverages.append(sum_cov)
				print(coverage_info, coverages)
				# print(vcf_infos, coverages)
				##filter vars for coverage, qual and repeats
				if min(coverages) >= 10 and max(coverages) <= 100 and qual >= 50 and rmsk == '.':
					##get male parent is heterozygous (with high confidence) and female is homozygous for ref
					if male_parent_genotype == 'het' and female_parent_genotype == '.':
						mout_fh.write(delim.join(line) + '\n')
						f_gt_coverages = female_parent_vcf_info.split(':')[4].split(',')
						# print(f_gt_coverages)
						if len(f_gt_coverages) == 2 and f_gt_coverages[1] == '0':
							m2out_fh.write(delim.join(line) + '\n')
						else:
							if len(f_gt_coverages) != 2:
								print('multiallele:', line)
					if female_parent_genotype == 'het' and male_parent_genotype == '.':
						fout_fh.write(delim.join(line) + '\n')
						m_gt_coverages = male_parent_vcf_info.split(':')[4].split(',')
						if len(m_gt_coverages) == 2 and m_gt_coverages[1] == '0':
							f2out_fh.write(delim.join(line) + '\n')
						else:
							if len(m_gt_coverages) != 2:
								print('multiallele:', line)


def make_bed_from_filtered_file(infile):
	outbed = infile.rsplit('.', 1)[0] + '.bed'
	with open(outbed, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				genotypes = line[9:11]
				allele_info = [i.split(':')[4].split(',') for i in genotypes]
				nafs = []
				for ai in allele_info:
					alelles = [int(i) for i in ai]
					naf = str(float(alelles[1]) / sum(alelles))
					nafs.append(naf)
				print(line,allele_info, nafs)

				line_out = [line[0], str(int(line[1]) - 1), line[2]] + nafs
				out_fh.write(delim.join(line_out) + '\n')

def make_files_for_graphing(snp_beds, window_size, step_size, genome_fai, working_dir):
	for snp_bed in snp_beds:
		for ws in window_size:
			#make bed file with windows and returns genome name and window size variable
			genome_and_window = homozygosity_mapping_cybertron.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]
			print(genome_and_window)
			window_bed = genome_and_window + '.bed'
			##all shared snps
			out_file_naf = snp_bed.rsplit('.', 1)[0] + '.' + genome_and_window + '.aaf_for_r.txt'
			out_file_combined = snp_bed.rsplit('.', 1)[0] + '.' + genome_and_window + '.combined_for_r.txt'
			##bedtools intersect 
			with open('temp.bed', "w") as naf_fh: 
				# hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', window_bed, '-b', snp_bed, '-c'], stdout=naf_fh)
				# hom_bt_intersect.wait()
				hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', window_bed, '-b', snp_bed, '-wa', '-wb'], stdout=naf_fh)
				hom_bt_intersect.wait()
			with open('temp.bed', "r") as in_fh:
				aaf_dict = {}
				for line in in_fh:
					line = line.rstrip().split(delim)
					chr_start_end = delim.join(line[:3])
					maaf = float(line[6])
					faaf = float(line[7])
					if chr_start_end in aaf_dict:
						aaf_dict[chr_start_end][0].append(maaf)
						aaf_dict[chr_start_end][1].append(faaf)
					else:
						aaf_dict[chr_start_end] = [[maaf],[faaf]]
			##write outfile from dict
			with open(out_file_naf, "w") as outn_fh, open(out_file_combined, "w") as outc_fh:
				outn_fh.write(delim.join(['chr', 'start', 'end', 'test', 'average_aaf']) + '\n')
				outc_fh.write(delim.join(['chr', 'start', 'end', 'test', 'value']) + '\n')

				for window in aaf_dict:
					ave_maaf = sum(aaf_dict[window][0]) / len(aaf_dict[window][0])
					ave_faaf = sum(aaf_dict[window][1]) / len(aaf_dict[window][1])
					aaf_diff = ave_maaf - ave_faaf
					snp_count = len(aaf_dict[window][0])
					log_snp_count = math.log(snp_count, 2)
					outn_fh.write(window + delim + 'male_aaf' + delim + str(ave_maaf) + '\n')
					outn_fh.write(window + delim + 'female_aaf' + delim + str(ave_faaf) + '\n')
					outc_fh.write(window + delim + 'aaf_male' + delim + str(ave_maaf) + '\n')
					outc_fh.write(window + delim + 'aaf_female' + delim + str(ave_faaf) + '\n')
					outc_fh.write(window + delim + 'aaf_difference' + delim + str(aaf_diff) + '\n')
					# outc_fh.write(window + delim + 'snp_count' + delim + str(snp_count) + '\n')
					outc_fh.write(window + delim + 'log2_snp_count' + delim + str(log_snp_count) + '\n')



def filter_parent_snps(in_file, female_snps, male_snps):
	female_snp_noalt_outfile = female_snps + '.q50_cov10.no_alt.txt'
	male_snp_noalt_outfile = male_snps + '.q50_cov10.no_alt.txt'
	with open(in_file, 'r') as in_fh, open(female_snp_noalt_outfile, 'w') as f2out_fh, open(male_snp_noalt_outfile, 'w') as m2out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc == 1:
				m2out_fh.write(delim.join(line) + '\n')
				f2out_fh.write(delim.join(line) + '\n')
			else:
				qual = float(line[5])
				rmsk = line[21]
				male_parent_vcf_info = line[7]
				female_parent_vcf_info = line[8]
				male_pool_vcf_info = line[9]
				female_pool_vcf_info = line[10]
				male_parent_genotype = line[22]
				female_parent_genotype = line[23]
				male_pool_genotype = line[24]
				female_pool_genotype = line[25]
				vcf_infos = line[7:11]
				func_refgene = line[11]
				func_ensgene = line[16]
				# coverages = [int(i.split(':')[2]) for i in vcf_infos]
				coverage_info = [i.split(':')[4].split(',') for i in vcf_infos]
				coverages =[]
				for c in coverage_info:
					cov = [int(i) for i in c]
					sum_cov = sum(cov)
					coverages.append(sum_cov)
				# print(coverage_info, coverages)
				# print(vcf_infos, coverages)
				##filter from genic variants
				if func_refgene == 'exonic' or func_ensgene == 'exonic':
					##filter vars for coverage, qual and repeats
					if min(coverages) >= 10 and max(coverages) <= 100 and qual >= 50 and rmsk == '.':
						##get male parent is heterozygous (with high confidence) and female is homozygous for ref
						if male_parent_genotype == 'het' and female_parent_genotype == '.':
							f_gt_coverages = female_parent_vcf_info.split(':')[4].split(',')
							# print(f_gt_coverages)
							if len(f_gt_coverages) == 2 and f_gt_coverages[1] == '0':
								m2out_fh.write(delim.join(line) + '\n')
							else:
								if len(f_gt_coverages) != 2:
									print('multiallele:', line)
						if female_parent_genotype == 'het' and male_parent_genotype == '.':
							m_gt_coverages = male_parent_vcf_info.split(':')[4].split(',')
							if len(m_gt_coverages) == 2 and m_gt_coverages[1] == '0':
								f2out_fh.write(delim.join(line) + '\n')
							else:
								if len(m_gt_coverages) != 2:
									print('multiallele:', line)







##run methods

##redoinf analysis from 0920, just using genic SNPs

work_dir = '/home/atimms/ngs_data/enu_mapping/dave_zf_trd_0821'
os.chdir(work_dir)

##params
project_name = 'dave_zf_sex_0920'
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
bamlist = 'bams.list'
st_vcf_suffix = '.st.vcf.gz'
result_vcf = project_name + st_vcf_suffix
fq_dict = {'SH_1_A5': ['SH_1_A5_CSFP200005686-1a_HJMKJDSXY_L1_1.fq.gz', 'SH_1_A5_CSFP200005686-1a_HJMKJDSXY_L1_2.fq.gz'],
	'SH_2_B3':['SH_2_B3_CSFP200005687-1a_HJMJNDSXY_L1_1.fq.gz', 'SH_2_B3_CSFP200005687-1a_HJMJNDSXY_L1_2.fq.gz'],
	'SH_3_m':['SH_3_m_CSFP200005688-1a_HJMKJDSXY_L1_1.fq.gz', 'SH_3_m_CSFP200005688-1a_HJMKJDSXY_L1_2.fq.gz'],
	'SH_4_f':['SH_4_f_CSFP200005689-1a_HJLM3DSXY_L1_1.fq.gz', 'SH_4_f_CSFP200005689-1a_HJLM3DSXY_L1_2.fq.gz']}

##filter vars for graphing - for parent being het and pool being hom reference 
annotated_file = project_name + '.combined.annotated.txt'
mom_het_dad_hom_prefix = project_name + '.mom_het_dad_wt.genic'
dad_het_mom_hom_prefix = project_name + '.dad_het_mom_wt.genic'
##filter 2 ways for both test
filter_parent_snps(annotated_file, mom_het_dad_hom_prefix, dad_het_mom_hom_prefix)

##make naf bed from filtered files
filtered_files = ['dave_zf_sex_0920.dad_het_mom_wt.genic.q50_cov10.no_alt.txt', 'dave_zf_sex_0920.mom_het_dad_wt.genic.q50_cov10.no_alt.txt']
# filtered_files = ['temp.txt'] ##testing
##make file into bed for garphing
for filtered_file in filtered_files:
	make_bed_from_filtered_file(filtered_file)



##make files for graphing
filtered_beds = ['dave_zf_sex_0920.dad_het_mom_wt.genic.q50_cov10.no_alt.bed', 'dave_zf_sex_0920.mom_het_dad_wt.genic.q50_cov10.no_alt.bed']
window_sizes = [10000000]
step_size = 1000000
make_files_for_graphing(filtered_beds, window_sizes, step_size, fasta_fai, work_dir)





