#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
# import dobyns_gemini_pipeline_cybertron_v9
import statistics

'''
load these:
module load java/1.8.0_121 
module load biobuilds/2016.11
module load vt/0.5772-60f436c3 
module load mono/5.10.1.47
module load Pisces/5.1.6.54
module load local_python/3.6.4 
'''

##set input variables and parameters
delim = '\t'
##working directory
working_dir = '/home/atimms/ngs_data/misc/kim_sumarize_exomes_0818'
os.chdir(working_dir)

samtools = '/cm/shared/apps/biobuilds/biobuilds-2016.11/bin/samtools'

def get_coverage_data(peds, infile_suffix, outfile):
	with open(outfile, "w") as out_fh:
		ped_count = 0
		for ped in peds:
			ped_count += 1	
			infile = ped + infile_suffix
			with open(infile, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line = line.split(delim)
					line_count += 1
					if line_count == 1:
						if ped_count == 1:
							out_fh.write(delim.join(line))
					else:
						##get all non parental and unaffected sib info
						ind_id = line[0].split('.')[0]
						if ind_id[0] == 'L':
							# print(len(ind_id), ind_id)
							if len(ind_id) <= 8:
								print(ped, ind_id)
								out_fh.write(delim.join(line))
							else:
								if ind_id[8] != 'f' and ind_id[8] != 'm' and ind_id[8] != 's':
									print(ped, ind_id, ind_id[8])
									out_fh.write(delim.join(line))
						else:
							if ind_id[-1] != 'F' and ind_id[-1] != 'M':
								print(ped, ind_id)
								out_fh.write(delim.join(line))



def make_dict_from_expression_files(brainspan_file, gtex_file):
	expression_dict = {}
	with open(brainspan_file, "U") as bs_fh:
		line_count = 0
		for line in bs_fh:
			line = line.rstrip().split(',')
			line_count += 1
			if line_count > 1:
				gene = line[3].split('.')[0]
				expressed = int(line[19])
				if gene in expression_dict:
					expression_dict[gene][0].append(expressed)
				else:
					expression_dict[gene] = [[expressed], []]
	with open(gtex_file, "U") as gt_fh:
		line_count = 0
		for line in gt_fh:
			line = line.rstrip().split(',')
			line_count += 1
			if line_count > 1:
				gene = line[1].split('.')[0]
				expressed = int(line[6])
				if gene in expression_dict:
					expression_dict[gene][1].append(expressed)
				else:
					expression_dict[gene] = [[], [expressed]]
	return(expression_dict)
	# for g in expression_dict:
	# 	print(g, expression_dict[g])



def get_peds_add_expression_data(peds, expression_dict, infile, outfile):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				header = line + ['expressed in brainspan', 'expressed in gtex', '\n']
				out_fh.write(delim.join(header))
			else:
				gene = line[5]
				ped_id = line[20]
				if ped_id in peds:
					print(gene, ped_id)
					if gene in expression_dict:
						if len(expression_dict[gene][0]) == 0:
							bs_exp = 'na'
						else:
							if sum(expression_dict[gene][0]) == 0:
								bs_exp = 'no'
							else:
								bs_exp = 'yes'
						if len(expression_dict[gene][1]) == 0:
							gt_exp = 'na'
						else:
							if sum(expression_dict[gene][1]) == 0:
								gt_exp = 'no'
							else:
								gt_exp = 'yes'			
					line_out = line + [bs_exp, gt_exp, '\n']
					out_fh.write(delim.join(line_out))


def get_counts_from_var_file(infile, outfile):
	count_dict = {}
	with open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count > 1:
				ped_id = line[20]
				analysis = line[24]
				chrom = line[0]
				max_aaf = float(line[12])
				cadd = line[11]
				if cadd == 'None':
					cadd = 0
				cadd = float(cadd)
				bs_exp = line[25]
				gt_exp = line[26]
				##set up dict for each pedigree
				if ped_id not in count_dict:
					count_dict[ped_id] = [0,0,0,0,0,0,0,0,0,0]
				##add to unfiltered numbers
				if analysis == 'de_novo':
					if chrom != 'chrX':
						count_dict[ped_id][0] += 1
					else:
						count_dict[ped_id][4] += 1
				elif analysis == 'x_linked_dn' or analysis == 'x_linked_de_novo':
					count_dict[ped_id][4] += 1	
				elif analysis == 'autosomal_recessive' or analysis == 'auto_rec':
					count_dict[ped_id][1] += 1					
				elif analysis == 'comp_hets' or analysis == 'compound_het':
					count_dict[ped_id][2] += 1	
				elif analysis == 'x_linked_recessive' or analysis == 'x_linked':
					count_dict[ped_id][3] += 1
				##filter the vars
				# if max_aaf <= 0.001 and cadd > 10 and (bs_exp == 'yes' or gt_exp == 'yes'):
				##not using exp data
				if max_aaf <= 0.001 and cadd > 10:
					# print(max_aaf, cadd, bs_exp, gt_exp, ped_id)
					if analysis == 'de_novo':
						if chrom != 'chrX':
							count_dict[ped_id][5] += 1
						else:
							count_dict[ped_id][9] += 1
					elif analysis == 'x_linked_dn' or analysis == 'x_linked_de_novo':
						count_dict[ped_id][9] += 1	
					elif analysis == 'autosomal_recessive' or analysis == 'auto_rec':
						count_dict[ped_id][6] += 1					
					elif analysis == 'comp_hets' or analysis == 'compound_het':
						count_dict[ped_id][7] += 1	
					elif analysis == 'x_linked_recessive' or analysis == 'x_linked':
						count_dict[ped_id][8] += 1

	# for p in count_dict:
	# 	print(p, count_dict[p])
	with open(outfile, "w") as out_fh:
		header = ['pedigree', 'autosomal_de_novo', 'autosomal_homozygous', 'compound_het', 'x_linked', 'x_linked_dn', 
				'filtered_autosomal_de_novo', 'filtered_autosomal_homozygous', 'filtered_compound_het', 
				'filtered_x_linked', 'filtered_x_linked_dn', '\n']
		out_fh.write(delim.join(header))
		for p in count_dict:
			count_data = [str(i) for i in count_dict[p]]
			line_out = [p] + count_data + ['\n']
			out_fh.write(delim.join(line_out))

def caculate_median_read_length(bam):
	##calc
	temp_sam = bam.split('.')[0] + '.read_length_temp.sam'
	with open(temp_sam, 'w') as temp_fh:
		# run_samtools = subprocess.Popen(['samtools view', bam, "| awk '{print length($10)}' | head -10000 | sort | uniq -c '])"])
		run_samtools = subprocess.Popen([samtools, 'view', bam],stdout=subprocess.PIPE)
		head_st = subprocess.Popen(['head', '-10000'], stdin=run_samtools.stdout,stdout=temp_fh)
		head_st.wait()
	read_lengths = []		
	with open(temp_sam, 'r') as temp_fh:
		for line in temp_fh:
			line = line.split(delim)
			read = line[9]
			# print(read)
			read_lengths.append(len(read))
	med_rl = statistics.median(read_lengths)
	# print(med_rl)
	return(med_rl)

def get_read_counts(bam):
	temp_txt = bam.split('.')[0] + '.read_count_temp.txt'
	with open(temp_txt, 'w') as temp_fh:
		run_samtools = subprocess.Popen([samtools, 'flagstat', bam], stdout=temp_fh)
		run_samtools.wait()
	with open(temp_txt, 'r') as temp_fh:
		lc = 0
		for line in temp_fh:
			lc +=1
			line = line.split(' ')
			if lc == 1:
				if line[4] == 'total':
					total = line[0]
				else:
					print('issue with total count in flagstat file', temp_txt)
			elif lc == 5:
				if line[3] == 'mapped':
					mapped = line[0]
				else:
					print('issue with mapped count in flagstat file', temp_txt)
	return(total, mapped)

def get_read_coverage_info(bams, outfile):
	with open(outfile, 'w') as out_fh:
		header = ['bam', 'median_read_length', 'total_reads', 'mapped_reads', '\n']
		out_fh.write(delim.join(header))
		for bam in bams:
			print(bam)
			median_read_length = caculate_median_read_length(bam)
			print(median_read_length)
			total_reads, mapped_reads = get_read_counts(bam)
			print(total_reads, mapped_reads)
			line_out = [bam, str(int(median_read_length)), total_reads, mapped_reads, '\n']
			out_fh.write(delim.join(line_out))




##run methods

##parameters
# project_name = 'kim_cblm_0818'
# pedigree_names = ['3C-4', 'LR11-330', 'LR11-331', 'LR05-160', 'LR05-203a1', 'LR10-102', 'LR11-033', 'LR04-185', 'LR04-414', 'LR06-105', 'LR14-098', 
# 		'LR03-305', 'LR04-106', 'LR04-186', 'LR05-203', 'LR05-354', 'LR08-323', 'LR08-390', 'LR08-396', 'LR10-230', 'LR11-169', 'LR12-032', 'LR12-115', 
# 		'LR12-434', 'LR12-443', 'LR12-463', 'LR12-464', 'LR13-002', 'LR13-003', 'LR13-037', 'LR13-085', 'LR13-153', 'LR13-199', 'LR13-200', 'LR13-315', 
# 		'LR03-055', 'LR03-206', 'LR03-274', 'LR03-332', 'LR04-020', 'LR04-371', 'LR06-085', 'LR09-227', 'LR09-280', 'LR11-241', 'LR03-077', 'LR03-169', 
# 		'LR03-223', 'LR03-278', 'LR03-298', 'LR04-017', 'LR04-084', 'LR04-208', 'LR04-233', 'LR05-118', 'LR05-396', 'LR12-475', 'LR01-079', 'LR02-263', 
# 		'LR03-039', 'LR03-304', 'LR03-340', 'LR04-022a1', 'LR04-222', 'LR04-399', 'LR05-162', 'LR05-265', 'LR05-398', 'LR08-056', 'LR10-016', 'LR12-313a2', 
# 		'LR12-439', 'DWM3', 'DWM10', 'DWM13', 'LR04-239', 'LR14-071', 'LR14-221', 'LR16-079', 'LR16-451', 'LR10-199', 'LR06-157', 'LR10-222', 'LR11-152', 
# 		'LR03-120', 'LR04-376', 'LR05-120', 'LR06-278', 'LR08-002', 'LR03-130', 'LR12-316', 'LR06-207', 'LR09-416', 'LR10-026', 'LR10-228a1', 'LR11-042', 
# 		'LR05-007', 'LR09-023', 'LR10-243', 'LR05-035']
project_name = 'kim_cblm_0219'
pedigree_names = ['LR05-203', 'LR05-250', 'LR06-137', 'LR09-055', 'LR09-227', 'LR09-416', 'LR10-228a1', 'LR11-169', 'LR11-330', 'LR11-331', 'LR14-221', 'LR15-267', 
		'3C-4', 'LR01-079', 'LR02-263', 'LR03-039', 'LR03-055', 'LR03-077', 'LR03-120', 'LR03-130', 'LR03-169', 'LR03-206', 'LR03-223', 'LR03-274', 'LR03-278', 
		'LR03-298', 'LR03-304', 'LR03-305', 'LR03-332', 'LR03-340', 'LR04-017', 'LR04-020', 'LR04-022a1', 'LR04-084', 'LR04-106', 'LR04-185', 'LR04-186', 'LR04-208', 
		'LR04-222', 'LR04-233', 'LR04-239', 'LR04-371', 'LR04-376', 'LR04-399', 'LR04-414', 'LR05-007', 'LR05-035', 'LR05-118', 'LR05-120', 'LR05-160', 'LR05-162', 
		'LR05-203a1', 'LR05-265', 'LR05-354', 'LR05-396', 'LR05-398', 'LR06-085', 'LR06-105', 'LR06-157', 'LR06-207', 'LR06-278', 'LR08-002', 'LR08-056', 'LR08-323', 
		'LR08-390', 'LR08-396', 'LR09-023', 'LR09-280', 'LR10-016', 'LR10-026', 'LR10-102', 'LR10-199', 'LR10-222', 'LR10-230', 'LR10-243', 'LR11-033', 'LR11-042', 
		'LR11-152', 'LR11-241', 'LR12-032', 'LR12-115', 'LR12-313a2', 'LR12-316', 'LR12-434', 'LR12-439', 'LR12-443', 'LR12-463', 'LR12-464', 'LR12-475', 'LR13-002', 
		'LR13-003', 'LR13-037', 'LR13-085', 'LR13-153', 'LR13-199', 'LR13-200', 'LR13-315', 'LR14-071', 'LR14-098', 'LR16-079', 'LR16-451']


coverage_file_suffix = '.coverage_data.txt'
combined_cov_file = project_name + coverage_file_suffix
brainspan_data = 'BrainSpan_prenatal_cb_rpkm.csv'
gtex_data = 'GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm_cbonly.csv'
##combined files manually, changed id and added individual files
all_variant_data = 'all_exome_data_std_pipeline.manual_0818.xls'
cohort_var_data = 'ped_specific_std_pipeline_expressed_0818.xls'
counts_file = project_name + '.var_counts.xls'
counts_file_no_exp = project_name + '.var_counts_no_exp.xls'
##for more coverage info
# bams_to_analyze = glob.glob('*.bam')
# bams_to_analyze = ['LR10-016.bwa_gatk.bam', 'LR10-016m.bwa_gatk.bam']
extra_cov_file = project_name + '.extra_coverage_info.txt'
new_bams = ['LR05-250.bwa_gatk.bam', 'LR06-137.bwa_gatk.bam', 'LR09-055.bwa_gatk.bam', 'LR15-267.bwa_gatk.bam']


##get coverage data
# get_coverage_data(pedigree_names, coverage_file_suffix, combined_cov_file)
##get extra coverge info - total reads, mapped reads and read length
# get_read_coverage_info(bams_to_analyze, extra_cov_file)
get_read_coverage_info(new_bams, extra_cov_file)


##get variants for all peds -- so all the same format
# exp_dict = make_dict_from_expression_files(brainspan_data, gtex_data)
# get_peds_add_expression_data(pedigree_names, exp_dict, all_variant_data, cohort_var_data)
##get counts from that variant file
# get_counts_from_var_file(cohort_var_data, counts_file)
# get_counts_from_var_file(cohort_var_data, counts_file_no_exp)


