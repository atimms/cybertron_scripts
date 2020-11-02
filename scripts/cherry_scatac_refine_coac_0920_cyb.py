#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parms
delim = '\t'


##methods
def make_dict_from_info_file(infile, file_type):
	peak_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				idx_num = line[0]
				if file_type == 'rna':
					peak_info = line[1:4] + [line[6]]
				elif file_type == 'peak':
					peak_info = line[1:4] + line[6:16]
				peak_dict[idx_num] = peak_info
	return(peak_dict)

def format_p2g_file(infile, outfile, p2g_dict, atac_dict):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim) 
			if lc == 1:
				header = line + ['atac_chr', 'atac_start', 'atac_end', 'score', 'replicateScoreQuantile', 'groupScoreQuantile', 
				'Reproducibility', 'GroupReplicate', 'distToGeneStart', 'nearestGene', 'peakType', 'distToTSS', 'nearestTSS', 
				'rna_chr', 'rna_start', 'rna_end', 'rna_gene']
				out_fh.write(delim.join(header) + '\n')
			else:
				idxATAC = line[0]
				idxRNA = line[1]
				atac_info = atac_dict[idxATAC]
				rna_info = p2g_dict[idxRNA]
				line_out = line + atac_info + rna_info
				out_fh.write(delim.join(line_out) + '\n')

def format_coac_file(infile, outfile, atac_dict):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim) 
			if lc == 1:
				header = line + ['query_chr', 'query_start', 'query_end', 'query_score', 'query_replicateScoreQuantile', 
				'query_groupScoreQuantile', 'query_Reproducibility', 'query_GroupReplicate', 'query_distToGeneStart', 
				'query_nearestGene', 'query_peakType', 'query_distToTSS', 'query_nearestTSS', 'subject_chr', 'subject_start', 
				'subject_end', 'subject_score', 'subject_replicateScoreQuantile', 'subject_groupScoreQuantile', 'subject_Reproducibility', 
				'subject_GroupReplicate', 'subject_distToGeneStart', 'subject_nearestGene', 'subject_peakType', 'subject_distToTSS', 
				'subject_nearestTSS']
				out_fh.write(delim.join(header) + '\n')
			else:
				idxquery = line[0]
				idxsubject = line[1]
				query_info = atac_dict[idxquery]
				subject_info = atac_dict[idxsubject]
				subject_peaktype = subject_info[10]
				# print(query_info, subject_info)
				if subject_peaktype == 'Promoter':
					line_out = line + query_info + subject_info
					out_fh.write(delim.join(line_out) + '\n')



def make_dict_from_formatted_file(infile, file_type):
	out_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				idx_num = line[0]
				if file_type == 'rna':
					idxRNA = line[1]
					correlation = float(line[2])
					fdr = line[3]
					pos = '_'.join(line[19:22])
					gene = line[22]
					if idx_num in out_dict:
						# print(correlation, out_dict[idx_num][1][2])
						if correlation > max(out_dict[idx_num][1][2]):
							out_dict[idx_num][0] = gene
						# print(idx_num, out_dict[idx_num])
						out_dict[idx_num][1][0].append(gene)
						out_dict[idx_num][1][1].append(pos)
						out_dict[idx_num][1][2].append(correlation)
						out_dict[idx_num][1][3].append(fdr)
						out_dict[idx_num][1][4].append(idxRNA)

					else:
						out_dict[idx_num] = [gene, [[gene], [pos], [correlation], [fdr], [idxRNA]]]
				elif file_type == 'peak':
					idxsubject = line[1]
					correlation = float(line[3])
					pos = '_'.join(line[17:20])
					gene = line[26]
					distToGeneStart = line[25]
					distToTSS = line[28]
					nearestTSS = line[29]
					GroupReplicate = line[24]
					score = line[20]
					if idx_num in out_dict:
						# print(correlation, out_dict[idx_num][1][2])
						if correlation > max(out_dict[idx_num][1][2]):
							out_dict[idx_num][0] = gene
						# print(idx_num, out_dict[idx_num])
						out_dict[idx_num][1][0].append(gene)
						out_dict[idx_num][1][1].append(pos)
						out_dict[idx_num][1][2].append(correlation)
						out_dict[idx_num][1][3].append(idxsubject)
						out_dict[idx_num][1][4].append(distToGeneStart)
						out_dict[idx_num][1][5].append(distToTSS)
						out_dict[idx_num][1][6].append(nearestTSS)
						out_dict[idx_num][1][7].append(GroupReplicate)
						out_dict[idx_num][1][8].append(score)
					else:
						out_dict[idx_num] = [gene, [[gene], [pos], [correlation], [idxsubject], [distToGeneStart], [distToTSS], [nearestTSS], [GroupReplicate], [score]]]
	return(out_dict)


def combine_info_per_peak(in_file, coac_file, p2g_file, out_file):
	peak_dict = make_dict_from_formatted_file(coac_file, 'peak')
	p2g_dict = make_dict_from_formatted_file(p2g_file, 'rna')
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc == 1:
				header = ['idx'] + line[:3] + line[5:15] + ['p2g_most_corr_gene', 'p2g_genes', 'p2g_positions', 'p2g_correlation_scores', 
				'p2g_fdr_scores', 'p2g_RNAidxs', 'coac_most_corr_gene', 'coac_genes', 'coac_positions', 'coac_correlation_scores', 'coac_idxs' , 
				'coac_dist_to_gene_starts', 'coac_dist_to_tss', 'coac_tss', 'coac_GroupReplicate', 'coac_score']
				out_fh.write(delim.join(header) + '\n')
			else:
				idx_num = line[0]
				if idx_num in p2g_dict:
					# print(p2g_dict[idx_num][1])
					p2g_gene = p2g_dict[idx_num][0]
					p2g_genes = ', '.join(p2g_dict[idx_num][1][0])
					p2g_pos = ', '.join(p2g_dict[idx_num][1][1])
					p2g_cors = [str(i) for i in p2g_dict[idx_num][1][2]]
					p2g_cors = ', '.join(p2g_cors)
					p2g_fdrs = ', '.join(p2g_dict[idx_num][1][3])
					p2g_rnaidxs = ', '.join(p2g_dict[idx_num][1][4])
					line_out = line[:4] + line[6:16] + [p2g_gene, p2g_genes, p2g_pos, p2g_cors, p2g_fdrs, p2g_rnaidxs]
				else:
					line_out = line[:4] + line[6:16] + ['', '', '', '', '', '']
				if idx_num in peak_dict:
					coac_gene = peak_dict[idx_num][0]
					coac_genes = ', '.join(peak_dict[idx_num][1][0])
					coac_pos = ', '.join(peak_dict[idx_num][1][1])
					coac_cors = [str(i) for i in peak_dict[idx_num][1][2]]
					coac_cors = ', '.join(coac_cors)
					coac_idxs = ', '.join(peak_dict[idx_num][1][3])
					coac_dist_to_genes = ', '.join(peak_dict[idx_num][1][4])
					coac_dist_to_tsss = ', '.join(peak_dict[idx_num][1][5])
					coac_tsss = ', '.join(peak_dict[idx_num][1][6])
					coac_group = ', '.join(peak_dict[idx_num][1][7])
					coac_score = ', '.join(peak_dict[idx_num][1][8])
					line_out.extend([coac_gene, coac_genes, coac_pos, coac_cors,coac_idxs,coac_dist_to_genes,coac_dist_to_tsss,coac_tsss,coac_group,coac_score])
				else:
					line_out.extend(['', '', '', '','','','','','',''])
				out_fh.write(delim.join(line_out) + '\n')

##run methods
work_dir = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all/co_acc_analysis_0920'
os.chdir(work_dir)
peak_file = 'human_all.min_cell20.peak_info.m2_q0.01.txt'
rna_file = 'human_all.min_cell20.m2_q0.01.rnaseq_info.txt'
coac_infile = 'human_all.min_cell20.co_acc.m2_q0.01.cor4.txt'
coac_formatted = 'human_all.co_acc.m2_q0.01.cor4.formatted.txt'
p2g_infile = 'human_all.min_cell20.p2g.m2_q0.01.cor4.txt'
p2g_formatted = 'human_all.p2g.m2_q0.01.cor4.formatted.txt'
combined_file = 'human_all.m2_q0.01.cor4.combined.xls'

##make dict from peak and rna info
peak_info_dict = make_dict_from_info_file(peak_file,'peak')
rna_info_dict =  make_dict_from_info_file(rna_file, 'rna')

##make new files with peak/rna info - also remove non promoter subject peaks
format_p2g_file(p2g_infile, p2g_formatted, rna_info_dict, peak_info_dict)
format_coac_file(coac_infile, coac_formatted, peak_info_dict)

##combined data into single file
combine_info_per_peak(peak_file, coac_formatted, p2g_formatted, combined_file)

##summarize data and subset and perform some test

