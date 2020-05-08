#!/usr/bin/env python
import sys
import subprocess
import os
import glob

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/exomes/exome_project_2019'
os.chdir(working_dir)

##file
# gof_lof_genes = 'gof_lof_genes.txt'

std_anal_header = ['chrom', 'start', 'end', 'ref', 'alt', 'gene', 'transcript', 'biotype', 'aa_change', 'impact', 'gerp_bp_score', 'cadd_scaled', 'max_aaf_all', 'annovar_ts_info', 
	'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG', 'cosmic88_coding', 'cosmic88_noncoding', 'dbdb inheritance', 'dbdb phenotype', 'dbdb syndrome', 'dbdb loe', 'mgi phenotype', 
	'omim phenotype', 'previous exomes', 'clinvar_gene_phenotype', 'rmsk', 'in_segdup', 'gnomMD_exome_AC', 'gnomMD_exome_AD', 'gnomMD_exome_AF', 'gnomMD_exome_nhom', 'gnomMD_exome_AC_non_neuro', 
	'gnomMD_exome_AD_non_neuro', 'gnomMD_exome_AF_non_neuro', 'gnomMD_exome_nhom_non_neuro', 'gnomMD_genome_AC', 'gnomMD_genome_AD', 'gnomMD_genome_AF', 'gnomMD_genome_nhom', 
	'gnomMD_genome_AC_non_neuro', 'gnomMD_genome_AD_non_neuro', 'gnomMD_genome_AF_non_neuro', 'gnomMD_genome_nhom_non_neuro', 'analysis']

mosaic_anal_header = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 
		'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 
		'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 
		'DANN_score', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred', 'integrated_fitCons_score', 'integrated_confidence_value', 
		'GERP++_RS', 'phyloP7way_vertebrate', 'phyloP20way_mammalian', 'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'PopFreqMax', '1000G_ALL', '1000G_AFR', '1000G_AMR', 
		'1000G_EAS', '1000G_EUR', '1000G_SAS', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'ESP6500siv2_ALL', 'ESP6500siv2_AA', 'ESP6500siv2_EA', 
		'CG46', 'avsnp147', 'CLINSIG', 'CLNDBN', 'CLNACC', 'CLNDSDB', 'CLNDSDBID', 'cosmic83_noncoding', 'cosmic83_coding','q', 'coverage', 'filter', 'format', 'info', 'var_combined', 'sample_alt_reads']

##methods


def check_for_var(pedigree_dict, solved_outfile, issue_outfile):
	var_count_dict = {}
	for pedigree in pedigree_dict:
		# print(pedigree)
		# print(pedigree_dict[pedigree])
		solved_status = pedigree_dict[pedigree][17]
		solved_gene = pedigree_dict[pedigree][16]
		solved_var = pedigree_dict[pedigree][18]
		ped_type = pedigree_dict[pedigree][0]
		dx_group = pedigree_dict[pedigree][3]
		# print(pedigree, solved_gene, solved_status)
		if solved_gene != '':
			var_count_dict[pedigree] = [0,[], [solved_gene, solved_status, solved_var, ped_type, dx_group]]
			# print(solved_gene, solved_status)
			std_file = pedigree + '.std_analysis.xls'
			with open(std_file, "r") as std_fh: 
				for line in std_fh:
					line = line.strip('\n').split(delim)
					gene = line[5]
					if len(solved_gene.split('/')) >1:
						print(gene, solved_gene.split('/'))
					if gene in solved_gene.split('/'):
						var_count_dict[pedigree][0] += 1
						var_count_dict[pedigree][1].append(line[:47] + [line[-2]])

	with open(solved_outfile, "w") as sout_fh, open(issue_outfile, "w") as iout_fh:
		iout_fh.write(delim.join(['ped', 'gene', 'status', 'var_info', 'ped_type', 'dxgroup1']) + '\n')
		sout_fh.write(delim.join(['ped', 'gene', 'status', 'var_info', 'ped_type', 'dxgroup1'] + std_anal_header) + '\n')
		for vc in var_count_dict:
			print(vc, var_count_dict[vc][0], var_count_dict[vc][2])
			##if can't find a var
			if var_count_dict[vc][0] == 0:
				iout_fh.write(delim.join([vc] + var_count_dict[vc][2]) + '\n')
			else:
				for var in var_count_dict[vc][1]:
					sout_fh.write(delim.join([vc] + var_count_dict[vc][2] + var) + '\n')

	print('all var count:', len(var_count_dict))


def read_ped_info_file_and_analyze(in_file, prefix_name):
	ped_dict = {}
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.strip('\n').split(delim)
			##get solved
			if lc >1 and line[5] == 'include':
				ped_id = line[0]
				ped_info = line[1:]
				if ped_id not in ped_dict:
					ped_dict[ped_id] = ped_info
				else:
					print(ped_id, 'seen multiple times')
	print(len(ped_dict))
	##what to summarize? solved% counts per dx, type etc
	# summarize_data(ped_dict)
	##check for variant in analysis file
	solved_var_file = prefix_name + '.solved_vars.xls'
	peds_no_var = prefix_name + '.no_vars.xls'
	check_for_var(ped_dict, solved_var_file, peds_no_var)


def check_de_novo_filters(gene_dict, dn_dict):
	kept, removed = 0,0
	for ped in dn_dict:
		for var_info in dn_dict[ped]:
			ped_name = var_info[0]
			gene = var_info[11]
			status = var_info[2]
			ped_type = var_info[4]
			dx_group = var_info[5]
			max_aaf = var_info[18]
			gerp = var_info[16]
			cadd = var_info[17]
			gnomad_ex_af = var_info[39]
			rmsk = var_info[35]
			seg_dup = var_info[36]
			biotype = var_info[13]
			if gene in gene_dict:
				genelist = gene_dict[gene]
			else:
				genelist = 'not seen'
			if cadd == 'None':
				cadd = '20'
			if gnomad_ex_af == 'na':
				gnomad_ex_af = '0'
			# if float(cadd) > 10 and float(gnomad_ex_af) <= 0.001 and genelist != 'not seen':
			# if float(cadd) > 10 and float(gnomad_ex_af) <= 0.001 and float(gerp) >= 0:
			if float(cadd) > 10 and float(gnomad_ex_af) <= 0.001 and biotype == 'protein_coding':
				kept += 1
				print(dx_group, status, gene, genelist)
			else:
				removed += 1
				# print(ped_name, dx_group, status, gene, cadd, gnomad_ex_af,dx_group, gerp)
				# print(dx_group, status, gene, ped_name)
	print(kept, removed)





def check_filter_on_solved(in_file, gene_file):
	genelist_dict, dn_var_dict  = {}, {}
	##make dict from genelist
	with open(gene_file, "r") as gf_fh:
		lc = 0
		for line in gf_fh:
			lc += 1
			if lc >1:
				line = line.strip('\n').split(delim)
				gene = line[0]
				group_sub = line[1] + '_' + line[2]
				if gene in genelist_dict:
					genelist_dict[gene].append(group_sub)
				else:
					genelist_dict[gene] = [group_sub]
	# print(genelist_dict)
	# for g in genelist_dict:
	# 	print(g, genelist_dict[g])
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc >1:
				line = line.strip('\n').split(delim)
				ped_name = line[0]
				analysis_type = line[53]
				if analysis_type == 'de_novo' or analysis_type == 'x_linked_de_novo':
					if ped_name in dn_var_dict:
						dn_var_dict[ped_name].append(line)
						print(ped_name + ' has multiple de novos')
					else:
						dn_var_dict[ped_name] = [line]
	##check dict
	# for var in dn_var_dict:
	# 	print(dn_var_dict[var])
	##check de novos
	check_de_novo_filters(genelist_dict, dn_var_dict)
	return(genelist_dict)


def get_denovo_vars_and_filter(in_file, prefix_name, solved_dn_genes, good_solved_genes, genelist_dict):
	ped_dict = {}
	dn_file = prefix_name + '.denovos_all.xls'
	dn_file_filtered = prefix_name + '.denovos_filtered.xls'
	dn_file_solved = prefix_name + '.denovos_solved.xls'
	dn_file_12solved = prefix_name + '.denovos_12solved.xls'
	dn_file_gl = prefix_name + '.denovos_genelists.xls'
	##get all denovo vars in one file
	with open(in_file, "r") as in_fh, open(dn_file, "w") as out_fh, open(dn_file_filtered, "w") as outf_fh, open(dn_file_solved, "w") as outs_fh, open(dn_file_12solved, "w") as outs12_fh, open(dn_file_gl, "w") as outgl_fh:
		##add header
		header = ['ped', 'gene', 'status', 'var_info', 'dxgroup1'] + std_anal_header
		out_fh.write(delim.join(header) + '\n')
		outf_fh.write(delim.join(header) + '\n')
		outs_fh.write(delim.join(header) + '\n')
		outs12_fh.write(delim.join(header) + '\n')
		outgl_fh.write(delim.join(header + ['genelist']) + '\n')
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.strip('\n').split(delim)
			##get peds we want and are trios
			if lc >1 and line[5] == 'include' and line[1] == 'trio':
				ped_id = line[0]
				solved_status = line[18]
				solved_gene = line[17]
				solved_var = line[19]
				ped_type = line[1]
				dx_group = line[4]
				std_file = ped_id + '.std_analysis.xls'
				# print(ped_id, line[5], line[1])
				with open(std_file, "r") as std_fh: 
					for line in std_fh:
						line = line.strip('\n').split(delim)
						anal = line[-2]
						cadd = line[11]
						gnomad_ex_af = line[33]
						gene = line[5]
						biotype = line[7]
						# print(ped_id, anal, line)

						if anal == 'de_novo' or anal == 'x_linked_de_novo':
							line_out = [ped_id, solved_gene, solved_status, solved_var, dx_group] + line[:47] + [line[-2]]
							out_fh.write(delim.join(line_out) + '\n')
							if cadd == 'None':
								cadd = '20'
							if gnomad_ex_af == 'na':
								gnomad_ex_af = '0'
							if float(cadd) > 10 and float(gnomad_ex_af) <= 0.001 and biotype == 'protein_coding':
								outf_fh.write(delim.join(line_out) + '\n')
							if gene in solved_dn_genes:
								outs_fh.write(delim.join(line_out) + '\n')
							if gene in good_solved_genes:
								outs12_fh.write(delim.join(line_out) + '\n')
							if gene in genelist_dict:
								line_out = line_out + [' '.join(genelist_dict[gene])]
								outgl_fh.write(delim.join(line_out) + '\n')
						

def get_recessive_vars_and_filter(in_file, prefix_name, solved_dn_genes, good_solved_genes, genelist_dict):
	ped_dict = {}
	rec_file = prefix_name + '.recessive_all.xls'
	rec_file_filtered = prefix_name + '.recessive_filtered.xls'
	rec_file_solved = prefix_name + '.recessive_solved.xls'
	rec_file_12solved = prefix_name + '.recessive_12solved.xls'
	rec_file_gl = prefix_name + '.recessive_genelists.xls'
	##get all denovo vars in one file
	with open(in_file, "r") as in_fh, open(rec_file, "w") as out_fh, open(rec_file_filtered, "w") as outf_fh, open(rec_file_solved, "w") as outs_fh, open(rec_file_12solved, "w") as outs12_fh, open(rec_file_gl, "w") as outgl_fh:
		##add header
		header = ['ped', 'gene', 'status', 'var_info', 'dxgroup1'] + std_anal_header
		out_fh.write(delim.join(header) + '\n')
		outf_fh.write(delim.join(header) + '\n')
		outs_fh.write(delim.join(header) + '\n')
		outs12_fh.write(delim.join(header) + '\n')
		outgl_fh.write(delim.join(header + ['genelist']) + '\n')
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.strip('\n').split(delim)
			##get peds we want and are trios
			if lc >1 and line[5] == 'include' and line[1] == 'trio':
				ped_id = line[0]
				solved_status = line[18]
				solved_gene = line[17]
				solved_var = line[19]
				ped_type = line[1]
				dx_group = line[4]
				std_file = ped_id + '.std_analysis.xls'
				# print(ped_id, line[5], line[1])
				with open(std_file, "r") as std_fh: 
					for line in std_fh:
						line = line.strip('\n').split(delim)
						anal = line[-2]
						cadd = line[11]
						gnomad_ex_af = line[33]
						gene = line[5]
						biotype = line[7]
						# print(ped_id, anal, line)
						# print(cadd, gnomad_ex_af)
						if anal == 'autosomal_recessive' or anal == 'comp_hets':
							line_out = [ped_id, solved_gene, solved_status, solved_var, dx_group] + line[:47] + [line[-2]]
							out_fh.write(delim.join(line_out) + '\n')
							if cadd == 'None':
								cadd = '20'
							if gnomad_ex_af == 'na' or gnomad_ex_af == '.':
								gnomad_ex_af = '0'
							if float(cadd) > 10 and float(gnomad_ex_af) <= 0.001 and biotype == 'protein_coding':
								outf_fh.write(delim.join(line_out) + '\n')
							if gene in solved_dn_genes:
								outs_fh.write(delim.join(line_out) + '\n')
							if gene in good_solved_genes:
								outs12_fh.write(delim.join(line_out) + '\n')
							if gene in genelist_dict:
								line_out = line_out + [' '.join(genelist_dict[gene])]
								outgl_fh.write(delim.join(line_out) + '\n')

def get_xlinked_vars_and_filter(in_file, prefix_name, solved_dn_genes, good_solved_genes, genelist_dict):
	ped_dict = {}
	xlinked_file = prefix_name + '.xlinked_all.xls'
	xlinked_file_filtered = prefix_name + '.xlinked_filtered.xls'
	xlinked_file_solved = prefix_name + '.xlinked_solved.xls'
	xlinked_file_12solved = prefix_name + '.xlinked_12solved.xls'
	xlinked_file_gl = prefix_name + '.xlinked_genelists.xls'
	##get all denovo vars in one file
	with open(in_file, "r") as in_fh, open(xlinked_file, "w") as out_fh, open(xlinked_file_filtered, "w") as outf_fh, open(xlinked_file_solved, "w") as outs_fh, open(xlinked_file_12solved, "w") as outs12_fh, open(xlinked_file_gl, "w") as outgl_fh:
		##add header
		header = ['ped', 'gene', 'status', 'var_info', 'dxgroup1'] + std_anal_header
		out_fh.write(delim.join(header) + '\n')
		outf_fh.write(delim.join(header) + '\n')
		outs_fh.write(delim.join(header) + '\n')
		outs12_fh.write(delim.join(header) + '\n')
		outgl_fh.write(delim.join(header + ['genelist']) + '\n')
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.strip('\n').split(delim)
			##get peds we want and are trios
			if lc >1 and line[5] == 'include' and line[1] == 'trio':
				ped_id = line[0]
				solved_status = line[18]
				solved_gene = line[17]
				solved_var = line[19]
				ped_type = line[1]
				dx_group = line[4]
				std_file = ped_id + '.std_analysis.xls'
				# print(ped_id, line[5], line[1])
				with open(std_file, "r") as std_fh: 
					for line in std_fh:
						line = line.strip('\n').split(delim)
						anal = line[-2]
						cadd = line[11]
						gnomad_ex_af = line[33]
						gene = line[5]
						biotype = line[7]
						# print(ped_id, anal, line)
						if anal == 'x_linked_recessive':
							line_out = [ped_id, solved_gene, solved_status, solved_var, dx_group] + line[:47] + [line[-2]]
							out_fh.write(delim.join(line_out) + '\n')
							if cadd == 'None':
								cadd = '20'
							if gnomad_ex_af == 'na':
								gnomad_ex_af = '0'
							if float(cadd) > 10 and float(gnomad_ex_af) <= 0.001 and biotype == 'protein_coding':
								outf_fh.write(delim.join(line_out) + '\n')
							if gene in solved_dn_genes:
								outs_fh.write(delim.join(line_out) + '\n')
							if gene in good_solved_genes:
								outs12_fh.write(delim.join(line_out) + '\n')
							if gene in genelist_dict:
								line_out = line_out + [' '.join(genelist_dict[gene])]
								outgl_fh.write(delim.join(line_out) + '\n')

def get_mosaic_vars_and_filter(in_file, prefix_name, solved_dn_genes, good_solved_genes, genelist_dict, meg_genes):
	ped_dict = {}
	mosaic_file = prefix_name + '.mosaic_all.xls'
	mosaic_file_filtered = prefix_name + '.mosaic_filtered.xls'
	mosaic_file_solved = prefix_name + '.mosaic_solved.xls'
	mosaic_file_12solved = prefix_name + '.mosaic_12solved.xls'
	mosaic_file_gl = prefix_name + '.mosaic_genelists.xls'
	mosaic_file_meg = prefix_name + '.mosaic_meg.xls'
	##get all denovo vars in one file
	with open(in_file, "r") as in_fh, open(mosaic_file, "w") as out_fh, open(mosaic_file_filtered, "w") as outf_fh, open(mosaic_file_solved, "w") as outs_fh, open(mosaic_file_12solved, "w") as outs12_fh, open(mosaic_file_gl, "w") as outgl_fh, open(mosaic_file_meg, "w") as outmeg_fh:
		##add header
		header = ['ped', 'gene', 'status', 'var_info', 'dxgroup1'] + mosaic_anal_header
		out_fh.write(delim.join(header) + '\n')
		outf_fh.write(delim.join(header) + '\n')
		outs_fh.write(delim.join(header) + '\n')
		outs12_fh.write(delim.join(header) + '\n')
		outgl_fh.write(delim.join(header + ['genelist']) + '\n')
		outmeg_fh.write(delim.join(header) + '\n')
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.strip('\n').split(delim)
			##get peds we want and are trios
			if lc >1 and line[5] == 'include' and line[1] == 'trio':
				ped_id = line[0]
				solved_status = line[18]
				solved_gene = line[17]
				solved_var = line[19]
				ped_type = line[1]
				dx_group = line[4]
				mosaic_infiles = glob.glob(ped_id + '*pisces.xls')
				# print(ped_id, mosaic_infiles)
				for mosaic_infile in mosaic_infiles:
					with open(mosaic_infile, "r") as mosin_fh:
						lc2 = 0
						for line in mosin_fh:
							lc2 += 1
							if lc2 == 1:
								line = line.strip('\n').split(delim)
								print(ped_id,mosaic_infile, len(line))
							if lc2 > 1:
								line = line.strip('\n').split(delim)
								cadd = line[30]
								exac_af = line[53]
								gene = line[6]
								line_out = [ped_id, solved_gene, solved_status, solved_var, dx_group] + line
								out_fh.write(delim.join(line_out) + '\n')
								if cadd == '.':
									cadd = '20'
								if exac_af == '.':
									exac_af = '0'
								if float(cadd) > 10 and float(exac_af) <= 0.001:
									outf_fh.write(delim.join(line_out) + '\n')
								if gene in solved_dn_genes:
									outs_fh.write(delim.join(line_out) + '\n')
								if gene in good_solved_genes:
									outs12_fh.write(delim.join(line_out) + '\n')
								if gene in meg_genes:
									outmeg_fh.write(delim.join(line_out) + '\n')
								if gene in genelist_dict:
									line_out = line_out + [' '.join(genelist_dict[gene])]
									outgl_fh.write(delim.join(line_out) + '\n')

def get_single_duo_vars_and_filter(in_file, prefix_name):
	ped_dict = {}
	sd_all_file = prefix_name + '.single_duo_all.xls'
	##get all singles vars in one file
	with open(in_file, "r") as in_fh, open(sd_all_file, "w") as out_fh:
		##add header
		header = ['ped', 'gene', 'status', 'var_info', 'dxgroup1'] + std_anal_header
		out_fh.write(delim.join(header) + '\n')
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.strip('\n').split(delim)
			##get peds we want and are trios
			if lc >1 and line[5] == 'include' and (line[1] == 'singleton' or line[1] == 'duo'):
				print(line[0], line[1], line[5])
				# ped_id = line[0]
				# solved_status = line[18]
				# solved_gene = line[17]
				# solved_var = line[19]
				# ped_type = line[1]
				# dx_group = line[4]
				# std_file = ped_id + '.std_analysis.xls'
				# # print(ped_id, line[5], line[1])
				# with open(std_file, "r") as std_fh: 
				# 	for line in std_fh:
				# 		line = line.strip('\n').split(delim)
				# 		anal = line[-2]
				# 		cadd = line[11]
				# 		gnomad_ex_af = line[33]
				# 		gene = line[5]
				# 		biotype = line[7]
				# 		# print(ped_id, anal, line)
				# 		if anal == 'x_linked_recessive':
				# 			line_out = [ped_id, solved_gene, solved_status, solved_var, dx_group] + line[:47] + [line[-2]]
				# 			out_fh.write(delim.join(line_out) + '\n')
				# 			if cadd == 'None':
				# 				cadd = '20'
				# 			if gnomad_ex_af == 'na':
				# 				gnomad_ex_af = '0'
				# 			if float(cadd) > 10 and float(gnomad_ex_af) <= 0.001 and biotype == 'protein_coding':
				# 				outf_fh.write(delim.join(line_out) + '\n')
				# 			if gene in solved_dn_genes:
				# 				outs_fh.write(delim.join(line_out) + '\n')
				# 			if gene in good_solved_genes:
				# 				outs12_fh.write(delim.join(line_out) + '\n')
				# 			if gene in genelist_dict:
				# 				line_out = line_out + [' '.join(genelist_dict[gene])]
				# 				outgl_fh.write(delim.join(line_out) + '\n')




##run methods
# analysis_file_dir = '/active/mirzaa_g/Projects/Projects-Active/ML_WES/Tier1/analysis_files'
project_name = 'solve_exomes_1119'
ped_info_file = 'exome_ped_data_1119.txt'


##step 1. check for variants in analysis files
##step 2. test filters, make sure solved are found
##step 3. find vars; apply filters, genelist, omim and clinvar, multiples within dx

##read original file and analyze:
##1. to get solved vars and problem vars
# read_ped_info_file_and_analyze(ped_info_file, project_name)

##2. apply some filters and learn (initially on de novos)
genelist_file = 'genelists_1119.txt'
solved_exomes = 'solve_exomes_1119.solved_vars_updated.xls'
genes_from_genelists_dict = check_filter_on_solved(solved_exomes, genelist_file)

##3. after manually checked the solved vars, apply filters (so get all and good i.e. 1 or 2a)
all_dn_solved_genes = ['DDX3X', 'TUBA1A', 'KIF1A', 'PIK3R2', 'DYNC1H1', 'ARID1B', 'MED13L', 'KIF11', 'TUBG1', 'CLCN7', 'CRNKL1', 'DDHD2', 'CHEK1', 'MAGEC1', 'SEMA6B', 
	'MACF1', 'FGFR1', 'BACH1', 'BRAF', 'DCST1', 'CACNA1G', 'AUTS2', 'TUBB2B', 'PDGFRB', 'STXBP1', 'TUBB2A', 'PPP1CB', 'ITSN1', 'KIAA1211', 'VPS13B', 'SLC1A2', 'BCL11A', 
	'SSH3', 'CREBBP', 'TUBB8', 'CDC42', 'CASK', 'FRMD1', 'GRIA2', 'ZEB2', 'SHOC2', 'ACTG1', 'KCNQ2', 'ATP1A3', 'BUB1', 'IGFBP2', 'RARG', 'TRIO', 'TAF15', 'SFI1', 'PIK3CA', 
	'MTOR', 'SETD2', 'FOXP1', 'SP4', 'SLC25A5', 'TRPC4', 'DNAJC13', 'KIDINS220', 'WDR37', 'KCTD1', 'KRT36', 'ADCY1', 'LRRN4CL', 'DLX6', 'PABPC3', 'ENAH', 'CTNNA2', 'PDHA1', 
	'CHD3', 'PPP2R5D', 'TBR1', 'DNM1', 'CCDC39', 'GRIN2B', 'SPTAN1', 'GXYLT1', 'PPARD', 'PCSK4', 'NCL', 'GALNT9', 'ZNF292']
good_dn_solved_genes = ['TUBA1A', 'DYNC1H1', 'KIF11', 'TUBG1', 'FGFR1', 'TUBB2B', 'TUBB2A', 'STXBP1', 'DDX3X', 'CASK', 'FRMD1', 'SHOC2', 'KCNQ2', 'PIK3CA', 'BRAF', 'CACNA1G', 
	'KIAA1211', 'CREBBP', 'ZEB2', 'ACTG1', 'ATP1A3', 'MTOR', 'SETD2', 'PDHA1', 'TBR1', 'DNM1', 'GRIN2B', 'CDC42', 'MACF1', 'DCST1', 'BCL11A', 'TUBB8', 'RARG', 'FOXP1', 'KCTD1']

all_rec_solved_genes = ['ASPM', 'TCTN3', 'GPSM2', 'PEX6', 'NDE1', 'CRADD', 'KATNAL1', 'CARS2', 'L1CAM', 'SLC5A5', 'WDR27', 'RBM28', 'FZD3', 'VEZF1', 'ZNF77', 'GLRA2', 'ASNS', 
		'QARS', 'RARS2', 'RPGRIP1L', 'C5orf42', 'RTTN', 'FIG4', 'VRK1', 'TOE1', 'KATNB1', 'CENPE', 'ANKLE2', 'SYNE1', 'PCLO', 'CEP290', 'RELN', 'MACF1', 'GLIS3', 'KIAA1033', 
		'PCDH1', 'SMPD4', 'C2CD3', 'DDX31', 'WDR91', 'DCLRE1C', 'PMM2', 'PTCHD2', 'ATP13A5', 'MKI67', 'AP1S3', 'PCSK6', 'MMP1', 'PIGQ', 'PI4KA', 'CIT', 'ANKLE1', 'CDK5RAP2', 
		'KRT38', 'PNPLA8', 'PUS3', 'ARMC9', 'PIBF1']
good_rec_solved_genes = ['ASPM', 'TCTN3', 'GPSM2', 'PEX6', 'NDE1', 'CRADD', 'KATNAL1', 'CARS2', 'L1CAM', 'ASNS', 'QARS', 'RARS2', 'RPGRIP1L', 'C5orf42', 'RTTN', 'FIG4', 'VRK1', 
		'TOE1', 'KATNB1', 'CENPE', 'ANKLE2', 'SYNE1', 'PCLO', 'CEP290', 'RELN', 'MACF1']

all_xlinked_solved_genes = ['ATRX', 'FLNA', 'KDM6A', 'SLC9A6', 'HDAC6', 'NHS', 'KIAA2022', 'AP1S2', 'P2RY4', 'KIF4A', 'TMLHE']
good_xlinked_solved_genes = ['ATRX', 'FLNA', 'KDM6A', 'SLC9A6', 'HDAC6']

meg_gene_list = [ 'AKT1', 'AKT3', 'CCND2', 'GNAQ', 'GNAS', 'MTOR', 'PIK3CA', 'PIK3R2', 'NPRL2', 'NPRL3', 'STRADA', 'TBC1D7', 'TSC1', 'TSC2', 'PTEN', 'HEPACAM', 'NFIA', 'NFIX', 'NSD1', 
		'DEPDC5', 'CDKN1C', 'CUL4B', 'DNMT3A', 'EED', 'EZH2', 'GLI3', 'GPC3', 'KCNJ8', 'MED12', 'MLC1', 'PTCH1', 'RAB39B', 'RIN2', 'RNF135', 'SETD2']

##de novo vars
# get_denovo_vars_and_filter(ped_info_file, project_name, all_dn_solved_genes, good_dn_solved_genes, genes_from_genelists_dict)

##recessive and comp het var
# get_recessive_vars_and_filter(ped_info_file, project_name, all_rec_solved_genes, good_rec_solved_genes, genes_from_genelists_dict)

##xlinked vars
# get_xlinked_vars_and_filter(ped_info_file, project_name, all_xlinked_solved_genes, good_xlinked_solved_genes, genes_from_genelists_dict)

##mosaic meg vars
# get_mosaic_vars_and_filter(ped_info_file, project_name, all_dn_solved_genes, good_dn_solved_genes, genes_from_genelists_dict, meg_gene_list)

##get singletons/duos and filter
get_single_duo_vars_and_filter(ped_info_file, project_name)


