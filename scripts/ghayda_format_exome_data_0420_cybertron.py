#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'


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
		'CG46', 'avsnp147', 'CLINSIG', 'CLNDBN', 'CLNACC', 'CLNDSDB', 'CLNDSDBID', 'cosmic83_noncoding', 'cosmic83_coding', 'sample_ref/alt', 'parent1_ref/alt', 'parent2_ref/alt']

##methods
def change_file_names(work_dir, ped_dict, file_suffix):
	os.chdir(work_dir)
	for ped in ped_dict:
		infile = ped + file_suffix
		outfile = ped_dict[ped] + file_suffix
		shutil.copy(infile, outfile)

def make_genelist_dict(gene_file):
	genelist_dict = {}
	with open(gene_file, "r") as gf_fh:
		lc = 0
		for line in gf_fh:
			lc += 1
			if lc >1:
				line = line.strip('\n').split(delim)
				gene = line[0]
				group = line[1]
				if gene in genelist_dict:
					if group not in genelist_dict[gene]:
						genelist_dict[gene].append(group)
				else:
					genelist_dict[gene] = [group]
	return(genelist_dict)

def filter_singles_duos(info_file, outfile_prefix, gene_dict, single_duo_definitions, analysis_file_dir, dx_wanted):
	##outfile names
	clinvar_file = outfile_prefix + '.single_duos.clinvar.xls'
	candidates_file = outfile_prefix + '.single_duos.candidates.xls'
	##get all denovo vars in one file
	with open(info_file, "U") as in_fh, open(clinvar_file, "w") as clinvar_fh, open(candidates_file, "w") as candidates_fh:
		##add header
		header = ['ped', 'solved_status', 'solved_gene', 'ped_type', 'dxgroup1'] + std_anal_header
		clinvar_fh.write(delim.join(header) + '\n')
		candidates_fh.write(delim.join(header + ['gene list']) + '\n')
		lc, ped_count, var_count, genelist_count, clinvar_count = 0, 0, 0, 0, 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				# print(line)
				# line = line.replace("\r","")
				# line = line.replace("\n","")
				line = line.rstrip('\n').rstrip('\r').split(delim)
				ped_id = line[0]
				# print(ped_id, len(line))
				solved_status = line[34].rstrip()
				solved_gene = line[35].rstrip()
				ped_type = line[12].rstrip()
				dx_group = line[22].rstrip()
				
				##get singles/duos
				if ped_type in single_duo_definitions and dx_group in dx_wanted:
					print(ped_id, solved_status, solved_gene, ped_type, dx_group, ped_count)
					ped_count += 1
					std_file = analysis_file_dir + ped_id + '.std_analysis.xls'
					with open(std_file, "r") as std_fh:
						lc2 = 0
						for line in std_fh:
							lc2 += 1
							if lc2 >1:
								line = line.strip('\n').split(delim)
								anal = line[-2]
								cadd = line[11]
								gnomad_ex_af = line[33]
								gene = line[5]
								biotype = line[7]
								clinsig = line[18]
								print(cadd, gnomad_ex_af, clinsig)
								if cadd == 'None':
									cadd = '20'
								if gnomad_ex_af == 'na' or gnomad_ex_af == '.':
									gnomad_ex_af = '0'
								if float(cadd) > 15 and float(gnomad_ex_af) <= 0.0001 and biotype == 'protein_coding':
									var_count += 1
									if gene in gene_dict:
										line_out = [ped_id, solved_status, solved_gene, ped_type, dx_group] + line[:47] + [line[-2]] + ['_'.join(gene_dict[gene])]
										genelist_count += 1
										candidates_fh.write(delim.join(line_out) + '\n')
									if clinsig != '.':
										if 'enign' not in clinsig:
											line_out = [ped_id, solved_status, solved_gene, ped_type, dx_group] + line[:47] + [line[-2]]
											clinvar_count += 1
											clinvar_fh.write(delim.join(line_out) + '\n')
	print(ped_count, var_count, genelist_count, clinvar_count)

def filter_mosaic(info_file, outfile_prefix, gene_dict, analysis_file_dir, dx_wanted):
	##outfile names
	clinvar_file = outfile_prefix + '.pisces.clinvar.xls'
	candidates_file = outfile_prefix + '.pisces.candidates.xls'
	##ped with no file
	no_mosiac_peds = []
	##get all denovo vars in one file
	with open(info_file, "U") as in_fh, open(clinvar_file, "w") as clinvar_fh, open(candidates_file, "w") as candidates_fh:
		##add header
		header = ['ped', 'solved_status', 'solved_gene', 'ped_type', 'dxgroup1'] + mosaic_anal_header
		clinvar_fh.write(delim.join(header) + '\n')
		candidates_fh.write(delim.join(header + ['gene list']) + '\n')
		lc, ped_count, var_count, genelist_count, clinvar_count = 0, 0, 0, 0, 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				# print(line)
				line = line.rstrip('\n').rstrip('\r').split(delim)
				ped_id = line[0]
				# print(ped_id, len(line))
				solved_status = line[34].rstrip()
				solved_gene = line[35].rstrip()
				ped_type = line[12].rstrip()
				dx_group = line[22].rstrip()
				##get peds we're interested in
				if dx_group in dx_wanted:
					# print(ped_id, solved_status, solved_gene, ped_type, dx_group, ped_count)
					ped_count += 1
					pisces_files = glob.glob(analysis_file_dir + ped_id + '*')
					# print(ped_id, pisces_files)
					if len(pisces_files) == 0:
						no_mosiac_peds.append(ped_id)
					else:
						for pisces_file in pisces_files:
							with open(pisces_file, "r") as pis_fh:
								lc2 = 0
								for line in pis_fh:
									lc2 += 1
									# if lc2 ==1:
									# 	line = line.strip('\n').split(delim)
									# 	print(pisces_file, len(line))
									if lc2 >1:
										line = line.strip('\n').split(delim)
										cadd = line[30]
										PopFreqMax = line[46]
										gene = line[6]
										clinsig = line[66]

										# print(cadd, PopFreqMax, clinsig)
										if cadd == '.':
											cadd = '20'
										if PopFreqMax == '.':
											PopFreqMax = '0'
										if float(cadd) > 15 and float(PopFreqMax) <= 0.0001:
											var_count += 1
											if gene in gene_dict:
												line_out = [ped_id, solved_status, solved_gene, ped_type, dx_group] + line[:73] + ['_'.join(gene_dict[gene])]
												genelist_count += 1
												candidates_fh.write(delim.join(line_out) + '\n')
											if clinsig != '.':
												if 'athogenic' in clinsig or 'Uncertain' in clinsig:
													if len(line) == 74:
														extra = ['.', '.', '.']
													elif len(line) == 84:
														# print(line[77], line[80], line[81])
														extra = [line[77].split(':')[2], line[80].split(':')[1], line[81].split(':')[1] ]
													elif len(line) == 83:
														extra = [line[77].split(':')[2], line[80].split(':')[1], '.']
													elif len(line) == 81:
														extra = [line[77].split(':')[2], '.', '.']
													line_out = [ped_id, solved_status, solved_gene, ped_type, dx_group] + line[:73] + extra
													clinvar_count += 1
													clinvar_fh.write(delim.join(line_out) + '\n')
										
		print(ped_count, var_count, genelist_count, clinvar_count)
		print(no_mosiac_peds)


def get_var_info(info_file, var_file, std_analysis_file_dir, pisces_analysis_file_dir, dx_wanted):
	##get all denovo vars in one file
	with open(info_file, "U") as in_fh, open(var_file, "w") as out_fh:
		##add header
		header = ['ped', 'solved_status', 'solved_gene', 'inheritance', 'ped_type', 'dxgroup1', 'var_inheritence', 'variant_info', 'genomic_coordinates', 'VAF']
		out_fh.write(delim.join(header) + '\n')
		lc, ped_count, var_count = 0, 0, 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip('\n').rstrip('\r').split(delim)
				ped_id = line[0]
				# print(ped_id, len(line))
				solved_status = line[34].rstrip()
				solved_gene = line[35].rstrip()
				inheritance = line[37].rstrip()
				ped_type = line[12].rstrip()
				dx_group = line[22].rstrip()
				##get solved/candidates
				if solved_gene != '':
					##tidy the solved gene
					# solved_gene_formated = solved_gene.rstrip().replace(' ', '').replace('[', '').replace(']', '').split('/')
					ped_count += 1
					# print(ped_id, solved_status, solved_gene, inheritance, ped_type, dx_group)
					##check std analysis file
					std_file = std_analysis_file_dir + ped_id + '.std_analysis.xls'
					with open(std_file, "r") as std_fh:
						lc2 = 0
						for line in std_fh:
							lc2 += 1
							if lc2 >1:
								line = line.strip('\n').split(delim)
								gene = line[5]
								anal = line[-2]
								variant_info = line[13]
								genomic_coordinates = line[0] + ':g.[' + line[2] + line[3] + '>' + line[4] + ']'
								if gene in solved_gene:
									line_out = [ped_id, solved_status, solved_gene, inheritance, ped_type, dx_group, anal, variant_info, genomic_coordinates]
									out_fh.write(delim.join(line_out) + '\n')
									var_count += 1
					##check pisces analysis file
					pisces_files = glob.glob(pisces_analysis_file_dir + ped_id + '*')
					# print(ped_id, pisces_files)
					if len(pisces_files) > 0:
						for pisces_file in pisces_files:
							with open(pisces_file, "r") as pis_fh:
								lc3 = 0
								for line in pis_fh:
									lc3 += 1
									if lc3 >1:
										line = line.strip('\n').split(delim)
										gene = line[6]
										clinsig = line[66]
										variant_info = line[9]
										if variant_info == '.':
											variant_info = line[7]
										genomic_coordinates = 'chr' + line[0] + ':g.[' + line[2] + line[3] + '>' + line[4] + ']'
										if len(line) > 74:
											vaf = line[77].split(':')[2]
										if gene in solved_gene:
											print(pisces_file, len(line))
											if len(line) > 76:
												# print(pisces_file, len(line))
												vaf = line[77].split(':')[2]
											else:
												vaf = ''
											line_out = [ped_id, solved_status, solved_gene, inheritance, ped_type, dx_group, 'mosaic', variant_info, genomic_coordinates, vaf]
											out_fh.write(delim.join(line_out) + '\n')
											var_count += 1

	print(ped_count,var_count)


def get_single_var_ar_genes(infiles, outfile):
	with open(outfile, "w") as out_fh:
		fc = 0
		for infile in infiles:
			ped = infile.rsplit('/', 1)[1].split('.')[0]
			print(infile, ped)
			with open(infile, "r") as in_fh:
				lc = 0
				fc += 1
				for line in in_fh:
					lc += 1
					line = line.rstrip().split(delim)
					if lc == 1:
						if fc == 1:
							header = ['ped'] + line[:47]
							out_fh.write(delim.join(header) + '\n')
					else:
						gnomad_exome_af = line[33]
						gnomad_genome_af = line[41]
						biotype = line[7]
						clinsig = line[18]
						dbdb_inh = line[21]

						if gnomad_exome_af == 'na' or gnomad_exome_af == '.':
							gnomad_exome_af = '0'
						if gnomad_genome_af == 'na' or gnomad_genome_af == '.':
							gnomad_genome_af = '0'
						if float(gnomad_exome_af) <= 0.01 and float(gnomad_genome_af) <= 0.01 and biotype == 'protein_coding':
							if 'Pathogenic' in clinsig or 'Likely_pathogenic' in clinsig or 'AR' in dbdb_inh:
								line_out = [ped] + line[:47]
								out_fh.write(delim.join(line_out) + '\n')


##run methods
# working_dir = '/home/atimms/ngs_data/exomes/working/exome_project_20'
working_dir = '/home/atimms/ngs_data/exomes/working/exome_project_20/exome_project_0420'
os.chdir(working_dir)
# project_name = 'exomes_0420'
project_name = 'exomes_0420_updated'
ped_info_file = 'master_ped_tracking_41720.txt'
candidate_genes = 'genelists_0420.txt'
single_duo_types = ['Singleton', 'Sibship', 'Duo', 'Parent-sibs', 'Singleton*', 'Duo*']
std_anal_file_dir = '/home/atimms/ngs_data/exomes/working/exome_project_20/src_files/std_analysis/'
# pisces_file_dir = '/home/atimms/ngs_data/exomes/working/exome_project_20/src_files/all_pisces/'
pisces_file_dir = '/home/atimms/ngs_data/exomes/working/exome_project_20/exome_project_0420/src_files/all_pisces/'
ped_name_change_dict = {'DWM10':'LR16-461', 'DWM13':'LR16-462', 'DWM3':'LR16-463' }
dx_groups_wanted = ['MEG', 'MIC', 'MCD', 'MHM', 'DEVN', 'Other']
var_info_file = 'variant_info_0420.xls'

##step1. manually format and get files
## copy src_files.zip to cybertron
## use clean() function to copy data from master_ped_tracking_41720 to new excel file and then copy to text file
##change names of some files (been given lr numbers)
# change_file_names(std_anal_file_dir, ped_name_change_dict, '.std_analysis.xls')
# change_file_names(pisces_file_dir, ped_name_change_dict, '.pisces.xls')


##step2. filter singles/duos and pisces - rare in gnomad, cadd 15 and either in candidate gene or clivar 
##make dict from gene candidats
candidate_gene_dict = make_genelist_dict(candidate_genes)
##filter singles and duos
# filter_singles_duos(ped_info_file, project_name, candidate_gene_dict, single_duo_types, std_anal_file_dir, dx_groups_wanted)
##filter mosaic
filter_mosaic(ped_info_file, project_name, candidate_gene_dict, pisces_file_dir, dx_groups_wanted)

##step2b. get all AR single vars - use auto dom files and filter on clinvar/dbdb
# auto_dom_files = glob.glob('/home/atimms/ngs_data/exomes/working/exome_project_20/src_files/auto_dom/*xls')
# get_single_var_ar_genes(auto_dom_files, project_name + '.auto_dom.clinvar_dbdb.xls')



##step3. get var info including solved by other methods
# get_var_info(ped_info_file, var_info_file, std_anal_file_dir, pisces_file_dir, dx_groups_wanted)


