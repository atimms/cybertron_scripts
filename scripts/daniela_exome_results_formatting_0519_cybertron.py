#!/usr/bin/env python
import sys
import subprocess
import os
import glob



'''
##load modules required for analysis
'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/exomes/working/daniela_0519'
os.chdir(working_dir)

def combine_results_add_gene_family_count(infiles, outfile, outfile_nor):
	print 'combining files:', infiles
	## get dict with ped counts for genenames
	gene_dict, gene_nor_dict = {}, {}
	for infile in infiles:
		with open(infile, "r") as in_fh:
			line_count = 0
			for line in in_fh:
				line = line.rstrip().split(delim)
				gene = line[5]
				ped = line[36]
				gene_type = line[8]
				rmsk = line[33]
				segdup = line[34]
				if gene_type == 'protein_coding' and rmsk == 'None' and segdup == '0':
					if gene in gene_nor_dict:
						if ped not in gene_nor_dict[gene]:
							gene_nor_dict[gene].append(ped)
					else:
						gene_nor_dict[gene] = [ped]
				if gene_type == 'protein_coding':
					if gene in gene_dict:
						if ped not in gene_dict[gene]:
							gene_dict[gene].append(ped)
					else:
						gene_dict[gene] = [ped]
					# print gene, ped, gene_type, rmsk, segdup, 'PASSSSSSED'

	print gene_dict
	for g in gene_dict:
		if len(gene_dict[g]) >1:
			print g, gene_dict[g]


	##combine all files, making sure col lengths are all the same
	with open(outfile, "w") as out_fh, open(outfile_nor, "w") as out_nor_fh:
		file_count, total_line_count = 0, 0
		for infile in infiles:
			# file_count += 1
			##add header from first file
			with open(infile, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line_count +=1
					# total_line_count += 1
					line = line.rstrip().split(delim)
					if line_count == 1:
						file_count += 1
						if file_count == 1:
							header = delim.join(line + ['peds_with_gene_count', 'peds_with_gene']) + '\n'
							out_fh.write(header)
							out_nor_fh.write(header)
					else:
						
						gene = line[5]
						gene_type = line[8]
						rmsk = line[33]
						segdup = line[34]
						if gene_type == 'protein_coding' and rmsk == 'None' and segdup == '0':
							line_out = delim.join(line + [str(len(gene_nor_dict[gene])), ','.join(gene_nor_dict[gene])]) + '\n'
							out_nor_fh.write(line_out)
						
						if gene_type == 'protein_coding':
							line_out = delim.join(line + [str(len(gene_dict[gene])), ','.join(gene_dict[gene])]) + '\n'
							out_fh.write(line_out)

def add_microtia_exome_genes(infiles, gene_dict, out_suffix):
	print 'quering files:', infiles
	for infile in infiles:
		outfile = infile.rsplit('.', 1)[0] + out_suffix 
		with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
			line_count = 0
			for line in in_fh:
				line = line.rstrip().split(delim)
				line_count += 1
				if line_count ==1:
					line_out = delim.join(line + ['microtia_exomes'] ) + '\n'
				else:
					gene = line[5]
					if gene in gene_dict:
						line_out = delim.join(line + [gene_dict[gene]] ) + '\n'
					else:
						line_out = delim.join(line + ['-'] ) + '\n'
				out_fh.write(line_out)				



##run methods
project = 'CFMCLP_exomes_0519'
std_results_files = glob.glob('*std_analysis.xls')
combined_std_file = project + '.standard_analysis.xls'
combined_std_norpts_file = project + '.standard_analysis.no_rpts.xls'
auto_dom_files = glob.glob('*auto_dominant_analysis.xls')
combined_ad_file = project + '.auto_dom_analysis.xls'
combined_ad_norpts_file = project + '.auto_dom_analysis.no_rpts.xls'
microtia_exome_genes = ['SLC25A5', 'ASMTL', 'DYNC1LI2', 'TMEM44', 'ABCA13', 'ABCA7', 'ABCD1', 'AEBP1', 'AKAP9', 'ANK2', 
		'ANKLE2', 'APC', 'APEX2', 'ARHGAP36', 'ARHGAP4', 'ARHGEF11', 'ARHGEF6', 'ARRB2', 'ASAH1', 'ATP10B', 'ATP11C', 'ATP13A5', 
		'ATP9B', 'BAI1', 'BICD1', 'BSX', 'BTN2A2', 'C7', 'CCDC138', 'CCDC19', 'CD109', 'CDH23', 'CDK10', 'CDK16', 'CDKL2', 'CEP170B', 
		'CEP350', 'CFTR', 'CLCN2', 'CMYA5', 'COG1', 'COL11A1', 'COL4A6', 'COQ10A', 'CORO7-PAM16', 'CPAMD8', 'CTSE', 'CYP1A2', 'CYP46A1', 
		'DCTN4', 'DHX35', 'DLC1', 'DLG5', 'DMD', 'DNAH11', 'DOCK11', 'DPH6', 'DSPP', 'DST', 'EFCAB4A', 'ENTPD2', 'EPHA4', 'ESF1', 'ESPL1', 
		'FAM13C', 'FAM69C', 'FAT4', 'FGGY', 'FLNA', 'FRG1', 'FZR1', 'GIPC1', 'GLG1', 'GNL3L', 'GPR112', 'GPR98', 'GRIA3', 'HEPH', 'HIVEP3', 
		'HK1', 'HMCN1', 'ITFG3', 'ITM2A', 'JADE3', 'KDM4B', 'KIAA1257', 'KIAA1522', 'KIAA2022', 'KIF4A', 'KLHL4', 'KYNU', 'LAMA1', 'LAMA2', 
		'LAS1L', 'LLGL1', 'LMO7', 'LRRC8D', 'LYN', 'MACF1', 'MEGF8', 'MFF', 'MOSPD2', 'MPP3', 'MTMR8', 'MTMR9', 'MUM1L1', 'MYO15A', 'MYO18A', 
		'MYO1C', 'NDST2', 'NHSL1', 'NOD1', 'NOP9', 'NOTCH1', 'NRAP', 'NYAP2', 'OPN1SW', 'PARP12', 'PCDH19', 'PCDHA6', 'PCDHGA11', 'PCM1', 
		'PFKP', 'PKHD1L1', 'PLCG2', 'PLXNA1', 'PLXNA3', 'PLXNA4', 'POLA1', 'PRIM1', 'PTMS', 'QPCTL', 'RGAG1', 'RLIM', 'RNASEL', 'RPL36A', 
		'RTEL1', 'SALL3', 'SATL1', 'SH2D1A', 'SHROOM4', 'SLC22A16', 'SLC2A10', 'SLC6A14', 'SPARC', 'SRGAP1', 'STAB2', 'STARD3', 'STRADA', 
		'SYAP1', 'SYNE1', 'TAAR6', 'THAP9', 'THEG', 'TMEM242', 'TMEM95', 'TPSD1', 'TRAPPC8', 'TRIM60', 'TRPC7', 'TULP4', 'UGGT2', 'URB1', 
		'USP31', 'UTP14A', 'VCPKMT', 'VSIG4', 'WAS', 'WNK3', 'ZFHX3', 'ZMYM4', 'ZNF521']
microtia_exome_type = ['xlinked_recessive, xlinked_de_novo,comp_het', 'xlinked_recessive,comp_het', 'de_novo,comp_het', 'homozygous,comp_het', 
		'de_novo', 'comp_het', 'xlinked_de_novo', 'comp_het', 'comp_het', 'comp_het', 'comp_het', 'homozygous', 'xlinked_recessive', 
		'xlinked_recessive', 'xlinked_recessive', 'comp_het', 'xlinked_recessive', 'de_novo', 'de_novo', 'comp_het', 'xlinked_recessive', 
		'comp_het', 'comp_het', 'comp_het', 'comp_het', 'de_novo', 'comp_het', 'comp_het', 'comp_het', 'comp_het', 'comp_het', 'de_novo', 
		'homozygous', 'xlinked_recessive', 'homozygous', 'de_novo', 'comp_het', 'comp_het', 'comp_het', 'comp_het', 'comp_het', 'comp_het', 
		'xlinked_recessive', 'homozygous', 'comp_het', 'comp_het', 'homozygous', 'homozygous', 'de_novo', 'de_novo', 'comp_het', 'homozygous', 
		'homozygous', 'xlinked_recessive', 'de_novo', 'xlinked_recessive', 'comp_het', 'homozygous', 'comp_het', 'homozygous', 'comp_het', 
		'de_novo', 'homozygous', 'homozygous', 'de_novo', 'homozygous', 'comp_het', 'homozygous', 'xlinked_recessive', 'de_novo', 'de_novo', 
		'de_novo', 'comp_het', 'xlinked_recessive', 'xlinked_recessive', 'comp_het', 'xlinked_recessive', 'xlinked_recessive', 'comp_het', 
		'de_novo', 'comp_het', 'homozygous', 'xlinked_recessive', 'xlinked_recessive', 'de_novo', 'de_novo', 'comp_het', 'xlinked_recessive', 
		'xlinked_recessive', 'xlinked_recessive', 'de_novo', 'comp_het', 'comp_het', 'xlinked_recessive', 'comp_het', 'comp_het', 'de_novo', 
		'homozygous', 'comp_het', 'comp_het', 'comp_het', 'xlinked_recessive', 'homozygous', 'xlinked_recessive', 'homozygous', 
		'xlinked_recessive', 'comp_het', 'homozygous', 'comp_het', 'de_novo', 'comp_het', 'comp_het', 'comp_het', 'comp_het', 'comp_het', 
		'homozygous', 'comp_het', 'de_novo', 'xlinked_recessive', 'homozygous', 'comp_het', 'comp_het', 'de_novo', 'comp_het', 'homozygous', 
		'comp_het', 'xlinked_recessive', 'comp_het', 'xlinked_recessive', 'homozygous', 'de_novo', 'de_novo', 'xlinked_recessive', 'xlinked_de_novo', 'comp_het', 'xlinked_recessive', 'homozygous', 'homozygous', 'xlinked_recessive', 'xlinked_recessive', 'xlinked_recessive', 'homozygous', 'homozygous', 'xlinked_recessive', 'homozygous', 'de_novo', 'comp_het', 'homozygous', 'comp_het', 'xlinked_recessive', 'comp_het', 'comp_het', 'homozygous', 'comp_het', 'homozygous', 'homozygous', 'comp_het', 'comp_het', 'de_novo', 'de_novo', 'comp_het', 'homozygous', 'comp_het', 'homozygous', 'xlinked_recessive', 'comp_het', 'xlinked_recessive', 'xlinked_recessive', 'xlinked_recessive', 'comp_het', 'comp_het', 'de_novo']
microtia_exome_dict = dict(zip(microtia_exome_genes, microtia_exome_type))
files_to_add_mic_exomes = std_results_files + auto_dom_files
microtia_out_suffix = '.microtia_exome_genes.xls'
##combine the 4 std analysis files and get gene counts i.e. how many peds have a var in thta ped
# combine_results_add_gene_family_count(std_results_files, combined_std_file, combined_std_norpts_file)

##manually make autodom files for peds 1 and 2 i.e. get the autdom analysis rows and remove analysis column
##combine the auto dom results and get gene counts
# combine_results_add_gene_family_count(auto_dom_files, combined_ad_file, combined_ad_norpts_file)
add_microtia_exome_genes(files_to_add_mic_exomes, microtia_exome_dict, microtia_out_suffix)


