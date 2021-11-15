#!/usr/bin/python
import os
import subprocess
import glob
import shutil
import math

##set input variables and parameters
delim = '\t'

'''
methods etc from:
cunn_rnaseq_0621_heatmaps_cyb.py

and graphing using
cunn_flexcell_rnaseq_0321_heatmaps.R

'''
##methods

def get_genenames_from_file(infile):
	genes = []
	with open(infile, "r") as in_fh:
		for line in in_fh:
			gene = line.rstrip()
			genes.append(gene)
	return(genes)

def filter_ratio_file_genelist(ratio_file, genelists, out_prefix):
	for genelist in genelists:
		gl_genes = get_genenames_from_file(genelist)
		out_file = ratio_file.split('.', 1)[0] + '.' + genelist.split('.', 1)[0] + out_prefix
		with open(ratio_file, "r") as in_fh, open(out_file, "w") as out_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				if lc == 1:
					out_fh.write(line)
				else:
					gene = line.split(delim)[0]
					if gene in gl_genes:
						out_fh.write(line)

def filter_order_ratio_files_by_sample(ratio_files, samplelist_dict):
	for ratio_file in ratio_files:
		for samplelist in samplelist_dict:
			samples_wanted = samplelist_dict[samplelist]
			out_file = ratio_file.rsplit('.', 2)[0] + '.' + samplelist + '.heatmap.txt'
			with open(ratio_file, "r") as in_fh, open(out_file, "w") as out_fh:
				out_fh.write(delim.join(['gene'] + samples_wanted) + '\n')
				lc = 0
				for line in in_fh:
					line = line.rstrip().split(delim)
					lc += 1
					if lc == 1:
						header = line
					else:
						gene = line[0]
						ratios = []
						for sample_wanted in samples_wanted:
							ratio = line[header.index(sample_wanted)]
							ratios.append(ratio)
						out_fh.write(delim.join([gene] + ratios) + '\n')


##run methods
working_dir = '/home/atimms/ngs_data/rnaseq/deema_heatmaps_1021'
os.chdir(working_dir)

##flexcell data
norm_count_ratio_file = 'cunn_flexcell_rnaseq_0321_ctl_flna.norm_counts.ratios.txt'
genelist_files = ['Bone_Morphogenic_Protein_Signalling_and_Regulation.txt', 'CELL_CELL_ADHESION_MEDIATED_BY_INTEGRIN.txt', 'Chondrocyte_development.txt', 'Chondrocyte_differentiation.txt', 'Chondroyte_proliferation.txt', 'focal_adhesion.txt', 'NEGATIVE_REGULATION_OF_BMP_SIGNALING_PATHWAY.txt', 'Osteoblast_Development.txt', 'Osteoblast_Differentiation.txt', 'Osteoblast_Proliferation.txt', 'POSITIVE_REGULATION_OF_BMP_SIGNALING_PATHWAY.txt', 'POSITIVE_REGULATION_OF_BONE_MINERALIZATION.txt', 'POSITIVE_REGULATION_OF_CELL_ADHESION_MEDIATED_BY_INTEGRIN.txt', 'POSITIVE_REGULATION_OF_CELL_ADHESION.txt', 'Reactome_RUNX2_regulates_osteoblast_diffferentiation.txt', 'REGULATION_OF_BMP_SIGNALING_PATHWAY.txt', 'REGULATION_OF_CHONDROCYTE_DIFFERENTIATION.txt', 'Regulation_of_osteoblast_differentiation.txt', 'Regulation_of_osteoblast_proliferation.txt', 'TGF-beta_Signaling_Pathway.txt']
##filter by genelist
# filter_ratio_file_genelist(norm_count_ratio_file, genelist_files, '.heatmap.txt')

##matrix data
norm_count_ratio_file = 'matrix_rnaseq_0220.norm_counts.ratios.txt'
genelist_files = ['Bone_Morphogenic_Protein_Signalling_and_Regulation.txt', 'CELL_CELL_ADHESION_MEDIATED_BY_INTEGRIN.txt', 'Chondrocyte_development.txt', 'Chondrocyte_differentiation.txt', 'Chondroyte_proliferation.txt', 'focal_adhesion.txt', 'NEGATIVE_REGULATION_OF_BMP_SIGNALING_PATHWAY.txt', 'Osteoblast_Development.txt', 'Osteoblast_Differentiation.txt', 'Osteoblast_Proliferation.txt', 'POSITIVE_REGULATION_OF_BMP_SIGNALING_PATHWAY.txt', 'POSITIVE_REGULATION_OF_BONE_MINERALIZATION.txt', 'POSITIVE_REGULATION_OF_CELL_ADHESION_MEDIATED_BY_INTEGRIN.txt', 'POSITIVE_REGULATION_OF_CELL_ADHESION.txt', 'Reactome_RUNX2_regulates_osteoblast_diffferentiation.txt', 'REGULATION_OF_BMP_SIGNALING_PATHWAY.txt', 'REGULATION_OF_CHONDROCYTE_DIFFERENTIATION.txt', 'Regulation_of_osteoblast_differentiation.txt', 'Regulation_of_osteoblast_proliferation.txt', 'TGF-beta_Signaling_Pathway.txt']
genelist_ratio_files = ['matrix_rnaseq_0220.Bone_Morphogenic_Protein_Signalling_and_Regulation.all_samples.txt', 'matrix_rnaseq_0220.CELL_CELL_ADHESION_MEDIATED_BY_INTEGRIN.all_samples.txt', 'matrix_rnaseq_0220.Chondrocyte_development.all_samples.txt', 'matrix_rnaseq_0220.Chondrocyte_differentiation.all_samples.txt', 'matrix_rnaseq_0220.Chondroyte_proliferation.all_samples.txt', 'matrix_rnaseq_0220.focal_adhesion.all_samples.txt', 'matrix_rnaseq_0220.NEGATIVE_REGULATION_OF_BMP_SIGNALING_PATHWAY.all_samples.txt', 'matrix_rnaseq_0220.Osteoblast_Development.all_samples.txt', 'matrix_rnaseq_0220.Osteoblast_Differentiation.all_samples.txt', 'matrix_rnaseq_0220.Osteoblast_Proliferation.all_samples.txt', 'matrix_rnaseq_0220.POSITIVE_REGULATION_OF_BMP_SIGNALING_PATHWAY.all_samples.txt', 'matrix_rnaseq_0220.POSITIVE_REGULATION_OF_BONE_MINERALIZATION.all_samples.txt', 'matrix_rnaseq_0220.POSITIVE_REGULATION_OF_CELL_ADHESION.all_samples.txt', 'matrix_rnaseq_0220.POSITIVE_REGULATION_OF_CELL_ADHESION_MEDIATED_BY_INTEGRIN.all_samples.txt', 'matrix_rnaseq_0220.Reactome_RUNX2_regulates_osteoblast_diffferentiation.all_samples.txt', 'matrix_rnaseq_0220.REGULATION_OF_BMP_SIGNALING_PATHWAY.all_samples.txt', 'matrix_rnaseq_0220.REGULATION_OF_CHONDROCYTE_DIFFERENTIATION.all_samples.txt', 'matrix_rnaseq_0220.Regulation_of_osteoblast_differentiation.all_samples.txt', 'matrix_rnaseq_0220.Regulation_of_osteoblast_proliferation.all_samples.txt', 'matrix_rnaseq_0220.TGF-beta_Signaling_Pathway.all_samples.txt']
sample_dict = {'flna_ctl': ['S_AUT25_1', 'S_AUT25_2', 'S_AUT25_3', 'S_OST85_1', 'S_OST85_2', 'S_OST85_3', 'S_STL19_1', 'S_STL19_2', 'S_STL19_3', 'S_C1625_1', 'S_C1625_2', 'S_C1625_3', 'S_C2084_1', 'S_C2084_2', 'S_C2084_3', 'S_C3049_1', 'S_C3049_2', 'S_C3049_3']}

#filter by genelist
filter_ratio_file_genelist(norm_count_ratio_file, genelist_files, '.all_samples.txt')
##filter/order by sample
filter_order_ratio_files_by_sample(genelist_ratio_files, sample_dict)


