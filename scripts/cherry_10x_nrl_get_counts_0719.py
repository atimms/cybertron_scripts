#!/usr/bin/python
import os
import subprocess
import glob
import shutil
import numpy as np

'''
module load local_python/3.6.4 
'''

##set input variables and parameters
delim = '\t'


def transpose_file(infile, outfile):
	with open(infile, 'r') as file, open(outfile, 'w') as final:
		#make list for each line
		a = [x.split(',') for x in file]
		#print a
		#turn into array
		b = np.array(a)
		#print b
		#transpose array
		c = b.T
		#convert array back to list of list
		d = c.tolist()
		#print d
		#convert each list back to 
		for e in d:
			line = delim.join(e) + '\n'
			final.write(line)

def make_dict_from_matrix(matrix):
	counts_dict = {}
	with open(matrix, 'r') as m_fh:
		line_count = 0
		for line in m_fh:
			line_count += 1
			if line_count == 1:
				header = line
			else:
				line = line.split(',')
				ens_id = line[0]
				counts = line[1:]
				if ens_id not in counts_dict:
					counts_dict[ens_id] = counts
				else:
					print('ens_id seen multiple times', ens_id)
	return(header, counts_dict)

def make_dict_from_features(infile):
	id_dict = {}
	with open(infile, 'r') as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			ens_id = line[0]
			gene_id = line[1].upper()
			if gene_id in id_dict:
				id_dict[gene_id].append(ens_id)
			else:
				id_dict[gene_id] = [ens_id]
	return(id_dict)


def get_gene_counts(genes, feature_file, matrix):
	header, counts_dict = make_dict_from_matrix(matrix)
	ens_gene_dict = make_dict_from_features(feature_file)
	print(header)
	for gene in genes:
		ens_ids = ens_gene_dict[gene]
		outfile = matrix.split('.')[0] + '.' + gene + '.csv'
		trans_outfile = matrix.split('.')[0] + '.' + gene + '.transposed.csv'
		with open(outfile, 'w') as out_fh:
			out_fh.write(header)
			for ens_id in ens_ids:
				line_out = ','.join([ens_id] + counts_dict[ens_id])
				out_fh.write(line_out)
		transpose_file(outfile, trans_outfile)









##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_10x_nrl_0319/10x_nrl_combined/outs'
genes_wanted = ['NR2E3', 'ESRRB', 'MEF2C', 'KDM5B', 'NONO', 'ARR3', 'PDE6H']
features_file = '/home/atimms/ngs_data/misc/cherry_10x_nrl_0319/10x_nrl_combined/outs/filtered_feature_bc_matrix/features.tsv'
matrix_file = '10x_nrl_combined.csv'

os.chdir(working_dir)
get_gene_counts(genes_wanted, features_file, matrix_file)






