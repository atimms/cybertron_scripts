#!/usr/bin/env python
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

'''
info....

##scrublet
##go to worker node
qsub -Iq cdbrmq -l mem=100gb,ncpus=5,walltime=10:00:00 -P 19833a08-f6fb-4bea-8526-8a79069da87

##to install scrublet
module load local_python/3.7.6
conda create --name scrublet
source activate scrublet
pip install scrublet

##to load
module load local_python/3.7.6
source activate scrublet

'''

def run_scrublet_rna(input_dir):
	counts_matrix = scipy.io.mmread(input_dir + 'matrix.mtx').T.tocsc()
	genes = np.array(scr.load_genes(input_dir + 'features.tsv', delimiter='\t', column=1))
	print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
	print('Number of genes in gene list: {}'.format(len(genes)))
	scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
	doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
	np.savetxt(input_dir + 'predicted_doublet_mask.txt', scrub.predicted_doublets_, fmt='%s')
	np.savetxt(input_dir + 'doublet_scores.txt', scrub.doublet_scores_obs_, fmt='%.4f')

def run_scrublet_atac(input_dir):
	counts_matrix = scipy.io.mmread(input_dir + 'matrix.mtx').T.tocsc()
	print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
	scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
	doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
	np.savetxt(input_dir + 'predicted_doublet_mask.txt', scrub.predicted_doublets_, fmt='%s')
	np.savetxt(input_dir + 'doublet_scores.txt', scrub.doublet_scores_obs_, fmt='%.4f')


##parameters
delim = '\t'


##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_scrublet_test_0820'
os.chdir(working_dir)

##data dirs 
hu5_rna_dir = '/home/atimms/ngs_data/misc/cherry_scrublet_test_0820/Hu5_scRNA/'
hu7_rna_dir = '/home/atimms/ngs_data/misc/cherry_scrublet_test_0820/Hu7_scRNA/'
hu37_rna_dir = '/home/atimms/ngs_data/misc/cherry_scrublet_test_0820/Hu37_scRNA/'
rna_data_dirs = [hu5_rna_dir, hu7_rna_dir, hu37_rna_dir]
hu5_atac_dir = '/home/atimms/ngs_data/misc/cherry_scrublet_test_0820/Hu5_scATAC/'
hu7_atac_dir = '/home/atimms/ngs_data/misc/cherry_scrublet_test_0820/Hu7_scATAC/'
hu8_atac_dir = '/home/atimms/ngs_data/misc/cherry_scrublet_test_0820/Hu8_scATAC/'
atac_data_dirs = [hu5_atac_dir, hu7_atac_dir, hu8_atac_dir]


##rna data
# for data_dir in rna_data_dirs:
# 	run_scrublet_rna(data_dir)
##atac data
for data_dir in atac_data_dirs:
	run_scrublet_atac(data_dir)