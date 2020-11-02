#!/usr/bin/env python
import subprocess
import os
import glob

'''
install/use mageck via conda:

install:
conda create --name mageck_env
conda activate mageck_env

from now on:
conda activate mageck_env

'''



##parameters
delim = '\t'

##methods
def get_fold_change_per_guide_and_gene(infile, guide_outfile, gene_outfile):
	with open(infile, "r") as in_fh, open(guide_outfile, "w") as guide_fh, open(gene_outfile, "w") as gene_fh:
		##make dict of guides in a gene
		gene_count_dict = {}
		lc = 0
		guide_fh.write(delim.join(['guide', 'gene', 'LH1 (static)', 'LH2 (top 10% sorted)', 'LH3 (bottom 25% sorted)', 'LH2 fold change', 'LH3 fold change']) + '\n')
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				guide_id = line[0]
				lh1 = line[2]
				lh2 = line[3]
				lh3 = line[4]
				gene = line[7]
				if lh1 != '0':
					if lh3 == '0':
						lh2_fc = 'inf'
					else:
						lh2_fc = float(lh2) / float(lh3)
					if lh2 == '0':
						lh3_fc = 'inf'
					else:			
						lh3_fc = float(lh3) / float(lh2)
					guide_fh.write(delim.join([guide_id, gene, lh1, lh2, lh3, str(lh2_fc), str(lh3_fc) ]) +'\n')
					if gene in gene_count_dict:
						gene_count_dict[gene][0].append(int(lh1))
						gene_count_dict[gene][1].append(int(lh1))
						gene_count_dict[gene][2].append(int(lh1))
					else:
						gene_count_dict[gene] = [[float(lh1)], [float(lh2)], [float(lh3)]]
		gene_fh.write(delim.join(['gene', 'LH1 mean', 'LH2 mean', 'LH3 mean', 'LH2 fold change', 'LH3 fold change']) + '\n')
		for g in gene_count_dict:
			sum_lh1 = sum(gene_count_dict[g][0])
			sum_lh2 = sum(gene_count_dict[g][1])
			sum_lh3 = sum(gene_count_dict[g][2])
			lh1_av = sum_lh1 / len(gene_count_dict[g][0])
			lh2_av = sum_lh2 / len(gene_count_dict[g][0])
			lh3_av = sum_lh3 / len(gene_count_dict[g][0])
			lh2_comb_fc = sum_lh2/sum_lh3
			lh3_comb_fc = sum_lh3/sum_lh2
			gene_fh.write(delim.join([g, str(lh1_av), str(lh2_av), str(lh3_av), str(lh2_comb_fc), str(lh3_comb_fc) ]) + '\n')

def run_mageck_test(infile, lh2_prefix, lh3_prefix):
	#mageck test -k sample.txt -t HL60.final,KBM7.final -c HL60.initial,KBM7.initial  -n demo
	run_megeck = subprocess.Popen(['mageck', 'test', '-k', infile, '-t', 'LH2', '-c', 'LH3', '-n', lh2_prefix])
	run_megeck.wait()
	run_megeck = subprocess.Popen(['mageck', 'test', '-k', infile, '-t', 'LH3', '-c', 'LH2', '-n', lh3_prefix])
	run_megeck.wait()

def run_mageck_mle(infile, matrixfile, out_prefix):
	#mageck mle -k leukemia.new.csv -d designmat.txt -n beta_leukemia
	run_megeck = subprocess.Popen(['mageck', 'mle', '-k', infile, '-d', matrixfile , '-n', out_prefix])
	run_megeck.wait()


##run methods
working_dir = '/home/atimms/ngs_data/misc/vishal_crispr_screen_1020'
os.chdir(working_dir)

##just fold change fifferences
data_file = 'vishal_crispr_data_s6612_1020.txt'
fc_per_guide = 'vishal_crispr_data_1020.fc_per_guide.txt'
fc_per_gene = 'vishal_crispr_data_1020.fc_per_gene.txt'
# get_fold_change_per_guide_and_gene(data_file, fc_per_guide, fc_per_gene)

##run mageck
counts_file = 'vishal_crispr_screen_counts.txt'
counts_noctls_file = 'vishal_crispr_screen_counts.no_ctls.txt'
lh2_mageck_test_prefix = 'vishal_crispr_1020.mageck_test_lh2'
lh3_mageck_test_prefix =  'vishal_crispr_1020.mageck_test_lh3'
design_matrix_file = 'vishal_crispr_screen_matrix.txt'
mle_results_prefix = 'vishal_crispr_1020.mageck_mle'

##run mageck test command
# run_mageck_test(counts_file, lh2_mageck_test_prefix, lh3_mageck_test_prefix)

##run mageck mle
run_mageck_mle(counts_noctls_file, design_matrix_file, mle_results_prefix)




