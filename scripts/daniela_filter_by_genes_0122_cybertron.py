#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
requires:
qsub -Iq cdbrmq -l mem=30gb,ncpus=2 -P 9aa67182-09d6-405c-b7f7-93b0320828c1
qsub -Iq cdbrmq -l mem=100gb,ncpus=10 -P 9aa67182-09d6-405c-b7f7-93b0320828c1
or use daniela_genomes.pbs i.e.
qsub -q cdbrmq daniela_genomes.pbs
'''

##set input variables and parameters
delim = '\t'


def filter_by_gene_lists(results_files, out_prefix, genelist_files):
	for genelist_file in genelist_files:
		gl_name = genelist_file.split('.')[0]
		with open(genelist_file) as file:
			genes = file.readlines()
			genes = [gene.rstrip() for gene in genes]
		# print(gl_name, genes)
		rfc = 0
		outfile = out_prefix + '.' + gl_name + '.xls'
		with open(outfile, 'w') as out_fh:
			for results_file in results_files:
				batch = results_file.split('.')[0]
				rfc += 1
				lc = 0
				with open(results_file, 'r') as in_fh:
					for line in in_fh:
						line = line.split(delim)
						lc += 1
						if lc == 1:
							if rfc == 1:
								out_fh.write(delim.join(line[:3] + ['batch'] + line[3:]))	
						else:
							gene = line[8]
							if gene in genes:
								out_fh.write(delim.join(line[:3] + [batch] + line[3:]))




##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/daniela_analysis_0122'
os.chdir(working_dir)


slivar_files_multiplex = ['exomes_0221.multiplex.std_analysis.xls', 'genome_0221.multiplex.std_analysis.xls']
slivar_files_single_duo = ['exomes_0221.single_duos.std_analysis.xls', 'exomes_1221.single_duos.std_analysis.xls', 
		'genome_0219.single_duos.std_analysis.xls', 'genome_0221.single_duos.std_analysis.xls', 'genome_0620.single_duos.std_analysis.xls']
slivar_files_trio = ['exomes_0221.trios.std_analysis.xls', 'exomes_1221.trios.std_analysis.xls', 
		'genome_0219.trios.std_analysis.xls', 'genome_0221.trios.std_analysis.xls', 'genome_0620.trios.std_analysis.xls']
gene_lists = ['human_syndrome.txt', 'jacs_absent_pinna.txt', 'jacs_small_pinna.txt', 
		'SF3_family_members.txt', 'acmg_secondary_genes.txt']

##filter by the gene lists
filter_by_gene_lists(slivar_files_trio, 'combined_' + slivar_files_trio[0].split('.')[1], gene_lists)
filter_by_gene_lists(slivar_files_multiplex, 'combined_' + slivar_files_multiplex[0].split('.')[1], gene_lists)
filter_by_gene_lists(slivar_files_single_duo, 'combined_' + slivar_files_single_duo[0].split('.')[1], gene_lists)

