#!/usr/bin/python
import os
import subprocess
import glob
import shutil

##set input variables and parameters
delim = '\t'


##methods
def get_all_de_genes(de_files, work_dir):
	os.chdir(work_dir)
	for de_file in de_files:
		print(de_file)
		outfile = de_file.rsplit('.',1)[0] + '.sig_genes.xls'
		with open(de_file, "r") as in_fh, open(outfile, "w") as out_fh:
			line_count = 0
			for line in in_fh:
				line_count += 1
				line = line.rstrip().split(',')
				if line_count ==1:
					out_fh.write(delim.join(line) + '\n')
				else:
					p_adj = line[6]
					if p_adj == 'NA':
						p_adj = '1'
					# print(line)
					if float(p_adj) <= 0.05:
						out_fh.write(delim.join(line) + '\n')

def get_shared_de_genes(work_dir, in_files, number_needed, up_outfile, down_outfile):
	os.chdir(work_dir)
	gene_dict = {}
	for in_file in in_files:
		sample = in_file.split('.')[1]
		with open(in_file, "r") as in_fh:
			line_count = 0
			for line in in_fh:
				line_count += 1
				if line_count >1:
					line = line.rstrip().split(',')
					gene = line[0].strip('"')
					logFC = line[2]
					p_adj = line[6]
					if p_adj == 'NA':
						p_adj = '1'
					# print(line)
					if float(p_adj) <= 0.05:
						##up genes i.e. overexpressed in first thing (+ve value)
						if float(logFC) > 0:
							if gene in gene_dict:
								gene_dict[gene][0].append(sample)
							else:
								gene_dict[gene] = [[sample], []]
						##sig down genes
						else:
							if gene in gene_dict:
								gene_dict[gene][1].append(sample)
							else:
								gene_dict[gene] = [[], [sample]]
	print(len(gene_dict))					
	with open(up_outfile, "w") as up_out_fh, open(down_outfile, "w") as down_out_fh:
		up_out_fh.write(delim.join(['gene', 'count', 'sample_names']) + '\n')
		down_out_fh.write(delim.join(['gene', 'count', 'sample_names']) + '\n')
		for g in gene_dict:
			if len(gene_dict[g][0]) >= number_needed:
				up_out_fh.write(delim.join([g, str(len(gene_dict[g][0])), '_'.join(gene_dict[g][0])]) + '\n')
			if len(gene_dict[g][1]) >= number_needed:
				down_out_fh.write(delim.join([g, str(len(gene_dict[g][1])), '_'.join(gene_dict[g][1])]) + '\n')





##run method

##get all DE genes
working_dir = '/archive/cunningham_m/strain_rnaseq_20/de_results_0320/all_results'
os.chdir(working_dir)
all_de_files = glob.glob('*.csv')
# print(len(working_dir))
get_all_de_genes(all_de_files, working_dir)


##get list of gene shared in all samples
working_dir = '/home/atimms/ngs_data/rnaseq/cunn_rnaseq_0220'
sample_gene_file = 'sample_gene.txt'

##in all samples
'''
for strain_type in ['strain0_vs_strain100', 'strain0_vs_strain15', 'strain15_vs_strain100']:
	infiles = glob.glob('*' + strain_type + '.csv')
	print(len(infiles))
	# print(infiles)
	number_req = 15
	outfile_up = strain_type + '.all_samples.up.xls'
	outfile_down = strain_type + '.all_samples.down.xls'
	get_shared_de_genes(working_dir, infiles, number_req, outfile_up, outfile_down)
'''

##in gene groups
gene_group_dict = {'AXL':['1017', 'C1653', 'C1856', '3007'],
		'FLNA':['C2084', 'C1625', 'C3049'], 'FLNB':['C3066', '1032', '4032'],
		'FLNC':['3035', 'C2038', '4025'], 'CTRL':['STL19', 'AUT25', 'OST85'],
		'PIEZO1':['C1671', 'C1859', 'C1926', 'C1957']}
##comparing strain types
'''
for strain_type in ['strain0_vs_strain100', 'strain0_vs_strain15', 'strain15_vs_strain100']:
	for gene in gene_group_dict:
		infiles = []
		for sample in gene_group_dict[gene]:
			infile = 'strain_rnaseq_0320.' + sample + '.' + strain_type + '.csv'
			infiles.append(infile)
		print(gene,len(infiles))
		# print(infiles)
		number_req = 2
		outfile_up = strain_type + '.' + gene + '.up.xls'
		outfile_down = strain_type + '.' + gene + '.down.xls'
		get_shared_de_genes(working_dir, infiles, number_req, outfile_up, outfile_down)
'''
##comparing case control
'''
for strain_type in ['strain0', 'strain15', 'strain100']:
	for gene in gene_group_dict:
		if gene != 'CTRL':
			infiles = []
			for sample in gene_group_dict[gene]:
				infile = 'strain_rnaseq_0320.' + sample + '_vs_control.' + strain_type + '.csv'
				infiles.append(infile)
			print(gene,len(infiles))
			# print(infiles)
			number_req = 2
			outfile_up = strain_type + '_vs_control.' + gene + '.up.xls'
			outfile_down = strain_type + '_vs_control.' + gene + '.down.xls'
			get_shared_de_genes(working_dir, infiles, number_req, outfile_up, outfile_down)
'''







