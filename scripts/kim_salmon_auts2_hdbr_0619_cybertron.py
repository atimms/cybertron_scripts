#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil
import numpy

##parameters
delim = '\t'
threads = '10'

##programs
salmon = '/home/atimms/programs/salmon-latest_linux_x86_64/bin/salmon'

##ref files etc
##standard version
# salmon_index = '/home/atimms/ngs_data/references/salmon/human/human_ens_index'
# salmon_human_ens_fa = '/home/atimms/ngs_data/references/salmon/human/gentrome.fa'
salmon_index = '/home/atimms/ngs_data/references/salmon/human/human_ens_JQ670867_index'
salmon_human_ens_fa = '/home/atimms/ngs_data/references/salmon/human/gentrome_JQ670867.fa'
salmon_human_ens_decoy = '/home/atimms/ngs_data/references/salmon/human/decoys.txt'
# ens_gtf = '/home/atimms/ngs_data/references/ensembl/hg19/Homo_sapiens.GRCh37.75.gtf'
ens_gtf = '/home/atimms/ngs_data/references/ensembl/hg38/Homo_sapiens.GRCh38.99.gtf'

auts2_nms = ['ENST00000342771', 'ENST00000406775', 'ENST00000403018', 'ENST00000611706', 'ENST00000615871', 
		'ENST00000443672', 'ENST00000449547', 'ENST00000439256', 'ENST00000476695', 'ENST00000416482', 'ENST00000418686', 
		'ENST00000483297', 'ENST00000475660', 'ENST00000489774', 'ENST00000464768', 'ENST00000481994', 'ENST00000498384', 
		'ENST00000465899', 'JQ670867']

auts2_nms = ['ENST00000403625', 'MAST4_alt']

def make_salmon_index(fa_file, decoy_file, index):
	#./bin/salmon index -t transcripts.fa -i transcripts_index -decoys decoys.txt -k 31
	st_n_sort = subprocess.Popen([salmon, 'index', '-t', fa_file, '-i', index, '-d', decoy_file, '-k', '31', '-p', threads])
	st_n_sort.wait()

def make_salmon_index_no_decoy(fa_file, index):
	#./bin/salmon index -t transcripts.fa -i transcripts_index -decoys decoys.txt -k 31
	st_n_sort = subprocess.Popen([salmon, 'index', '-t', fa_file, '-i', index, '-k', '31', '-p', threads])
	st_n_sort.wait()

def salmon_quant(sample_id, fq_files, index):
	if len(fq_files) == 2:
		#./bin/salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq --validateMappings -o transcripts_quant
		st_n_sort = subprocess.Popen([salmon, 'quant', '-i', index, '-l', 'A', '-1', fq_files[0], '-2', fq_files[1], 
				'--validateMappings', '-o', sample_id + '.trans_quant' ,  '-p', threads])
		st_n_sort.wait()
	elif len(fq_files) == 1:
		#./bin/salmon quant -i transcripts_index -l <LIBTYPE> -r reads.fq --validateMappings -o transcripts_quant
		print('single ended not setup yet')

def get_tpms(infile):
	with open(infile, "r") as infh:
		nm_list, tpm_list = [], []
		for line in infh:
			line = line.rstrip().split(delim)
			nm = line[0]
			tpm = line[3]
			if nm in auts2_nms:
				nm_list.append(nm)
				tpm_list.append(tpm)
	return(nm_list, tpm_list)

def call_salmon_on_samples_file_and_combine(work_dir, info_file, salmon_index_file, outfile):
	os.chdir(work_dir)
	##get info on project
	with open(info_file, "r") as infh:
		sample_dict = {}
		line_count = 0
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print conditions
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				sample_dict[sample] = condition_by_sample
				fq1 = line[1].split(',')
				fq2 = line[2].split(',')
				##if there were multiple fqs then I combined previously
				if len(fq1) == 2:
					fqs = [sample + '.r1.fq.gz', sample + '.r2.fq.gz']
				else:
					fqs = [fq1[0], fq2[0]]
				print(sample, fqs, salmon_index_file)
				##run salmon on each sample
				salmon_quant(sample, fqs, salmon_index_file)
				nms, tpms = get_tpms(sample + '.trans_quant/quant.sf')
				sample_dict[sample].extend(tpms)
	with open(outfile, "w") as outfh:
		header = ['sample', 'developmental_stage', 'organism_part'] + nms
		outfh.write(delim.join(header) + '\n')
		for s in sample_dict:
			line_out = [s] + sample_dict[s]
			outfh.write(delim.join(line_out) + '\n')

def summarise_data(infile, t_out, ta_out):
	with open(infile, "r") as infh:
		t_dict, ta_dict = {}, {}
		lc = 0
		for line in infh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc == 1:
				t_header = [line[2], 'sample #'] + line[3:]
				ta_header =  [line[2], line[1], 'sample #'] + line[3:]
			else:
				tissue = line[2]
				tissue_age = line[2] + ',' + line[1]
				tpms = [float(i) for i in line[3:]]
				##populate tissue dict
				if tissue in t_dict:
					t_dict[tissue].append(tpms)
				else:
					t_dict[tissue] = [tpms]
				##and tissue/age dict
				if tissue_age in ta_dict:
					ta_dict[tissue_age].append(tpms)
				else:
					ta_dict[tissue_age] = [tpms]
	##write tissue summary
	with open(t_out, "w") as toutfh:
		toutfh.write(delim.join(t_header) + '\n')
		for t in t_dict:
			sample_no = len(t_dict[t])
			a = numpy.array(t_dict[t])
			tpm_means = numpy.mean(a, axis=0)
			tpm_means = [str(i) for i in list(tpm_means)]
			print(t,sample_no, t_dict[t])
			print(tpm_means)
			line_out = [t, str(sample_no)] + tpm_means
			print(line_out)
			toutfh.write(delim.join(line_out) + '\n')			
	##write tissue/age summary
	with open(ta_out, "w") as taoutfh:
		taoutfh.write(delim.join(ta_header) + '\n')
		for ta in ta_dict:
			sample_no = len(ta_dict[ta])
			a = numpy.array(ta_dict[ta])
			tpm_means = numpy.mean(a, axis=0)
			tpm_means = [str(i) for i in list(tpm_means)]
			print(ta,sample_no, ta_dict[ta])
			print(tpm_means)
			ta_split = ta.split(',')
			line_out = ta_split + [str(sample_no)] + tpm_means
			print(line_out)
			taoutfh.write(delim.join(line_out) + '\n')	


def make_tpm_matrix(infile_suffix, outfile):
	infiles = glob.glob('*' + infile_suffix)
	tpm_dict = {}
	sample_list = []

	for infile in infiles:
		sample = infile.split('.')[0]
		sample_list.append(sample)
		with open(infile, "r") as infh:
			for line in infh:
				line = line.rstrip().split(delim)
				nm = line[0]
				tpm = line[3]
				if nm in tpm_dict:
					tpm_dict[nm].append(tpm)
				else:
					tpm_dict[nm] = [tpm]
	with open(outfile, "w") as outfh:
		outfh.write(delim.join(['transcript'] + sample_list) + '\n')
		for n in tpm_dict:
			outfh.write(delim.join([n] + tpm_dict[n]) + '\n')

def add_gene_name_tpms(infile, outfile, ensemble_gtf):
	ens_dict = {}
	with open(ensemble_gtf, "r") as ens_fh:
		for line in ens_fh:
			if line[0] != '#':
				info = line.split(delim)[8].split(';')
				# print(info)
				trans_id = info[2].split('"')[1]
				if 'gene_name' in info[5]:
					gene_id = info[5].split('"')[1]
					# print(trans_id, gene_id)
					ens_dict[trans_id] = gene_id




	with open(infile, "r") as infh, open(outfile, "w") as outfh:
		lc = 0
		for line in infh:
			line = line.split(delim)
			lc += 1
			if lc == 1:
				outfh.write(delim.join(['genename'] + line))
			else:
				transcript = line[0]
				if transcript in ens_dict:
					genename = ens_dict[transcript]
				else:
					genename = 'na'
				outfh.write(delim.join([genename] + line))




##run methods
data_dir = '/home/atimms/ngs_data/rnaseq/kim_hdbr_rnaseq_0619'
sample_file = 'hdbr_all.txt'
##test version
# sample_file = 'hdbr_test.txt'
results_file = 'hdbr_auts2_transcripts.by_sample.xls'
sum_tissue_file = 'hdbr_auts2_transcripts.by_tissue.xls'
sum_tissue_age_file = 'hdbr_auts2_transcripts.by_tissue_age.xls'

##make index for human ens transcript 
# make_salmon_index(salmon_human_ens_fa, salmon_human_ens_decoy, salmon_index)
##test
# sample = 'HDBR872'
# fqs = ['ERR1473601_1.fastq.gz', 'ERR1473601_2.fastq.gz']
# salmon_quant(data_dir, sample, fqs, salmon_index)
##process text file to run salmon and collate results
# call_salmon_on_samples_file_and_combine(data_dir, sample_file, salmon_index, results_file)
##sumrrize data
# summarise_data(results_file, sum_tissue_file, sum_tissue_age_file)


##repeat add auts2 tanscript
##make index for human ens transcript 
# make_salmon_index_no_decoy(salmon_human_ens_fa, salmon_index)
##process text file to run salmon and collate results
# call_salmon_on_samples_file_and_combine(data_dir, sample_file, salmon_index, results_file)
##sumarize data
# summarise_data(results_file, sum_tissue_file, sum_tissue_age_file)

##get all tpms
# make_tpm_matrix('.trans_quant/quant.sf', 'hdbr_tpms_1219.xls')
##then add genenames
# add_gene_name_tpms('hdbr_tpms_1219.xls', 'hdbr_tpms_genename_0420.xls', ens_gtf)


##repeat again for MAST4 alt transcript
##manuallt add scotts transcript
salmon_index = '/home/atimms/ngs_data/references/salmon/human/gentrome_mast4_index'
salmon_human_ens_fa = '/home/atimms/ngs_data/references/salmon/human/gentrome_mast4.fa'
results_file = 'hdbr_mast4_transcripts.by_sample.xls'
##make index for human ens transcript 
make_salmon_index_no_decoy(salmon_human_ens_fa, salmon_index)
##process text file to run salmon and collate results
call_salmon_on_samples_file_and_combine(data_dir, sample_file, salmon_index, results_file)

