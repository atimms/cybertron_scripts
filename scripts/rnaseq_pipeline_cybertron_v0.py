#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'
star_threads = '16'

##programs
star = '/home/atimms/programs/STAR-2.5.3a/bin/Linux_x86_64_static/STAR'
# feature_counts = '/home/atimms/programs/subread-1.6.0-Linux-x86_64/bin/featureCounts'
# feature_counts = '/home/atimms/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts'
feature_counts = '/home/atimms/programs/subread-2.0.0-Linux-x86_64/bin/featureCounts'
picard = '/home/atimms/programs/picard.jar'
qorts = '/home/atimms/programs/hartleys-QoRTs-39cd1fc/QoRTs.jar'
##ref files etc
star_ref_dir = '/home/atimms/ngs_data/references/star/'
fa_file_path = ['/home/atimms/ngs_data/references/igenomes/', '/genome.fa']
gtf_file_path = ['/home/atimms/ngs_data/references/igenomes/', '/genes.gtf']
# gtf_file_path = ['/home/atimms/ngs_data/references/igenomes/', '/genes_temp.gtf']
ens_hg19_gtf_file = '/home/atimms/ngs_data/references/ensembl/hg19/Homo_sapiens.GRCh37.75.gtf'
deseq_r_template = ''
# star_bam_suffix = '.Aligned.out.bam'
star_bam_suffix = 'Aligned.sortedByCoord.out.bam'
sorted_bam_suffix = '.sorted.bam'

##methods

def combine_fqs(fqs_to_combine, fq_name):
	with open(fq_name, 'w') as cfq_fh:
		print "combining %s files named %s to make file %s"%(len(fqs_to_combine), fqs_to_combine, fq_name)
		cat_files = subprocess.Popen(['cat'] + fqs_to_combine, stdout=cfq_fh)
		cat_files.wait()

def star_align_paired_end(sample_name, fq_files, star_genome_dir):
	r1_fq = fq_files[0]
	r2_fq = fq_files[1]
	print sample_name, r1_fq, r2_fq, star_genome_dir
	# star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align.wait()

def star_align_single_end(sample_name, fq_file, star_genome_dir):
	print sample_name, fq_file
	# star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', fq_file, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', fq_file, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align.wait()

def feature_count_paired_end(bam_suffix, genome_gtf, outfile, samples):
	bam_files = []
	for s in samples:
		bam = s + '.' + bam_suffix
		bam_files.append(bam)
	print bam_files, len(bam_files)
	feat_count = subprocess.Popen([feature_counts,'-p', '-t', 'exon', '-g', 'gene_id', '-a', genome_gtf, '-o', outfile, '-T', '10', '-B'] + bam_files)
	feat_count.wait()

def feature_count_paired_end_allow_overlap(bam_suffix, genome_gtf, outfile, samples):
	bam_files = []
	for s in samples:
		bam = s + '.' + bam_suffix
		bam_files.append(bam)
	print bam_files, len(bam_files)
	feat_count = subprocess.Popen([feature_counts,'-p', '-t', 'exon', '-g', 'gene_id', '-O', '-a', genome_gtf, '-o', outfile, '-T', '10', '-B'] + bam_files)
	feat_count.wait()

def feature_count_paired_end_saf(bam_suffix, genome_gtf, outfile, samples):
	bam_files = []
	for s in samples:
		bam = s + '.' + bam_suffix
		bam_files.append(bam)
	print bam_files, len(bam_files)
	feat_count = subprocess.Popen([feature_counts,'-p', '-t', 'exon', '-F', 'SAF', '-a', genome_gtf, '-o', outfile, '-T', '10', '-B'] + bam_files)
	feat_count.wait()

def feature_count_single_end(bam_suffix, genome_gtf, outfile, samples):
	bam_files = []
	for s in samples:
		bam = s + '.' + bam_suffix
		bam_files.append(bam)
	print bam_files, len(bam_files)
	feat_count = subprocess.Popen([feature_counts,'-t', 'exon', '-g', 'gene_id', '-a', genome_gtf, '-o', outfile, '-T', '10'] + bam_files)
	feat_count.wait()

def feature_count_single_end_allow_overlap(bam_suffix, genome_gtf, outfile, samples):
	bam_files = []
	for s in samples:
		bam = s + '.' + bam_suffix
		bam_files.append(bam)
	print bam_files, len(bam_files)
	feat_count = subprocess.Popen([feature_counts,'-t', 'exon', '-g', 'gene_id', '-O', '-a', genome_gtf, '-o', outfile, '-T', '10'] + bam_files)
	feat_count.wait()

def format_feature_counts_file(infile, outfile):
	with open(outfile, "w") as outph, open(infile, "r") as inph:
		line_count = 0
		for line in inph:
			if line[0] != '#':
				line_count += 1
				line = line.strip('\n').split(delim)
				if line_count == 1:
					print line
					samples = line[6:]
					# print samples
					samples = [s.split('.')[0] for s in samples]
					# print samples
					# samples = [sample_dict_cond[s][0] for s in samples]
					print samples
					header = ['gene'] + samples + ['\n']
					outph.write(delim.join(header))
				else:
					gene = line[0]
					if gene != '':
						lineout = [gene] + line[6:] + ['\n']
						# print lineout
						outph.write(delim.join(lineout))


def make_metadata_file_format_samplenames_order_as_counts(sample_dict, outfile, counts_file, header):
	with open(outfile, "w") as outph, open(counts_file, "r") as cou_ph:
		outph.write(delim.join(header + ['\n']))
		line_count = 0
		for line in cou_ph:
			line_count += 1
			if line_count == 2:
				print line
				line = line.strip('\n').split(delim)
				samples = line[6:]
				samples = [s.split('.')[0] for s in samples]
				for sample in samples:
					out_info = sample_dict[sample]
					print sample, sample_dict[sample], out_info
					outph.write(delim.join([sample] + out_info + ['\n']))


def call_rnaseq_methods_for_deseq_analysis(work_dir, info_file, genome, covariates):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	print 'working on project:', project_name
	##get info on project
	with open(info_file, "U") as infh:
		metadata_dict = {}
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print conditions
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				metadata_dict[sample] = condition_by_sample
				sample_list.append(sample)
				fq1_info = line[1].split(',')
				fq2_info = line[2].split(',')
				# print fq1_info, len(fq1_info)
				# print fq2_info, len(fq2_info)
				##do we need to combine fqs
				if len(fq1_info) == 1:
					fq1 = fq1_info[0]
				elif len(fq1_info) >= 1:
					fq1 = sample + '.r1.fq.gz'
					combine_fqs(fq1_info, fq1)
				##single or paired sequencing
				if fq2_info == ['']:
					seq_type = 'single_end'
				else:
					seq_type = 'paired_end'
					##if same fqs for r1 and r2
					if len(fq1_info) == len(fq2_info):
						##combine files
						if len(fq2_info) == 1:
							fq2 = fq2_info[0]
						elif len(fq2_info) > 1:
							fq2 = sample + '.r2.fq.gz'
							combine_fqs(fq2_info, fq2)
					else:
						print 'issue with sample %s, there are different number of fq files for read 1 and 2'%sample
				##run star
				star_index_dir = star_ref_dir + genome
				print 'seq_type = ', seq_type
				if seq_type == 'single_end':
					star_align_single_end(sample, fq1, star_index_dir)
				elif seq_type == 'paired_end':
					star_align_paired_end(sample, [fq1, fq2], star_index_dir)

	##run feature count on bams
	gtf_file = gtf_file_path[0] + genome + gtf_file_path[1]
	feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
	deseq_feature_count_file = project_name + '.star_fc.counts.txt'
	if seq_type == 'single_end':
		feature_count_single_end(star_bam_suffix, gtf_file, feature_count_results_file, sample_list)
	elif seq_type == 'paired_end':
		feature_count_paired_end(star_bam_suffix, gtf_file, feature_count_results_file, sample_list)
	format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)
	##make metadata
	deseq_metadata_file = project_name + '.star_fc.metadata.txt'
	make_metadata_file_format_samplenames_order_as_counts(metadata_dict, deseq_metadata_file, feature_count_results_file, ['sample'] + conditions)

def count_and_make_files_for_deseq(work_dir, info_file, genome, covariates):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	print 'working on project:', project_name
	##get info on project
	with open(info_file, "U") as infh:
		metadata_dict = {}
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print conditions
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				metadata_dict[sample] = condition_by_sample
				sample_list.append(sample)
				fq1_info = line[1].split(',')
				fq2_info = line[2].split(',')
				##single or paired sequencing
				if fq2_info == ['']:
					seq_type = 'single_end'
				else:
					seq_type = 'paired_end'
					##if same fqs for r1 and r2
					if len(fq1_info) == len(fq2_info):
						pass
					else:
						print 'issue with sample %s, there are different number of fq files for read 1 and 2'%sample


	##run feature count on bams
	gtf_file = gtf_file_path[0] + genome + gtf_file_path[1]
	feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
	deseq_feature_count_file = project_name + '.star_fc.counts.txt'
	if seq_type == 'single_end':
		feature_count_single_end(star_bam_suffix, gtf_file, feature_count_results_file, sample_list)
	elif seq_type == 'paired_end':
		feature_count_paired_end(star_bam_suffix, gtf_file, feature_count_results_file, sample_list)
	format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)
	##make metadata
	deseq_metadata_file = project_name + '.star_fc.metadata.txt'
	make_metadata_file_format_samplenames_order_as_counts(metadata_dict, deseq_metadata_file, feature_count_results_file, ['sample'] + conditions)


def sort_index_bam(sample_name, genome):
	in_bam = sample_name + star_bam_suffix
	out_bam = sample_name + sorted_bam_suffix
	fa_file = fa_file_path[0] + genome + fa_file_path[1]
	# picard_rs = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'ReorderSam', 'INPUT=' + in_bam, 'OUTPUT=' + out_bam, 'CREATE_INDEX=true', 'REFERENCE=' + fa_file])
	picard_rs = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'SortSam', 'INPUT=' + in_bam, 'OUTPUT=' + out_bam, 'CREATE_INDEX=true', 'SORT_ORDER=coordinate'])
	picard_rs.wait()

def run_quorts_qc(sample_name, genome, out_dir):
	bam = sample_name + sorted_bam_suffix
	gtf_file = gtf_file_path[0] + genome + gtf_file_path[1]
	out_dir = out_dir + '/' + sample_name
	#java -jar ~/programs/hartleys-QoRTs-39cd1fc/QoRTs.jar QC HDBR525.sorted.bam /home/atimms/ngs_data/references/igenomes/hg19/genes.gtf qorts_test
	qourts_qc = subprocess.Popen(['java', '-Xmx100g', '-jar', qorts, 'QC', '--stranded', '--runFunctions', 'writeKnownSplices,writeNovelSplices,writeSpliceExon', bam, gtf_file, out_dir])
	qourts_qc.wait()

def run_quorts_mns(sample_file, genome, out_dir):
	gtf_file = gtf_file_path[0] + genome + gtf_file_path[1]
	qourts_mns = subprocess.Popen(['java', '-Xmx100g', '-jar', qorts, 'mergeNovelSplices', '--stranded', '--minCount', '6', out_dir, sample_file, gtf_file, out_dir])
	qourts_mns.wait()

def make_quorts_sample_file(s_dict, out_file):
	with open(out_file, "w") as outph:
		outph.write(delim.join(['sample.ID', 'group.ID', '\n']))
		for s in s_dict:
			print s, s_dict[s]
			outph.write(delim.join([s, s_dict[s][1], '\n']))

def map_with_star_and_sort_index_bams(work_dir, info_file, genome, covariates):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	quorts_out_folder = 'qorts_counts'

	print 'working on project:', project_name
	##get info on project
	with open(info_file, "U") as infh:
		metadata_dict = {}
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print conditions
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				metadata_dict[sample] = condition_by_sample
				sample_list.append(sample)
				fq1_info = line[1].split(',')
				fq2_info = line[2].split(',')
				# print fq1_info, len(fq1_info)
				# print fq2_info, len(fq2_info)
				##do we need to combine fqs
				if len(fq1_info) == 1:
					fq1 = fq1_info[0]
				elif len(fq1_info) >= 1:
					fq1 = sample + '.r1.fq.gz'
					combine_fqs(fq1_info, fq1)
				##single or paired sequencing
				if fq2_info == ['']:
					seq_type = 'single_end'
				else:
					seq_type = 'paired_end'
					##if same fqs for r1 and r2
					if len(fq1_info) == len(fq2_info):
						##combine files
						if len(fq2_info) == 1:
							fq2 = fq2_info[0]
						elif len(fq2_info) >= 1:
							fq2 = sample + '.r2.fq.gz'
							combine_fqs(fq2_info, fq2)
					else:
						print 'issue with sample %s, there are different number of fq files for read 1 and 2'%sample
				##run star

				# '''
				star_index_dir = star_ref_dir + genome
				print 'seq_type = ', seq_type
				if seq_type == 'single_end':
					star_align_single_end(sample, fq1, star_index_dir)
				elif seq_type == 'paired_end':
					star_align_paired_end(sample, [fq1, fq2], star_index_dir)
				##sort and index bam file from star
				sort_index_bam(sample, genome)
				##run qorts in prep for junctionseq -- arrayexpress data
				# run_quorts_qc(sample, genome, quorts_out_folder)
				# '''
	##make sample file
	# sample_info = 'decoder.bySample.txt'
	# make_quorts_sample_file(metadata_dict, sample_info)
	# ##run qorts
	# run_quorts_mns(sample_info, genome, quorts_out_folder)


def map_with_star_and_index_bams(work_dir, info_file, genome, covariates):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	quorts_out_folder = 'qorts_counts'

	print 'working on project:', project_name
	##get info on project
	with open(info_file, "U") as infh:
		metadata_dict = {}
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print conditions
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				metadata_dict[sample] = condition_by_sample
				sample_list.append(sample)
				fq1_info = line[1].split(',')
				fq2_info = line[2].split(',')
				# print fq1_info, len(fq1_info)
				# print fq2_info, len(fq2_info)
				##do we need to combine fqs
				if len(fq1_info) == 1:
					fq1 = fq1_info[0]
				elif len(fq1_info) >= 1:
					fq1 = sample + '.r1.fq.gz'
					combine_fqs(fq1_info, fq1)
				##single or paired sequencing
				if fq2_info == ['']:
					seq_type = 'single_end'
				else:
					seq_type = 'paired_end'
					##if same fqs for r1 and r2
					if len(fq1_info) == len(fq2_info):
						##combine files
						if len(fq2_info) == 1:
							fq2 = fq2_info[0]
						elif len(fq2_info) >= 1:
							fq2 = sample + '.r2.fq.gz'
							combine_fqs(fq2_info, fq2)
					else:
						print 'issue with sample %s, there are different number of fq files for read 1 and 2'%sample
				##run star

				# '''
				star_index_dir = star_ref_dir + genome
				print 'seq_type = ', seq_type
				if seq_type == 'single_end':
					star_align_single_end(sample, fq1, star_index_dir)
				elif seq_type == 'paired_end':
					star_align_paired_end(sample, [fq1, fq2], star_index_dir)
				##index bam file from star
				# sort_index_bam(sample, genome)
				bam_file = sample + '.' + star_bam_suffix
				samtools_index = subprocess.Popen(['samtools', 'index', bam_file])




def sort_index_bams_from_txt_file(work_dir, info_file, genome):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	print 'working on project:', project_name
	##get info on project
	with open(info_file, "U") as infh:
		line_count = 0
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print conditions
			else:
				sample = line[0]
				sort_index_bam(sample, genome)



def count_and_make_files_for_deseq_custom_annotation(work_dir, info_file, genome, covariates, ann_file, ann_type):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	print 'working on project:', project_name
	##get info on project
	with open(info_file, "U") as infh:
		metadata_dict = {}
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print conditions
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				metadata_dict[sample] = condition_by_sample
				sample_list.append(sample)
				fq1_info = line[1].split(',')
				fq2_info = line[2].split(',')
				##single or paired sequencing
				if fq2_info == ['']:
					seq_type = 'single_end'
				else:
					seq_type = 'paired_end'
					##if same fqs for r1 and r2
					if len(fq1_info) == len(fq2_info):
						pass
					else:
						print 'issue with sample %s, there are different number of fq files for read 1 and 2'%sample


	##run feature count on bams
	# gtf_file = gtf_file_path[0] + genome + gtf_file_path[1]
	feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
	deseq_feature_count_file = project_name + '.star_fc.counts.txt'
	if seq_type == 'single_end':
		if ann_type == 'gtf':
			feature_count_single_end(star_bam_suffix, ann_file, feature_count_results_file, sample_list)
	elif seq_type == 'paired_end':
		if ann_type == 'gtf':
			feature_count_paired_end(star_bam_suffix, ann_file, feature_count_results_file, sample_list)
		elif ann_type == 'saf':
			feature_count_paired_end_saf(star_bam_suffix, ann_file, feature_count_results_file, sample_list)
	format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

	##make metadata
	deseq_metadata_file = project_name + '.star_fc.metadata.txt'
	make_metadata_file_format_samplenames_order_as_counts(metadata_dict, deseq_metadata_file, feature_count_results_file, ['sample'] + conditions)

def count_and_make_files_for_deseq_allow_overlap(work_dir, info_file, genome, covariates):
	os.chdir(work_dir)
	project_name = info_file.split(".")[0]
	print 'working on project:', project_name
	##get info on project
	with open(info_file, "U") as infh:
		metadata_dict = {}
		line_count = 0
		sample_list = []
		for line in infh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				conditions = line[3:]
				print conditions
			else:
				sample = line[0]
				condition_by_sample = line[3:]
				metadata_dict[sample] = condition_by_sample
				sample_list.append(sample)
				fq1_info = line[1].split(',')
				fq2_info = line[2].split(',')
				##single or paired sequencing
				if fq2_info == ['']:
					seq_type = 'single_end'
				else:
					seq_type = 'paired_end'
					##if same fqs for r1 and r2
					if len(fq1_info) == len(fq2_info):
						pass
					else:
						print 'issue with sample %s, there are different number of fq files for read 1 and 2'%sample


	##run feature count on bams
	gtf_file = gtf_file_path[0] + genome + gtf_file_path[1]
	feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
	deseq_feature_count_file = project_name + '.star_fc.counts.txt'
	if seq_type == 'single_end':
		feature_count_single_end_allow_overlap(star_bam_suffix, gtf_file, feature_count_results_file, sample_list)
	elif seq_type == 'paired_end':
		feature_count_paired_end_allow_overlap(star_bam_suffix, gtf_file, feature_count_results_file, sample_list)
	format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)
	##make metadata
	deseq_metadata_file = project_name + '.star_fc.metadata.txt'
	make_metadata_file_format_samplenames_order_as_counts(metadata_dict, deseq_metadata_file, feature_count_results_file, ['sample'] + conditions)


##call methods


##laure 1117 --- example from60
'''
working_dir = '/data/atimms/laura_rnaseq_1117'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq
call_all_rnaseq_methods(working_dir, 'laura_1117_dys.txt', 'mm10', ['sample'])
'''

##kim array express data
working_dir = '/home/atimms/ngs_data/kim_array_express_1117'
##call all methods, so align, count and make filees for deseq
##just pons samples
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'ae_pons.txt', 'hg19', ['sample'])
# map_with_star_and_sort_index_bams(working_dir, 'ae_pons.txt', 'hg19', ['sample'])
##all samples
# map_with_star_and_sort_index_bams(working_dir, 'ae_all.txt', 'hg19', ['sample'])

##milena rnaseq  0818
working_dir = '/home/atimms/ngs_data/rnaseq/milena_rnaseq_0818'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'milena_rnaseq_0818.all_samples.txt', 'hg19', ['sample'])
##then make files for each time point - and rpt all, used . in wrong place for file names
# count_and_make_files_for_deseq(working_dir, 'milena_rnaseq_0818_all_samples.txt', 'hg19', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'milena_rnaseq_0818_day22.txt', 'hg19', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'milena_rnaseq_0818_day4.txt', 'hg19', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'milena_rnaseq_0818_day6.txt', 'hg19', ['sample'])
##remove lig4_1 wt sample
# count_and_make_files_for_deseq(working_dir, 'milena_rnaseq_0818_all_trimmed.txt', 'hg19', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'milena_rnaseq_0818_day22_trimmed.txt', 'hg19', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'milena_rnaseq_0818_day4_trimmed.txt', 'hg19', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'milena_rnaseq_0818_day6_trimmed.txt', 'hg19', ['sample'])
# sort_index_bams_from_txt_file(working_dir, 'milena_rnaseq_0818_all_samples.txt', 'hg19')

##tim cherry f1 data 0918
working_dir = '/home/atimms/ngs_data/misc/cherry_f1_hybrid/spreb_pwb_rnaseq_0918'
conditions_to_use = []
##
# count_and_make_files_for_deseq(working_dir, 'milena_rnaseq_0818_all_samples.txt', 'hg19', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'milena_rnaseq_0818_all_samples.txt', 'hg19', ['sample'])

##kim rnaseq  1018
working_dir = '/home/atimms/ngs_data/rnaseq/kim_rnaseq_hb_1018'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'kim_rnaseq_hb_1018.txt', 'hg19', ['sample'])


##greb1l rnaseq 1118
working_dir = '/home/atimms/ngs_data/rnaseq/dave_greb1l_1118'
conditions_to_use = []
##call all methods, so align, count and make files for deseq
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'greb1l_all_1118.txt', 'mm10', ['genotype'])
# count_and_make_files_for_deseq_custom_annotation(working_dir, 'greb1l_all_1118.txt', 'mm10', ['sample'], 'mm10-tRNAs.saf', 'saf')
# sort_index_bams_from_txt_file(working_dir, 'greb1l_all_1118.txt', 'hg19')

##greb1l hek293 rnaseq 1118
working_dir = '/home/atimms/ngs_data/rnaseq/dave_hek293_greb1l_1118'
conditions_to_use = []
##call all methods, so align, count and make files for deseq
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'hek293_greb1l_1118.txt', 'hg19', ['genotype'])

##kim/paul rnaseq and acomy
working_dir = '/home/atimms/ngs_data/rnaseq/paul_kim_rnaseq_1218'
conditions_to_use = []
##call all methods, so align, count and make files for deseq
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'paul_rnaseq_1218.txt', 'mm10', ['genotype'])
##had mislabelled 
# count_and_make_files_for_deseq(working_dir, 'paul_rnaseq_1218.txt', 'mm10', ['genotype'])

##kim/paul rnaseq and acomy
working_dir = '/home/atimms/ngs_data/rnaseq/kim_rett_rnaseq_0119'
conditions_to_use = []
##call all methods, so align, count and make files for deseq
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'kim_rett_rnaseq_0119.txt', 'hg19', ['region', 'dx'])
##sort and index bams
# sort_index_bams_from_txt_file(working_dir, 'kim_rett_rnaseq_0119.txt', 'hg19')
##redo so have counts for MECP2 transcript ie custom gtf (handmade)
# count_and_make_files_for_deseq_allow_overlap(working_dir, 'kim_rett_rnaseq_mecp1_0319.txt', 'hg19', ['sample'])

##kim/paul rnaseq and acomy
working_dir = '/home/atimms/ngs_data/misc/jimmy_sc_pik3ca_0119'
# map_with_star_and_sort_index_bams(working_dir, 'jimmy_sc_pikc3a_0119.txt', 'hg19', ['sample'])

##sophie rnaseq 0319
working_dir = '/home/atimms/ngs_data/rnaseq/sophie_carm1_0319'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'sophie_carm1_0319.txt', 'mm10', ['sample'])

##amrie rnaseq 0319
working_dir = '/home/atimms/ngs_data/rnaseq/amrei_jck_rnaseq_0319'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'amrei_jck_0319.txt', 'mm10', ['sample'])
##had mislabelled 
# count_and_make_files_for_deseq(working_dir, 'amrei_jck_0319_p10.txt', 'mm10', ['genotype'])
# count_and_make_files_for_deseq(working_dir, 'amrei_jck_0319_p15.txt', 'mm10', ['genotype'])
working_dir = '/home/atimms/ngs_data/rnaseq/amrei_pkd_rnaseq_0319'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'amrei_pkd_0319.txt', 'mm10', ['sample'])

##sophie rnaseq 0319 with extra timepoint
working_dir = '/home/atimms/ngs_data/rnaseq/sophie_carm1_2_0319'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'sophie_carm1_2_0319.txt', 'mm10', ['sample'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'sophie_carm1_e10_0319.txt', 'mm10', ['sample'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'sophie_carm1_e12_0319.txt', 'mm10', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'sophie_carm1_e10_0619.txt', 'mm10', ['sample'])

##amrie rnaseq 0419
working_dir = '/home/atimms/ngs_data/rnaseq/amrei_inv_rnaseq_0419'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'amrei_inv_0419.txt', 'mm10', ['sample'])

working_dir = '/home/atimms/ngs_data/rnaseq/amrei_jck_rnaseq_0419'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'amrei_jck_0419.txt', 'mm10', ['sample'])
##individual timepoints
# count_and_make_files_for_deseq(working_dir, 'amrei_jck_p5_0419.txt', 'mm10', ['genotype'])
# count_and_make_files_for_deseq(working_dir, 'amrei_jck_p10_0419.txt', 'mm10', ['genotype'])
# count_and_make_files_for_deseq(working_dir, 'amrei_jck_p15_0419.txt', 'mm10', ['genotype'])

working_dir = '/home/atimms/ngs_data/rnaseq/cherry_mactel_0519'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'cherry_mactel_0519.txt', 'GRCh38', ['sample'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'cherry_mactel_0519_nohu3.txt', 'GRCh38', ['sample'])

##combined scc nad displasia rnaseq for laura 0619`
working_dir = '/home/atimms/ngs_data/rnaseq/laura_rnaseq_0619'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'laura_rnaseq_0619.txt', 'mm10', ['sample'])


##amrie rnaseq 0619 
working_dir = '/home/atimms/ngs_data/rnaseq/amrei_inv_rnaseq_0619'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'amrei_inv_0619.txt', 'mm10', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'amrei_inv_p0_0619.txt', 'mm10', ['genotype'])
# count_and_make_files_for_deseq(working_dir, 'amrei_inv_p7_0619.txt', 'mm10', ['genotype'])

working_dir = '/home/atimms/ngs_data/rnaseq/amrei_pkd_rnaseq_0619'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'amrei_pkd_0619.txt', 'mm10', ['sample'])
##individual timepoints
# count_and_make_files_for_deseq(working_dir, 'amrei_pkd_p10_0619.txt', 'mm10', ['genotype'])
# count_and_make_files_for_deseq(working_dir, 'amrei_pkd_p20_0619.txt', 'mm10', ['genotype'])
# count_and_make_files_for_deseq(working_dir, 'amrei_pkd_p30_0619.txt', 'mm10', ['genotype'])
# count_and_make_files_for_deseq(working_dir, 'amrei_pkd_p30f_0619.txt', 'mm10', ['genotype'])

##kim array express data
working_dir = '/home/atimms/ngs_data/rnaseq/kim_hdbr_rnaseq_0619'
# map_with_star_and_index_bams(working_dir, 'hdbr_test.txt', 'hg19', ['sample'])
# map_with_star_and_index_bams(working_dir, 'hdbr_all.txt', 'hg19', ['sample'])


##kim/paul rnaseq and acomy
working_dir = '/home/atimms/ngs_data/rnaseq/kim_rett_rnaseq_0619'
conditions_to_use = []
##call all methods, so align, count and make files for deseq
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'kim_rett_rnaseq_0619.txt', 'hg19', ['region', 'dx'])


##tim/kevin rnaseq 0619
working_dir = '/home/atimms/ngs_data/rnaseq/kevin_rnaseq_0619/kevin_organoid_tc_0619'
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'kevin_organiod_tc_0619.txt', 'GRCh38', ['experiment'])
# count_and_make_files_for_deseq(working_dir, 'kevin_organiod_tc_0619.txt', 'GRCh38', ['experiment'])
working_dir = '/home/atimms/ngs_data/rnaseq/kevin_rnaseq_0619/kevin_rpe_library_0619'
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'kevin_rpe_ribo_depleted_0619.txt', 'GRCh38', ['experiment'])
# count_and_make_files_for_deseq(working_dir, 'kevin_rpe_ribo_depleted_0619.txt', 'GRCh38', ['experiment'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'kevin_rpe_polya_0619.txt', 'GRCh38', ['experiment'])
##use just 'affected' 
# count_and_make_files_for_deseq(working_dir, 'kevin_rpe_polya_aff_0619.txt', 'GRCh38', ['experiment'])
# count_and_make_files_for_deseq(working_dir, 'kevin_rpe_ribo_depleted_aff_0619.txt', 'GRCh38', ['experiment'])
# count_and_make_files_for_deseq(working_dir, 'kevin_rpe_polya_heta_homu_0619.txt', 'GRCh38', ['experiment'])
# count_and_make_files_for_deseq(working_dir, 'kevin_rpe_ribo_depleted_heta_homu_0619.txt', 'GRCh38', ['experiment'])


##vishal/lan rnaseq 0719
working_dir = '/home/atimms/ngs_data/rnaseq/lan_pig_cpb_rnaseq_0719'
##pig data: call all methods, so align, count and make filees for deseq 
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'lan_pig_cpb_0719.txt', 'Sscrofa11.1', ['treatment'])
##mouse data: call all methods, so align, count and make filees for deseq then analyze the 2 litters sepratley
working_dir = '/home/atimms/ngs_data/rnaseq/lan_mouse_hypoxia_rnaseq_0719'
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'lan_mouse_p15_hypoxia_0719_all.txt', 'mm10', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'lan_mouse_p15_hypoxia_0719_all.txt', 'mm10', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'lan_mouse_p15_hypoxia_0719_l1.txt', 'mm10', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'lan_mouse_p15_hypoxia_0719_l2.txt', 'mm10', ['sample'])

##tim/kevin rnaseq 0819
working_dir = '/home/atimms/ngs_data/rnaseq/kevin_rnaseq_0819'
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'sorted_cell_0819_all.txt', 'GRCh38', ['experiment'])
# count_and_make_files_for_deseq(working_dir, 'sorted_cell_0819_c4.txt', 'GRCh38', ['experiment'])
# count_and_make_files_for_deseq(working_dir, 'sorted_cell_0819_c5.txt', 'GRCh38', ['experiment'])
# count_and_make_files_for_deseq(working_dir, 'sorted_cell_0819_with_organoid.txt', 'GRCh38', ['experiment'])
# count_and_make_files_for_deseq(working_dir, 'sorted_cell_0819_with_organoid28.txt', 'GRCh38', ['experiment'])

##tim/kevin rnaseq 0819
working_dir = '/home/atimms/ngs_data/rnaseq/kevin_rnaseq_1019'
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'kevin_rnaseq_1019.txt', 'GRCh38', ['experiment'])


##tim rnaseq 1019
working_dir = '/home/atimms/ngs_data/rnaseq/cherry_rnaseq_1019'
##run on all sampeles
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'cherry_rnaseq_1019.all.txt', 'GRCh38', ['experiment'])
##then just analyze retina and rpe
# count_and_make_files_for_deseq(working_dir, 'cherry_rnaseq_1019.ret_rpe.txt', 'GRCh38', ['experiment'])
##had to change name of txt file
# count_and_make_files_for_deseq(working_dir, 'cherry_rnaseq_1019_ret_rpe.txt', 'GRCh38', ['experiment'])
# count_and_make_files_for_deseq(working_dir, 'cherry_rnaseq_1019_all.txt', 'GRCh38', ['experiment'])


##sophie rnaseq 1119
working_dir = '/home/atimms/ngs_data/rnaseq/sophie_carm1_1119'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'sophie_carm1_1119.txt', 'mm10', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'sophie_carm1_1119_e10.txt', 'mm10', ['sample'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'sophie_carm1_1119_e10all.txt', 'mm10', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'sophie_carm1_1119_e10split.txt', 'mm10', ['sample'])
##added biad metric from picard, so make new files
# count_and_make_files_for_deseq(working_dir, 'sophie_carm1_1119_e10_bias.txt', 'mm10', ['sample'])


##sophie rnaseq 1119
working_dir = '/home/atimms/ngs_data/rnaseq/lisa_1119'
conditions_to_use = []
##call all methods, so align, count and make filees for deseq -- for all time points
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'lisa_zf_1119.txt', 'danRer10', ['sample'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'lisa_zf_1119_group.txt', 'danRer10', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'lisa_zf_1119_g1.txt', 'danRer10', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'lisa_zf_1119_g2.txt', 'danRer10', ['sample'])

##kim array express data
working_dir = '/home/atimms/ngs_data/rnaseq/kim_hdbr_rnaseq_0619'
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'kim_hdbr_rnaseq_all_1219.txt', 'GRCh38', ['sample'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'kim_hdbr_rnaseq_test_1219.txt', 'GRCh38', ['sample'])

##michael strain data 0120
working_dir = '/home/atimms/ngs_data/rnaseq/cunn_rnaseq_0120'
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'strain_rnaseq_0120_all.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0120_all.txt', 'GRCh38', ['sample'])

##michael strain data 0120 adn 0220 combined
working_dir = '/home/atimms/ngs_data/rnaseq/cunn_rnaseq_0220'
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'strain_rnaseq_0220_all.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_all.txt', 'GRCh38', ['sample'])
##by individual sample
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_1017.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_1032.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_3007.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_3035.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_4025.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_4032.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_AUT25.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C1625.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C1653.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C1671.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C1856.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C1859.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C1926.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C1957.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C2038.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C2084.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C3049.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_C3066.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_OST85.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_STL19.txt', 'GRCh38', ['sample'])
##do by individual sample 2 ways, by itself and then including controls
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_1017.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_1017_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_1032.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_1032_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_3007.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_3007_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_3035.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_3035_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_4025.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_4025_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_4032.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_4032_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_AUT25.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1625.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1625_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1653.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1653_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1671.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1671_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1856.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1856_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1859.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1859_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1926.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1926_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1957.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C1957_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C2038.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C2038_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C2084.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C2084_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C3049.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C3049_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C3066.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_C3066_ctrls.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_OST85.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_STL19.txt', 'GRCh38', ['sample'])
# ##rpt this one
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_all.txt', 'GRCh38', ['sample'])
##by gene group
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_AXL.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_CTRL.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_FLNA.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_FLNB.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_FLNC.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'strain_rnaseq_0220_PIEZO1.txt', 'GRCh38', ['sample'])

##kenny 2x experiments 0620
working_dir = '/home/atimms/ngs_data/misc/kenny_mouse_0620/rnaseq'
conditions_to_use = []
##call all methods, so align, count and make files for deseq
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'wegner_0620.txt', 'mm10', ['sample'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'eucomm_0620.txt', 'mm10', ['sample'])


##vishal data 0620
working_dir = '/home/atimms/ngs_data/rnaseq/vishal_rnaseq_0620'
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'dHL60_0620.txt', 'GRCh38', ['sample'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'U118_0620.txt', 'GRCh38', ['sample'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'dSY5Y_0620.txt', 'GRCh38', ['sample'])
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'WTC11_0620.txt', 'GRCh38', ['sample'])

##kim new bunch of rnaseq 0720 
working_dir = '/home/atimms/ngs_data/rnaseq/kim_rnaseq_0720'
# call_rnaseq_methods_for_deseq_analysis(working_dir, 'kim_rnaseq_0720_all.txt', 'GRCh38', ['sample'])
# count_and_make_files_for_deseq(working_dir, 'kim_rnaseq_0720_all.txt', 'GRCh38', ['sample'])
count_and_make_files_for_deseq(working_dir, 'kim_rnaseq_0720_DWM_Bulk_vs_CTL_Bulk.txt', 'GRCh38', ['sample'])
count_and_make_files_for_deseq(working_dir, 'kim_rnaseq_0720_DWM_PCL_vs_CTL_PCL.txt', 'GRCh38', ['sample'])
count_and_make_files_for_deseq(working_dir, 'kim_rnaseq_0720_DWM_RL_vs_CTL_RL.txt', 'GRCh38', ['sample'])
count_and_make_files_for_deseq(working_dir, 'kim_rnaseq_0720_DWM_vs_CTL.txt', 'GRCh38', ['sample'])
count_and_make_files_for_deseq(working_dir, 'kim_rnaseq_0720_DWM_EGL_vs_CTL_EGL.txt', 'GRCh38', ['sample'])






