#!/usr/bin/env python
import subprocess
import os
import glob

'''
info....

##for sinto part i.e. getting bam files 
##go to worker node and load biobuilds and python3 and open env 
qsub -Iq cdbrmq -l mem=200gb,ncpus=20 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load biobuilds
module load local_python/3.7.6
source activate sinto
##program is now at:
/home/atimms/programs/sinto/scripts/sinto

##for macs2
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load local_python/3.7.6 
source activate macs2

##homer
load modules:
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load biobuilds
module load homer/4.9.1

##for later homer, can't use homer module
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878
conda activate homer
module load R/4.0.3

code combied from:
cherry_scatac_archr_call_peaks_0920_cyb.py
cherry_scatac_collapse_corr_1020_cyb.py
cherry_scatac_corr_pairwise_0221.py

'''
##parameters
delim = '\t'
thread_number = '20'
##programs
sinto = '/home/atimms/programs/sinto/scripts/sinto'
genrich = '/home/atimms/programs/Genrich/Genrich'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
bg_to_bw = '/home/atimms/programs/bedGraphToBigWig'


##methods
def get_cell_info_org_human(infile, outfile_suffix):
	##make dict of all cell barcode by sample
	cell_info_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.replace('"', '').rstrip().split(',')
				# print(line)
				sample = line[1].split('#')[0]
				cell_barcode = line[1].split('#')[1]
				cell_class = line[2].replace(' ', '_').replace('/', '_')
				# print(sample, cell_barcode)
				if sample in cell_info_dict:
					cell_info_dict[sample].append(cell_barcode + ',' +  cell_class)
				else:
					cell_info_dict[sample] = [cell_barcode + ',' +  cell_class]
	##print out just for ipscs with ipsc as cell type
	out_files = []
	for s in cell_info_dict:
		print(s, len(cell_info_dict[s]))
		outfile = s + outfile_suffix
		out_files.append(outfile)
		with open(outfile, "w") as out_fh:
			for bc in cell_info_dict[s]:
				cell_info = bc.split(',')
				if len(cell_info) != 2:
					print(cell_info, 'looks odd')
				out_fh.write(delim.join(cell_info) + '\n')
	return(out_files)

def run_sinto_method(sample_name, bam_file, cells_file, w_dir):
	directory = './' + sample_name
	if not os.path.exists(sample_name):
		os.makedirs(sample_name)
	os.chdir(directory)
	run_sinto = subprocess.Popen([sinto, 'filterbarcodes', '-b', w_dir + bam_file,  '-c', w_dir + cells_file, '-p', thread_number])
	run_sinto.wait()
	os.chdir(w_dir)

def collapse_samples_by_cell_class(infile):
	coll_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				cell_class = line[2]
				sample = line[1]
				bam_file = sample + '/' + cell_class + '.bam'
				type_cc = line[0] + '.' + cell_class
				if type_cc in coll_dict:
					coll_dict[type_cc].append(bam_file)
				else:
					coll_dict[type_cc] = [bam_file]
	return(coll_dict)

def merge_bams(collapsed_dict, w_dir):
	for sample in collapsed_dict:
		merged_bam = sample + '.bam'
		# input_bams = ['I=' + b for b in collapsed_dict[sample]]
		input_bams = collapsed_dict[sample]
		print(merged_bam, input_bams)
		##the wk had slightly differeent reference, so we got an error, use samtools
		# picard_md = subprocess.Popen(['picard', 'MergeSamFiles', 'CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'] + input_bams + ['O=' + merged_bam, 'TMP_DIR=' + w_dir])
		# picard_md.wait()
		st_merge = subprocess.Popen(['samtools', 'merge', '-O', 'bam', '-@', '16', merged_bam] + input_bams)
		st_merge.wait()
		st_index = subprocess.Popen(['samtools', 'index', merged_bam])
		st_index.wait()

def make_tag_dirs(samples):
	for sample in samples:
		bam = sample + '.bam'
		outdir = sample + '.tag_dir'
		mk_tag_dir = subprocess.Popen(['makeTagDirectory', outdir, bam])
		mk_tag_dir.wait()

def run_macs2(samples, qvalues):
	##run macs2
	for sample in samples:
		for qvalue in qvalues:
			bam = sample + '.bam'
			outname = sample + '.macs2_q' + qvalue
			print(sample, bam, outname)
			##run macs2
			run_macs2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', outname, '-g', 'hs', '-q', qvalue, '--keep-dup', 'all', '--tempdir', '.'])
			run_macs2.wait()

def merge_peaks(out_prefix, d_value, peak_files, out_suffix):
	out_file = out_prefix + d_value + 'd' + out_suffix
	print(len(peak_files), out_file)
	with open(out_file, 'w') as out_fh:
		# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
		run_merge_peaks = subprocess.Popen(['mergePeaks', '-d', d_value] + peak_files, stdout=out_fh)
		run_merge_peaks.wait()

def annotate_peaks(samples, peak_file, d_value, tag_suffix, out_prefix, sizes_req, peak_suffix):
	tag_dirs = []
	##get list of tag_dirs
	for sample in samples:
		tag_dir = sample + tag_suffix
		tag_dirs.append(tag_dir)
	print(len(tag_dirs), peak_file)
	##annotate
	for size_req in sizes_req:
		if size_req == 'none':
			out_file_no_size = out_prefix + d_value + 'd.' + size_req + 'size' + peak_suffix
			with open(out_file_no_size, 'w') as out_ns_fh:
				mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-d'] + tag_dirs, stdout=out_ns_fh)
				mk_tag_dir.wait()
		else:
			out_file_using_size = out_prefix + d_value + 'd.' + size_req + 'size' + peak_suffix
			with open(out_file_using_size, 'w') as out_fh:
				mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-size', size_req, '-d'] + tag_dirs, stdout=out_fh)
				mk_tag_dir.wait()

def get_correlation_info_homer(samples, d_values, size_values, bed_suffices, merge_prefix, annotate_prefix):

	##combine bed files
	# '''
	comb_beds = []
	for s in samples:
		bed = s + bed_suffix
		comb_beds.append(bed)
	for d in d_values:
		merge_peaks(merge_prefix, d, comb_beds, bed_suffix)
	
	##annotate files
	tag_dir_suffix = '.tag_dir'
	# correlation_prefix = 'correlation.'
	for d in d_values:
		merged_bed = merge_prefix + d + 'd' + bed_suffix
		outfile_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
		##annotate those peaks, and then make a correlation file
		annotate_peaks(samples, merged_bed, d, tag_dir_suffix, annotate_prefix, size_values, outfile_suffix)
	# '''

def format_correlation_all_samples(annotate_prefix, annotate_suffix, d_values, sizes_req, out_prefix):
	##format file by file
	for size_req in sizes_req:
		for d_value in d_values:
			infile = annotate_prefix + d_value + 'd.' + size_req + 'size' + annotate_suffix
			out_file = out_prefix  + infile.split('.', 1)[1]
			print(infile, out_file)
			with open(out_file, 'w') as out_fh, open(infile, 'r') as in_fh:
				lc = 0
				for line in in_fh:
					lc += 1
					line = line.strip('\n').split(delim)
					if lc == 1:
						sample_info = line[19:]
						sample_info = [s.split('.')[0] + '.' + s.split('.')[1] for s in sample_info]
						out_fh.write(delim.join(['peak_id'] + sample_info + ['\n']))
						# print(sample_info)
					else:
						line_out = [line[0]] + line[19:]
						out_fh.write(delim.join(line_out + ['\n']))

def get_correlation_info_homer_pairwise(cc_dict, d_values, size_values, bed_suffixes, tag_dir_suffix):
	for bed_suffix in bed_suffixes:
		for cell_class in cc_dict:
			samples = cc_dict[cell_class]
			##combine bed files
			merge_prefix = cell_class + '_merged.'
			annotate_prefix = cell_class + '_annotated.'
			annotate_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
			corr_prefix = cell_class + '_correlation.'
			# '''
			comb_beds = []
			for s in samples:
				bed = s + bed_suffix
				comb_beds.append(bed)
			for d in d_values:
				merge_peaks(merge_prefix, d, comb_beds, bed_suffix)
			##annotate files
			tag_dir_suffix = '.tag_dir'
			# correlation_prefix = 'correlation.'
			for d in d_values:
				merged_bed = merge_prefix + d + 'd' + bed_suffix
				##annotate those peaks, and then make a correlation file
				annotate_peaks(samples, merged_bed, d, tag_dir_suffix, annotate_prefix, size_values, annotate_suffix)
			# '''
			##format annotation files
			for size_req in size_values:
				for d_value in d_values:
					infile = annotate_prefix + d_value + 'd.' + size_req + 'size' + annotate_suffix
					out_file = corr_prefix + infile.split('.', 1)[1]
					# print(infile, out_file)
					with open(out_file, 'w') as out_fh, open(infile, 'r') as in_fh:
						lc = 0
						for line in in_fh:
							lc += 1
							line = line.strip('\n').split(delim)
							if lc == 1:
								sample_info = line[19:]
								sample_info = [s.split('.')[0] + '.' + s.split('.')[1] for s in sample_info]
								out_fh.write(delim.join(['peak_id'] + sample_info) + '\n')
								# print(sample_info)
							else:
								chrom = line[1]
								if '_' not in chrom and chrom != 'chrM' and chrom != 'chrY':
									human_count = float(line[19])
									org_count = float(line[20])
									# if min([human_count, org_count]) > 1:
									line_out = [line[0]] + line[19:]
									out_fh.write(delim.join(line_out) + '\n')

def make_summary_files_for_scatterplots(infiles, outfile_suffix):
	for infile in infiles:
		out_file = infile.rsplit('.', 1)[0] + outfile_suffix
		# print(infile, out_file)
		with open(out_file, 'w') as out_fh, open(infile, 'r') as in_fh:
			lc = 0
			out_fh.write(delim.join(['peak', 'human', 'organoid']) + '\n')
			for line in in_fh:
				lc += 1
				line = line.strip().split(delim)
				if lc > 1:
					peak = line[0]
					chrom = line[1]
					if '_' not in chrom and 'KI' not in chrom and 'GL' not in chrom and chrom != 'chrM' and chrom != 'chrY':
						human_counts = [float(i) for i in line[19:36]]
						org_counts = [float(i) for i in line[36:]]
						human_mean = sum(human_counts) / len(human_counts)
						org_mean = sum(org_counts) / len(org_counts)
						# if min([human_count, org_count]) > 1:
						line_out = [peak, str(human_mean), str(org_mean)]
						out_fh.write(delim.join(line_out) + '\n')

def get_interesting_peaks(infiles, outfile_suffix):
	for infile in infiles:
		out_file = infile.rsplit('.', 1)[0] + outfile_suffix
		# print(infile, out_file)
		with open(out_file, 'w') as out_fh, open(infile, 'r') as in_fh:
			lc = 0
			
			for line in in_fh:
				lc += 1
				line = line.strip().split(delim)
				if lc == 1:
					out_fh.write(delim.join(line[:19] + ['human mean', 'organoid mean', 'overall mean', 'organoid mean / human mean', 'human mean / organoid mean'] + line[19:]) + '\n')
				else:
					peak = line[0]
					chrom = line[1]
					if '_' not in chrom and 'KI' not in chrom and 'GL' not in chrom and chrom != 'chrM' and chrom != 'chrY':
						human_counts = [float(i) for i in line[19:36]]
						org_counts = [float(i) for i in line[36:]]
						human_mean = sum(human_counts) / len(human_counts)
						org_mean = sum(org_counts) / len(org_counts)
						all_mean = sum(human_counts + org_counts) / len(human_counts + org_counts)
						# print(org_counts)
						# if org_mean == 0.0:
						# 	org_mean = 0.01
						# if human_mean == 0.0:
						# 	human_mean = 0.01
						# if human_mean > 400 and org_mean <10:
						# 	print(human_mean, org_mean, all_mean, human_mean/org_mean)
						if all_mean >= 20:
							if org_mean/human_mean >= 5 or human_mean/org_mean >= 5:
								line_out = line[:19] + [str(human_mean), str(org_mean), str(all_mean), str(org_mean/human_mean), str(human_mean/org_mean)] + line[19:]
								out_fh.write(delim.join(line_out) + '\n')


def combine_tag_dirs(samples_to_combine, combined_name):
	outdir = combined_name + '.tag_dir'
	tag_dir_cmds = []
	for sample in samples_to_combine:
		tag_dir_cmd = [ '-d', sample + '.tag_dir']
		tag_dir_cmds.extend(tag_dir_cmd)
	print(tag_dir_cmds)
	mk_tag_dir = subprocess.Popen(['makeTagDirectory', outdir] + tag_dir_cmds)
	mk_tag_dir.wait()

def get_diff_peaks_diff_sizes(out_prefix1, out_prefix2, d_values, q_values, sizes_req, human_name, org_name):
	human_tag_dir = human_name  + '.tag_dir'
	org_tag_dir = org_name  + '.tag_dir'
	for d_value in d_values:
		for q_value in q_values:
			for size_req in sizes_req:
				merged_file = 'merged.' + d_value + 'd.macs2.bampe_q' + q_value + '_keepdups_summits.txt'
				out_file1 = out_prefix1 + 'size' + size_req + '.' + merged_file.split('.', 1)[1]
				out_file2 = out_prefix2 + 'size' + size_req + '.' + merged_file.split('.', 1)[1]
				with open(out_file1, 'w') as out_fh:
					# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
					run_diff_peaks1 = subprocess.Popen(['getDifferentialPeaks', merged_file, human_tag_dir, org_tag_dir, '-size', size_req], stdout=out_fh)
					run_diff_peaks1.wait()
				with open(out_file2, 'w') as out_fh:
					# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
					run_diff_peaks2 = subprocess.Popen(['getDifferentialPeaks', merged_file, org_tag_dir, human_tag_dir, '-size', size_req], stdout=out_fh)
					run_diff_peaks2.wait()

def get_diff_peaks_diff_sig2(out_prefix1, out_prefix2, d_values, q_values, p_values, fc_values, sizes_req, human_name, org_name):
	human_tag_dir = human_name  + '.tag_dir'
	org_tag_dir = org_name  + '.tag_dir'
	for d_value in d_values:
		for q_value in q_values:
			for p_value in p_values:
				for fc_value in fc_values:
					for size_req in sizes_req:
						merged_file = 'merged.' + d_value + 'd.macs2.bampe_q' + q_value + '_keepdups_summits.txt'
						out_file1 = out_prefix1 + 'size' + size_req + '.p_' + p_value + '.fc_' + fc_value + '.'  + merged_file.split('.', 1)[1]
						out_file2 = out_prefix2 + 'size' + size_req + '.p_' + p_value + '.fc_' + fc_value + '.'  + merged_file.split('.', 1)[1]
						with open(out_file1, 'w') as out_fh:
							# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
							run_diff_peaks1 = subprocess.Popen(['getDifferentialPeaks', merged_file, human_tag_dir, org_tag_dir, '-F', fc_value, '-P', p_value, '-size', size_req], stdout=out_fh)
							run_diff_peaks1.wait()
						with open(out_file2, 'w') as out_fh:
							# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
							run_diff_peaks2 = subprocess.Popen(['getDifferentialPeaks', merged_file, org_tag_dir, human_tag_dir, '-F', fc_value, '-P', p_value, '-size', size_req], stdout=out_fh)
							run_diff_peaks2.wait()


def get_diff_peaks_diff_sig(out_prefix1, out_prefix2, d_values, q_values, p_values, fc_values, human_name, org_name):
	human_tag_dir = human_name  + '.tag_dir'
	org_tag_dir = org_name  + '.tag_dir'
	for d_value in d_values:
		for q_value in q_values:
			for p_value in p_values:
				for fc_value in fc_values:
					merged_file = 'merged.' + d_value + 'd.macs2.bampe_q' + q_value + '_keepdups_summits.txt'
					out_file1 = out_prefix1 + 'p_' + p_value + '.fc_' + fc_value + '.'  + merged_file.split('.', 1)[1]
					out_file2 = out_prefix2 + 'p_' + p_value + '.fc_' + fc_value + '.'  + merged_file.split('.', 1)[1]
					with open(out_file1, 'w') as out_fh:
						# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
						run_diff_peaks1 = subprocess.Popen(['getDifferentialPeaks', merged_file, human_tag_dir, org_tag_dir, '-F', fc_value, '-P', p_value], stdout=out_fh)
						run_diff_peaks1.wait()
					with open(out_file2, 'w') as out_fh:
						# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
						run_diff_peaks2 = subprocess.Popen(['getDifferentialPeaks', merged_file, org_tag_dir, human_tag_dir, '-F', fc_value, '-P', p_value], stdout=out_fh)
						run_diff_peaks2.wait()


def get_diff_peaks_replicates(out_prefix1, out_prefix2, human_names, org_names, all_peaks, ctl_type):
	human_tag_dirs = [i  + '.tag_dir' for i in human_names]
	org_tag_dirs = [i  + '.tag_dir' for i in org_names]

	out_file1 = out_prefix1 + 'balanced.txt'
	out_file2 = out_prefix2 + 'balanced.txt'
	if ctl_type == 'input':
		if all_peaks == 'no':
			with open(out_file1, 'w') as out_fh:
				run_diff_peaks1 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + human_tag_dirs + ['-i'] + org_tag_dirs + ['-balanced', '-genome', 'hg38'], stdout=out_fh)
				run_diff_peaks1.wait()
			with open(out_file2, 'w') as out_fh:
				run_diff_peaks2 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + org_tag_dirs + ['-i'] + human_tag_dirs + ['-balanced', '-genome', 'hg38'], stdout=out_fh)
				run_diff_peaks2.wait()
		elif all_peaks == 'yes':
			with open(out_file1, 'w') as out_fh:
				run_diff_peaks1 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + human_tag_dirs + ['-i'] + org_tag_dirs + ['-balanced', '-genome', 'hg38', '-all'], stdout=out_fh)
				run_diff_peaks1.wait()
			with open(out_file2, 'w') as out_fh:
				run_diff_peaks2 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + org_tag_dirs + ['-i'] + human_tag_dirs + ['-balanced', '-genome', 'hg38', '-all'], stdout=out_fh)
				run_diff_peaks2.wait()
	elif ctl_type == 'background':
		if all_peaks == 'no':
			with open(out_file1, 'w') as out_fh:
				run_diff_peaks1 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + human_tag_dirs + ['-b'] + org_tag_dirs + ['-balanced', '-genome', 'hg38'], stdout=out_fh)
				run_diff_peaks1.wait()
			with open(out_file2, 'w') as out_fh:
				run_diff_peaks2 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + org_tag_dirs + ['-b'] + human_tag_dirs + ['-balanced', '-genome', 'hg38'], stdout=out_fh)
				run_diff_peaks2.wait()
		elif all_peaks == 'yes':
			with open(out_file1, 'w') as out_fh:
				run_diff_peaks1 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + human_tag_dirs + ['-b'] + org_tag_dirs + ['-balanced', '-genome', 'hg38', '-all'], stdout=out_fh)
				run_diff_peaks1.wait()
			with open(out_file2, 'w') as out_fh:
				run_diff_peaks2 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + org_tag_dirs + ['-b'] + human_tag_dirs + ['-balanced', '-genome', 'hg38', '-all'], stdout=out_fh)
				run_diff_peaks2.wait()


def find_peaks_cell_classes(names, peak_suffix):
	for name in names:
		tag_dir = name + '.tag_dir'
		peak_file = name + peak_suffix
		mk_tag_dir = subprocess.Popen(['findPeaks', tag_dir, '-o', peak_file])
		mk_tag_dir.wait()



def make_beds_from_diff_peaks_file(homer_dp_files):
	for homer_dp_file in homer_dp_files:
		out_bed = homer_dp_file.rsplit('.', 1)[0] + '.bed'
		with open(out_bed, 'w') as out_fh, open(homer_dp_file, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.strip('\n').split(delim)
				if lc > 1:
					chrom = line[1]
					if '.' not in chrom and '_' not in chrom:
						out_fh.write(delim.join(line[1:4]) + '\n')

def make_beds_from_find_peaks_files(samples, infile_suffix):
	for sample in samples:
		infile = sample + infile_suffix
		out_bed = infile.rsplit('.', 1)[0] + '.bed'
		with open(out_bed, 'w') as out_fh, open(infile, 'r') as in_fh:
			for line in in_fh:
				if line[0] != '#':
					line = line.strip('\n').split(delim)
					chrom = line[1]
					if '.' not in chrom and '_' not in chrom:
						out_fh.write(delim.join(line[1:4]) + '\n')

def bt_int_homer_beds_vs_cc_beds(homer_result_beds, cc_samples, cc_bed_suffix):
	cc_beds = [i + cc_bed_suffix for i in cc_samples]
	for homer_result_bed in homer_result_beds:
		out_bti = homer_result_bed.rsplit('.', 1)[0] + '.bt_int.mature_cc.txt'
		##bedtools intersect 
		with open(out_bti, "w") as out_fh: 
			hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', homer_result_bed, '-b'] + cc_beds + ['-C', '-filenames'], stdout=out_fh)
			hom_bt_intersect.wait()

def filter_and_make_beds_from_diff_peaks_file(homer_dp_files, cut_offs_list):
	for homer_dp_file in homer_dp_files:
		for cut_offs in cut_offs_list:
			fc_wanted = cut_offs[0]
			p_wanted = cut_offs[1]
			out_bed = homer_dp_file.rsplit('.', 1)[0] +'.lfc' + str(fc_wanted) + '_p' + str(p_wanted) + '.bed'
			with open(out_bed, 'w') as out_fh, open(homer_dp_file, 'r') as in_fh:
				lc = 0
				for line in in_fh:
					lc += 1
					line = line.strip('\n').split(delim)
					if lc > 1:
						chrom = line[1]
						logfc = float(line[24])
						adjp = float(line[26])
						if '.' not in chrom and '_' not in chrom:
							if logfc >= fc_wanted and adjp <= p_wanted:
								out_fh.write(delim.join(line[1:4]) + '\n')

##run methods

##params etc
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scatac_corr_rpts_1021/'
os.chdir(working_dir)


##just for new organoid data....

##step1. make cell info files for sinto
org_raw_info = 'org_all.cell_class_info.101221.csv'
cell_info_file_suffix = '.cell_types.txt'
##format files
# org_cell_info_files = get_cell_info_org_human(org_raw_info, cell_info_file_suffix)
# print(org_cell_info_files)

##step2. run sinto
bam_dict = {'20wk1': '20wk.possorted.bam', '20wk2': '20wk_c5_1.possorted.bam', '28wk1': '28-1.possorted.bam', 
		'28wk2': '28-2.possorted.bam', '5wk1': '5wk.possorted.bam', '5wk2': '5wk_c5_1.possorted.bam', 
		'12wk1': '12wk1.possorted.bam', '12wk2': '12wk2.possorted.bam', '12wk3': '12wk3.possorted.bam'}
'''
for cell_info_file in org_cell_info_files:
	sample = cell_info_file.split('.')[0]
	bam = bam_dict[sample]
	# print(sample, bam, cell_info_file)
	run_sinto_method(sample, bam, cell_info_file, working_dir)
'''

##step3. merge bams
##cell class info file made from 'cell_types.txt' files and remove any cell class with less than 20 cells per sample
sample_cell_class_info = 'org_samples_cell_classes_gt20_1021.txt'
##get info
# coll_cell_class_dict = collapse_samples_by_cell_class(sample_cell_class_info)
# print(coll_cell_class_dict)
##merge bams
# merge_bams(coll_cell_class_dict, working_dir)

##step4. make tag dir
tag_dir_samples = ['organoid.AC_HC_Precursors', 'organoid.Cones', 'organoid.Early_RPCs', 'organoid.Late_RPCs', 'organoid.PR_BC_Precursors', 
		'organoid.RGCs', 'organoid.Rods', 'organoid.Amacrine_Horizontal_Cells', 'organoid.Bipolar_Cells', 'organoid.Muller_Glia', 'organoid.Developing_RGCs']
# make_tag_dirs(tag_dir_samples)


##for human and organoid data...

##mv tag dir and cell class bams for human data to working dir
# mv ../cherry_scatac_corr_pairwise_0221/human*dir .
# mv ../cherry_scatac_archr_collapse_corr_1020/human.*bam* .

##step5. call macs2 q0.01 and q0.000001 on all cell classes
human_org_ccs = ['human.AC_HC_GC_Precursors', 'human.Developing_Amacrines', 'human.Developing_Bipolars', 'human.Developing_Cones', 
	'human.Developing_Ganglions', 'human.Developing_Horizontals', 'human.Developing_Rods', 'human.Early_Progenitors', 
	'human.Ganglion_Precursors', 'human.Late_Progenitors', 'human.Mature_Amacrines', 'human.Mature_Bipolars', 
	'human.Mature_Cones', 'human.Mature_Horizontals', 'human.Mature_Mullers', 'human.Mature_Rods', 'human.Photoreceptor_Bipolar_Precursors',
	'organoid.AC_HC_Precursors', 'organoid.Cones', 'organoid.Early_RPCs', 'organoid.Late_RPCs', 'organoid.PR_BC_Precursors', 
	'organoid.RGCs', 'organoid.Rods', 'organoid.Amacrine_Horizontal_Cells', 'organoid.Bipolar_Cells', 'organoid.Muller_Glia', 'organoid.Developing_RGCs']

q_values = ['0.01', '0.000001']
# run_macs2(human_org_ccs, q_values)

##step6. homer merge peaks/annotate... pretty sure I didn't need this
## 3. run homer i.e. make tag dirs, combine beds and annotate
d_values_wanted = ['100', '500']
size_values_wanted = ['500', '2000'] ##does none work for homer
bed_suffix = '.macs2_q0.000001_summits.bed'
txt_suffix =  '.macs2_q0.000001_summits.txt'
merged_peak_prefix = 'merged.collapsed_cc.'
annotated_peak_prefix = 'annotated.collapsed_cc.'
correlation_peak_prefix = 'correlation.collapsed_cc.'
## a. run homer using q0.000001
# get_correlation_info_homer(human_org_ccs, d_values_wanted, size_values_wanted, bed_suffix, merged_peak_prefix, annotated_peak_prefix)
## b. format corrlation file for graphing using q0.000001
# format_correlation_all_samples(annotated_peak_prefix, txt_suffix, d_values_wanted, size_values_wanted, correlation_peak_prefix)

##repeat for q0.01
bed_suffix = '.macs2_q0.01_summits.bed'
txt_suffix =  '.macs2_q0.01_summits.txt'
# get_correlation_info_homer(human_org_ccs, d_values_wanted, size_values_wanted, bed_suffix, merged_peak_prefix, annotated_peak_prefix)
# format_correlation_all_samples(annotated_peak_prefix, txt_suffix, d_values_wanted, size_values_wanted, correlation_peak_prefix)


##step7. pairwise comparisons

cell_class_dict = {'cones':['human.Mature_Cones', 'organoid.Cones'], 'rods':['human.Mature_Rods', 'organoid.Rods'],
		'bipolars':['human.Mature_Bipolars', 'organoid.Bipolar_Cells'], 'mullers':['human.Mature_Mullers', 'organoid.Muller_Glia'],
		'horizontals':['human.Mature_Horizontals', 'organoid.Amacrine_Horizontal_Cells'], 'amacrines':['human.Mature_Amacrines', 'organoid.Amacrine_Horizontal_Cells'],
		'early_progenitor':['human.Early_Progenitors', 'organoid.Early_RPCs'], 'late_progenitor':['human.Late_Progenitors', 'organoid.Late_RPCs']}

# cell_class_dict = {'cones':['human.Mature_Cones', 'organoid.Cones']}
d_values_wanted = ['100', '500']
size_values_wanted = ['500', '2000'] ##does none work for homer
bed_suffixes = ['.macs2_q0.01_summits.bed','.macs2_q0.000001_summits.bed']
tag_dir_suffix = '.tag_dir'
# get_correlation_info_homer_pairwise(cell_class_dict, d_values_wanted, size_values_wanted, bed_suffixes, tag_dir_suffix)
##rpt for none size
d_values_wanted = ['100']
size_values_wanted = ['none']
# get_correlation_info_homer_pairwise(cell_class_dict, d_values_wanted, size_values_wanted, bed_suffixes, tag_dir_suffix)


##step8. formatting homer files for human/organoid comparision -- not used
annotated_files = ['annotated.collapsed_cc.100d.2000size.macs2_q0.000001_summits.txt', 'annotated.collapsed_cc.100d.2000size.macs2_q0.01_summits.txt', 
		'annotated.collapsed_cc.100d.500size.macs2_q0.000001_summits.txt', 'annotated.collapsed_cc.100d.500size.macs2_q0.01_summits.txt', 
		'annotated.collapsed_cc.500d.2000size.macs2_q0.000001_summits.txt', 'annotated.collapsed_cc.500d.2000size.macs2_q0.01_summits.txt', 
		'annotated.collapsed_cc.500d.500size.macs2_q0.000001_summits.txt', 'annotated.collapsed_cc.500d.500size.macs2_q0.01_summits.txt']
ann_files_for_peak_summary= ['annotated.collapsed_cc.100d.2000size.macs2_q0.000001_summits.txt', 'annotated.collapsed_cc.100d.500size.macs2_q0.000001_summits.txt', 
		'annotated.collapsed_cc.500d.2000size.macs2_q0.000001_summits.txt', 'annotated.collapsed_cc.500d.500size.macs2_q0.000001_summits.txt']
summary_suffix = '.average_for_scatterplot.txt'
summary_suffix_two = '.diff_peaks.txt'
##get peaks on regular chromosomes, and get average in human/org
# make_summary_files_for_scatterplots(annotated_files, summary_suffix)
# get_interesting_peaks(ann_files_for_peak_summary, summary_suffix_two)




##step 9. compare human/organoid data using homer getDifferentialPeaks
human_ccs = ['human.AC_HC_GC_Precursors', 'human.Developing_Amacrines', 'human.Developing_Bipolars', 'human.Developing_Cones', 
	'human.Developing_Ganglions', 'human.Developing_Horizontals', 'human.Developing_Rods', 'human.Early_Progenitors', 
	'human.Ganglion_Precursors', 'human.Late_Progenitors', 'human.Mature_Amacrines', 'human.Mature_Bipolars', 
	'human.Mature_Cones', 'human.Mature_Horizontals', 'human.Mature_Mullers', 'human.Mature_Rods', 'human.Photoreceptor_Bipolar_Precursors']
org_ccs =['organoid.AC_HC_Precursors', 'organoid.Cones', 'organoid.Early_RPCs', 'organoid.Late_RPCs', 'organoid.PR_BC_Precursors', 
	'organoid.RGCs', 'organoid.Rods', 'organoid.Amacrine_Horizontal_Cells', 'organoid.Bipolar_Cells', 'organoid.Muller_Glia', 'organoid.Developing_RGCs']
human_combined = 'human_combined'
org_combined = 'organoid_combined'
q_values_req = ['0.01', '0.000001']
d_values_req = ['100', '500']
size_values_wanted = ['500', '2000']
p_values_req = ['0.01', '0.001']
fc_values_req = ['3', '2']
diff_peak_prefix_human_org = 'diff_peaks.human_organoid.'
diff_peak_prefix_org_human = 'diff_peaks.organoid_human.'
##combine individul tag dir for human/org
# combine_tag_dirs(human_ccs, human_combined)
# combine_tag_dirs(human_ccs, org_combined)

##run getDifferentialPeaks
# get_diff_peaks_diff_sizes(diff_peak_prefix_human_org, diff_peak_prefix_org_human, d_values_req, q_values_req, size_values_wanted, human_combined, org_combined)
# get_diff_peaks_diff_sig(diff_peak_prefix_human_org, diff_peak_prefix_org_human, d_values_req, q_values_req, p_values_req, fc_values_req, human_combined, org_combined)


##try again but add sizes for get diff peaks
q_values_req = ['0.000001']
d_values_req = ['100']
size_values_wanted = ['500', '2000']
p_values_req = ['0.05']
fc_values_req = ['1.5']
# get_diff_peaks_diff_sig2(diff_peak_prefix_human_org, diff_peak_prefix_org_human, d_values_req, q_values_req, p_values_req, fc_values_req, size_values_wanted, human_combined, org_combined)



##step 10. compare adult human/organoid data using homer getDifferentialPeaksReplicates
## homer doesn't work with the module? so need to install independantley including a conda env
human_samples = ['hu5_bulk', 'hu7_bulk', 'hu8_bulk']
org_samples = ['28-1_bulk', '28-2_bulk']
diff_peak_rep_prefix_human_org = 'diff_peaks_rep.adult_human_28wk_organoid.'
diff_peak_rep_prefix_org_human = 'diff_peaks_rep.28wk_organoid_adult_human.'
diff_peak_rep_prefix_human_org_all = 'diff_peaks_rep.adult_human_28wk_organoid.all_peaks.'
diff_peak_rep_prefix_org_human_all = 'diff_peaks_rep.28wk_organoid_adult_human.all_peaks.'
diff_peak_rep_bg_prefix_human_org = 'diff_peaks_rep.bg.adult_human_28wk_organoid.'
diff_peak_rep_bg_prefix_org_human = 'diff_peaks_rep.bg.28wk_organoid_adult_human.'
diff_peak_rep_bg_prefix_human_org_all = 'diff_peaks_rep.bg.adult_human_28wk_organoid.all_peaks.'
diff_peak_rep_bg_prefix_org_human_all = 'diff_peaks_rep.bg.28wk_organoid_adult_human.all_peaks.'
# get_diff_peaks_replicates(diff_peak_rep_prefix_human_org, diff_peak_rep_prefix_org_human, human_samples, org_samples, 'no', 'input')
# get_diff_peaks_replicates(diff_peak_rep_prefix_human_org_all, diff_peak_rep_prefix_org_human_all, human_samples, org_samples, 'yes', 'input')
# get_diff_peaks_replicates(diff_peak_rep_bg_prefix_human_org, diff_peak_rep_bg_prefix_org_human, human_samples, org_samples, 'no', 'background')
# get_diff_peaks_replicates(diff_peak_rep_bg_prefix_human_org_all, diff_peak_rep_bg_prefix_org_human_all, human_samples, org_samples, 'yes', 'background')


##step 11 compare results from 10 against cell classes
mature_ccs = ['human.Mature_Amacrines', 'human.Mature_Bipolars', 'human.Mature_Cones', 'human.Mature_Horizontals', 
	'human.Mature_Mullers', 'human.Mature_Rods', 'organoid.Cones', 'organoid.RGCs', 'organoid.Rods', 
	'organoid.Amacrine_Horizontal_Cells', 'organoid.Bipolar_Cells', 'organoid.Muller_Glia']
homer_peak_suffix = '.homer_peaks.default.txt'
homer_peak_bed_suffix = '.homer_peaks.default.bed'
diff_peak_rep_files = ['diff_peaks_rep.28wk_organoid_adult_human.all_peaks.balanced.txt', 'diff_peaks_rep.28wk_organoid_adult_human.balanced.txt', 
		'diff_peaks_rep.adult_human_28wk_organoid.all_peaks.balanced.txt', 'diff_peaks_rep.adult_human_28wk_organoid.balanced.txt', 
		'diff_peaks_rep.bg.28wk_organoid_adult_human.all_peaks.balanced.txt', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.txt', 
		'diff_peaks_rep.bg.adult_human_28wk_organoid.all_peaks.balanced.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.txt']
diff_peak_rep_beds = [i.rsplit('.', 1)[0] + '.bed' for i in diff_peak_rep_files]

##get homer peaks for cell classes
# find_peaks_cell_classes(mature_ccs, homer_peak_suffix)

##make bed files from homer getDifferentialPeaksReplicates results and cell class findpeaks
# make_beds_from_diff_peaks_file(diff_peak_rep_files)
# make_beds_from_find_peaks_files(mature_ccs, homer_peak_suffix)

##compare cell class peaks with getDifferentialPeaksReplicates results
# bt_int_homer_beds_vs_cc_beds(diff_peak_rep_beds, mature_ccs, homer_peak_bed_suffix)

def format_int_files_for_upset(infiles):
	for infile in infiles:
		outfile = infile.rsplit('.', 1)[0] + '.upset.txt'
		with open(outfile, 'w') as out_fh, open(infile, 'r') as in_fh:
			out_fh.write(delim.join(['peak', 'cell_class', 'count']) + '\n')
			for line in in_fh:
				line = line.strip('\n').split(delim)
				peak = '_'.join(line[:3])
				cc = line[3].rsplit('.',3)[0]
				count = line[4]
				if int(count) > 1:
					count = '1'
				out_fh.write(delim.join([peak, cc, count]) + '\n')


##filter background files by fc and adj pvalue and make bed
files_to_filter = ['diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.txt']
#logfc/p value pairs, first is default
params_to_filter = [[1.0, 0.05], [2.0, 0.001], [3.0, 0.001]]
# filter_and_make_beds_from_diff_peaks_file(files_to_filter, params_to_filter)

##redo bt intersect
diff_peak_rep_beds_cutoffs = ['diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc1_p0.05.bed', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc2_p0.001.bed', 
		'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc3_p0.001.bed', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc1_p0.05.bed', 
		'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc2_p0.001.bed', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc3_p0.001.bed']
# bt_int_homer_beds_vs_cc_beds(diff_peak_rep_beds_cutoffs, mature_ccs, homer_peak_bed_suffix)

##reformat for graphing
files_to_reformat = ['diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc1_p0.05.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc2_p0.001.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.28wk_organoid_adult_human.balanced.lfc3_p0.001.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc1_p0.05.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc2_p0.001.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.balanced.lfc3_p0.001.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.28wk_organoid_adult_human.all_peaks.balanced.bt_int.mature_cc.txt', 'diff_peaks_rep.bg.adult_human_28wk_organoid.all_peaks.balanced.bt_int.mature_cc.txt']
format_int_files_for_upset(files_to_reformat)








