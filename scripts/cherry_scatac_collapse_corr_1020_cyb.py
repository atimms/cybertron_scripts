#!/usr/bin/env python
import subprocess
import os
import glob


'''
info....

##for merge peaks
module load biobuilds

##for macs2
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load local_python/3.7.6 
source activate macs2

##homer
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da87
module load biobuilds
module load homer/4.9.1

'''

##parameters
delim = '\t'
thread_number = '20'
##programs
sinto = '/home/atimms/programs/sinto/scripts/sinto'
genrich = '/home/atimms/programs/Genrich/Genrich'
bedtools = '/home/atimms/programs/bedtools2.28/bin/bedtools'
bg_to_bw = '/home/atimms/programs/bedGraphToBigWig'


##methods
def collapse_samples_by_cell_class(infile, bam_dir):
	coll_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				cell_class = line[2]
				sample = line[1]
				bam_file = bam_dir + sample + '/' + cell_class + '.bam'
				type_cc = line[0] + '.' + cell_class
				if type_cc in coll_dict:
					coll_dict[type_cc].append(bam_file)
				else:
					coll_dict[type_cc] = [bam_file]
	return(coll_dict)

def collapse_samples_by_cell_class_tp(infile, bam_dir):
	##extras is just samples not included in the just cell class analysis
	complete_dict, extras_dict = {}, {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				cell_class = line[2]
				sample = line[1]
				bam_file = bam_dir + sample + '/' + cell_class + '.bam'
				type_tp = line[0]
				type_cc = type_tp + '.' + cell_class
				if type_cc in complete_dict:
					complete_dict[type_cc].append(bam_file)
				else:
					complete_dict[type_cc] = [bam_file]
				if '_' in type_tp:
					if type_cc in extras_dict:
						extras_dict[type_cc].append(bam_file)
					else:
						extras_dict[type_cc] = [bam_file]			
	return(complete_dict, extras_dict)

def merge_bams(collapsed_dict, w_dir):
	for sample in collapsed_dict:
		merged_bam = sample + '.bam'
		input_bams = ['I=' + b for b in collapsed_dict[sample]]
		print(merged_bam, input_bams)
		picard_md = subprocess.Popen(['picard', 'MergeSamFiles', 'CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'] + input_bams + ['O=' + merged_bam, 'TMP_DIR=' + w_dir])
		picard_md.wait()

def run_macs2_on_all(samples):
	##run macs2
	for sample in samples:
		bam = sample + '.bam'
		outname = sample + '.macs2_q0.01'
		print(sample, bam, outname)
		##run macs2
		run_macs2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', outname, '-g', 'hs', '-q', '0.01', '--keep-dup', 'all'])
		run_macs2.wait()

def make_tag_dirs(name, bam):
	outdir = name + '.tag_dir'
	mk_tag_dir = subprocess.Popen(['makeTagDirectory', outdir, bam])
	mk_tag_dir.wait()

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
		out_file_using_size = out_prefix + d_value + 'd.' + size_req + 'size' + peak_suffix
		with open(out_file_using_size, 'w') as out_fh:
			mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-size', size_req, '-d'] + tag_dirs, stdout=out_fh)
			mk_tag_dir.wait()

def get_correlation_info_homer(samples, d_values, size_values, bed_suffix, merge_prefix, annotate_prefix):
	##make tag dirs
	'''
	for sample in samples:
		bam = sample + '.bam'
		print(sample, bam)
		##make tag dirs so can analyze
		make_tag_dirs(sample, bam)
	'''
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



##run methods

##params etc
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scatac_archr_collapse_corr_1020/'
os.chdir(working_dir)

##sample, cell class and bam info
ind_bam_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scatac_archr_call_peaks_0920/'
sample_cell_class_info = 'samples_cell_classes_gt20.txt'
sample_timepoint_cell_class_info = 'samples_timepint_cell_classes_gt20.txt'

##just collapse by cell class and if human or organiod

##get info
coll_cell_class_dict = collapse_samples_by_cell_class(sample_cell_class_info, ind_bam_dir)
## 1. merge bams
# merge_bams(coll_cell_class_dict, working_dir)
## 2. call peaks using macs2
samples = coll_cell_class_dict.keys()
# run_macs2_on_all(samples)
## 3. run homer i.e. make tag dirs, combine beds and annotate
d_values_wanted = ['100', '500']
size_values_wanted = ['500', '2000'] ##does none work for homer
bed_suffix = '.macs2_q0.01_summits.bed'
txt_suffix =  '.macs2_q0.01_summits.txt'
merged_peak_prefix = 'merged.collapsed_cc_only.'
annotated_peak_prefix = 'annotated.collapsed_cc_only.'
correlation_peak_prefix = 'correlation.'
# get_correlation_info_homer(samples, d_values_wanted, size_values_wanted, bed_suffix, merged_peak_prefix, annotated_peak_prefix)
## 4. format corrlation file for graphing
# format_correlation_all_samples(annotated_peak_prefix, txt_suffix, d_values_wanted, size_values_wanted, correlation_peak_prefix)

##just collapse by cell class and by human or organiod and organiod time point

##get info
coll_cell_class_tp_dict, coll_cell_class_tp_extras_dict = collapse_samples_by_cell_class_tp(sample_timepoint_cell_class_info, ind_bam_dir)
## 1. merge bams (just the extra samples)
# merge_bams(coll_cell_class_tp_extras_dict, working_dir)
## 2. call peaks using macs2 (just the extra samples)
extra_samples = coll_cell_class_tp_extras_dict.keys()
# run_macs2_on_all(extra_samples)
## 3. run homer i.e. make tag dirs, combine beds and annotate (all the samples)
all_samples = coll_cell_class_tp_dict.keys()
d_values_wanted = ['100', '500']
size_values_wanted = ['500', '2000'] ##does none work for homer
bed_suffix = '.macs2_q0.01_summits.bed'
txt_suffix =  '.macs2_q0.01_summits.txt'
merged_peak_prefix = 'merged.collapsed_cc_tp.'
annotated_peak_prefix = 'annotated.collapsed_cc_tp.'
correlation_peak_prefix = 'correlation.'
##extra samples just for make tag dirs, all for everything else
# get_correlation_info_homer(extra_samples, d_values_wanted, size_values_wanted, bed_suffix, merged_peak_prefix, annotated_peak_prefix)
# get_correlation_info_homer(all_samples, d_values_wanted, size_values_wanted, bed_suffix, merged_peak_prefix, annotated_peak_prefix)
## 4. format corrlation file for graphing
format_correlation_all_samples(annotated_peak_prefix, txt_suffix, d_values_wanted, size_values_wanted, correlation_peak_prefix)





