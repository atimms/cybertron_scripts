#!/usr/bin/env python
import subprocess
import os
import glob

'''
info....

##for sinto part i.e. getting bam files 
##go to worker node and load biobuilds and python3 and open env 
qsub -Iq cdbrmq -l mem=200gb,ncpus=20,walltime=10:00:00 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load biobuilds
module load local_python/3.7.6
source activate sinto
##program is now at:
/home/atimms/programs/sinto/scripts/sinto

##for macs2
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load local_python/3.7.6 
source activate macs2

##idr 
can just load module i.e.
module load idr/2.0.4
https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html

##homer
load modules:
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

##files
##copied from /active/cherry_t/NGS_Data/ChIP_input_for_bkgrnd/Hu1-ret-input-peakcalling-R1_trimmed-fixchr.bed.gz
control_bed = 'Hu1-ret-input-peakcalling-R1_trimmed-fixchr.bed'
chrom_sizes = '/home/atimms/ngs_data/references/10x/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.chrom_sizes'

##methods

def call_macs2_on_original_bams(sample_list):
	print(sample_list)
	# '''
	##run macs2
	for sample in sample_list:
		bam = sample + '.possorted.bam'
		name = sample + '_bulk'
		##outfile names
		out_bampe_keepdups_1 =  name + '.macs2.bampe_p1e-2_keepdups'
		out_bampe_keepdups_2 =  name + '.macs2.bampe_p1e-5_keepdups'
		out_bampe_keepdups_3 =  name + '.macs2.bampe_p1e-10_keepdups'
		out_bampe_keepdups_4 =  name + '.macs2.bampe_q0.01_keepdups'
		out_bampe_keepdups_5 =  name + '.macs2.bampe_q0.000001_keepdups'
		##run macs2 from bam files
		# run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_2, '-g', 'hs', '-p', '1e-5', '--keep-dup', 'all'])
		# run_macs2_2.wait()		
		# run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_3, '-g', 'hs', '-p', '1e-10', '--keep-dup', 'all', '--tempdir', '.'])
		# run_macs2_2.wait()
		# run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_4, '-g', 'hs', '-q', '0.01', '--keep-dup', 'all', '--tempdir', '.'])
		# run_macs2_2.wait()
		run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_5, '-g', 'hs', '-q', '0.000001', '--keep-dup', 'all', '--tempdir', '.'])
		run_macs2_2.wait()
	# '''

def make_tag_dirs(name, bam):
	outdir = name + '.tag_dir'
	mk_tag_dir = subprocess.Popen(['makeTagDirectory', outdir, bam])
	mk_tag_dir.wait()

def merge_peaks_original_bams(out_prefix, d_value, peak_files, out_suffix):
	out_file = out_prefix + d_value + 'd' + out_suffix
	print(len(peak_files), out_file)
	with open(out_file, 'w') as out_fh:
		# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
		mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value] + peak_files, stdout=out_fh)
		mk_tag_dir.wait()

def annotate_peaks_original_bams(samples, peak_prefix, peak_suffix, d_value, tag_suffix, out_prefix, sizes_req):
	peak_file = peak_prefix + d_value + 'd' + peak_suffix
	tag_dirs = []
	##get list of tag_dirs
	for sample in samples:
		tag_dir = sample + '_bulk' + tag_suffix
		tag_dirs.append(tag_dir)
	out_file_no_size = out_prefix + d_value + 'd' + peak_suffix
	print(len(tag_dirs), peak_file, out_file_no_size)
	with open(out_file_no_size, 'w') as out_ns_fh:
		mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-d'] + tag_dirs, stdout=out_ns_fh)
		mk_tag_dir.wait()
	for size_req in sizes_req:
		out_file_using_size = out_prefix + d_value + 'd.' + size_req + 'size' + peak_suffix
		with open(out_file_using_size, 'w') as out_fh:
			#annotatePeaks.pl pu1peaks.txt mm8 -size 400 -d Macrophage-PU.1/ Bcell-PU.1/ > output.txt
			##alternatives
			# mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg19', '-size', size_req], stdout=out_fh)
			# mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg19', '-noann', '-size', size_req, '-norm', '10000000', '-d'] + tag_dirs, stdout=out_fh)
			mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-size', size_req, '-d'] + tag_dirs, stdout=out_fh)
			mk_tag_dir.wait()

def compute_correlation_original_bams(annotate_prefix, annotate_suffix, d_value, sizes_req, out_prefix):
	##get list of infile
	infiles = []
	no_size_file = annotate_prefix + d_value + 'd' + annotate_suffix
	infiles.append(no_size_file)
	for size_req in sizes_req:
		out_file_using_size = annotate_prefix + d_value + 'd.' + size_req + 'size' + annotate_suffix
		infiles.append(out_file_using_size)
	# print(len(infiles), infiles)
	##format file by file
	for infile in infiles:
		out_file = out_prefix  + infile.split('.', 1)[1]
		print(infile, out_file)
		with open(out_file, 'w') as out_fh, open(infile, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.strip('\n').split(delim)
				if lc == 1:
					sample_info = line[19:]
					sample_info = [s.split('.')[0] for s in sample_info]
					out_fh.write(delim.join(['peak_id'] + sample_info + ['\n']))
					print(sample_info)
				else:
					line_out = [line[0]] + line[19:]
					out_fh.write(delim.join(line_out + ['\n']))

def homer_correlation_on_original_bams(samples, d_values, size_values, bed_suffixes, tag_dirs_needed):
	##make tag dirs
	if tag_dirs_needed == 'yes':
		for sample in samples:
			bam = sample + '.possorted.bam'
			name = sample + '_bulk'
			##make tag dirs so can analyze
			make_tag_dirs(name, bam)
	##for different parameters merge peaks and then annoate
	tag_dir_suffix = '.tag_dir'
	merge_prefix = 'merged.'
	annotate_prefix = 'annotated.'
	correlation_prefix = 'correlation.'
	for d in d_values:
		for bed_suffix in bed_suffixes:
			##merge peak files for annotating the peaks
			peak_beds = [s + '_bulk' + bed_suffix  for s in samples]
			merge_file_suffix = bed_suffix.rsplit('.', 1)[0] + '.txt'
			merge_peaks_original_bams(merge_prefix, d, peak_beds, merge_file_suffix)
			##annotate those peaks, and then make a correlation file
			annotate_peaks_original_bams(samples, merge_prefix, merge_file_suffix, d, tag_dir_suffix, annotate_prefix, size_values)
			##make file to use in r for heatmaps etc
			compute_correlation_original_bams(annotate_prefix, merge_file_suffix, d, size_values, correlation_prefix)

##run methods
# working_dir = '/home/atimms/ngs_data/misc/cherry_scatac_call_peaks_0820'
# os.chdir(working_dir)

##experiment one - correlation on original bams i.e. not split by cell class
##sample names
combined_sample_list = ['d53', 'd59', 'd74', 'd78', 'd113', 'd132', 'hu5', 'hu7', 'hu8', 
		'IPSC_c4', 'IPSC_c5_1', '5wk', '5wk_c5_1', '20wk', '20wk_c5_1', '28-1', '28-2']
human_sample_list = ['d53', 'd59', 'd74', 'd78', 'd113', 'd132', 'hu5', 'hu7', 'hu8']
org_sample_list = ['IPSC_c4', 'IPSC_c5_1', '5wk', '5wk_c5_1', '20wk', '20wk_c5_1', '28-1', '28-2']
##homer values
d_values = ['100', '200', '500']
size_values = ['500', '2000']
size_values_plus = ['500', '2000', 'none']
peak_bed_suffices = ['.macs2.bampe_p1e-10_keepdups_summits.bed', '.macs2.bampe_q0.01_keepdups_summits.bed']

# combined_sample_list = ['5wk_c5_1']

##call peaks using macs2
# call_macs2_on_original_bams(combined_sample_list)
##split 4 ways for speed
# call_macs2_on_original_bams(combined_sample_list[:3])
# call_macs2_on_original_bams(combined_sample_list[3:6])
# call_macs2_on_original_bams(combined_sample_list[6:9])
# call_macs2_on_original_bams(combined_sample_list[9:12])
# call_macs2_on_original_bams(combined_sample_list[12:15])
# call_macs2_on_original_bams(combined_sample_list[15:])

##use homer to look for correlation
# homer_correlation_on_original_bams(combined_sample_list, d_values, size_values, peak_bed_suffices, 'yes')
# homer_correlation_on_original_bams(combined_sample_list, d_values, size_values, peak_bed_suffices, 'no')

##reanalysis 1021 with additional samples
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scatac_corr_rpts_1021'
os.chdir(working_dir)
new_samples = ['12wk1', '12wk2', '12wk3']
combined_sample_list = ['d53', 'd59', 'd74', 'd78', 'd113', 'd132', 'hu5', 'hu7', 'hu8', 
		'IPSC_c4', 'IPSC_c5_1', '5wk', '5wk_c5_1', '20wk', '20wk_c5_1', '28-1', '28-2',
		'12wk1', '12wk2', '12wk3']
# peak_bed_suffices = ['.macs2.bampe_q0.01_keepdups_summits.bed', '.macs2.bampe_q0.000001_keepdups_summits.bed']
# peak_bed_suffices = ['.macs2.bampe_q0.01_keepdups_summits.bed']
peak_bed_suffices = ['.macs2.bampe_q0.000001_keepdups_summits.bed']

##just call macs2 on 3 new samples, just using q0.01 a
# call_macs2_on_original_bams(new_samples)
##and on all bams using q0.000001
# call_macs2_on_original_bams(combined_sample_list)

##use homer to look for correlation
# homer_correlation_on_original_bams(combined_sample_list, d_values, size_values, peak_bed_suffices, 'yes')
homer_correlation_on_original_bams(combined_sample_list, d_values, size_values, peak_bed_suffices, 'no')


