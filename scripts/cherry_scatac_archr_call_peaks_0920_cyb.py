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
qsub -Iq cdbrmq -l mem=100gb,ncpus=5,walltime=10:00:00 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load local_python/3.7.6 
source activate macs2

##idr 
can just load module i.e.
module load idr/2.0.4
https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html

##homer
load modules:
qsub -Iq cdbrmq -l mem=100gb,ncpus=5,walltime=10:00:00 -P 19833a08-f6fb-4bea-8526-8a79069da87
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



##methods
def run_sinto_method(sample_name, bam_file, cells_file, w_dir):
	directory = './' + sample_name
	if not os.path.exists(sample_name):
		os.makedirs(sample_name)
	os.chdir(directory)
	run_sinto = subprocess.Popen([sinto, 'filterbarcodes', '-b', w_dir + bam_file,  '-c', w_dir + cells_file, '-p', thread_number])
	run_sinto.wait()
	os.chdir(w_dir)

def get_cell_info_ipscs(infile, outfile_suffix):
	##make dict of all cell barcode by sample
	cell_info_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.replace('"', '').rstrip().split(',')
				sample = line[1].split('#')[0]
				cell_barcode = line[1].split('#')[1]
				# print(sample, cell_barcode)
				if sample in cell_info_dict:
					cell_info_dict[sample].append(cell_barcode)
				else:
					cell_info_dict[sample] = [cell_barcode]
	##print out just for ipscs with ipsc as cell type
	out_files = []
	for s in cell_info_dict:
		print(s, len(cell_info_dict[s]))
		if s.startswith('i'):
			outfile = s + outfile_suffix
			out_files.append(outfile)
			with open(outfile, "w") as out_fh:
				for bc in cell_info_dict[s]:
					out_fh.write(delim.join([bc, 'IPSCs']) + '\n')
	return(out_files)

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

def run_macs2_on_all(samples):
	# '''
	##run macs2
	for sample in samples:
		bams = glob.glob(sample + '/' + '*bam')
		for bam in bams:
			outname = sample + '.' + bam.split('/')[1].split('.')[0] + '.macs2_q0.01'
			print(sample, bam, outname)
			##run macs2
			run_macs2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', outname, '-g', 'hs', '-q', '0.01', '--keep-dup', 'all'])
			run_macs2.wait()


def run_macs2_cell_class_collapsed(infile):
	macs2_dict = {}
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				name = '.'.join([line[0], line[2]]) + '.macs2_q0.01'
				bam = line[3]
				if name in macs2_dict:
					macs2_dict[name].append(bam)
				else:
					macs2_dict[name] = [bam]
	for outname in macs2_dict:
		# bams = ' '.join(macs2_dict[outname])
		bams = macs2_dict[outname]
		print(outname, bams)
		##run macs2
		# run_macs2 = subprocess.Popen(['macs2', 'callpeak', '-t', bams, '-f', 'BAMPE', '-n', outname, '-g', 'hs', '-q', '0.01', '--keep-dup', 'all'])
		run_macs2 = subprocess.Popen(['macs2', 'callpeak', '-t'] + bams + ['-f', 'BAMPE', '-n', outname, '-g', 'hs', '-q', '0.01', '--keep-dup', 'all'])
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

def get_correlation_info_homer(infile, d_values, size_values, bed_suffix):
	##read in info file
	bam_dict = {}
	ind_beds, comb_beds = [], []
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				name = '.'.join([line[1], line[2]])
				combined_bed = '.'.join([line[0], line[2]]) + bed_suffix
				ind_beds.append(name + bed_suffix)
				if combined_bed not in comb_beds:
					comb_beds.append(combined_bed)
				bam = line[3]
				if name in bam_dict:
					print(name, bam, 'name seen mutiple times')
				else:
					bam_dict[name] = bam
	##make tag dirs
	'''
	for td_name in bam_dict:
		bam = bam_dict[td_name]
		print(td_name, bam)
		##make tag dirs so can analyze
		make_tag_dirs(td_name, bam)
	'''
	##combine bed files
	'''
	print(len(ind_beds), len(comb_beds))
	for d in d_values:
		ind_merge_prefix = 'merged_ind_samples.'
		comb_merge_prefix = 'merged_combined.'
		merge_peaks(ind_merge_prefix, d, ind_beds, bed_suffix)
		merge_peaks(comb_merge_prefix, d, comb_beds, bed_suffix)
	'''
	##annotate files
	samples = bam_dict.keys()
	tag_dir_suffix = '.tag_dir'
	ind_annotate_prefix = 'annotated.ind_sample_beds.'
	comb_annotate_prefix = 'annotated.combined_beds.'
	# correlation_prefix = 'correlation.'
	for d in d_values:
		ind_merge_bed = 'merged_ind_samples.' + d + 'd' + bed_suffix
		comb_merge_bed = 'merged_combined.' + d + 'd' + bed_suffix
		outfile_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
		##annotate those peaks, and then make a correlation file
		annotate_peaks(samples, ind_merge_bed , d, tag_dir_suffix, ind_annotate_prefix, size_values, outfile_suffix)
		annotate_peaks(samples, comb_merge_bed , d, tag_dir_suffix, comb_annotate_prefix, size_values, outfile_suffix)


def format_correlation_files(samples_required, annotate_prefixes, annotate_suffix, d_values, sizes_req):
	##get list of infile
	infiles = []
	for d_value in d_values:
		for annotate_prefix in annotate_prefixes:
			for size_req in sizes_req:
				out_file_using_size = annotate_prefix + d_value + 'd.' + size_req + 'size' + annotate_suffix
				infiles.append(out_file_using_size)
	# print(len(infiles), infiles)
	##get list of samples required
	samples_req = []
	with open(samples_required, 'r') as samp_fh:
		lc1 = 0
		for line1 in samp_fh:
			lc1 += 1
			line1 = line1.rstrip().split(delim)	
			if lc1 > 1:
				sample = line1[0]
				samples_req.append(sample)
	print(samples_req)
	##format file by file
	for infile in infiles:
		out_prefix = 'correlation.' + samples_required.split('.')[0]
		out_file = out_prefix  + '.' + infile.split('.', 1)[1]
		print(infile, out_file)
		with open(out_file, 'w') as out_fh, open(infile, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.rstrip().split(delim)
				if lc == 1:
					sample_info = line[19:]
					# sample_info = [s.split('.')[0] for s in sample_info]
					sample_info = [s.split('.')[0] + '.' + s.split('.')[1] for s in sample_info]
					##get sample we want
					sample_index_req = []
					for s_req in samples_req:
						sample_index_req.append(sample_info.index(s_req))
					sample_info_wanted = [sample_info[i] for i in sample_index_req]
					out_fh.write(delim.join(['peak_id'] + sample_info_wanted) + '\n')

				else:
					counts = line[19:]
					counts_wanted = [counts[i] for i in sample_index_req]
					line_out = [line[0]] + counts_wanted
					out_fh.write(delim.join(line_out) + '\n')


def annotate_subset_homer(samples, d_values, size_values, bed_suffix):
	##annotate files
	tag_dir_suffix = '.tag_dir'
	ind_annotate_prefix = 'annotated.samples_gt20.ind_sample_beds.'
	comb_annotate_prefix = 'annotated.samples_gt20.combined_beds.'
	# correlation_prefix = 'correlation.'
	for d in d_values:
		ind_merge_bed = 'merged_ind_samples.' + d + 'd' + bed_suffix
		comb_merge_bed = 'merged_combined.' + d + 'd' + bed_suffix
		outfile_suffix = bed_suffix.rsplit('.',1)[0] + '.txt'
		##annotate those peaks, and then make a correlation file
		annotate_peaks(samples, ind_merge_bed , d, tag_dir_suffix, ind_annotate_prefix, size_values, outfile_suffix)
		annotate_peaks(samples, comb_merge_bed , d, tag_dir_suffix, comb_annotate_prefix, size_values, outfile_suffix)


def combine_corr_by_cell_class(infiles):
	for infile in infiles:
		outfile1 = infile.split('.')[0] + '.cell_class_collapsed_sum.' + infile.split('.',3)[3]
		outfile2 = infile.split('.')[0] + '.cell_class_collapsed_ave.' + infile.split('.',3)[3]
		print(infile, outfile1)
		##make dict with all cell classes and their index location
		cell_class_dict = {}
		with open(infile, 'r') as in_fh:
			first_line = in_fh.readline()
			first_line = first_line.rstrip().split(delim)
			sample_info = first_line[1:]
			##get list of all cell classes with origin outlined
			cc_list = []
			for s in sample_info:
				if s.startswith('H') or s.startswith('d'):
					cc = 'human.' + s.split('.')[1]
				else:
					cc = 'organoid.' + s.split('.')[1]
				cc_list.append(cc)
			cc_set = list(set(cc_list))
			# print(cc_list, cc_set)
			for cc in cc_set:
				##get all indices not just the first
				cc_indices = [i for i, x in enumerate(cc_list) if x == cc]
				if cc in cell_class_dict:
					print(cc, 'seen multiple times')
				else:
					cell_class_dict[cc] = cc_indices
		##print index info
		# for cell_class in cell_class_dict:
		# 	print(cell_class, cell_class_dict[cell_class])
		##combine values for each cell class and write outfile
		with open(outfile1, 'w') as out1_fh,open(outfile2, 'w') as out2_fh, open(infile, 'r') as in_fh:
			lc = 0
			header = ['peak_id'] + list(cell_class_dict.keys())
			out1_fh.write(delim.join(header) + '\n')
			out2_fh.write(delim.join(header) + '\n')
			for line in in_fh:
				lc += 1
				line = line.rstrip().split(delim)
				if lc > 1:
					peak_id = line[0]
					counts = line[1:]
					sums, aves = [], []
					for cell_class in cell_class_dict:
						counts_wanted = cell_class_dict[cell_class]
						counts_list = [float(counts[i]) for i in counts_wanted]
						sum_counts = sum(counts_list)
						ave_counts = sum_counts / len(counts_list)
						# print(peak_id, cell_class, counts_wanted, counts_list, sum_counts, ave_counts)
						sums.append(str(sum_counts))
						aves.append(str(ave_counts))
					out1_fh.write(delim.join([peak_id] + sums) + '\n')
					out2_fh.write(delim.join([peak_id] + aves) + '\n')






##run methods

##params etc
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scatac_archr_call_peaks_0920/'
os.chdir(working_dir)

##step1. make cell info files for sinto
ipsc_raw_info = 'org_all_w_ipscs.092420.csv'
org_raw_info = 'org_all.cell_class_info.092420.csv'
human_raw_info = 'human_all.cell_class_info.092520.csv'
cell_info_file_suffix = '.cell_types.txt'
##get files x3
'''
ipsc_cell_info_files = get_cell_info_ipscs(ipsc_raw_info, cell_info_file_suffix)
print(ipsc_cell_info_files)
org_cell_info_files = get_cell_info_org_human(org_raw_info, cell_info_file_suffix)
print(org_cell_info_files)
human_cell_info_files = get_cell_info_org_human(human_raw_info, cell_info_file_suffix)
print(human_cell_info_files)
all_cell_info_files = ipsc_cell_info_files + org_cell_info_files + human_cell_info_files
'''

##step2. run sinto
bam_dict = {'20wk1': '20wk.possorted.bam', '20wk2': '20wk_c5_1.possorted.bam', '28wk1': '28-1.possorted.bam', 
		'28wk2': '28-2.possorted.bam', '5wk1': '5wk.possorted.bam', '5wk2': '5wk_c5_1.possorted.bam', 
		'd113': 'd113.possorted.bam', 'd132': 'd132.possorted.bam', 'd53': 'd53.possorted.bam', 
		'd59': 'd59.possorted.bam', 'd74': 'd74.possorted.bam', 'd78': 'd78.possorted.bam', 
		'Hu5': 'hu5.possorted.bam', 'Hu7': 'hu7.possorted.bam', 'Hu8': 'hu8.possorted.bam', 
		'ipsc1': 'IPSC_c4.possorted.bam', 'ipsc2': 'IPSC_c5_1.possorted.bam'}
# all_cell_info_files = ['ipsc2.cell_types.txt'] ##for testing
'''
for cell_info_file in all_cell_info_files:
	sample = cell_info_file.split('.')[0]
	bam = bam_dict[sample]
	# print(sample, bam, cell_info_file)
	run_sinto_method(sample, bam, cell_info_file, working_dir)
'''

##step3. run macs2 many ways
samples = bam_dict.keys()
bam_file_combo_info = 'bam_file_and_cell_class_info.txt'
##on all samples and cell classes individually
# run_macs2_on_all(samples)
##on all cell classes
# run_macs2_cell_class_collapsed(bam_file_combo_info)

##step4. run homer i.e. make tag dirs, combine beds and annotate
d_values_wanted = ['100', '500']
size_values_wanted = ['500', '2000'] ##does none work for homer
bed_suffix = '.macs2_q0.01_summits.bed'
# get_correlation_info_homer(bam_file_combo_info, d_values_wanted, size_values_wanted, bed_suffix)

##step5. format correlation files
##make file with samples wanted in first column (txt before first . used in final files)
sample_file_all = 'all_samples_cell_classes.txt'
sample_file_gt20 = 'gt20_samples_cell_classes.txt'
sample_file_gt100 = 'gt100_samples_cell_classes.txt'
sample_files = [sample_file_all, sample_file_gt20, sample_file_gt100]
annotated_suffix = '.macs2_q0.01_summits.txt'
annotated_prefixes = ['annotated.combined_beds.', 'annotated.ind_sample_beds.']
##for testing
# sample_file = 'test.txt'
# d_values_wanted = ['100']
# size_values_wanted = ['500']
# annotated_prefixes = ['annotated.combined_beds.']
# for sample_file in sample_files:
# 	format_correlation_files(sample_file, annotated_prefixes, annotated_suffix, d_values_wanted, size_values_wanted)


##step6. repeat annotation and formatting correlation on file with > 20 cells
samples_gt_20_cells = ['d132.AC_HC_GC_Precursors', 'Hu8.Mature_Cones', 'd74.Developing_Rods', 'd59.Developing_Rods', 
		'd53.Developing_Rods', 'd78.Developing_Rods', 'd113.AC_HC_GC_Precursors', 'd132.Developing_Ganglions', 'd78.Early_Progenitors', 
		'd53.Developing_Ganglions', 'd53.Mature_Rods', '20wk2.Early_Progenitors', '5wk1.Amacrine_Horizontal_Precursors', 
		'28wk1.Early_Progenitors', 'Hu7.Mature_Horizontals', '28wk1.Amacrine_Cells', 'Hu5.Mature_Horizontals', 
		'5wk1.Bipolar_Photoreceptor_Precursors', 'd74.Early_Progenitors', '5wk2.Amacrine_Horizontal_Precursors', 'd53.AC_HC_GC_Precursors', 
		'28wk2.Early_Progenitors', 'd132.Developing_Horizontals', '20wk1.Early_Progenitors', 'Hu5.Mature_Cones', 'd113.Developing_Ganglions', 
		'Hu7.Mature_Cones', 'Hu8.Mature_Horizontals', '28wk1.Horizontal_Cells', '5wk2.Bipolar_Photoreceptor_Precursors', 'd59.Early_Progenitors', 
		'Hu5.Mature_Amacrines', '28wk2.Amacrine_Cells', 'd132.Developing_Cones', '28wk2.Horizontal_Cells', 'd113.Developing_Horizontals', 
		'28wk1.Late_Progenitors', 'd59.Photoreceptor_Bipolar_Precursors', 'Hu8.Mature_Amacrines', '20wk1.Bipolar_Cells', '20wk2.Late_Progenitors', 
		'5wk1.RPE_Early_Progenitors', 'd132.Developing_Bipolars', 'Hu7.Mature_Amacrines', '5wk2.Rods', 'd53.Late_Progenitors', 'd113.Developing_Cones', 
		'20wk1.Horizontal_Cells', '20wk2.Horizontal_Cells', '20wk2.Bipolar_Cells', '20wk1.Late_Progenitors', 'Hu7.Mature_Mullers', 
		'28wk2.Late_Progenitors', 'd74.Photoreceptor_Bipolar_Precursors', '28wk1.Bipolar_Cells', '5wk2.RPE_Early_Progenitors', 'Hu8.Mature_Mullers', 
		'Hu5.Mature_Mullers', '28wk2.Bipolar_Cells', 'd113.Developing_Bipolars', '5wk1.Retinal_Ganglion_Cells', '20wk1.Amacrine_Cells', 
		'5wk1.Amacrine_Horizontal_Ganglion_Precursors', 'Hu5.Mature_Bipolars', '5wk2.Amacrine_Horizontal_Ganglion_Precursors', '28wk1.Muller_Glia', 
		'5wk2.Retinal_Ganglion_Cells', '5wk1.Rods', '20wk2.Amacrine_Cells', 'Hu7.Mature_Bipolars', '28wk2.Muller_Glia', 'Hu8.Mature_Bipolars', 
		'd53.Ganglion_Precursors', 'd132.Developing_Amacrines', 'd78.Photoreceptor_Bipolar_Precursors', '20wk2.Cones', '28wk2.Cones', 
		'd74.AC_HC_GC_Precursors', '20wk2.Muller_Glia', '20wk1.Muller_Glia', 'd59.AC_HC_GC_Precursors', 'd113.Developing_Amacrines', 
		'28wk1.Cones', 'd74.Developing_Ganglions', '20wk1.Cones', 'd78.Developing_Ganglions', 'd59.Developing_Ganglions', 'd78.AC_HC_GC_Precursors', 
		'd74.Late_Progenitors', 'd132.Late_Progenitors', 'ipsc2.IPSCs', 'd113.Late_Progenitors', 'd59.Late_Progenitors', 'Hu5.Mature_Rods', 
		'd113.Developing_Rods', 'Hu7.Mature_Rods', '20wk1.Rods', '28wk2.Rods', '28wk1.Rods', '20wk2.Rods', 'd132.Developing_Rods', 
		'd78.Late_Progenitors', 'Hu8.Mature_Rods', 'd53.Early_Progenitors', '5wk1.Early_Progenitors', 'ipsc1.IPSCs', '5wk2.Early_Progenitors']
##run homer using subest
d_values_wanted = ['100', '500']
size_values_wanted = ['500', '2000']
# annotate_subset_homer(samples_gt_20_cells, d_values_wanted, size_values_wanted, bed_suffix)
sample_file_gt20 = 'gt20_samples_cell_classes.txt'
sample_file_gt100 = 'gt100_samples_cell_classes.txt'
sample_files = [sample_file_gt20, sample_file_gt100]
annotated_suffix = '.macs2_q0.01_summits.txt'
annotated_prefixes = ['annotated.samples_gt20.ind_sample_beds.', 'annotated.samples_gt20.combined_beds.']
annotated_prefixes = ['annotated.samples_gt20.combined_beds.']
# for sample_file in sample_files:
# 	format_correlation_files(sample_file, annotated_prefixes, annotated_suffix, d_values_wanted, size_values_wanted)


##step7. combine counts and look for correlation in group
ind_correlation_files = ['correlation.all_samples_cell_classes.combined_beds.100d.2000size.macs2_q0.01_summits.txt', 
		'correlation.all_samples_cell_classes.combined_beds.100d.500size.macs2_q0.01_summits.txt', 
		'correlation.all_samples_cell_classes.combined_beds.500d.2000size.macs2_q0.01_summits.txt', 
		'correlation.all_samples_cell_classes.combined_beds.500d.500size.macs2_q0.01_summits.txt']
# ind_correlation_files = ['test.test.test.100d.2000size.macs2_q0.01_summits.txt']
combine_corr_by_cell_class(ind_correlation_files)
