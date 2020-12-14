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



##for getting stranded use bamCoverage from deepTools:
##to setup:
module load local_python/3.7.6
conda create --name deeptools
source activate deeptools
conda install -c bioconda deeptools

##from now on:
qsub -Iq cdbrmq -l mem=200gb,ncpus=20 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load local_python/3.7.6
source activate deeptools
module load biobuilds
'''

##parameters
delim = '\t'
thread_number = '20'
##programs
sinto = '/home/atimms/programs/sinto/scripts/sinto'
bedtools = '/home/atimms/programs/bedtools2.28/bin/bedtools'
bg_to_bw = '/home/atimms/programs/bedGraphToBigWig'

##methods
def run_sinto_method(in_file, work_dir):
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip('\n').rstrip('\r').split(delim)
			sample_name = line[0]
			bam_file = line[1]
			cells_file = line[2]
			directory = './' + sample_name
			if not os.path.exists(sample_name):
				os.makedirs(sample_name)
			os.chdir(directory)
			run_sinto = subprocess.Popen(['sinto', 'filterbarcodes', '-b', bam_file,  '-c', cells_file, '-p', thread_number])
			run_sinto.wait()
			os.chdir(work_dir)

def merge_bams(collapsed_dict, w_dir):
	for sample in collapsed_dict:
		merged_bam = sample + '.bam'
		input_bams = ['I=' + b for b in collapsed_dict[sample]]
		print(merged_bam, input_bams)
		picard_md = subprocess.Popen(['picard', 'MergeSamFiles', 'CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'] + input_bams + ['O=' + merged_bam, 'TMP_DIR=' + w_dir])
		picard_md.wait()


def make_bigwig_files_rpkm(in_file):
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip('\n').rstrip('\r').split(delim)
			sample = line[0]
			bam = line[1]
			out_for_bw = sample + '.fwd.bigwig'
			out_rev_bw = sample + '.rev.bigwig'
			print(sample, bam)
			if os.path.exists(bam):
				##index file
				samtools_index = subprocess.Popen(['samtools', 'index', bam])
				samtools_index.wait()
				##get bigwigs
				bc_f = subprocess.Popen(['bamCoverage', '-b', bam, '-o', out_for_bw, '--filterRNAstrand', 'forward', '-p', '15', '--normalizeUsing', 'RPKM'])
				bc_f.wait()
				bc_r = subprocess.Popen(['bamCoverage', '-b', bam, '-o', out_rev_bw, '--filterRNAstrand', 'reverse', '-p', '15', '--normalizeUsing', 'RPKM'])
				bc_r.wait()
			else:
				print("bam file doesn't exist", bam)

def split_barcode_file(in_file):
	bc_dict = {}
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.replace('"', '')
				line = line.rstrip('\n').rstrip('\r').split(',')
				sample = line[0].split('_')[0]
				barcode = line[0].split('_')[1]
				cell_class = line[1].replace(' ', '_').replace('/', '_')
				bc_cc = barcode + ',' + cell_class
				if sample in bc_dict:
					bc_dict[sample].append(bc_cc)
				else:
					bc_dict[sample] = [bc_cc]
	for s in bc_dict:
		print(s, len(bc_dict[s]))
		outfile = 'org_scRNA.' + s + '.cell_types.txt'
		with open(outfile, "w") as out_fh:
			for info in bc_dict[s]:
				line_out = info.split(',')
				out_fh.write(delim.join(line_out) + '\n')


def merge_bam_make_bigwigs(in_file, work_dir):
	bam_dict = {}
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			sample = line[0]
			bam = line[1]
			if sample in bam_dict:
				bam_dict[sample].append(bam)
			else:
				bam_dict[sample] = [bam]
	##merge bams
	merge_bams(bam_dict, work_dir)
	##make file with merged bams/sample to make bigwigs
	file_for_making_bigwig = in_file.rsplit('.', 1)[0] + '.make_bigwigs.temp.txt'
	with open(file_for_making_bigwig, "w") as bw_fh:
		for s in bam_dict:
			line_out = delim.join([s, s + '.bam'])
			bw_fh.write(line_out + '\n')
	##make bigwigs
	make_bigwig_files_rpkm(file_for_making_bigwig)








##run methods

##info/files wtc
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_scRNA_make_bigwigs_1220'
os.chdir(working_dir)
pb_info = 'all_pseudobulk_info.txt'
organoid_barcode_info = 'org_barcodes_per_cellclass.csv'
organoid_sinto_info = 'org_scRNA.sinto_info.txt'
# organoid_sinto_info = 'test.txt'
human_merge_tp_info = 'human_scRNA.timepoint_cellclass.merge_info.txt'
human_merge_cc_info = 'human_scRNA.cellclass.merge_info.txt'
org_merge_tp_info = 'org_scRNA.timepoint_cellclass.merge_info.txt'
org_merge_cc_info = 'org_scRNA.cellclass.merge_info.txt'

##step1. make psedubulk bigwigs
# make_bigwig_files_rpkm(pb_info)


##step2. make individual bams for each sample/cell class
##copied human files from cherry_scRNA_make_bigwigs_1220

##split eric's file into barcodes per sample
# split_barcode_file(organoid_barcode_info)

##run sinto on all organoid samples
# run_sinto_method(organoid_sinto_info, working_dir)


##step3. merge bams and make bigs for different versions
##manually make files with 2 cols: name and bam file, will combined all bams with same name
# merge_bam_make_bigwigs(human_merge_tp_info, working_dir)
merge_bam_make_bigwigs(human_merge_cc_info, working_dir)
# merge_bam_make_bigwigs(org_merge_tp_info, working_dir)
# merge_bam_make_bigwigs(org_merge_cc_info, working_dir)




