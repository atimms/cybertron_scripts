#!/usr/bin/env python
import subprocess
import os
import glob

'''
info....

##for getting stranded use bamCoverage from deepTools:
##to setup:
module load local_python/3.7.6
conda create --name deeptools
source activate deeptools
conda install -c bioconda deeptools

##from now on:
qsub -Iq cdbrmq -l mem=200gb,ncpus=20,walltime=10:00:00 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load biobuilds
module load local_python/3.7.6
source activate deeptools

'''

##parameters
delim = '\t'
thread_number = '20'
##programs
sinto = '/home/atimms/programs/sinto/scripts/sinto'
bedtools = '/home/atimms/programs/bedtools2.28/bin/bedtools'
bg_to_bw = '/home/atimms/programs/bedGraphToBigWig'


def make_bigwig_files(in_file):
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip('\n').rstrip('\r').split(delim)
			sample = line[0]
			bam = line[1]
			out_bw = sample + '.bigwig'
			print(sample, bam)
			if os.path.exists(bam):
				##index file
				samtools_index = subprocess.Popen(['samtools', 'index', bam])
				samtools_index.wait()
				##get bigwigs
				bc_f = subprocess.Popen(['bamCoverage', '-b', bam, '-o', out_bw, '-p', '18', '--normalizeUsing', 'RPKM'])
				bc_f.wait()

			else:
				print("bam file doesn't exist", bam)


##run methods, 3x sets of data
pb_working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_make_archr_bigwigs_1020/pseudobulk'
pb_info = 'pseudobulk.info.txt'
sc_working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_make_archr_bigwigs_1020/sample_cellclass'
sc_info = 'sample_cellclass.info.txt'
tc_working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_make_archr_bigwigs_1020/timepoint_cellclass'
tc_info = 'timepoint_cellclass.info.txt'

##test
# os.chdir('/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cherry_make_archr_bigwigs_1020')
# make_bigwig_files('test.txt')


##pseudobulk
# os.chdir(pb_working_dir)
# make_bigwig_files(pb_info)
##sample_cellclass
# os.chdir(sc_working_dir)
# make_bigwig_files(sc_info)
##timepoint_cellclass
os.chdir(tc_working_dir)
make_bigwig_files(tc_info)






