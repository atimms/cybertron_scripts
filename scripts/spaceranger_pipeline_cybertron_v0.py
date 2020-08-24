#!/usr/bin/env python
import subprocess
import os

##parameters
delim = '\t'
spaceranger = '/home/atimms/programs/spaceranger-1.1.0/spaceranger'
grc38_ref = '/home/atimms/ngs_data/references/10x/refdata-gex-GRCh38-2020-A'

##methods
def make_dict_from_info_file(in_file):
	idict = {}
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				sample = line[0]
				fqs = line[1]
				image = line[2]
				slide = line[3]
				area = line[4]
				idict[sample] = [fqs, image, slide, area]
	return(idict)

def run_spaceranger_count(infodict, refdir):
	for sample in infodict:
		#cellranger count --id=Mut2-6 --transcriptome=/gpfs/home/atimms/ngs_data/references/10x/refdata-cellranger-GRCh38-3.0.0 --fastqs=/home/atimms/ngs_data/misc/cherry_10x_organoid_0120/MacTel_organoid_wk5_rerun_done/outs/fastq_path/HHT75BGXC/Mut2-6,/home/atimms/ngs_data/misc/cherry_10x_organoid_0120/MacTel_Organoid_5wk_done/outs/fastq_path/HFLF2BGXC/Mut2-6 --sample=Mut2-6 --expect-cells=1000
		fqs = infodict[sample][0]
		image = infodict[sample][1]
		slide = infodict[sample][2]
		area = infodict[sample][3]
		sr_count = subprocess.Popen([spaceranger, 'count', '--id=' + sample, '--transcriptome=' + refdir, '--fastqs=' + fqs, '--sample=' + sample, '--image=' + image, '--slide=' + slide, '--area=' + area, '--localcores=10', '--localmem=100'])
		sr_count.wait()


def spaceranger_master_no_aggr(work_dir, infile, ref_dir):
	os.chdir(work_dir)
	info_dict = make_dict_from_info_file(infile)
	for s in info_dict:
		print(s, info_dict[s])
	run_spaceranger_count(info_dict, ref_dir)


##run methods

##kim data 0820
working_dir = '/home/atimms/ngs_data/cellranger/kim_visium_0820'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
info_file = 'kim_visium_0820.txt'
spaceranger_master_no_aggr(working_dir, info_file, grc38_ref)


