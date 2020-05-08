#!/usr/bin/env python
import sys
import subprocess
import os
import glob


##note
'''
for mscentipede....
need python 2
load modules:
module load biobuilds/2017.11
module load local_python/2.7.14
install packages:
conda install -c conda-forge cvxopt
to plot profile need to alter the python scripts plot_accessibility_profile.py
add the lines before import matplotlib.pyplot as plot:
import matplotlib
matplotlib.use('agg')

##for defcom
module load local_python/2.7.14
##to install
conda install scikit-learn
##in defcon
python setup.py test
python setup.py install


'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/misc/cherry_atac_analysis_1018'
os.chdir(working_dir)

##program
call_binding = '/home/atimms/programs/msCentipede/call_binding.py'
plt_accesability = '/home/atimms/programs/msCentipede/plot_accessibility_profile_adj.py'

def split_and_format_motif_bed(in_bed, out_file_suffix):
	motif_dict = {}
	with open(in_bed, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			tf = line[3].split("(")[0]
			info = line[:3] + [line[5], line[4]]
			if tf in motif_dict:
				motif_dict[tf].append(info)
			else:
				motif_dict[tf] = [info]
	for trans in motif_dict:
		print trans
		outfile = trans + out_file_suffix
		with open(outfile, "w") as out_fh:
			out_fh.write(delim.join(['Chr', 'Start', 'Stop', 'Strand', 'PwmScore']) + '\n')
			for inf in motif_dict[trans]:
				out_fh.write(delim.join(inf) + '\n')
		run_gzip = subprocess.Popen(['gzip', outfile])
		run_gzip.wait()

def run_mscentipede(motif_files, bam_files):
	##python call_binding.py --task learn test/CTCF_chr10_motifs.txt.gz test/Gm12878_Rep1.bam test/Gm12878_Rep2.bam
	for motif_file in motif_files:
		print(motif_file, bam_files)
		run_msc = subprocess.Popen(['python', call_binding, '--task', 'learn', '--protocol=ATAC_seq', '--batch', '500', motif_file] + bam_files)
		run_msc = subprocess.Popen(['python', call_binding, '--task', 'learn', '--protocol=ATAC_seq', motif_file] + bam_files)
		run_msc.wait()
		run_msc_infer = subprocess.Popen(['python', call_binding, '--task', 'infer', '--protocol=ATAC_seq', motif_file] + bam_files)
		run_msc_infer.wait()
		run_msc_plot = subprocess.Popen(['python', plt_accesability, '--protocol=ATAC_seq', motif_file])
		run_msc_plot.wait()

def get_motifs_in_peaks_bed(in_bed, out_bed, peaks_bed):
	with open(out_bed, "w") as out_fh: 
		hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', in_bed, '-b', peaks_bed, '-wa'], stdout=out_fh)
		hom_bt_intersect.wait()



def filter_bam_files(in_bams, final_bam_suffix):
	for in_bam in in_bams:
		sample = in_bam.split('.')[0]
		outbam = sample + final_bam_suffix
		##bwa and convert to bam
		st_sam_bam_pe = subprocess.Popen(['samtools', 'view', '-b', '-@', '15', '-q', '1', '-o', outbam, in_bam, 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX'])
		st_sam_bam_pe.wait()
		st_index = subprocess.Popen(['samtools', 'index', outbam])
		st_index.wait()


def graph_mscentipede(motif_files, plot_file):
	##python call_binding.py --task learn test/CTCF_chr10_motifs.txt.gz test/Gm12878_Rep1.bam test/Gm12878_Rep2.bam
	for motif_file in motif_files:
		run_msc_plot = subprocess.Popen(['python', plot_file, '--protocol=ATAC_seq', motif_file])
		run_msc_plot.wait()


##run methods
peak_bed = 'Hu-ret-ATAC_pooled_peaks_conservative_postidr.bed'
homer_motifs = 'homer_motifs_hg38_1018.bed'
homer_motifs_in_peaks = 'homer_motifs_hg38_in_peak_bed_1018.bed'
motif_bed_suffix = '.homer_hg38_1018.txt'
motif_bed_suffix = '_all.homer_hg38_1018.txt'
mscent_motifs = glob.glob('*' + motif_bed_suffix + '.gz')
bams_to_filter = glob.glob('*.bwa_mkdup.bam')
# bams_to_filter = ['Hu1_ret_ATAC.bwa_mkdup.bam']
filtered_bam_suffix = '.bwa_mkdup_filtered.bam'
bams_to_analyze = glob.glob('*.bwa_mkdup_filtered.bam')
# bams_to_analyze = ['Hu1_ret_ATAC.bwa_mkdup_filtered.bam']
##get motifs that intersect with peaks bed, better input for mscentipede?
# get_motifs_in_peaks_bed(homer_motifs, homer_motifs_in_peaks, peak_bed)

##split in 7 ts factors and reformat for mscent
# split_and_format_motif_bed(homer_motifs_in_peaks, motif_bed_suffix)
##and do that for all motifs
# split_and_format_motif_bed(homer_motifs, motif_bed_suffix)

##filter bam, rm q0 reads and non std chr
# filter_bam_files(bams_to_filter, filtered_bam_suffix)
##std way, all the tfs
# mscent_motifs = ['BORIS.homer_hg38_1018.txt.gz', 'CRE.homer_hg38_1018.txt.gz', 'CRX.homer_hg38_1018.txt.gz', 'MafA.homer_hg38_1018.txt.gz', 'Mef2d.homer_hg38_1018.txt.gz', 'Otx2.homer_hg38_1018.txt.gz', 'RORgt.homer_hg38_1018.txt.gz']
##using all motifs not just those in peaks
# mscent_motifs = ['BORIS_all.homer_hg38_1018.txt.gz']
##temp
# mscent_motifs = ['BORIS_test.homer_hg38_1018.txt.gz']
mscent_motifs = ['BORIS.homer_hg38_1018.txt.gz', 'CRX.homer_hg38_1018.txt.gz', 'MafA.homer_hg38_1018.txt.gz', 'Mef2d.homer_hg38_1018.txt.gz', 'Otx2.homer_hg38_1018.txt.gz', 'RORgt.homer_hg38_1018.txt.gz', 'BORIS_all.homer_hg38_1018.txt.gz']

##run mscentipede
# run_mscentipede(mscent_motifs, bams_to_analyze)

##test graphing
plot_acc_file = '/home/atimms/programs/msCentipede/plot_accessibility_profile_test.py'
graph_mscentipede(mscent_motifs, plot_acc_file)




