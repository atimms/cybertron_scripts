#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/exomes/working/kim_0220'
os.chdir(working_dir)

'''
load these:
module load java/1.8.0_121 
module load biobuilds/2016.11
'''


##programs and program files
laser = '/home/atimms/programs/LASER-2.03/laser'
pileup2seq = '/home/atimms/programs/LASER-2.03/pileup2seq/pileup2seq.py'
hgdp_site = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.site'
hgdp_bed = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.bed'
hgdp_geno = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.geno'
hgdp_coord = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.RefPC.coord'
fasta = '/home/atimms/ngs_data/references/hg19/human_g1k_v37.fasta'



hgdp_continent_dict = {'Adygei' : 'Europe', 'Balochi' : 'Central and South Asia', 'BantuKenya' : 'Africa', 
		'BantuSouthAfrica' : 'Africa', 'Basque' : 'Europe', 'Bedouin' : 'Middle East', 
		'BiakaPygmy' : 'Africa', 'Brahui' : 'Central and South Asia', 
		'Burusho' : 'Central and South Asia', 'Cambodian' : 'East Asia', 'Colombian' : 'America', 
		'Dai' : 'East Asia', 'Daur' : 'East Asia', 'Druze' : 'Middle East', 'French' : 'Europe', 
		'Han' : 'East Asia', 'Han-NChina' : 'East Asia', 'Hazara' : 'Central and South Asia', 
		'Hezhen' : 'East Asia', 'Italian' : 'Europe', 'Japanese' : 'East Asia', 
		'Kalash' : 'Central and South Asia', 'Karitiana' : 'America', 'Lahu' : 'East Asia', 
		'Makrani' : 'Central and South Asia', 'Mandenka' : 'Africa', 'Maya' : 'America', 
		'MbutiPygmy' : 'Africa', 'Melanesian' : 'Oceania', 'Miao' : 'East Asia', 
		'Mongola' : 'East Asia', 'Mozabite' : 'Middle East', 'Naxi' : 'East Asia', 
		'Orcadian' : 'Europe', 'Oroqen' : 'East Asia', 'Palestinian' : 'Middle East', 
		'Papuan' : 'Oceania', 'Pathan' : 'Central and South Asia', 'Pima' : 'America', 
		'Russian' : 'Europe', 'San' : 'Africa', 'Sardinian' : 'Europe', 'She' : 'East Asia', 
		'Sindhi' : 'Central and South Asia', 'Surui' : 'America', 'Tu' : 'East Asia', 
		'Tujia' : 'East Asia', 'Tuscan' : 'Europe', 'Uygur' : 'Central and South Asia', 
		'Xibo' : 'East Asia', 'Yakut' : 'East Asia', 'Yi' : 'East Asia', 'Yoruba' : 'Africa'}
# study_site_dict = {'101':'Site_Bogota', '103':'Site_Cali', '106':'Site_Cali', '107':'Site_Pereira'}
# select_continents = ['Europe', 'America', 'Africa']


##methods
def samtools_pileup(bam_files):
	'samtools mpileup -q 30 -Q 20 -f ../../LASER-resource/reference/hs37d5.fa -l HGDP_938.bed exampleBAM/NA12878.chrom22.recal.bam > NA12878.chrom22.pileup'
	##iterate over bamfiles
	for bam in bam_files:
		sample = bam.split('.')[0]
		pileup_file = sample + '.pileup'
		with open(pileup_file, 'w') as pu_fh:
			bcf_index = subprocess.Popen(['samtools', 'mpileup', '-q', '30', '-Q', '20', '-f', fasta, '-l', hgdp_bed, bam], stdout=pu_fh)
			bcf_index.wait()

def pileup_names_from_bams(bam_files):
	pileup_names = []
	for bam in bam_files:
		pileup_file = bam.split('.')[0] + '.pileup'
		pileup_names.append(pileup_file)
	return pileup_names

def combine_mpileup_to_seq_file(mpileup_files, out_prefix):
	'python pileup2seq.py -f hs37d5.fa -m ref.site -b target.bed -i example.id -o output A.pileup B.pileup C.pileup'
	run_pileup2seq = subprocess.Popen(['python', pileup2seq, '-f', fasta, '-m', hgdp_site, '-o', out_prefix] + mpileup_files)
	run_pileup2seq.wait()

def run_laser(geno_file, coord_file, file_prefix, pcs_req, rpts_req):
	run_laser_cmd = subprocess.Popen([laser, '-g', geno_file, '-c', coord_file, '-k', pcs_req, '-r', rpts_req, '-s', file_prefix.split('.')[0] + '.seq', '-o', file_prefix])
	run_laser_cmd.wait()

def combine_format_coord_files_for_r(sample_file_prefix, coord_file, output_file, populations_req):
	sample_file = sample_file_prefix + '.SeqPC.coord'
	with open(coord_file, "r") as ctl_fh, open(sample_file, "r") as samp_fh, open(output_file, 'w') as out_fh:
		##write header
		out_fh.write(delim.join(['continent', 'population', 'sample', 'pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'pc7', 'pc8', '\n']))
		line_count = 0
		for line in ctl_fh:
			line = line.strip('\n').split(delim)
			line_count += 1
			if line_count > 1:
				popID = line[0]
				contID = hgdp_continent_dict[popID]
				# print line[:4]
				if populations_req == 'na':
					out_fh.write(delim.join([contID] + line[:10] + ['\n']))
				else:
					if contID in populations_req:
						out_fh.write(delim.join([contID] + line[:10] + ['\n']))

		line_count = 0
		for line in samp_fh:
			line = line.strip('\n').split(delim)
			line_count += 1
			if line_count > 1:
				sample_id = line[0]
				pcs = line[6:14]
				# print sample_id, pcs
				out_fh.write(delim.join(['Sample', 'cblm', sample_id] + pcs + ['\n']))

##run methods

##file names etc
bams_for_analysis = glob.glob('*.bam')
project_prefix = 'exomes_0320'
hgdp_all_for_r = project_prefix + '.hgdp_all.txt'
##bam files
bams_for_analysis = ['B03_12.bwa_gatk.bam', 'B05_13.bwa_gatk.bam', 'B28_17.bwa_gatk.bam', 'B36_17.bwa_gatk.bam', 'B40_17.bwa_gatk.bam', 'M23_17.bwa_gatk.bam', 'N04_15.bwa_gatk.bam', 'NBB_1347.bwa_gatk.bam', 'NBB_13646.bwa_gatk.bam', 'NBB_1714.bwa_gatk.bam', 'NBB_4599.bwa_gatk.bam', 'NBB_4787.bwa_gatk.bam', 'NBB_4916.bwa_gatk.bam', 'NBB_5242.bwa_gatk.bam', 'NBB_5294.bwa_gatk.bam', 'NBB_5334.bwa_gatk.bam', 'NBB_5419.bwa_gatk.bam', 'NBB_5531.bwa_gatk.bam', 'NBB_5565.bwa_gatk.bam', 'NBB_5841.bwa_gatk.bam', 'NBB_5864.bwa_gatk.bam', 'NBB_5936.bwa_gatk.bam', 'NBB_5978.bwa_gatk.bam', 'NBB_6033.bwa_gatk.bam', 'U30_17.bwa_gatk.bam', 'U38_18.bwa_gatk.bam']

##cblm exomes
##make mpileup files
# samtools_pileup(bams_for_analysis)

##covert pileup and combine in seq file
##get list of pileup from bams and then combine into seq file
# pileup_files = pileup_names_from_bams(bams_for_analysis)
# print pileup_files
# combine_mpileup_to_seq_file(pileup_files, project_prefix)
##run laser
run_laser(hgdp_geno, hgdp_coord, project_prefix + '.k10_r5', '10', '5')
run_laser(hgdp_geno, hgdp_coord, project_prefix + '.k6_r5', '6', '5')
run_laser(hgdp_geno, hgdp_coord, project_prefix + '.k10_r10', '10', '10')
##combine sample and control coord files for graphing
combine_format_coord_files_for_r(project_prefix + '.k10_r5', hgdp_coord, project_prefix + '.k10_r5.hgdp_all.txt', 'na')
combine_format_coord_files_for_r(project_prefix + '.k6_r5', hgdp_coord, project_prefix + '.k6_r5.hgdp_all.txt', 'na')
combine_format_coord_files_for_r(project_prefix + '.k10_r10', hgdp_coord, project_prefix + '.k10_r10.hgdp_all.txt', 'na')

##move results files to ratchet and graph in R

