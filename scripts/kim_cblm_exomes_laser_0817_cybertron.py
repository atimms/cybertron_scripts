#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/exomes/cblm_laser_0817'
os.chdir(working_dir)

##programs and program files
laser = '/home/atimms/programs/LASER-2.03/laser'
pileup2seq = '/home/atimms/programs/LASER-2.03/pileup2seq/pileup2seq.py'
hgdp_site = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.site'
hgdp_bed = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.bed'
hgdp_geno = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.geno'
hgdp_coord = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.RefPC.coord'
fasta = '/home/atimms/ngs_data/references/hg19/human_g1k_v37.fasta'
##file names
b14_bams = ['3C-2P.bwa_gatk.bam', '3C-4P.bwa_gatk.bam', '3C-6P.bwa_gatk.bam', 'DWM10.bwa_gatk.bam', 'DWM13.bwa_gatk.bam', 'DWM3.bwa_gatk.bam', 'LR01-079.bwa_gatk.bam', 'LR02-263.bwa_gatk.bam', 'LR03-039.bwa_gatk.bam', 'LR03-055.bwa_gatk.bam', 'LR03-077.bwa_gatk.bam', 'LR03-120.bwa_gatk.bam', 'LR03-130.bwa_gatk.bam', 'LR03-169.bwa_gatk.bam', 'LR03-206.bwa_gatk.bam', 'LR03-223.bwa_gatk.bam', 'LR03-274.bwa_gatk.bam', 'LR03-278.bwa_gatk.bam', 'LR03-298.bwa_gatk.bam', 'LR03-304.bwa_gatk.bam', 'LR03-305.bwa_gatk.bam', 'LR03-332.bwa_gatk.bam', 'LR03-340.bwa_gatk.bam', 'LR04-017.bwa_gatk.bam', 'LR04-020.bwa_gatk.bam', 'LR04-022a1.bwa_gatk.bam', 'LR04-084.bwa_gatk.bam', 'LR04-106.bwa_gatk.bam', 'LR04-185.bwa_gatk.bam', 'LR04-186.bwa_gatk.bam', 'LR04-208.bwa_gatk.bam', 'LR04-222.bwa_gatk.bam', 'LR04-233.bwa_gatk.bam', 'LR04-341.bwa_gatk.bam', 'LR04-371.bwa_gatk.bam', 'LR04-376.bwa_gatk.bam', 'LR04-399.bwa_gatk.bam', 'LR04-414.bwa_gatk.bam', 'LR05-007_sample1.bwa_gatk.bam', 'LR05-035.bwa_gatk.bam', 'LR05-118.bwa_gatk.bam', 'LR05-120.bwa_gatk.bam', 'LR05-160.bwa_gatk.bam', 'LR05-162.bwa_gatk.bam', 'LR05-203a1.bwa_gatk.bam', 'LR05-203a2_sample2.bwa_gatk.bam', 'LR05-265.bwa_gatk.bam', 'LR05-354.bwa_gatk.bam', 'LR05-396.bwa_gatk.bam', 'LR05-398.bwa_gatk.bam', 'LR06-085.bwa_gatk.bam', 'LR06-105.bwa_gatk.bam', 'LR06-157.bwa_gatk.bam', 'LR06-207.bwa_gatk.bam', 'LR06-278.bwa_gatk.bam', 'LR08-002.bwa_gatk.bam', 'LR08-056.bwa_gatk.bam', 'LR08-323.bwa_gatk.bam', 'LR08-390.bwa_gatk.bam', 'LR08-396.bwa_gatk.bam', 'LR09-023_sample1.bwa_gatk.bam', 'LR09-227a1.bwa_gatk.bam', 'LR09-280.bwa_gatk.bam', 'LR09-416a1.bwa_gatk.bam', 'LR10-016.bwa_gatk.bam', 'LR10-026.bwa_gatk.bam', 'LR10-102.bwa_gatk.bam', 'LR10-199.bwa_gatk.bam', 'LR10-222.bwa_gatk.bam', 'LR10-228a1.bwa_gatk.bam', 'LR10-230.bwa_gatk.bam', 'LR10-243.bwa_gatk.bam', 'LR11-033.bwa_gatk.bam', 'LR11-042.bwa_gatk.bam', 'LR11-152.bwa_gatk.bam', 'LR11-169.bwa_gatk.bam', 'LR11-241.bwa_gatk.bam', 'LR11-330a1.bwa_gatk.bam', 'LR11-331a4.bwa_gatk.bam', 'LR12-032.bwa_gatk.bam', 'LR12-115.bwa_gatk.bam', 'LR12-313a2.bwa_gatk.bam', 'LR12-316.bwa_gatk.bam', 'LR12-434.bwa_gatk.bam', 'LR12-439.bwa_gatk.bam', 'LR12-443.bwa_gatk.bam', 'LR12-463.bwa_gatk.bam', 'LR12-464.bwa_gatk.bam', 'LR12-475.bwa_gatk.bam', 'LR13-002.bwa_gatk.bam', 'LR13-037.bwa_gatk.bam', 'LR13-085.bwa_gatk.bam', 'LR13-153.bwa_gatk.bam', 'LR13-199.bwa_gatk.bam', 'LR13-200.bwa_gatk.bam', 'LR13-315.bwa_gatk.bam', 'LR14-071.bwa_gatk.bam', 'LR14-098.bwa_gatk.bam', 'LR14-221.bwa_gatk.bam', 'LR16-079.bwa_gatk.bam']
# b14_bams = ['3C-2P.bwa_gatk.bam', '3C-4P.bwa_gatk.bam']
b14_prefix = 'cblm_probands'
b14_hgdp_all_for_r = b14_prefix + '.hgdp_all.txt'
b14_hgdp_select_for_r = b14_prefix + '.hgdp_select.txt'
##testing with danileas data
test_bams = ['101000201.bwa_gatk.bam', '101000601.bwa_gatk.bam', '103000101.bwa_gatk.bam']
test_prefix = 'dan_test'
test_hgdp_all_for_r = b14_prefix + '.hgdp_all.txt'
test_hgdp_select_for_r = b14_prefix + '.hgdp_select.txt'


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

##cblm exomes
##make mpileup files
# samtools_pileup(b14_bams)
##covert pileup and combine in seq file
##get list of pileup from bams and then combine into seq file
# pileup_files = pileup_names_from_bams(b14_bams)
# print pileup_files
# combine_mpileup_to_seq_file(pileup_files, b14_prefix)
##run laser
run_laser(hgdp_geno, hgdp_coord, b14_prefix + '.k10_r5', '10', '5')
run_laser(hgdp_geno, hgdp_coord, b14_prefix + '.k6_r5', '6', '5')
run_laser(hgdp_geno, hgdp_coord, b14_prefix + '.k10_r10', '10', '10')
##combine sample and control coord files for graphing
combine_format_coord_files_for_r(b14_prefix + '.k10_r5', hgdp_coord, b14_prefix + '.k10_r5.hgdp_all.txt', 'na')
combine_format_coord_files_for_r(b14_prefix + '.k6_r5', hgdp_coord, b14_prefix + '.k6_r5.hgdp_all.txt', 'na')
combine_format_coord_files_for_r(b14_prefix + '.k10_r10', hgdp_coord, b14_prefix + '.k10_r10.hgdp_all.txt', 'na')

##cblm exomes
##make mpileup files
# samtools_pileup(test_bams)
##covert pileup and combine in seq file
##get list of pileup from bams and then combine into seq file
# pileup_files = pileup_names_from_bams(test_bams)
# print pileup_files
# combine_mpileup_to_seq_file(pileup_files, test_prefix)
##run laser
# run_laser(hgdp_geno, hgdp_coord, test_prefix, '10', '3')
# combine_format_coord_files_for_r(test_prefix, hgdp_coord, test_prefix + '.hgdp_all.txt' , 'na')

