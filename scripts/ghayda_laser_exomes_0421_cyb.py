#!/usr/bin/python
import os
import subprocess
import glob
import shutil

##set input variables and parameters
delim = '\t'

'''
qsub -Iq cdbrmq -l mem=30gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
load these:
module load java/1.8.0_121 
module load biobuilds
'''

##programs and program files
laser = '/home/atimms/programs/LASER-2.03/laser'
pileup2seq = '/home/atimms/programs/LASER-2.03/pileup2seq/pileup2seq.py'
hgdp_site = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.site'
hgdp_bed = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.bed'
hgdp_geno = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.geno'
hgdp_coord = '/home/atimms/ngs_data/references/laser/HGDP/HGDP_938.RefPC.coord'
fasta = '/home/atimms/ngs_data/references/hg19/human_g1k_v37.fasta'
rss_dir = '/archive/.snapshot/202005281800-Dobyns_W-Archive/dobyns_w/exome_data/all_exome_files/'

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

##methods
def samtools_pileup(bam_files):
	'samtools mpileup -q 30 -Q 20 -f ../../LASER-resource/reference/hs37d5.fa -l HGDP_938.bed exampleBAM/NA12878.chrom22.recal.bam > NA12878.chrom22.pileup'
	##iterate over bamfiles
	for bam in bam_files:
		sample = bam.rsplit('/', 1)[1].split('.')[0]
		print(sample)
		pileup_file = sample + '.pileup'
		with open(pileup_file, 'w') as pu_fh:
			bcf_index = subprocess.Popen(['samtools', 'mpileup', '-q', '30', '-Q', '20', '-f', fasta, '-l', hgdp_bed, bam], stdout=pu_fh)
			bcf_index.wait()

def pileup_names_from_bams(bam_files):
	pileup_names = []
	for bam in bam_files:
		pileup_file = bam.rsplit('/', 1)[1].split('.')[0] + '.pileup'
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

def make_dict_from_info_file(in_file):
	ped_bam_dict = {}
	with open(in_file, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line = line.strip('\n').split(delim)
			line_count += 1
			if line_count > 1:
				ped_id = line[2]
				sample_id = line[0]
				bam = line[7] + line[6]
				if ped_id in ped_bam_dict:
					ped_bam_dict[ped_id].append(bam)
				else:
					ped_bam_dict[ped_id] = [bam]
	return(ped_bam_dict)

def laser_master(work_dir, infile):
	os.chdir(work_dir)
	ped_dict = make_dict_from_info_file(infile)
	for p in ped_dict:
		bams = ped_dict[p]
		print(p, bams)
		##get pileupfiles 
		samtools_pileup(bams)
		##get names and then comibine and convery
		pileup_files = pileup_names_from_bams(bams)
		print(pileup_files)
		combine_mpileup_to_seq_file(pileup_files, p)
		##run laser 3 ways
		# run_laser(hgdp_geno, hgdp_coord, p + '.k10_r5', '10', '5')
		run_laser(hgdp_geno, hgdp_coord, p + '.k6_r5', '6', '5')
		run_laser(hgdp_geno, hgdp_coord, p + '.k10_r10', '10', '10')
		##combine sample and control coord files for graphing
		# combine_format_coord_files_for_r(p + '.k10_r5', hgdp_coord, p + '.k10_r5.hgdp_all.txt', 'na')
		combine_format_coord_files_for_r(p + '.k6_r5', hgdp_coord, p + '.k6_r5.hgdp_all.txt', 'na')
		combine_format_coord_files_for_r(p + '.k10_r10', hgdp_coord, p + '.k10_r10.hgdp_all.txt', 'na')



##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_laser_exomes_0421'

# ped_info_file = 'ped_bam_info.test.txt'
# ped_info_file = 'ped_bam_info.1.txt'
# ped_info_file = 'ped_bam_info.2.txt'
# ped_info_file = 'ped_bam_info.3.txt'
# ped_info_file = 'ped_bam_info.4.txt'
# ped_info_file = 'ped_bam_info.5.txt'
ped_info_file = 'ped_bam_info.rpt.txt'

laser_master(working_dir, ped_info_file)
