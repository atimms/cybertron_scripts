#!/usr/bin/env python
import sys
import subprocess
import os
import glob


"""
to set up TPMCalculator:
conda create --name tpm_calculator
conda activate tpm_calculator
conda install tpmcalculator

to run:
start interactive session i.e.
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878
conda activate tpm_calculator

"""

##parameters
delim = '\t'

##ref files etc
genome_name = 'GRCh38'
fa_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name + '/genome.fa'
gtf_file = '/home/atimms/ngs_data/references/igenomes/' + genome_name +  '/genes.gtf'


def calculate_tpm_paired_end(sample_names, bamfile_suffix, out_file):
	#Usage: ./bin/TPMCalculator -g GTF_file [-d BAM_files_directory|-b BAM_file]
	tpm_dict = {}
	for sample_name in sample_names:
		bam = sample_name + bamfile_suffix
		tpm_file = bam.rsplit('.', 1)[0] + '_genes.out'
		# tpmcalc_cmd = subprocess.Popen(['TPMCalculator', '-g', gtf_file, '-b', bam, '-p'])
		# tpmcalc_cmd.wait()
		##open tpm files and add to dict
		lc = 0
		with open(tpm_file, 'r') as in_tpm_fh:
			for line in in_tpm_fh:
				lc += 1
				if lc > 1:
					line = line.split(delim)
					gene = line[0]
					tpm = line[6]
					if gene in tpm_dict:
						tpm_dict[gene][sample_names.index(sample_name)] = tpm
					else:
						tpm_dict[gene] = ['0'] * len(sample_names)
						tpm_dict[gene][sample_names.index(sample_name)] = tpm
	##print to outfile
	with open(out_file, 'w') as out_fh:
		out_fh.write(delim.join(['gene'] + sample_names) + '\n')
		for g in tpm_dict:
			out_fh.write(delim.join([g] + tpm_dict[g]) + '\n')



def calculate_tpm_paired_end_from_file(sample_file, out_file):
	#Usage: ./bin/TPMCalculator -g GTF_file [-d BAM_files_directory|-b BAM_file]
	##get info file per bam/sample and makes list of sample
	sample_names = []
	with open(sample_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)	
			sample_name = line[0]
			sample_names.append(sample_name)
			bam = line[1]	
			print(sample_name, bam)
			##do the calculations
			# tpmcalc_cmd = subprocess.Popen(['TPMCalculator', '-g', gtf_file, '-b', bam, '-p'])
			# tpmcalc_cmd.wait()
	##add tpm values to dict
	tpm_dict = {}
	with open(sample_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)	
			sample_name = line[0]
			bam = line[1]
			tpm_file = bam.rsplit('/', 1)[1].rsplit('.', 1)[0] + '_genes.out'
			##open tpm files and add to dict
			lc = 0
			with open(tpm_file, 'r') as in_tpm_fh:
				for line in in_tpm_fh:
					lc += 1
					if lc > 1:
						line = line.split(delim)
						gene = line[0]
						tpm = line[6]
						if gene in tpm_dict:
							tpm_dict[gene][sample_names.index(sample_name)] = tpm
						else:
							tpm_dict[gene] = ['0'] * len(sample_names)
							tpm_dict[gene][sample_names.index(sample_name)] = tpm
	##print to outfile
	with open(out_file, 'w') as out_fh:
		out_fh.write(delim.join(['gene'] + sample_names) + '\n')
		for g in tpm_dict:
			out_fh.write(delim.join([g] + tpm_dict[g]) + '\n')


##kim rnaseq 0720 just RL data
'''
working_dir = '/home/atimms/ngs_data/rnaseq/kim_rnaseq_0720'
os.chdir(working_dir)
samples = ['Ctrl_13086_RL', 'Ctrl_H26074_RL', 'DW22840EE_RL', 'DW05950_RL']
bam_suffix = '.Aligned.sortedByCoord.out.bam'
combined_tpm_file = 'kim_rnaseq_0720_b6_rl.tpm.txt'
calculate_tpm_paired_end(samples, bam_suffix, combined_tpm_file)
'''

##kim lots of rnaseq samples, from 0720 and 0218
working_dir = '/home/atimms/ngs_data/rnaseq/kim_calc_tpm_0321'
os.chdir(working_dir)
sample_bam_info_file = 'kim_rnaseq_0321.info.txt'
combined_tpm_file = 'kim_rnaseq_0321.tpms.txt'
##info files split x3
# sample_bam_info_file = 'temp_infoaa'
# sample_bam_info_file = 'temp_infoab'
# sample_bam_info_file = 'temp_infoac'
calculate_tpm_paired_end_from_file(sample_bam_info_file, combined_tpm_file)

