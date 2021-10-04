#!/usr/bin/env python
import os
import subprocess
import shutil

'''
working from 10 threads...

for testing
qsub -Iq cdbrmq -l mem=120gb,ncpus=10 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load java/1.8.0_202 ##for picard/gatk

'''


##parameters
delim = '\t'
threads = '10'
##programs
samtools = '/home/atimms/programs/samtools-1.11/bin/samtools'
picard = '/home/atimms/programs/picard_2.25.6/picard.jar'
sarek_dir = '/home/atimms/programs/sarek'
hg38_exome_capture = '/home/atimms/ngs_data/references/exome_beds_hg38/hg38_targets_combined_padded_0721.bed'
hg38_refseq_exons = '/home/atimms/ngs_data/references/hg38/hg38_RefSeq_exons.bed'

##directory names
fq_dir = 'fq_files/'
original_file_dir = 'original_files/'


##methods
def make_ped_files(input_dict, file_prefix):
	combined_ped = file_prefix + '.ped'
	header = ['#family_id', 'name', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
	##combined ped
	with open(combined_ped, "w") as cpout_fh:
		cpout_fh.write(delim.join(header) + '\n')
		for ped in input_dict:
			# print ped, input_dict[ped]
			outfile = ped + '.ped'
			with open(outfile, "w") as out_fh:
				out_fh.write(delim.join(header) + '\n')
				for outline in input_dict[ped]:
					cpout_fh.write(delim.join(outline) + '\n')
					out_fh.write(delim.join(outline) + '\n')


def make_ped_files_and_dict_from_info(input_file, file_prefix):
	ped_file_dict, analysis_dict = {}, {}
	with open(input_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				ped_name = line[0]
				sample_name = line[1]
				ped_type = line[6]
				sequenced = line[7]
				filetype = line[8]
				filename = line[9]
				# cram_fa = line[10]
				mosaic = line[10]
				sex = line[4]
				# print ped_name, line
				##add data to ped file dict
				ped_file_info = line[:6]
				# print ped_file_info
				if ped_name in ped_file_dict:
					ped_file_dict[ped_name].append(ped_file_info)
				else:
					ped_file_dict[ped_name] = [ped_file_info]
				##add data to analysis_dict
				if sequenced == 'yes':
					if ped_name in analysis_dict:
						analysis_dict[ped_name][0].append(sample_name)
						analysis_dict[ped_name][1].append(filetype)
						analysis_dict[ped_name][2].append(filename)
						# analysis_dict[ped_name][3].append(cram_fa)
						analysis_dict[ped_name][3].append(ped_type)
						analysis_dict[ped_name][4].append(mosaic)
						analysis_dict[ped_name][5].append(sex)
					else:
						analysis_dict[ped_name] = [[sample_name], [filetype], [filename], [ped_type], [mosaic], [sex]]
	##make the ped files
	make_ped_files(ped_file_dict, file_prefix)
	##retrun analysis info
	return analysis_dict		


##does this work, or need to go cram>bam>fq
def convert_cram_fastq(cramfile, r1_fastq, r2_fastq):
	# print(cramfile, r1_fastq, r2_fastq)
	sample_name = cramfile.split('.')[0]
	temp_cram = sample_name + 'temp.cram'
	##sort with/without threading
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-n', '-O', 'cram', '-T', sample_name, '-o', temp_cram, '-@', threads, '-m', '10G', cramfile])
	# st_sort_pe = subprocess.Popen([samtools, 'sort', '-n', '-O', 'cram', '-T', sample_name, '-o', temp_cram, '-m', '10G', cramfile])
	st_sort_pe.wait()
	# st_index = subprocess.Popen(['samtools', 'index', '-@', '20', temp_cram])
	# st_index.wait()
	##make fqs with/without threading
	st_fq = subprocess.Popen([samtools, 'fastq', '-@', threads, '-c', '6', '-1', r1_fastq, '-2', r2_fastq, temp_cram])
	# st_fq = subprocess.Popen([samtools, 'fastq', '-c', '6', '-1', r1_fastq, '-2', r2_fastq, temp_cram])
	st_fq.wait()
	##mv cram to origanl file dir
	shutil.move(cramfile, original_file_dir)

def convert_bam_fastq_picard(bamfile, r1_fastq, r2_fastq):
	# picard_sam_fq = subprocess.Popen(['java', '-Xmx8g', '-jar', picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
	picard_sam_fq = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'SamToFastq', 'INPUT=' + bamfile, 'FASTQ=' + r1_fastq, 'SECOND_END_FASTQ=' + r2_fastq, 'VALIDATION_STRINGENCY=SILENT'])
	picard_sam_fq.wait()
	##mv to origanl file dir
	shutil.move(bamfile, original_file_dir)


def rename_fastq_files(infiles, r1_fastq, r2_fastq):
	##cp to std file name
	shutil.copy(infiles[0], r1_fastq)
	shutil.copy(infiles[1], r2_fastq)
	##mv to origanl file dir
	shutil.move(infiles[0], original_file_dir)
	shutil.move(infiles[1], original_file_dir)


def prepare_seq_files(ped_dict):
	fastq_dict = {}
	for ped in ped_dict:
		print(ped, ped_dict[ped])
		samples = ped_dict[ped][0]
		for sample in samples:
			sample_pos = samples.index(sample)
			filetype = ped_dict[ped][1][sample_pos]
			filename = ped_dict[ped][2][sample_pos]
			sex = ped_dict[ped][5][sample_pos]
			mosaic = ped_dict[ped][4][sample_pos]
			##make dir for fq files and original files
			if not os.path.isdir(fq_dir):
				os.mkdir(fq_dir)
			if not os.path.isdir(original_file_dir):
				os.mkdir(original_file_dir)
			read1_fastq = fq_dir + sample + '.r1.fastq.gz'
			read2_fastq = fq_dir + sample + '.r2.fastq.gz'
			if filetype == 'fastq':
				original_fastqs = filename.split(',')
				rename_fastq_files(original_fastqs, read1_fastq, read2_fastq)
			elif filetype == 'bam':
				convert_bam_fastq_picard(filename, read1_fastq, read2_fastq)
			elif filetype == 'cram':
				convert_cram_fastq(filename, read1_fastq, read2_fastq)


##master method
def make_exome_files(working_dir, info_file):
	os.chdir(working_dir)
	project_name = info_file.rsplit('.',1)[0]
	##make ped files per ped and return analysis dict
	analysis_dict = make_ped_files_and_dict_from_info(info_file, project_name)
	##prep files i.e, convery cram/bams to fastq
	prepare_seq_files(analysis_dict)





##run methods
# work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_genedx_0821/'
# exome_info_file = 'ghayda_genedx_0821.txt'
# make_exome_files(work_dir, exome_info_file)

# work_dir = '/home/atimms/ngs_data/exomes/working/kim_exomes_0621/'
# exome_info_file = 'kim_exomes_0621.txt'
# exome_info_file = 'kim_exomes_0621_temp.txt'
# make_exome_files(work_dir, exome_info_file)

work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_genedx_0921/'
exome_info_file = 'ghayda_genedx_0921.txt'
make_exome_files(work_dir, exome_info_file)

