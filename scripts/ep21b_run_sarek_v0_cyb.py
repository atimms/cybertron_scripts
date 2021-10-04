#!/usr/bin/env python
import os
import subprocess

'''

conda activate nextflow

'''


##parameters
delim = '\t'
##programs
samtools = '/home/atimms/programs/samtools-1.11/bin/samtools'
picard = '/home/atimms/programs/picard_2.25.6/picard.jar'
sarek_dir = '/home/atimms/programs/sarek'
hg38_exome_capture = '/home/atimms/ngs_data/references/exome_beds_hg38/hg38_targets_combined_padded_0721.bed'
hg38_refseq_exons = '/home/atimms/ngs_data/references/hg38/hg38_RefSeq_exons.bed'

##params
fq_dir = 'fq_files/'


##methods

def make_sarek_tsv_by_cohort(infile, outfile, w_dict):
	with open(infile, "r") as in_fh, open(outfile, "w") as outfh:
		lc = 0
		ped_list = []
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				ped_id = line[0]
				sample = line[1]
				sequenced = line[7]
				sex_number = line[4]
				if sex_number == '1':
					sex = 'M'
				elif sex_number == '2':
					sex = 'F'
				else:
					print('sex not recognized:', ped_dict[sample][1])
				##set up for lane number, has to be unique per ped (issue with sarek)
				ped_list.append(ped_id)
				lane_number = ped_list.count(ped_id)
				fq1 = w_dict + fq_dir + sample + '.r1.fastq.gz'
				fq2 = w_dict + fq_dir + sample + '.r2.fastq.gz'
				##using lc instead of true lane id as issue withe read id
				line_out = [sample, sex, '1', sample, str(lane_number), fq1, fq2]
				outfh.write(delim.join(line_out) + '\n')




def run_sarek_by_cohort(working_dir, info_file):
	os.chdir(working_dir)
	project_name = info_file.rsplit('.',1)[0]
	sarek_tsv = project_name + '_fastq.tsv'
	##make tsv file
	make_sarek_tsv_by_cohort(info_file, sarek_tsv, working_dir)
	##run sarek
	# '''
	#nextflow run nf-core/sarek -bg -profile bennett,conda --input ./resources/LR20-210_fastq.tsv --tools HaplotypeCaller --generate_gvcf --genome GRCh38 --outdir /home/atimms/ngs_data/exomes/working/ghayda_new_pipeline_0721
	os.chdir(sarek_dir)
	##add -bg, wothwhile?
	# sarek_cmd = subprocess.Popen(['nextflow', 'run', 'nf-core/sarek', '-bg', '-r', '2.7.1', '-profile', 'bennett,conda', '--input', working_dir + sarek_tsv, '--tools', 'HaplotypeCaller,Mutect2,Strelka,FreeBayes', '--generate_gvcf', '--genome', 'GRCh38', '--target_bed', hg38_exome_capture, '--outdir', working_dir])
	# sarek_cmd = subprocess.Popen(['nextflow', 'run', 'nf-core/sarek', '-r', '2.7.1', '-profile', 'bennett,conda', '--input', working_dir + sarek_tsv, '--tools', 'HaplotypeCaller,Mutect2,Strelka,FreeBayes', '--generate_gvcf', '--genome', 'GRCh38', '--target_bed', hg38_exome_capture, '--outdir', working_dir])
	##without any genotyping
	sarek_cmd = subprocess.Popen(['nextflow', 'run', 'nf-core/sarek', '-r', '2.7.1', '-profile', 'bennett,conda', '--input', working_dir + sarek_tsv, '--genome', 'GRCh38', '--target_bed', hg38_exome_capture, '--outdir', working_dir])
	sarek_cmd.wait()
	os.chdir(working_dir)
	# '''







##run methods
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_genedx_0821/'
exome_info_file = 'ghayda_genedx_0821.txt'
# exome_info_file = 'LR10-064_0821.txt' ##repeat on just LR10-064 (issue with rg)
# run_sarek_by_cohort(work_dir, exome_info_file)

work_dir = '/home/atimms/ngs_data/exomes/working/kim_exomes_0621/'
# exome_info_file = 'kim_exomes_0621.txt'
##split into 2 files as issue with sarek prepocessing
# exome_info_file = 'kim_exomes_0621_1.txt'
# exome_info_file = 'kim_exomes_0621_2.txt'
# exome_info_file = 'kim_exomes_0621_3.txt'
# run_sarek_by_cohort(work_dir, exome_info_file)

work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_genedx_0921/'
exome_info_file = 'ghayda_genedx_0921.txt'
run_sarek_by_cohort(work_dir, exome_info_file)

