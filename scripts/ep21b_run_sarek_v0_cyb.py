#!/usr/bin/env python
import os
import subprocess

'''
working from 10 threads?

for testing
qsub -Iq cdbrmq -l mem=120gb,ncpus=10 -P 19833a08-f6fb-4bea-8526-8a79069da878
qsub -Iq cdbrmq -l mem=20gb,ncpus=1 -P 19833a08-f6fb-4bea-8526-8a79069da878
#module load java/1.8.0_202 ##for picard/gatk, ?needed if using nextflow
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

##methods

def make_sarek_tsv_by_cohort(infile, outfile, w_dict):
	with open(infile, "r") as in_fh, open(outfile, "w") as outfh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				sample = line[1]
				sequenced = line[7]
				sex_number = line[4]
				if sex_number == '1':
					sex = 'M'
				elif sex_number == '2':
					sex = 'F'
				else:
					print('sex not recognized:', ped_dict[sample][1])
				fq1 = w_dict + sample + '.r1.fastq.gz'
				fq2 = w_dict + sample + '.r2.fastq.gz'
				line_out = [sample, sex, '1', sample, '1', fq1, fq2]
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
run_sarek_by_cohort(work_dir, exome_info_file)
