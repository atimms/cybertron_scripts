#!/usr/bin/env python
import os
import subprocess
import glob

'''
qsub -Iq cdbrmq -l mem=10gb,ncpus=1 -P 0f241fcd-0dd3-491c-acac-15c6cff1c880
'''
##programs etc
delim = '\t'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
json_tsv = '/home/atimms/programs/nirvana_json_2_tsv_1221.py'
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,gnomad211_exome,gnomad211_genome,dbnsfp41a,cosmic90_coding,cosmic90_noncoding,bed,vcf', '-bedfile', 'minivan_hotspot_1121.bed', '-vcfdbfile', 'Illumina_somatic_hotspots_GRCh37.vcf']
av_operation = ['-operation', 'g,f,f,f,f,f,r,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']

##file names
fasta = '/home/atimms/ngs_data/references/hg19/genome.fa'
tsv_suffix = '.tsv'
norm_vcf_suffix = '.bcf_norm.vcf.gz'
multi_suffix = '.hg19_multianno.txt'
anno_suffix = '.all_variants.xls'
filtered_anno_suffix = '.0.05af_variants.xls'
mapping_metrics_suffix = '.mapping_metrics.csv'
bed_metrics_suffix = '.target_bed_coverage_metrics.csv'
qc_stats_suffix = '.qc_metrics.xls'

##methods
def run_annovar_vcf(input_vcf, normalized_vcf, file_prefix):
	##normalize
	bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', normalized_vcf, '-O', 'z', input_vcf])
	bcftools_norm.wait()
	##annotate vcfs with annovar i.e. run_table_annovar
	command = [table_annovar] + av_buildver + [normalized_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', file_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def convert_json_tsv(infile, outfile):
	json_to_tsv = subprocess.Popen(['python', json_tsv, '-i', infile, '-o', outfile])
	json_to_tsv.wait()	


def combine_files_filter_format(tsv_file, multianno_file, anno_file, filt_anno_file):
	##make dict from multianno file
	multi_dict = {}
	with open(multianno_file, 'r') as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.strip().split(delim)
			if lc ==1:
				multi_header = line[5:10] + ['gnomad211_genome_AF', 'gnomad211_genome_AF_popmax'] + ['gnomad211_exome_AF', 'gnomad211_exome_AF_popmax'] + [line[69]] + [line[84]] + line[95:97] + ['sc_hotspot', 'illumina_hotspot']
			else:
				##change name from annovar input
				hotspot_info1 = line[97]
				if hotspot_info1 == 'Name=NA':
					hotspot_info1 = 'flag_review'
				hotspot_info2 = line[98]
				if hotspot_info2 == 'NA':
					hotspot_info2 = 'flag_review'
				var = '_'.join(line[102:104] + line[105:107])
				# print(var)
				info = line[5:12] + line[27:29] + [line[69]] + [line[84]] + line[95:97] + [hotspot_info1,hotspot_info2]
				genome_popmax = line[11]
				exome_popmax = line[28]
				if genome_popmax == '.':
					genome_popmax = '0'
				if exome_popmax == '.':
					exome_popmax = '0'				
				multi_dict[var] = [info, [float(genome_popmax), float(exome_popmax)]]
	##open tsv and add multianno info
	with open(tsv_file, 'r') as in_fh, open(anno_file, 'w') as out_fh, open(filt_anno_file, 'w') as fout_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.strip().split(delim)
			if lc ==1:
				header = line + multi_header
				out_fh.write(delim.join(header) + '\n')
				fout_fh.write(delim.join(header) + '\n')
			else:
				var2 = '_'.join(line[:4])
				if var2 in multi_dict:
					line_out = line + multi_dict[var2][0]
					##print to outfile
					out_fh.write(delim.join(line_out) + '\n')
					if max(multi_dict[var2][1]) <= 0.05:
						##print to filtered outfile
						fout_fh.write(delim.join(line_out) + '\n')
				else:
					out_fh.write(delim.join(line) + '\n')



def remove_int_files(suffixes_to_rm):
	for suffix_to_rm in suffixes_to_rm:
		files_rm = glob.glob('*' + suffix_to_rm)
		print(files_rm)
		for file_rm in files_rm:
			os.remove(file_rm)
	
def make_tsv_add_extra_annotation(file_name_dict):
	for file_name in file_name_dict:
		json = file_name_dict[file_name][0]
		tsv = file_name + tsv_suffix
		in_vcf = file_name_dict[file_name][1]
		norm_vcf = file_name + norm_vcf_suffix
		multianno = file_name + multi_suffix
		anno = file_name + anno_suffix
		filtered_anno = file_name + filtered_anno_suffix
		##convert gzip json to tsv
		convert_json_tsv(json, tsv)
		##bcftools on original vcf andrun annovar
		run_annovar_vcf(in_vcf, norm_vcf, file_name)
		##combine and format multianno/tsv
		combine_files_filter_format(tsv, multianno, anno, filtered_anno)
	##rm intermediate files
	remove_int_files([tsv_suffix, norm_vcf_suffix, multi_suffix, 'put'])


def get_mapping_stats(samples, mapping_info_wanted, bed_info_wanted):
	for sample in samples:
		outfile = sample + qc_stats_suffix
		mapping_file = sample + mapping_metrics_suffix
		bed_file = sample + bed_metrics_suffix
		with open(outfile, 'w') as out_fh:
			with open(mapping_file, 'r') as in_fh:
				for line in in_fh:
					line = line.strip().split(',')
					# print(line)
					mapping_stat = line[2]
					info = line[3:]
					if line[1] == '':
						if mapping_stat in mapping_info_wanted:
							out_fh.write(delim.join([mapping_stat] + info) + '\n')
			with open(bed_file, 'r') as in_fh:
				for line in in_fh:
					print(line)
					line = line.strip().split(',')
					bed_stat = line[2]
					info = line[3:]
					if bed_stat in bed_info_wanted:
						out_fh.write(delim.join([mapping_stat] + info) + '\n')					





##run methods

##params etc
working_dir = '/active/bennett_j/miniVAN/VANseq_BSSH_Pipeline_Development'
os.chdir(working_dir)

# file_prefix_names = ['21H-164L0001.hard-filtered', '21H-267L0004.hard-filtered', '21H-301L0003.hard-filtered']
##21H-277L0003 309L0001 didn't work
file_prefix_dict = {'21H-164L0001.hard-filtered': ['21H-164L0001.hard-filtered.annotations.json.gz', '21H-164L0001.hard-filtered.vcf.gz'],
		'21H-277L0003.hard-filtered': ['21H-277L0003.hard-filtered.annotations.json.gz','21H-277L0003.hard-filtered.vcf.gz'],
		'21H-301L0003.hard-filtered': ['21H-301L0003.hard-filtered.annotations.json.gz','21H-301L0003.hard-filtered.vcf.gz'],
		'21H-306L0006.hard-filtered': ['21H-306L0006.hard-filtered.annotations.json.gz','21H-306L0006.hard-filtered.vcf.gz'],
		'21H-309L0001.hard-filtered': ['21H-309L0001.hard-filtered.annotations.json.gz','21H-309L0001.hard-filtered.vcf.gz']}
# file_prefix_dict = {'21H-164L0001.hard-filtered': ['21H-164L0001.hard-filtered.annotations.json.gz', '21H-164L0001.hard-filtered.vcf.gz']}
##changing json to annotated file
# make_tsv_add_extra_annotation(file_prefix_dict)


##get mappng stats
mapping_stats_wanted = ['Total input reads', 'Number of duplicate marked reads', 'Mapped reads', 'Number of unique & mapped reads (excl. duplicate marked reads)',
			'Reads with MAPQ [40:inf)', 'Reads with MAPQ [ 0:10)', 'Estimated read length', 'Bases in target bed [% of genome]', 
			'Average sequenced coverage over target region', 'Insert length: mean', 'Insert length: standard deviation']
bed_stats_wanted = ['Uniformity of coverage (PCT > 0.2*mean) over QC coverage region', 'PCT of QC coverage region with coverage [100x: inf)', 
			'PCT of QC coverage region with coverage [ 50x: inf)', 'PCT of QC coverage region with coverage [ 20x: inf)', 
			'PCT of QC coverage region with coverage [ 10x: inf)', 'PCT of QC coverage region with coverage [  3x: 10x)', 
			'PCT of QC coverage region with coverage [  1x:  3x)', 'PCT of QC coverage region with coverage [  0x:  1x)']
samples_stats_wanted = ['21H-277L0003', '21H-301L0003', '21H-306L0006', '21H-309L0001']



get_mapping_stats(samples_stats_wanted, mapping_stats_wanted, bed_stats_wanted)





