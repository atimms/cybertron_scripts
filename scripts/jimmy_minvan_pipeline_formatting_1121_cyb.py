#!/usr/bin/env python
import os
import subprocess
import glob

'''
qsub -Iq cdbrmq -l mem=10gb,ncpus=1 -P 0f241fcd-0dd3-491c-acac-15c6cff1c880
'''

delim = '\t'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
json_tsv = '/home/atimms/programs/nirvana_json_2_tsv_1121.py'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,gnomad211_exome,gnomad211_genome,dbnsfp41a,cosmic90_coding,cosmic90_noncoding,bed,vcf', '-bedfile', 'minivan_hotspot_1121.bed', '-vcfdbfile', 'Illumina_somatic_hotspots_GRCh37.vcf']
av_operation = ['-operation', 'g,f,f,f,f,f,r,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.']


##methods
def run_annovar(infile_prefixes, in_suffix):
	for infile_prefix in infile_prefixes:
		avinput_file = infile_prefix + in_suffix
		##annotate vcfs with annovar i.e. run_table_annovar
		command = [table_annovar] + av_buildver + [avinput_file] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', infile_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()



def convert_json_tsv(infile_prefixes, in_suffix, out_suffix):
	for infile_prefix in infile_prefixes:
		infile = infile_prefix + in_suffix
		outfile = infile_prefix + out_suffix
		json_to_tsv = subprocess.Popen(['python', json_tsv, '-i', infile, '-o', outfile])
		json_to_tsv.wait()	

def convert_tsv_avinput(infile_prefixes, in_suffix, out_suffix):
	for infile_prefix in infile_prefixes:
		infile = infile_prefix + in_suffix
		outfile = infile_prefix + out_suffix
		with open(outfile, 'w') as out_fh, open(infile, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.strip().split(delim)
				if lc ==1:
					header = line
				else:	
					chrom = line[0]
					pos = line[1]
					ref = line[2]
					alt = line[3]
					if len(ref) == len(alt):
						start = pos
						end = str(int(start) + (len(ref) -1))
						# end = pos
						# ref_out = ref
						# alt_out = alt
					elif len(alt) > len(ref):
						start = pos
						end = pos
						# ref_out = ref
						# alt_out = alt
					elif len(ref) > len(alt):
						start = pos
						end = str(int(pos) + (len(ref) - 1))
						# ref_out = ref
						# alt_out = alt
					# out_fh.write(delim.join([chrom, start, end, ref_out, alt_out] + line) + '\n')
					out_fh.write(delim.join([chrom, start, end, ref, alt] + line) + '\n')
	return(header)

def format_multianno(infile_prefixes, in_suffix, out_suffix, filt_out_suffix, original_header):
	for infile_prefix in infile_prefixes:
		infile = infile_prefix + in_suffix
		outfile = infile_prefix + out_suffix
		filtoutfile = infile_prefix + filt_out_suffix
		with open(outfile, 'w') as out_fh, open(filtoutfile, 'w') as fout_fh, open(infile, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.strip().split(delim)
				if lc ==1:
					header = original_header + line[5:10] + ['gnomad211_genome_AF', 'gnomad211_genome_AF_popmax'] + ['gnomad211_exome_AF', 'gnomad211_exome_AF_popmax'] + [line[69]] + [line[84]] + line[95:97] + ['sc_hotspot', 'illumina_hotspot']
					out_fh.write(delim.join(header) + '\n')
					fout_fh.write(delim.join(header) + '\n')
				else:
					##change name from annovar input
					hotspot_info1 = line[97]
					if hotspot_info1 == 'Name=NA':
						hotspot_info1 = 'flag_review'
					hotspot_info2 = line[98]
					if hotspot_info2 == 'NA':
						hotspot_info2 = 'flag_review'
					out_fh.write(delim.join(line[99:] + line[5:12] + line[27:29] + [line[69]] + [line[84]] + line[95:97] + [hotspot_info1,hotspot_info2]) + '\n')
					genome_popmax = line[11]
					exome_popmax = line[28]
					if genome_popmax == '.':
						genome_popmax = '0'
					if exome_popmax == '.':
						exome_popmax = '0'
					if float(genome_popmax) <= 0.05 and float(exome_popmax) <= 0.05:
						fout_fh.write(delim.join(line[99:] + line[5:12] + line[27:29] + [line[69]] + [line[84]] + line[95:97]+ [hotspot_info1,hotspot_info2]) + '\n')

def remove_int_files(suffixes_to_rm):
	for suffix_to_rm in suffixes_to_rm:
		files_rm = glob.glob('*' + suffix_to_rm)
		print(files_rm)
		for file_rm in files_rm:
			os.remove(file_rm)
	



##run methods

##params etc
working_dir = '/active/bennett_j/miniVAN/VANseq_BSSH_Pipeline_Development'
os.chdir(working_dir)

file_prefix_names = ['21H-164L0001.hard-filtered.annotations', '21H-267L0004.hard-filtered.annotations', 
		'21H-301L0003.hard-filtered.annotations']
##21H-277L0003 309L0001 didn't work
# file_prefix_names = ['21H-164L0001.hard-filtered.annotations']
json_suffix = '.json.gz'
tsv_suffix = '.tsv'
av_suffix = '.avinput'
multi_suffix = '.hg19_multianno.txt'
anno_suffix = '.all_variants.xls'
filtered_anno_suffix = '.0.05af_variants.xls'

##changing json to annotated file
##convert gzip json to tsv
# convert_json_tsv(file_prefix_names, json_suffix, tsv_suffix)
##convert tsv to avinput
# tsv_header = convert_tsv_avinput(file_prefix_names, tsv_suffix, av_suffix)
##run annovar
# run_annovar_vcf(file_prefix_names, av_suffix)
##format multianno
# format_multianno(file_prefix_names, multi_suffix, anno_suffix, filtered_anno_suffix, tsv_header)
##rm intermediate files
remove_int_files([tsv_suffix, av_suffix, multi_suffix])
