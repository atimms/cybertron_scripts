#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import gzip

##modules 
'''
module load biobuilds/2017.11
'''

##parameters
delim = '\t'
fa_file = '/home/atimms/ngs_data/references/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta'

##programs
table_annovar = '/home/atimms/programs/annovar_0618/table_annovar.pl'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'

#annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp35a,dbnsfp31a_interpro,avsnp150,gnomad211_genome,gnomad211_exome,clinvar_20190305']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.','-vcfinput']


def annotate_vars_from_vcf(out_prefix, in_vcf, norm_vcf):
	'''
	vcf_temp1 = 'temp111.vcf'
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-o', vcf_temp1, in_vcf])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fa_file, '-O', 'z', '-o', norm_vcf, vcf_temp1])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', norm_vcf])
	bcf_index.wait()
	
	##annotate vcf file
	command = [table_annovar] + av_buildver + [norm_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()
	'''
	##change multianno to annotated
	multi = out_prefix + '.hg38_multianno.txt'
	outfile = out_prefix + '.annotated.txt'
	with open(multi, "r") as av, open(outfile, "w") as final:
		lc = 0
		for line in av:
			lc += 1
			line = line.strip('\n').split(delim)
			if lc == 1:
				header = line[:123] + ['filter', 'format', 'info']
				final.write(delim.join(header) + '\n')
			else:
				line_out = line[:123] + [line[132], line[134], line[135]]
				final.write(delim.join(line_out) + '\n')
	# '''

def filter_anntxt_files(in_file, out_file):
	##filter anntxt files
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc, pc = 0, 0
		for line in in_fh:
			line = line.split(delim)
			lc += 1
			if lc == 1:
				out_fh.write(delim.join(line))
			else:
				gatk_filter = line[123]
				rmsk = line[10]
				segdup = line[11]
				rg_func = line[5]
				exonic_func = line[8]
				genome_af = line[85]
				exome_af = line[102]
				if genome_af == '.':
					genome_af = 0
				if exome_af == '.':
					exome_af = 0
				afs = [float(genome_af), float(exome_af)]
				if gatk_filter == 'PASS' and max(afs) <= 0.01:
					if rg_func == 'exonic' or rg_func == 'splicing':
						if exonic_func != 'synonymous SNV':
							pc +=1
							out_fh.write(delim.join(line))
	print(lc, pc)



##run methods
working_dir = '/home/atimms/ngs_data/genomes/milena_wtc_genome_0320'
os.chdir(working_dir)
project_prefix = 'WTC_genome_0320'
start_vcf = 'AH77TTBBXX_DS-229105_GCCAAT_recalibrated.vcf.gz'
norm_vcf = project_prefix + '.norm.vcf.gz'
# annotate_vars_from_vcf(project_prefix, start_vcf, norm_vcf)
filter_anntxt_files(project_prefix + '.annotated.txt', project_prefix + '.rare_exonic.xls')


