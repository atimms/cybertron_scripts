#!/usr/bin/env python
import sys
import subprocess
import os
import glob


##note
'''

load modules:
module load biobuilds/2017.11
module load vt/0.5772-60f436c3 
'''

##parameters
delim = '\t'
##programs
mem_get_me = '/home/atimms/programs/MEM/ME_UPD_Contamination_Detection_From_TrioVCF.V.0.0.3.pl'
mem_analysis = '/home/atimms/programs/MEM/MEM_windowAnalysis.V.0.0.4.pl'
snpeff_jar = '/home/atimms/programs/snpEff/snpEff.jar'
##file
supdup_bed = '/home/atimms/ngs_data/references/hg19/hg19_genomicSuperDups.bed'
genome_file = '/home/atimms/programs/MEM/b37.genome.bed'
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'


def make_vcf_for_mem(ped, input_vcf):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	temp_vcf = 'temp1.vcf'
	normalized_vcf = ped + '.int.norm.vcf'
	with open (temp_vcf, 'w') as tvcf_fh:
		zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		vt_decompose = subprocess.Popen(['vt', 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
		vt_normalize = subprocess.Popen(['vt', 'normalize', '-r', fasta, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
		vt_normalize.wait()
	##annotate with snpeff and compress and index
	with open (normalized_vcf, 'w') as nvcf_fh:
		snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'GRCh37.75', '-v', '-formatEff', '-classic', temp_vcf], stdout=nvcf_fh)
		snpeff_vcf.wait()
		bgzip_vcf = subprocess.Popen(['bgzip', normalized_vcf])
		bgzip_vcf.wait()
		tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', normalized_vcf + '.gz'])
		tabix_vcf.wait()

def make_case_ped(in_file, out_file_suffix):
	outfiles = []
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				aff = line[5]
				sample = line[1]
				if aff == '2':
					out_file = sample + out_file_suffix
					outfiles.append(out_file)
					with open(out_file, "w") as out_fh:
						out_fh.write(delim.join(line) + '\n')
	return(outfiles)

def run_extract_mendelian_errors(ped, vcf, output_name):
	##perl ME_UPD_Contamination_Detection_From_TrioVCF.V.0.0.3.pl -i VCF -p pedigreeFile -o output_name
	me_extract = subprocess.Popen(['perl', mem_get_me, '-i', vcf, '-p', ped, '-o', output_name])
	me_extract.wait()

def filter_me_snps(in_file, out_file):
	temp_bed = in_file.split('.')[0] + '.temp.bed'
	out_bed = in_file.split('.')[0] + '.supdup.temp.bed'
	##make bed from me snps
	with open(in_file, "r") as in_fh, open(temp_bed, "w") as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc == 1:
				header = line
			else:
				line = line.rstrip().split(delim)
				chrom = line[0]
				start = line[1]
				end = str(int(line[2]) + 1)
				out_fh.write(delim.join([chrom, start, end] + line) + '\n')
	##insect with sup dups
	with open(out_bed, "w") as out_fh: 
		hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', temp_bed, '-b', supdup_bed, '-v'], stdout=out_fh)
		hom_bt_intersect.wait()
	##get final file
	with open(out_bed, "r") as in_fh, open(out_file, "w") as out_fh:
		out_fh.write(header)
		for line in in_fh:
			print(line)
			line = line.split(delim)
			line_out = line[3:]
			out_fh.write(delim.join(line_out))

def make_list_of_files(file_list, list_file):
	with open(list_file, "w") as outfh:
		for file_name in file_list:
			outfh.write(file_name + '\n')


def run_mem_analysis(prefix, infile):
	##perl MEM_windowAnalysis.pl -a cases -c caseSample.list_of_files -b backgroundSample.list_of_files -w 2000000 -s 100000
	list_of_files = prefix + '.list_of_files'
	make_list_of_files(infile, list_of_files)
	me_extract = subprocess.Popen(['perl', mem_analysis, '-a', 'cases', '-g', genome_file, '-c', list_of_files, '-n', prefix, '-M', '1'])
	me_extract.wait()

def run_mem_ind_pipeline(work_dir, ped_name, make_norm_vcf):
	os.chdir(work_dir)
	ped_file = ped_name + '.ped'
	case_ped_suffix = '.case.ped'
	vcf_file = ped_name + '.int.norm.vcf.gz'

	if make_norm_vcf == 'yes':
		in_vcf = pedigree + '.intersected_vcfs/0002.vcf.gz'
		make_vcf_for_mem(ped_name, in_vcf)
	##take the cases ffrom ped file
	case_peds = make_case_ped(ped_file, case_ped_suffix)
	print(case_peds)
	snp_files_for_mem = []
	for case_ped in case_peds:
		sample_name = case_ped.split('.')[0]
		me_prefix = sample_name + '.me_output'
		me_snp_file = me_prefix + '.ME.snps.txt'
		me_filtered_snp_file = me_prefix + '.filtered.snps.txt'
		##run step one of mem
		run_extract_mendelian_errors(case_ped, vcf_file, me_prefix)
		##filter me snps
		filter_me_snps(me_snp_file, me_filtered_snp_file)
		snp_files_for_mem.append(me_filtered_snp_file)
	##run analysis
	run_mem_analysis(ped_name, snp_files_for_mem)







##run methods
'''
##testing
project_name = 'mem_test_1217'
pedigrees = ['LR18-470', 'LR18-493']
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_mem_1218/test_1218'
# pedigrees = ['LR06-005', 'LR06-207']
pedigrees = ['LR07-201']
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_mem_1218/test2_1218'
##starting from interescted vcf dir (yes) or .int.norm.vcf.gz (no)





##exomes after 1018
# working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_mem_1218/exomes_after_1018'
# pedigrees = ['LR18-416', 'LR18-470', 'LR18-492', 'LR18-493']
# normalise_vcf = 'no'
##exomes after 1217
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_mem_1218/exomes_after_1217'
# pedigrees = ['LR16-030', 'LR16-031', 'LR16-172']
pedigrees = ['LR16-030']
normalise_vcf = 'no'
##get me's for each pedigree and then filter
# for pedigree in pedigrees:
# 	run_mem_ind_pipeline(working_dir, pedigree, normalise_vcf)
'''

##second test
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_mem_1218/test_0119'
##starting from interescted vcf dir (yes) or .int.norm.vcf.gz (no)
normalise_vcf = 'yes'
pedigrees = ['LR11-343']
##get me's for each pedigree and then filter
for pedigree in pedigrees:
	run_mem_ind_pipeline(working_dir, pedigree, normalise_vcf)






