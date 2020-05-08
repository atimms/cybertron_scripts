#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys


'''
##load modules required for analysis
module load java/1.8.0_121 
module load biobuilds/2016.11
module load local_python/3.6.4 
'''

##parameters
delim = '\t'
threads = '16'
working_dir = '/home/atimms/ngs_data/targetted/jimmy_fetalseq_amplicons_1018'
os.chdir(working_dir)
##programs
bwa = 'bwa'
samtools = 'samtools'
bcftools = 'bcftools'
freebayes = '/home/atimms/programs/freebayes_0718/bin/freebayes'
bgzip = 'bgzip'
gatk = '/cm/shared/apps/GATK/3.7/GenomeAnalysisTK.jar'


##filenames
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
bwa_bam_suffix = '.bwa.bam'
var_bed = 'fetal_vars_1018.bed'
bamlist = 'bams.list'

def align_with_bwa_pe(sample_list, r1_suffix, r2_suffix):
	#bwa mem -R "@RG\tID:${SAMPLENAME}_RG\tSM:${SAMPLENAME}" $refdir/Homo_sapiens_assembly19.fasta ${SAMPLENAME}.indexed.fq > ${SAMPLENAME}.indexed.sam
	for sample in sample_list:
		r1_fq = sample + r1_suffix
		r2_fq = sample + r2_suffix
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample
		pe_bam = sample + bwa_bam_suffix
		with open(pe_bam, 'w') as pe_out:
			bwa_pe = subprocess.Popen([bwa, 'mem', '-t', threads, '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
			st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-u', '-'], stdin=bwa_pe.stdout, stdout=subprocess.PIPE)
			st_sort_pe = subprocess.Popen([samtools, 'sort','-O', 'bam', '-T', sample, '-'], stdin=st_sam_bam_pe.stdout, stdout=pe_out)
			st_sort_pe.wait()
		st_index = subprocess.Popen([samtools, 'index', pe_bam])


##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

def gatk_make_gvcfs(sample_dict, bed, bam_suffix):
	##get g.vcfs for each sample (add -L)
	for sample in sample_dict:
		gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', sample + bam_suffix,'-L', bed, '--emitRefConfidence', 'GVCF', '-o', sample + '.g.vcf'])
		# gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', sample + bam_suffix,'-L', bed, '--emitRefConfidence', 'BP_RESOLUTION', '--variant_index_type', 'LINEAR', '--variant_index_parameter', '128000', '-o', sample + '.g.vcf'])

		gatk_hc.wait()

def get_alleles_from_gvcfs(sample_list, bed, out_file):
	allele_dict = {}
	with open(bed, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			b_pos = '_'.join([line[0],line[2]])
			b_alleles = line[3:]
			allele_dict[b_pos] = [b_alleles, []]
	# print(len(allele_dict))
	# print(allele_dict.keys(), allele_dict.items())
	for s in sample_list:
		s_gvcf = s + '.g.vcf'
		with open(s_gvcf, "r") as sg_fh:
			for line in sg_fh:
				if line[0] != '#':
					line = line.rstrip().split(delim)
					v_pos = '_'.join(line[:2])
					v_ref = line[3][0]
					v_alt = line[4].split(',')
					v_info = line[9].split(':')
					# print(v_pos)
					if v_pos in allele_dict:
						# print(allele_dict[v_pos])
						##check ref alleles match
						if v_ref == allele_dict[v_pos][0][0]:
							pass
						else:
							print("ref alleles don't match", line)
						##the get ref allele info
						bed_allele = allele_dict[v_pos][0][1]
						##if no alt allele
						if len(v_alt) == 1:
							ref_count = v_info[1]
							alt_count = '0'
							# print(v_info, ref_count)
						##	
						else:
							
							v_alt_alt = [v[0] for v in v_alt]
							# print(v_alt, v_alt_alt)
							ref_count = v_info[1].split(',')[0]
							v_index = 0
							for v in v_alt_alt:
								v_index += 1
								if v == bed_allele:
									# print(v, v_index)
									if v_index == 1:
										alt_count = v_info[1].split(',')[1]
									elif v_index == 2:
										alt_count = v_info[1].split(',')[2]
									elif v_index == 3:
										alt_count = v_info[1].split(',')[3]
									else:
										print('how many alleles', line, bed_allele)
										print(v_alt, v_alt_alt)
						allele_dict[v_pos][1].extend([ref_count, alt_count])
					else:
						print(v_pos, 'not in dict !!!')
	
	for i in allele_dict:
		print(i, len(allele_dict[i][1]))
	with open(out_file, "w") as out_fh:
		header = ['chr', 'pos', 'ref', 'alt']
		for sample in sample_list:
			header.extend([sample + ' ref count', sample + ' alt count'])
		print(header)
		out_fh.write(delim.join(header) + '\n')
		for i in allele_dict:
			line_out = i.split('_') + allele_dict[i][0] + allele_dict[i][1]
			out_fh.write(delim.join(line_out) + '\n')

def make_list_of_bams(samples, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in samples:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

def genotype_vars_with_freebayes(bedfile, final_vcf):
	vcf_temp1 = 'temp_fb1.vcf'
	vcf_temp2 = 'temp_fb2.vcf'
	##included -i so 
	with open(vcf_temp1, 'w') as vcf_fh:
		run_freebayes = subprocess.Popen([freebayes, '-f', fasta, '-L', bamlist, '-t', bedfile, '-i', '--haplotype-length', '0', '--min-alternate-count', '1', '--min-alternate-fraction', '0', '--pooled-continuous', '--report-monomorphic'], stdout=vcf_fh)
		run_freebayes.wait()
	run_bgzip = subprocess.Popen(['bgzip', vcf_temp1])
	run_bgzip.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	os.remove( vcf_temp1 + '.gz')


def gatk_unified_genotyper(sample_dict, bed, bam_suffix):
	for sample in sample_dict:
		gatk_ug = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'UnifiedGenotyper', '-R', fasta, '-I', sample + bam_suffix,'-L', bed, '-dt', 'none', '--output_mode', 'EMIT_ALL_SITES', '-o', sample + '.ug.vcf'])
		gatk_ug.wait()

def get_alleles_from_ug_vcfs(sample_list, bed, out_file):
	allele_dict = {}
	with open(bed, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			b_pos = '_'.join([line[0],line[2]])
			b_alleles = line[3:]
			allele_dict[b_pos] = [b_alleles, []]
	# print(len(allele_dict))
	# print(allele_dict.keys(), allele_dict.items())
	for s in sample_list:
		s_gvcf = s + '.ug.vcf'
		with open(s_gvcf, "r") as sg_fh:
			for line in sg_fh:
				if line[0] != '#':
					line = line.rstrip().split(delim)
					v_pos = '_'.join(line[:2])
					v_ref = line[3][0]
					v_alt = line[4]
					v_info = line[9].split(':')
					print(v_info, line)
					# print(v_pos)
					if v_pos in allele_dict:
						# print(allele_dict[v_pos])
						##check ref alleles match
						if v_ref == allele_dict[v_pos][0][0]:
							pass
						else:
							print("ref alleles don't match", line)
						##the get ref allele info
						bed_allele = allele_dict[v_pos][0][1]
						##if no alt allele
						if v_alt == '.':
							if v_info == ['./.']:
								ref_count == '0'
							else:
								ref_count = v_info[1]
							alt_count = '0'
							# print(v_info, ref_count)
						##	
						else:
							ref_count = v_info[1].split(',')[0]
							alt_count = v_info[1].split(',')[1]
						allele_dict[v_pos][1].extend([ref_count, alt_count])
					else:
						print(v_pos, 'not in dict !!!')
	
	for i in allele_dict:
		print(i, len(allele_dict[i][1]))
	with open(out_file, "w") as out_fh:
		header = ['chr', 'pos', 'ref', 'alt']
		for sample in sample_list:
			header.extend([sample + ' ref count', sample + ' alt count'])
		print(header)
		out_fh.write(delim.join(header) + '\n')
		for i in allele_dict:
			line_out = i.split('_') + allele_dict[i][0] + allele_dict[i][1]
			out_fh.write(delim.join(line_out) + '\n')



##run methods
r1_fq_suffix = '_L001_R1_001.fastq.gz'
r2_fq_suffix = '_L001_R2_001.fastq.gz'
samples = ['Brain-a_S1', 'Brain-b_S4', 'Con1-a_S3', 'Con1-b_S6', 'Con2-a_S9', 'Con2-b_S12', 'Heart-a_S7', 'Heart-b_S10', 
		'Kidney-a_S13', 'Kidney-b_S16', 'Liver-a_S19', 'Liver-b_S21', 'Lung-a_S2', 'Lung-b_S5', 'Muscle-a_S8', 'Muscle-b_S11', 
		'NTC-a_S15', 'NTC-b_S18', 'Placenta-a_S20', 'Placenta-b_S22', 'Skin-a_S14', 'Skin-b_S17']
# samples = ['Brain-a_S1', 'Brain-b_S4']
project = 'fetal_amplicons_1018'
allele_file = project + '.allele_counts.xls'
bamlist = 'bams.list'
freebayes_vcf = project + '.fb.vcf'
allele_file2 = project + '.allele_counts_ug.xls'
##map
# align_with_bwa_pe(samples, r1_fq_suffix, r2_fq_suffix)
##made bed from list dana gave me then genotype using gatk
# gatk_make_gvcfs(samples, var_bed, bwa_bam_suffix)
# get_alleles_from_gvcfs(samples, var_bed, allele_file)
##try using freebays
# make_list_of_bams(samples,bwa_bam_suffix, bamlist)
# genotype_vars_with_freebayes(var_bed, freebayes_vcf)

##made bed from list dana gave me then genotype using unified genotype
# gatk_unified_genotyper(samples, var_bed, bwa_bam_suffix)
get_alleles_from_ug_vcfs(samples, var_bed, allele_file2)

