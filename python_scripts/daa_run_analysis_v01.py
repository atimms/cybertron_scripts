#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys
import pybedtools as pbt
import dobyns_amplicon_gemini_pipeline_cybertron_v0

##parameters
delim = '\t'
threads = '16'

##programs
stitcher = '/cm/shared/apps/Pisces/5.2.0.1/Stitcher/Stitcher.dll'
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
stitcher_bam_suffix = '.stitch.bam'
bam_list = 'bams.list'
freebayes_vcf_suffix = '.freebayes.vcf.gz'
freebayes_stitcher_vcf_suffix = '.stitcher_freebayes.vcf.gz'
gatk_vcf_suffix = '.gatk.vcf.gz'
gatk_stitcher_vcf_suffix = '.stitcher_gatk.vcf.gz'
bwa_cov_suffix = '.coverage.txt'
bwa_stitcher_cov_suffix = '.stitcher_coverage.txt'

def align_with_bwa_pe(sample_dict):
	#bwa mem -R "@RG\tID:${SAMPLENAME}_RG\tSM:${SAMPLENAME}" $refdir/Homo_sapiens_assembly19.fasta ${SAMPLENAME}.indexed.fq > ${SAMPLENAME}.indexed.sam
	for sample in sample_dict:
		r1_fq = sample_dict[sample][0]
		r2_fq = sample_dict[sample][1]
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample
		pe_bam = sample + bwa_bam_suffix
		# pe_out = open(pe_bam, 'w')
		# bwa_pe = subprocess.Popen([bwa, 'mem', '-t', threads, '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		# st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-u', '-'], stdin=bwa_pe.stdout, stdout=subprocess.PIPE)
		# st_sort_pe = subprocess.Popen([samtools, 'sort','-O', 'bam', '-T', sample, '-'], stdin=st_sam_bam_pe.stdout, stdout=subprocess.PIPE)
		# pe_bam_file = st_sort_pe.communicate()[0]
		# pe_out.write(pe_bam_file)
		# pe_out.close()
		with open(pe_bam, 'w') as pe_out:
			bwa_pe = subprocess.Popen([bwa, 'mem', '-t', threads, '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
			st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-u', '-'], stdin=bwa_pe.stdout, stdout=subprocess.PIPE)
			st_sort_pe = subprocess.Popen([samtools, 'sort','-O', 'bam', '-T', sample, '-'], stdin=st_sam_bam_pe.stdout, stdout=pe_out)
			st_sort_pe.wait()
		st_index = subprocess.Popen([samtools, 'index', pe_bam])
		
def run_stitcher_on_bam(sample_dict):
	bam_files = []
	for sample in sample_dict:
		in_bam = sample + bwa_bam_suffix
		out_bam = sample + '.bwa.stitched.bam'
		temp_bam = sample + '.temp.stitched.bam'
		sorted_bam = sample + stitcher_bam_suffix
		##run sticher, then sort and index bam
		run_stitcher = subprocess.Popen(['dotnet', stitcher, '-Bam', in_bam, '-t', '18', '-OutFolder', './'])
		run_stitcher.wait()
		# st_sort_pe = subprocess.Popen([samtools, 'sort','-O', 'bam', '-o',sorted_bam, '-T', sample, out_bam])
		# st_sort_pe.wait()
		# st_index = subprocess.Popen([samtools, 'index', sorted_bam])
		st_sort_pe = subprocess.Popen(['samtools', 'sort', '-O', 'bam', '-T', sample, '-o', temp_bam, '-@', '10', '-m', '10G', out_bam])
		st_sort_pe.wait()
		st_index = subprocess.Popen(['samtools', 'index', temp_bam])
		st_index.wait()
		##make new bam file with adjusted RG info
		picard_arrd = subprocess.Popen(['picard', 'AddOrReplaceReadGroups', 'INPUT=' + temp_bam, 'OUTPUT=' + sorted_bam, 'RGPL=illumina', 'RGLB=' + sample, 'RGPU=machine1', 'RGID=' + sample, 'RGSM=' + sample, 'TMP_DIR=./', 'CREATE_INDEX=true', 'VALIDATION_STRINGENCY=LENIENT'])
		picard_arrd.wait()


##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

##call freebayes in all sample listed in list of bamfiles
def variant_calling_freebayes(bamlist, name_prefix, final_vcf_prefix, bed_file):
	vcf_temp1 = name_prefix + '.temp_fb1.vcf'
	vcf_temp2 = name_prefix + '.temp_fb2.vcf'
	final_vcf = name_prefix + final_vcf_prefix
	vcf_fh = open(vcf_temp1, 'w')
	freebayes_run = subprocess.Popen([freebayes,'-f' , fasta, '-L', bamlist, '-t', bed_file], stdout=vcf_fh)
	freebayes_run.wait()
	vcf_fh.close()
	bgzip_run = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()




def variant_calling_gatk_ug(bamlist, name_prefix, final_vcf_prefix, bed_file):
	vcf_temp0 = name_prefix + '.temp_gatk0.vcf'
	vcf_temp1 = name_prefix + '.temp_gatk1.vcf'
	vcf_temp2 = name_prefix + '.temp_gatk2.vcf'
	final_vcf = name_prefix + final_vcf_prefix
	vcf_fh = open(vcf_temp1, 'w')
	gatk_ug = subprocess.Popen(['java', '-Xmx12g', '-jar', gatk, '-T', 'UnifiedGenotyper', '-R', fasta, '-rf', 'BadCigar', '-allowPotentiallyMisencodedQuals', '-L', bed_file, '-I', bamlist, '-o', vcf_temp0, '-dcov', '5000', '-dt', 'NONE', '-glm', 'both', '-A', 'AlleleBalance', '-minIndelFrac', '0.005', '-G', 'Standard' ])
	gatk_ug.wait()
	gatk_filter = subprocess.Popen(['java', '-Xmx12g', '-jar', gatk, '-T', 'VariantFiltration', '-R', fasta, '-rf', 'BadCigar', '-allowPotentiallyMisencodedQuals', '-V', vcf_temp0, '-o', vcf_temp1, '-window', '20', '-cluster', '5', '-filterName', 'ABFilter', '-filter' ,'ABHet>0.75', '-filterName', 'QDFilter', '-filter', 'QD<5.0', '-filterName', 'QUALFilter', '-filter', 'QUAL<30.0', '-filterName', 'LowCoverage', '-filter', 'DP<5'])
	gatk_filter.wait()
	vcf_fh.close()
	bgzip_run = subprocess.Popen([bgzip, vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()

def get_coverage_data(bam_list, output_file, bed_file, run_stitcher):
	b = pbt.BedTool(bed_file)
	if run_stitcher == 'yes':
		b_multi = b.multi_bam_coverage(bams=bam_list, f=0.9, output='coverage.temp')
		# pass
	elif run_stitcher == 'no':
		b_multi = b.multi_bam_coverage(bams=bam_list, output='coverage.temp')
		# pass
	with open('coverage.temp', "r") as cov_temp, open(output_file, 'w') as outfile:
		header_list = ['chr', 'start', 'end']
		for bam in bam_list:
			bam = bam.split('.')[0]
			header_list.append(bam)
		header = delim.join(header_list + ['\n'])
		outfile.write(header)
		for line in cov_temp:
			outfile.write(line)

def remove_intermediate_files():
	files_to_go = glob.glob('*temp*')
	print 'removing files:'
	for f in files_to_go:
		os.remove(f)
		print f

def intesect_two_vcf_files(vcf1, vcf2, output_dir):
	bcf_isec = subprocess.Popen([bcftools, 'isec', vcf1, vcf2, '-p', output_dir ])
	bcf_isec.wait()


def run_analysis(working_dir, analysis_dict, project):
	for pedigree in analysis_dict:
		analysis_info = analysis_dict[pedigree][0]
		sample_fq_dict = analysis_dict[pedigree][1]
		print pedigree
		print analysis_info
		print sample_fq_dict
		amp_bed = analysis_info[3]
		# '''
		##map and make pairs into single reads if required
		align_with_bwa_pe(sample_fq_dict)
		##call variants on non-stiched bams
		make_list_of_bams(sample_fq_dict, bwa_bam_suffix, bam_list)
		variant_calling_freebayes(bam_list, pedigree, freebayes_vcf_suffix, amp_bed)
		variant_calling_gatk_ug(bam_list, pedigree, gatk_vcf_suffix, amp_bed)
		remove_intermediate_files()
		intesect_two_vcf_files(pedigree + '.gatk.vcf.gz', pedigree + '.freebayes.vcf.gz', pedigree + '.intersected_vcfs')
		##get coverage info
		bams = [i + bwa_bam_suffix for i in sample_fq_dict]
		print bams
		cov_file = pedigree + bwa_cov_suffix
		get_coverage_data(bams, cov_file, amp_bed, 'no')
		
		##run gemini on non-stictched bams
		dobyns_amplicon_gemini_pipeline_cybertron_v0.standard_gemini_protocol(working_dir, pedigree, analysis_info[0])
		# '''



		##run stitcher and call variants etc on stitched bams....
		'''
		run_stitcher_on_bam(sample_fq_dict)
		make_list_of_bams(sample_fq_dict, stitcher_bam_suffix, bam_list)
		variant_calling_freebayes(bam_list, pedigree, freebayes_stitcher_vcf_suffix, amp_bed)
		variant_calling_gatk_ug(bam_list, pedigree, gatk_stitcher_vcf_suffix, amp_bed)
		remove_intermediate_files()

		##get coverage info
		bams = [i + stitcher_bam_suffix for i in sample_fq_dict]
		cov_file = pedigree + bwa_stitcher_cov_suffix
		get_coverage_data(bams, cov_file, amp_bed, 'yes')
		'''

def run_gemini(working_dir, analysis_dict, project):
	for pedigree in analysis_dict:
		analysis_info = analysis_dict[pedigree][0]
		sample_fq_dict = analysis_dict[pedigree][1]
		print pedigree
		print analysis_info
		print sample_fq_dict
		amp_bed = analysis_info[3]
		##run gemini on non-stictched bams
		dobyns_amplicon_gemini_pipeline_cybertron_v0.standard_gemini_protocol(working_dir, pedigree, analysis_info[0])


