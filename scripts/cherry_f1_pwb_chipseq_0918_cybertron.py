#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import numpy
import scipy
# import pytables

##note
'''
analysis of pwb chipseq K27ac

load modules:
module load local_python/2.7.14
module load R/3.3.3

must have environment set up i.e.
add local_python/2.7.14 and rm biobuilds
for lapels to run needed old version of pysam: conda install pysam=0.8.4
for find_intersected_snps needed pytables version 2 so:conda install pytables=2.4.0
to install/use macs use: conda install -c auto macs 
to use homer: module load homer/4.9.1

'''




##parameters
delim = '\t'
# working_dir = '/home/atimms/ngs_data/misc/cherry_f1_hybrid_0618'
working_dir = '/home/atimms/ngs_data/misc/cherry_f1_hybrid/pwb_H3K27ac_0918'
os.chdir(working_dir)
threads = '18'


##genomes - dl from http://www.csbio.unc.edu/CCstatus/index.py?run=Pseudo
psedogenome_dir = '/home/atimms/ngs_data/references/pseudogenomes/'
bl6_fa = psedogenome_dir + 'C57BL6J_b38.fa'
pwk_fa = psedogenome_dir + 'PWKPhJ_b38_f.fa'
bl6_name = bl6_fa.split('/')[-1].split('.')[0]
pwk_name = pwk_fa.split('/')[-1].split('.')[0]
GRCm38_fa = '/home/atimms/ngs_data/references/mm10/GRCm38_68.fa'
GRCm38_name = 'GRCm38'
GRCm38_name_and_dir = '/home/atimms/ngs_data/references/mm10/GRCm38_68'
GRCm38_info = '/home/atimms/ngs_data/references/mm10/GRCm38_68.chromInfo.txt.gz'
##snps made in rnaseq anaysis 0918
pwk_snp_vcf = psedogenome_dir + 'PWK_PhJ.mgp.v5.snps.dbSNP142.vcf.gz'
pwk_passed_snp_vcf = psedogenome_dir + 'PWK_PhJ.mgp.v5.snps.dbSNP142.passed.vcf.gz'
pwk_passed_bed = psedogenome_dir + 'PWK_PhJ.mgp.v5.snps.dbSNP142.passed.bed'
pwk_passed_aaseq = psedogenome_dir +'PWK_PhJ.mgp.v5.snps.dbSNP142.passed.txt'
pwk_snp_dir = psedogenome_dir + 'pwk_vcfs'
pwk_snp_wasp_dir = psedogenome_dir + 'pwk_wasp_snps'
mm10_fa = '/home/atimms/ngs_data/references/mm10/mm10.fa'
mm10_name = 'mm10'
mm10_name_and_dir = '/home/atimms/ngs_data/references/mm10/mm10'

##programs
bowtie2_index_builder = '/home/atimms/programs/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build'
bowtie2 = '/home/atimms/programs/bowtie2-2.3.4.1-linux-x86_64/bowtie2'
bwa = '/home/atimms/programs/bwa-0.7.17/bwa'
star = '/home/atimms/programs/STAR-2.5.3a/bin/Linux_x86_64_static/STAR'
samtools = '/home/atimms/programs/samtools-1.8/samtools'
snp2h5 = '/home/atimms/programs/WASP/snp2h5/snp2h5'
extract_snps = '/home/atimms/programs/WASP/mapping/extract_vcf_snps.sh'
find_intersecting_snps = '/home/atimms/programs/WASP/mapping/find_intersecting_snps.py'
filter_remapped_reads = '/home/atimms/programs/WASP/mapping/filter_remapped_reads.py'
rmdups_se = '/home/atimms/programs/WASP/mapping/rmdup.py'
rmdups_pe = '/home/atimms/programs/WASP/mapping/rmdup_pe.py'
vcf2bed = '/home/atimms/programs/bedops_linux_x86_64-v2.4.35/vcf2bed'
##set up asseq installed program and wrote this wrapper script
extractAsReads_wrapper = '/home/atimms/programs/aaseq_extractAsReads.R'
hisat2_build = '/home/atimms/programs/hisat2-2.1.0/hisat2-build'
hisat2 = '/home/atimms/programs/hisat2-2.1.0/hisat2'
tophat2 = '/home/atimms/programs/tophat-2.1.1.Linux_x86_64/tophat2'



##methods

def build_bowtie2_indexes(fa_file, name):
	bowtie2_index = subprocess.Popen([bowtie2_index_builder, fa_file, name])
	bowtie2_index.wait()

def make_hdf5_files_from_snps(out_prefix):
	# snp_chrs = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '3', '4', '5', '6', '7', '8', '9', 'X', 'Y']
	# for snp_chr in snp_chrs:
	in_vcf =  out_prefix + '.split.' + '*' + '.vcf'
	##  ./snp2h5/snp2h5 --chrom data/ucsc/hg19/chromInfo.txt.gz --format vcf --haplotype haplotypes.h5 --snp_index snp_index.h5 --snp_tab snp_tab.h5 data/1000G/ALL.chr*.vcf.gz
	hdf5_make = subprocess.Popen([snp2h5, '--chrom', GRCm38_info, '--format', 'vcf', '--haplotype', out_prefix + '.haplotypes.h5', '--snp_index', out_prefix + '.snp_index.h5', '--snp_tab', out_prefix + '.snp_tab.h5', in_vcf])
	hdf5_make.wait()

def get_passed_snps(in_vcf, out_vcf):
	##correct filtering??
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-o', out_vcf, '-O', 'z', in_vcf])
	bcftools_filter.wait()


def split_vcf_by_chr(vcf, out_prefix, snp_dir):
	snp_chrs = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '3', '4', '5', '6', '7', '8', '9', 'X', 'Y']
	tabix_vcf = subprocess.Popen(['tabix', vcf])
	tabix_vcf.wait()
	for snp_chr in snp_chrs:
		out_vcf = snp_dir + '/' + out_prefix + '.split.chr' + snp_chr + '.vcf'
		with open(out_vcf, 'w') as vcf_fh:
			tabix_split_vcf = subprocess.Popen(['tabix', '-h', vcf, snp_chr], stdout=vcf_fh)
			tabix_split_vcf.wait()

def make_txt_files_from_snps(in_dir, out_dir):
	vcfs = glob.glob(in_dir + '/*.vcf')
	for vcf in vcfs:
		bgzip_cmd = subprocess.Popen(['bgzip', vcf])
		bgzip_cmd.wait()
	run_extract_snps = subprocess.Popen([extract_snps, in_dir, out_dir])
	run_extract_snps.wait()

def map_using_bowtie2(name, fq_list, fa_name, out_bam):
	##bowtie2 -x /n/scratch2/el116/external/genomes/index/index/C57BL6J_b38 --dovetail -p 3 -1 ${describer}_R1_rep.fastq -2 ${describer}_R2_rep.fastq -S ${describer}_C57BL6.sam
	temp_bam = name + '.' + fa_name + '.temp.bam'
	index_prefix = psedogenome_dir + fa_name
	if len(fq_list) == 1:
		single_end_bowtie2 = subprocess.Popen([bowtie2, '-x', index_prefix, '--dovetail', '-p', threads, '-U', fq_list[0]], stdout=subprocess.PIPE)
		st_sam_bam = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-q', '10', '-o', temp_bam, '-'], stdin=single_end_bowtie2.stdout)
		st_sam_bam.wait()
	elif len(fq_list) == 2:
		print('paired end mapping')
		paired_end_bowtie2 = subprocess.Popen([bowtie2, '-x', index_prefix, '-p', threads, '-1', fq_list[0], '-2', fq_list[1]], stdout=subprocess.PIPE)
		st_sam_bam = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-q', '10', '-o', temp_bam, '-'], stdin=paired_end_bowtie2.stdout)
		# st_sam_bam = subprocess.Popen([samtools, 'view', '-b', '-@', '2', '-o', temp_bam, '-'], stdin=paired_end_bowtie2.stdout)
		st_sam_bam.wait()
	else:
		print("don't understand the fq list:", fq_list)
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', out_bam, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', out_bam])
	st_index.wait()


def align_bams_to_mm10_lapels(in_bam, out_bam, fa_name):
	mod_file = psedogenome_dir + fa_name + '.mod'
	# tabix_split_vcf = subprocess.Popen(['pylapels', '-n', '-o', out_bam, mod_file, in_bam])
	tabix_split_vcf = subprocess.Popen(['pylapels', '-p', '20', '-o', out_bam, mod_file, in_bam])
	# tabix_split_vcf = subprocess.Popen(['pylapels', '-o', out_bam, mod_file, in_bam])
	tabix_split_vcf.wait()

def run_samtools_fixmate_get_paired(name, fa_name, in_bam, out_bam):
	##samtools fixmate ${describer}_C57BL6.bam ${describer}_C57BL6_fixed.bam
	temp_bam = name + '.' + fa_name + 'temp2.bam'
	# temp2_bam = name + '.' + fa_name + 'temp2.bam'
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-n', '-O', 'bam', '-T', name, '-o', temp_bam, '-@', '10', '-m', '10G', in_bam])
	st_sort_pe.wait()
	st_fix = subprocess.Popen([samtools, 'fixmate', temp_bam, out_bam])
	st_fix.wait()

def get_allele_specific_reads(in_bam, out_prefix, snp_file):
	#Rscript /groups/neuroduo/emi/ChIP-Seq/tools/extractasReads.R $PWD ${describer}_fixed.bam /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/PWK_PhJ/mgp.v5.merged.snps.dbSNP142_PASS_PWK_PhJ_foraseq.txt ${describer}_fixed
	aaseq_asr = subprocess.Popen(['Rscript', extractAsReads_wrapper, working_dir, in_bam, snp_file, out_prefix])
	aaseq_asr.wait()

def find_intersecting_snps_single_end(name, fa_name, in_bam, out_dir, snp_dir):
	##"python /groups/neuroduo/emi/ChIP-Seq/tools/wasp/mapping/find_intersecting_snps.py --is_paired_end --snp_tab /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_snp_tab_het.h5 --snp_index /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_snp_index_het.h5 --haplotype /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_haplotypes_het.h5 --max_seqs 1048576 --max_snps 20 $PWD/${describer}.bam" 
	temp_bam = name + '.' + fa_name + '.bam'
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', temp_bam, '-@', '10', '-m', '10G', in_bam])
	st_sort_pe.wait()
	# int_snps = subprocess.Popen(['python', find_intersecting_snps, '--is_paired_end', '-s', '--snp_dir', snp_dir, '--max_seqs', '1048576', '--max_snps', '20', '--output_dir', out_dir, temp_bam])
	int_snps = subprocess.Popen(['python', find_intersecting_snps, '-s', '--snp_dir', snp_dir, '--output_dir', out_dir, temp_bam])
	int_snps.wait()

def find_intersecting_snps_paired_end(name, fa_name, in_bam, out_dir, snp_dir):
	##"python /groups/neuroduo/emi/ChIP-Seq/tools/wasp/mapping/find_intersecting_snps.py --is_paired_end --snp_tab /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_snp_tab_het.h5 --snp_index /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_snp_index_het.h5 --haplotype /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_haplotypes_het.h5 --max_seqs 1048576 --max_snps 20 $PWD/${describer}.bam" 
	temp_bam = name + '.' + fa_name + '.bam'
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', temp_bam, '-@', '10', '-m', '10G', in_bam])
	st_sort_pe.wait()
	# int_snps = subprocess.Popen(['python', find_intersecting_snps, '--is_paired_end', '-s', '--snp_dir', snp_dir, '--max_seqs', '1048576', '--max_snps', '20', '--output_dir', out_dir, temp_bam])
	int_snps = subprocess.Popen(['python', find_intersecting_snps, '--is_paired_end', '-s', '--snp_dir', snp_dir, '--output_dir', out_dir, temp_bam])
	int_snps.wait()

def filter_remapped_reads_sort(name, fa_name, fa_to_map_against, to_remap_bam, remap_bam, remap_keep_bam):
	##filter_remapped_reads.py [-h] to_remap_bam remap_bam keep_bam
	temp_bam = name + '.' + fa_name + '.temp3.bam'
	filt_remapped = subprocess.Popen(['python', filter_remapped_reads, to_remap_bam, remap_bam, temp_bam])
	filt_remapped.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', remap_keep_bam, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', remap_keep_bam])
	st_index.wait()


def remove_duplicates_paired(name, fa_name, in_bam, final_bam):
	temp_bam = name + '.' + fa_name + '.temp4.bam'
	rmdup_reads = subprocess.Popen(['python', rmdups_pe, in_bam, temp_bam])
	rmdup_reads.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', final_bam, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', final_bam])
	st_index.wait()

def remove_duplicates_single(name, fa_name, in_bam, final_bam):
	temp_bam = name + '.' + fa_name + '.temp6.bam'
	rmdup_reads = subprocess.Popen(['python', rmdups_se, in_bam, temp_bam])
	rmdup_reads.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', final_bam, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', final_bam])
	st_index.wait()

def merge_reads(name, fa_name, keep_bam, remap_keep_bam, merged_bam):
	temp_bam = name + '.' + fa_name + '.temp5.bam'
	st_merge = subprocess.Popen([samtools, 'merge', temp_bam, keep_bam, remap_keep_bam])
	st_merge.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', merged_bam, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', merged_bam])
	st_index.wait()


def make_tag_dirs(name, bam):
	outdir = name + '.tag_dir'
	mk_tag_dir = subprocess.Popen(['makeTagDirectory', outdir, bam])
	mk_tag_dir.wait()

def homer_getDifferentialPeaksReplicates(bl6_names, spret_names, outfile_prefix):
	bl6_tag_dirs = [i + '.tag_dir' for i in bl6_names]
	spret_tag_dirs = [i + '.tag_dir' for i in spret_names]
	bl6_outfile = outfile_prefix + 'bl6.outputPeaks.txt'
	spret_outfile = outfile_prefix + 'spret.outputPeaks.txt'
	with open(bl6_outfile, 'w') as bl6_fh:
		get_diff_peaks_bl6 = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + bl6_tag_dirs + ['-b'] + spret_tag_dirs + ['-genome', 'mm10'], stdout=bl6_fh)
		get_diff_peaks_bl6.wait()
	with open(spret_outfile, 'w') as spret_fh:
		get_diff_peaks_spret = subprocess.Popen(['getDifferentialPeaksReplicates.pl', '-t'] + spret_tag_dirs + ['-b'] + bl6_tag_dirs + ['-genome', 'mm10'], stdout=spret_fh)
		get_diff_peaks_spret.wait()
##run methods
# 
##stuff to do once


##make index files for all genomes
# '''
# build_bowtie2_indexes(bl6_fa, bl6_name)
# build_bowtie2_indexes(spret_fa, spret_name)
# build_bowtie2_indexes(pwk_fa, psedogenome_dir + pwk_name)
# '''

##make SNP files for PWK
##Convert SNP data to HDF5 format 
# get_passed_snps(pwk_snp_vcf, pwk_passed_snp_vcf)
# split_vcf_by_chr(pwk_passed_snp_vcf, pwk_name, pwk_snp_dir)
#make snp files for wasp (don't use hd5)
# make_txt_files_from_snps(pwk_snp_dir, pwk_snp_wasp_dir)
##get snps for allele specific reads
# convert_vcf_to_aaseq_snps(pwk_passed_snp_vcf, pwk_passed_bed, pwk_passed_aaseq



##run methods
pwb_fq_dict = {'PWB1-ret-input':['PWB1-ret-input.R1.fq.gz'], 'PWB1-ret-K27ac':['PWB1-ret-K27ac.R1.fq.gz'], 
		'PWB2-ret-K27ac':['PWB2-ret-K27ac.R1.fq.gz'], 'PWB3-ret-K27ac':['PWB3-ret-K27ac.R1.fq.gz']}

##for test
# pwb_fq_dict = {'test1':['test1.fastq.gz']}

for sample in pwb_fq_dict:
	##file and dir names
	##1 bams
	bam_bl6 = sample + '.' + bl6_name + '.bowtie2.bam'
	bam_pwk = sample + '.' + pwk_name + '.bowtie2.bam'
	##2 bams
	bam_bl6_mm10 = sample + '.' + bl6_name + '.mm10l.bowtie2.bam'
	bam_pwk_mm10 = sample + '.' + pwk_name + '.mm10l.bowtie2.bam'
	##3 bams
	bam_bl6_mm10_fixed = sample + '.' + bl6_name + '.mm10l.bowtie2.fixed.bam'
	bam_pwk_mm10_fixed = sample + '.' + pwk_name + '.mm10l.bowtie2.fixed.bam'
	##4 out_prefix
	bl6_allele = sample + '.' + bl6_name
	pwk_allele = sample + '.' + pwk_name
	##5 in bams and out_dir
	bl6_hap1_bam = bl6_allele + '_hap1.bam'
	pwk_hap2_bam = pwk_allele + '_hap2.bam'
	bl6_intsnp_dir = sample + '.' + bl6_name + '_hap1.int_snps/'
	pwk_intsnp_dir = sample + '.' + pwk_name + '_hap2.int_snps/'
	##6 in fqs and out bam
	bl6_remap_fqs = [bl6_intsnp_dir + sample + '.' + bl6_name + '.remap.single.fq.gz']
	pwk_remap_fqs = [pwk_intsnp_dir + sample + '.' + pwk_name + '.remap.single.fq.gz']
	bl6_remap_bam = sample + '.' + bl6_name + '.bowtie2_remap.bam'
	pwk_remap_bam = sample + '.' + pwk_name + '.bowtie2_remap.bam'
	##7 bams converted back to mm10
	bl6_remap_bam_mm10 = sample + '.' + bl6_name + '.bowtie2_remap.mm10l.bam'
	pwk_remap_bam_mm10 = sample + '.' + pwk_name + '.bowtie2_remap.mm10l.bam'
	##8 bam with original reads to map and the the out bams i.e ones to keep
	bl6_to_remap_bam = bl6_intsnp_dir + sample + '.' + bl6_name + '.to.remap.bam'
	pwk_to_remap_bam = pwk_intsnp_dir + sample + '.' + pwk_name + '.to.remap.bam'
	bl6_keep_bam = sample + '.' + bl6_name + '.keep.bam'
	pwk_keep_bam = sample + '.' + pwk_name + '.keep.bam'
	##9 dedup, so final names
	bl6_final_bam = sample + '.' + bl6_name + '.allele_wasp.bam'
	pwk_final_bam = sample + '.' + pwk_name + '.allele_wasp.bam'
	# '''
	# ##1. map reads to both genomes
	# map_using_bowtie2(sample, pwb_fq_dict[sample], bl6_name, bam_bl6)
	# map_using_bowtie2(sample, pwb_fq_dict[sample], pwk_name, bam_pwk)
	
	# ##2. convert to mm10 coordinates 
	# align_bams_to_mm10_lapels(bam_bl6, bam_bl6_mm10, bl6_name)
	# align_bams_to_mm10_lapels(bam_pwk, bam_pwk_mm10, pwk_name)

	# ##3. fix flags of orphaned reads with no mate
	# run_samtools_fixmate_get_paired(sample, bl6_name, bam_bl6_mm10, bam_bl6_mm10_fixed)
	# run_samtools_fixmate_get_paired(sample, pwk_name, bam_pwk_mm10, bam_pwk_mm10_fixed)

	# ##4. get allele specific reads
	# get_allele_specific_reads(bam_bl6_mm10_fixed, bl6_allele, pwk_passed_aaseq)
	# get_allele_specific_reads(bam_pwk_mm10_fixed, pwk_allele, pwk_passed_aaseq)

	# ##5. find_intersecting_snps ##only works on pwk snps
	# find_intersecting_snps_paired_end(sample, bl6_name, bl6_hap1_bam, bl6_intsnp_dir, pwk_snp_wasp_dir)
	# find_intersecting_snps_paired_end(sample, pwk_name, pwk_hap2_bam, pwk_intsnp_dir, pwk_snp_wasp_dir)

	# ##6. remap reads that overlap a snp to respective genome
	# map_using_bowtie2(sample, bl6_remap_fqs, bl6_name, bl6_remap_bam)
	# map_using_bowtie2(sample, pwk_remap_fqs, pwk_name, pwk_remap_bam)

	# ##7. convert remaped bam to mm10
	# align_bams_to_mm10_lapels(bl6_remap_bam, bl6_remap_bam_mm10, bl6_name)
	# align_bams_to_mm10_lapels(pwk_remap_bam, pwk_remap_bam_mm10, pwk_name)

	##8. filter_remapped_reads ##no need to merge as no 'kept reads'
	filter_remapped_reads_sort(sample, bl6_name, bl6_fa, bl6_to_remap_bam, bl6_remap_bam_mm10, bl6_final_bam)
	filter_remapped_reads_sort(sample, pwk_name, pwk_fa, pwk_to_remap_bam, pwk_remap_bam_mm10, pwk_final_bam)

	##9. rmdups --- doesn't work so didn't use
	# remove_duplicates_single(sample, bl6_name, bl6_keep_bam, bl6_final_bam)
	# remove_duplicates_single(sample, pwk_name, pwk_keep_bam, pwk_final_bam)
	
	# ##s9 make tag dirs for homer
	# make_tag_dirs(sample + '_' + bl6_name, bl6_final_bam)
	# make_tag_dirs(sample + '_' + spret_name, spret_final_bam)
	# '''

##run homer getDifferentialPeaksReplicates
# bl6_samples = ['ATAC-SPREB-1_C57BL6J_b38', 'ATAC-SPREB-2_C57BL6J_b38', 'ATAC-SPREB-3_C57BL6J_b38']
# spret_samples = ['ATAC-SPREB-1_SPRETEiJ_b38_f', 'ATAC-SPREB-2_SPRETEiJ_b38_f', 'ATAC-SPREB-3_SPRETEiJ_b38_f']
# results_file_prefix = 'f1_spreb.'
# homer_getDifferentialPeaksReplicates(bl6_samples, spret_samples, results_file_prefix)

