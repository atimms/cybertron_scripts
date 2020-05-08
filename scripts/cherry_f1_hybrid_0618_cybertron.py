#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import numpy
import scipy
# import pytables

##note
##must have environment set up i.e.
##add local_python/2.7.14 and rm biobuilds
##for lapels to run needed old version of pysam: conda install pysam=0.8.4
##for find_intersected_snps needed pytables version 2 so:conda install pytables=2.4.0
##to install/use macs use: conda install -c auto macs 
##to use homer: module load homer/4.9.1

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/misc/cherry_f1_hybrid_0618'
os.chdir(working_dir)
threads = '15'


##genomes - dl from http://www.csbio.unc.edu/CCstatus/index.py?run=Pseudo
bl6_fa = 'C57BL6J_b38.fa'
spret_fa = 'SPRETEiJ_b38_f.fa'
bl6_name = bl6_fa.split('.')[0]
spret_name = spret_fa.split('.')[0]
GRCm38_fa = '/home/atimms/ngs_data/references/mm10/GRCm38_68.fa'
GRCm38_name = 'GRCm38'
GRCm38_name_and_dir = '/home/atimms/ngs_data/references/mm10/GRCm38_68'
GRCm38_info = '/home/atimms/ngs_data/references/mm10/GRCm38_68.chromInfo.txt.gz'
spret_snp_vcf = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz'
spret_passed_snp_vcf = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.passed.vcf.gz'
spret_passed_bed = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.passed.bed'
spret_passed_aaseq = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.passed.txt'
spret_snp_dir = 'spret_vcfs'
spret_snp_wasp_dir = 'spret_wasp_snps'
mm10_fa = '/home/atimms/ngs_data/references/mm10/mm10.fa'
mm10_name = 'mm10'
mm10_name_and_dir = '/home/atimms/ngs_data/references/mm10/mm10'

##programs
bowtie2_index_builder = '/home/atimms/programs/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build'
bowtie2 = '/home/atimms/programs/bowtie2-2.3.4.1-linux-x86_64/bowtie2'
bwa = '/home/atimms/programs/bwa-0.7.17/bwa'
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


##methods

def combine_fqs(fqs_to_combine, fq_name):
	with open(fq_name, 'w') as cfq_fh:
		print "combining %s files named %s to make file %s"%(len(fqs_to_combine), fqs_to_combine, fq_name)
		cat_files = subprocess.Popen(['cat'] + fqs_to_combine, stdout=cfq_fh)
		cat_files.wait()

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
	if len(fq_list) == 1:
		single_end_bowtie2 = subprocess.Popen([bowtie2, '-x', fa_name, '--dovetail', '-p', threads, '-U', fq_list[0]], stdout=subprocess.PIPE)
		st_sam_bam = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-q', '10', '-o', temp_bam, '-'], stdin=single_end_bowtie2.stdout)
		st_sam_bam.wait()
	elif len(fq_list) == 2:
		print('paired end mapping')
		paired_end_bowtie2 = subprocess.Popen([bowtie2, '-x', fa_name, '-p', threads, '-1', fq_list[0], '-2', fq_list[1]], stdout=subprocess.PIPE)
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
	mod_file = fa_name + '.mod'
	# tabix_split_vcf = subprocess.Popen(['pylapels', '-n', '-o', out_bam, mod_file, in_bam])
	tabix_split_vcf = subprocess.Popen(['pylapels', '-p', '15', '-o', out_bam, mod_file, in_bam])
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

##combine fastq files for f1 atac seq
'''
for sample in ['ATAC-SPREB-1', 'ATAC-SPREB-2', 'ATAC-SPREB-3']:
	for read_no in ['R1', 'R2']:
		combine_fqs(sorted(glob.glob(sample + '*' + read_no + '*001.fastq.gz')), sample + '.' + read_no + '.fastq.gz')
'''

##make index files for all genomes
'''
build_bowtie2_indexes(bl6_fa, bl6_name)
build_bowtie2_indexes(spret_fa, spret_name)
#build_bowtie2_indexes(mm10_fa, mm10_name_and_dir)
'''

##run wasp and then map: from https://github.com/bmvdgeijn/WASP/tree/master/mapping
##Convert SNP data to HDF5 format 
# get_passed_snps(spret_snp_vcf, spret_passed_snp_vcf)
# split_vcf_by_chr(spret_passed_snp_vcf, spret_name, spret_snp_dir)
##doesn't work due to wildcard... run manually i.e. ~/programs/WASP/snp2h5/snp2h5 --chrom ~/ngs_data/references/mm10/GRCm38_68.chromInfo.txt.gz --format vcf --haplotype SPRETEiJ_b38.haplotypes.h5 --snp_index SPRETEiJ_b38.snp_index.h5 --snp_tab SPRETEiJ_b38.snp_tab.h5 SPRETEiJ_b38_f.split*vcf
# make_hdf5_files_from_snps(spret_name)
##as aren't phased use text based file moved to new dir
# make_txt_files_from_snps(spret_snp_dir, spret_snp_wasp_dir)

##get snps for allele specific reads
# convert_vcf_to_aaseq_snps(spret_passed_snp_vcf, spret_passed_bed, spret_passed_aaseq)



##for parental atac seq
parental_fq_dict =  {'ATAC_Bl6J_m1': ['ATAC-Bl6J_ret_adult_m1.fastq.gz', 'bl6'], 
		'ATAC_Bl6J_m2': ['ATAC-Bl6J_ret_adult_m2.fastq.gz', 'bl6'], 
		'ATAC_Spret_M1': ['ATAC-Spret_ret_adult_M1.fastq.gz', 'spret'], 
		'ATAC_Spret_M3': ['ATAC-Spret_ret_adult_M3.fastq.gz', 'spret']}
# parental_fq_dict =  {'test2': ['test2.fastq.gz', 'bl6'], 'test3': ['test3.fastq.gz', 'spret']}

for sample in parental_fq_dict:
	strain = parental_fq_dict[sample][1]
	if strain == 'bl6':
		ref_name = bl6_name
	elif strain == 'spret':
		ref_name = spret_name
	##1 fq file and fist bam
	fq_file = parental_fq_dict[sample][0]
	bt2_bam = sample + '.' + ref_name + '.bowtie2.bam'
	##2 mapped to mm10
	lapels_bam = sample + '.' + ref_name + '.bowtie2.mm10l.bam'
	##3 snps that intersect 
	snp_int_dir = sample + '.' + ref_name + 'int_snps/'
	##4 remap fq
	remap_fq = snp_int_dir + sample + '.' + ref_name + '.remap.fq.gz'
	remapped_bam = sample + '.' + ref_name + '.bowtie2_remap.bam'
	##5 convert remapped bam to mm10
	remapped_bam_mm10 = sample + '.' + ref_name + '.bowtie2_remap.mm10l.bam'
	##6 filter reads
	remap_bam = snp_int_dir + sample + '.' + ref_name + '.to.remap.bam'
	remap_and_keep_bam = snp_int_dir + sample + '.' + ref_name + '.to_keep.bam'
	##7 merge so original kept bam and the merged bam
	kept_bam = snp_int_dir + sample + '.' + ref_name + '.keep.bam'
	merge_bam = snp_int_dir + sample + '.' + ref_name + '.merged.bam'
	##8 rmdups 
	final_bam = sample + '.' + ref_name + '.wasp.bam'

	'''
	##1. map reads
	##using bowtie2
	map_using_bowtie2(sample, [fq_file], ref_name, bt2_bam)
	##2. convert to mm10 coordinates 
	align_bams_to_mm10_lapels(bt2_bam, lapels_bam, ref_name)
	##3. find_intersecting_snps ##only works on spret
	find_intersecting_snps_single_end(sample, ref_name, lapels_bam, snp_int_dir, spret_snp_wasp_dir)
	##s4. remap reads that overlap a snp 
	map_using_bowtie2(sample, [remap_fq], ref_name, remapped_bam)
	##s5. convert to mm10
	align_bams_to_mm10_lapels(remapped_bam, remapped_bam_mm10, ref_name)
	##s6. filter_remapped_reads
	filter_remapped_reads_sort(sample, ref_name, '', remap_bam, remapped_bam_mm10, remap_and_keep_bam)
	##s7. merge all reads
	merge_reads(sample, ref_name, kept_bam, remap_and_keep_bam, merge_bam)
	##s8. rmdups
	remove_duplicates_single(sample, ref_name, merge_bam, final_bam)
	##s9 make tag dirs for homer
	make_tag_dirs(sample, final_bam)
	'''
##run homer getDifferentialPeaksReplicates
bl6_samples = ['ATAC_Bl6J_m1', 'ATAC_Bl6J_m2']
spret_samples = ['ATAC_Spret_M1', 'ATAC_Spret_M3']
results_file_prefix = 'parental.'
homer_getDifferentialPeaksReplicates(bl6_samples, spret_samples, results_file_prefix)

##for f1 data

atac_f1_fq_dict = {'ATAC-SPREB-1':['ATAC-SPREB-1.R1.fastq.gz', 'ATAC-SPREB-1.R2.fastq.gz'], 
		'ATAC-SPREB-2':['ATAC-SPREB-2.R1.fastq.gz', 'ATAC-SPREB-2.R2.fastq.gz'], 
		'ATAC-SPREB-3':['ATAC-SPREB-3.R1.fastq.gz', 'ATAC-SPREB-3.R2.fastq.gz']}

##for test
# atac_f1_fq_dict = {'test1':['test1.fastq.gz', 'test2.fastq.gz']}

for sample in atac_f1_fq_dict:
	##file and dir names
	##1 bams
	bam_bl6 = sample + '.' + bl6_name + '.bowtie2.bam'
	bam_spret = sample + '.' + spret_name + '.bowtie2.bam'
	##2 bams
	bam_bl6_mm10 = sample + '.' + bl6_name + '.mm10l.bowtie2.bam'
	bam_spret_mm10 = sample + '.' + spret_name + '.mm10l.bowtie2.bam'
	##3 bams
	bam_bl6_mm10_fixed = sample + '.' + bl6_name + '.mm10l.bowtie2.fixed.bam'
	bam_spret_mm10_fixed = sample + '.' + spret_name + '.mm10l.bowtie2.fixed.bam'
	##4 out_prefix
	bl6_allele = sample + '.' + bl6_name
	spret_allele = sample + '.' + spret_name
	##5 in bams and out_dir
	bl6_hap1_bam = bl6_allele + '_hap1.bam'
	spret_hap2_bam = spret_allele + '_hap2.bam'
	bl6_intsnp_dir = sample + '.' + bl6_name + '_hap1.int_snps/'
	spret_intsnp_dir = sample + '.' + spret_name + '_hap2.int_snps/'
	##6 in fqs and out bam
	bl6_remap_fqs = [bl6_intsnp_dir + sample + '.' + bl6_name + '.remap.fq1.gz', bl6_intsnp_dir + sample + '.' + bl6_name + '.remap.fq2.gz']
	spret_remap_fqs = [spret_intsnp_dir + sample + '.' + spret_name + '.remap.fq1.gz', spret_intsnp_dir + sample + '.' + spret_name + '.remap.fq2.gz']
	bl6_remap_bam = sample + '.' + bl6_name + '.bowtie2_remap.bam'
	spret_remap_bam = sample + '.' + spret_name + '.bowtie2_remap.bam'
	##7 bams converted back to mm10
	bl6_remap_bam_mm10 = sample + '.' + bl6_name + '.bowtie2_remap.mm10l.bam'
	spret_remap_bam_mm10 = sample + '.' + spret_name + '.bowtie2_remap.mm10l.bam'
	##8 bam with original reads to map and the the out bams i.e ones to keep
	bl6_to_remap_bam = bl6_intsnp_dir + sample + '.' + bl6_name + '.to.remap.bam'
	spret_to_remap_bam = spret_intsnp_dir + sample + '.' + spret_name + '.to.remap.bam'
	bl6_keep_bam = sample + '.' + bl6_name + '.keep.bam'
	spret_keep_bam = sample + '.' + spret_name + '.keep.bam'
	##9 dedup, so final names
	bl6_final_bam = sample + '.' + bl6_name + '.allele_wasp.bam'
	spret_final_bam = sample + '.' + spret_name + '.allele_wasp.bam'
	'''
	##1. map reads to both genomes
	map_using_bowtie2(sample, atac_f1_fq_dict[sample], bl6_name, bam_bl6)
	map_using_bowtie2(sample, atac_f1_fq_dict[sample], spret_name, bam_spret)
	##2. convert to mm10 coordinates 
	align_bams_to_mm10_lapels(bam_bl6, bam_bl6_mm10, bl6_name)
	align_bams_to_mm10_lapels(bam_spret, bam_spret_mm10, spret_name)
	##3. fix flags of orphaned reads with no mate
	run_samtools_fixmate_get_paired(sample, bl6_name, bam_bl6_mm10, bam_bl6_mm10_fixed)
	run_samtools_fixmate_get_paired(sample, spret_name, bam_spret_mm10, bam_spret_mm10_fixed)
	##4. get allele specific reads
	get_allele_specific_reads(bam_bl6_mm10_fixed, bl6_allele, spret_passed_aaseq)
	get_allele_specific_reads(bam_spret_mm10_fixed, spret_allele, spret_passed_aaseq)
	##5. find_intersecting_snps ##only works on spret snps
	find_intersecting_snps_paired_end(sample, bl6_name, bl6_hap1_bam, bl6_intsnp_dir, spret_snp_wasp_dir)
	find_intersecting_snps_paired_end(sample, spret_name, spret_hap2_bam, spret_intsnp_dir, spret_snp_wasp_dir)
	##6. remap reads that overlap a snp to respective genome
	map_using_bowtie2(sample, bl6_remap_fqs, bl6_name, bl6_remap_bam)
	map_using_bowtie2(sample, spret_remap_fqs, spret_name, spret_remap_bam)
	##7. convert remaped bam to mm10
	align_bams_to_mm10_lapels(bl6_remap_bam, bl6_remap_bam_mm10, bl6_name)
	align_bams_to_mm10_lapels(spret_remap_bam, spret_remap_bam_mm10, spret_name)
	##8. filter_remapped_reads ##no need to merge as no 'kept reads'
	filter_remapped_reads_sort(sample, bl6_name, bl6_fa, bl6_to_remap_bam, bl6_remap_bam_mm10, bl6_keep_bam)
	filter_remapped_reads_sort(sample, spret_name, spret_fa, spret_to_remap_bam, spret_remap_bam_mm10, spret_keep_bam)
	##9. rmdups
	remove_duplicates_paired(sample, bl6_name, bl6_keep_bam, bl6_final_bam)
	remove_duplicates_paired(sample, spret_name, spret_keep_bam, spret_final_bam)
	'''
	##s9 make tag dirs for homer
	make_tag_dirs(sample + '_' + bl6_name, bl6_final_bam)
	make_tag_dirs(sample + '_' + spret_name, spret_final_bam)

##run homer getDifferentialPeaksReplicates
bl6_samples = ['ATAC-SPREB-1_C57BL6J_b38', 'ATAC-SPREB-2_C57BL6J_b38', 'ATAC-SPREB-3_C57BL6J_b38']
spret_samples = ['ATAC-SPREB-1_SPRETEiJ_b38_f', 'ATAC-SPREB-2_SPRETEiJ_b38_f', 'ATAC-SPREB-3_SPRETEiJ_b38_f']
results_file_prefix = 'f1_spreb.'
homer_getDifferentialPeaksReplicates(bl6_samples, spret_samples, results_file_prefix)

