#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import numpy
import scipy
# import pytables

bowtie2_index = subprocess.Popen(['which', 'python'])
bowtie2_index.wait()

##note
##must have environment set up i.e.
##add local_python/2.7.14 and rm biobuilds
##for lapels to run needed old version of pysam: conda install pysam=0.8.4
##for find_intersected_snps needed pytables version 2 so:conda install pytables=2.4.0

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/misc/cherry_f1_hybrid_0318'
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

def map_using_bowtie2(name, fq_list, fa_name):
	##bowtie2 -x /n/scratch2/el116/external/genomes/index/index/C57BL6J_b38 --dovetail -p 3 -1 ${describer}_R1_rep.fastq -2 ${describer}_R2_rep.fastq -S ${describer}_C57BL6.sam
	print name, fq_list, fa_name
	if 'mm10' in fa_name:
		final_name = 'mm10'
	else:
		final_name = fa_name
	temp_bam = name + '.' + final_name + '.temp.bam'
	out_bam_sorted = name + '.' + final_name + '.bowtie2.bam'
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
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', out_bam_sorted, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', out_bam_sorted])
	st_index.wait()

def map_using_bwa_mem(name, fq_list, fasta):
	fa_name = fasta.split('.')[0]
	out_bam = name + '.' + fa_name + '.bwa_mem.bam'
	if len(fq_list) == 1:
		bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '15', fasta, fq_list[0]], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-q', '10', '-o', out_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()
	if len(fq_list) == 2:
		bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '15', fasta, fq_list[0], fq_list[1]], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-q', '10', '-o', out_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()

def map_using_bwa_aln(name, fq_list, fasta):
	##bwa aln  -f "Hu1-ATAC-ret_20K_60min_R1_b1_p_trimmed_rep_hg19.sai"  -q 0  -t 4  -n 2  -k 2  -l 32  -e -1  -o 0 fa fq
	##bwa sampe  -f "Hu1-ATAC-ret_20K_60min_R1_R2_b1_p_trimmed_rep_hg19.sam"  -n 2  -r "@RG\tID:TJC_Retina\tSM:Hu_ATAC_Retina\tLB:Hu1_5-ATAC-ret_20K_60min\tDS:Hu1-ATAC-ret_20K_60min"  "/groups/neuroduo/david/Genomes/human_37_1/bwa/hg19_ref_GRCh37_chr.1-Y.fa"  "Hu1-ATAC-ret_20K_60min_R1_b1_p_trimmed_rep_hg19.sai" "Hu1-ATAC-ret_20K_60min_R2_b1_p_trimmed_rep_hg19.sai"  "Hu1-ATAC-ret_20K_60min_R1_b1_p_trimmed_rep.fastq" "Hu1-ATAC-ret_20K_60min_R2_b1_p_trimmed_rep.fastq" 
	fa_name = fasta.split('.')[0]
	out_sai = name + '.' + fa_name + '.sai'
	out_bam = name + '.' + fa_name + '.bwa_aln.bam'
	if len(fq_list) == 1:
		bwa_aln = subprocess.Popen([bwa, 'aln', '-t', '20', '-f', out_sai, fasta, fq_list[0]], stdout=subprocess.PIPE)
		bwa_aln.wait()
		bwa_samse = subprocess.Popen([bwa, 'samse', fasta, out_sai, fq_list[0]], stdout=subprocess.PIPE)
		st_sam_bam_se = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-q', '10', '-o', out_bam, '-'], stdin=bwa_samse.stdout)
		st_sam_bam_se.wait()

	if len(fq_list) == 2:
		out_sai_r1 = name + '.' + fa_name + '.r1.sai'
		out_sai_r2 = name + '.' + fa_name + '.r2.sai'
		bwa_aln_r1 = subprocess.Popen([bwa, 'aln', '-t', '20', '-f', out_sai_r1, fasta, fq_list[0]], stdout=subprocess.PIPE)
		bwa_aln_r1.wait()
		bwa_aln_r2 = subprocess.Popen([bwa, 'aln', '-t', '20', '-f', out_sai_r2, fasta, fq_list[1]], stdout=subprocess.PIPE)
		bwa_aln_r2.wait()
		bwa_sampe = subprocess.Popen([bwa, 'sampe', fasta, out_sai_r1, out_sai_r2, fq_list[0], fq_list[1]], stdout=subprocess.PIPE)
		st_sam_bam_se = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-q', '10', '-o', out_bam, '-'], stdin=bwa_samse.stdout)
		st_sam_bam_se.wait()

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

def align_bams_to_mm10_lapels(name, fa_name):
	in_bam = name + '.' + fa_name + '.bowtie2.bam'
	out_bam = name + '.' + fa_name + '_mm10lap.bowtie2.bam'
	# in_bam = name + '.test.bam'
	# out_bam = name + '.test_lap.bam'
	mod_file = fa_name + '.mod'
	# tabix_split_vcf = subprocess.Popen(['pylapels', '-n', '-o', out_bam, mod_file, in_bam])
	tabix_split_vcf = subprocess.Popen(['pylapels', '-p', '15', '-o', out_bam, mod_file, in_bam])
	tabix_split_vcf.wait()

def find_intersecting_snps_single_end(name, fa_name, snp_dir):
	##python /groups/neuroduo/emi/ChIP-Seq/tools/wasp/mapping/find_intersecting_snps.py --is_paired_end --snp_tab /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_snp_tab_het.h5 --snp_index /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_snp_index_het.h5 --haplotype /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_haplotypes_het.h5 --max_seqs 1048576 --max_snps 20 $PWD/${describer}.bam
	in_bam = name + '.' + fa_name + '_mm10lap.bowtie2.bam'
	out_dir = name + '.' + fa_name + '.int_snps'
	int_snps = subprocess.Popen(['python', find_intersecting_snps, '-s', '--snp_dir', snp_dir, '--max_seqs', '1048576', '--max_snps', '20', '--output_dir', out_dir, in_bam])
	int_snps.wait()

def find_intersecting_snps_paired_end(name, fa_name, bam_suffix, snp_dir):
	##"python /groups/neuroduo/emi/ChIP-Seq/tools/wasp/mapping/find_intersecting_snps.py --is_paired_end --snp_tab /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_snp_tab_het.h5 --snp_index /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_snp_index_het.h5 --haplotype /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/CAST_EiJ/CAST_EiJ_v5_haplotypes_het.h5 --max_seqs 1048576 --max_snps 20 $PWD/${describer}.bam" 
	in_bam = name + '.' + fa_name + bam_suffix + '.bam'
	temp_bam = name + '.' + fa_name + bam_suffix + '.sorted.bam'
	out_dir = name + '.' + fa_name + '.int_snps'
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', temp_bam, '-@', '10', '-m', '10G', in_bam])
	st_sort_pe.wait()
	# int_snps = subprocess.Popen(['python', find_intersecting_snps, '--is_paired_end', '-s', '--snp_dir', snp_dir, '--max_seqs', '1048576', '--max_snps', '20', '--output_dir', out_dir, temp_bam])
	int_snps = subprocess.Popen(['python', find_intersecting_snps, '--is_paired_end', '-s', '--snp_dir', snp_dir, '--output_dir', out_dir, temp_bam])
	int_snps.wait()

def remap_fq_single(name, fa_name, fa_to_map_against):
	fq  = name + '.' + fa_name + '.int_snps/' + name + '.' + fa_name + '_mm10lap.bowtie2.remap.fq.gz'
	new_name = name + '.' + fa_name + '.int_snps/' + name + fa_name + '_mm10lap.bowtie2.remapped'
	map_using_bowtie2(new_name, [fq], fa_to_map_against)

def remap_fq_paired(name, fa_name, bam_suffix, fa_to_map_against):
	fq1  = name + '.' + fa_name + '.int_snps/' + name + '.' + fa_name + bam_suffix + '.sorted.remap.fq1.gz'
	fq2  = name + '.' + fa_name + '.int_snps/' + name + '.' + fa_name + bam_suffix + '.sorted.remap.fq2.gz'
	new_name = name + '.' + fa_name + '.int_snps/' + name + fa_name + '.' + bam_suffix + '.bowtie2.remapped'
	print(fq1, fq2)
	map_using_bowtie2(new_name, [fq1, fq2], fa_to_map_against)

def filter_remapped_reads_sort(name, fa_name, fa_to_map_against):
	##filter_remapped_reads.py [-h] to_remap_bam remap_bam keep_bam
	to_remap_bam  = name + '.' + fa_name + '.int_snps/' + name + '.' + fa_name + '_mm10lap.bowtie2.to.remap.bam'
	remap_bam = name + '.' + fa_name + '.int_snps/' + name + fa_name + '_mm10lap.bowtie2.remapped.mm10.bowtie2.bam'
	temp_bam = 'temp2.bam'
	remap_keep_bam = name + '.' + fa_name + '.int_snps/' + name + '.' + fa_name + '_mm10lap.bowtie2.remapped_keep.bam'
	filt_remapped = subprocess.Popen(['python', filter_remapped_reads, to_remap_bam, remap_bam, temp_bam])
	filt_remapped.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', remap_keep_bam, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', remap_keep_bam])
	st_index.wait()


def filter_remapped_reads_sort_f1(name, fa_name, bam_suffix, fa_to_map_against):
	##filter_remapped_reads.py [-h] to_remap_bam remap_bam keep_bam
	to_remap_bam  = name + '.' + fa_name + '.int_snps/' + name + fa_name + '.' + bam_suffix + '.bowtie2.remapped.mm10.bowtie2.bam'
	remap_bam = name + '.' + fa_name + '.int_snps/' + name + fa_name + '_mm10lap.bowtie2.remapped.mm10.bowtie2.bam'
	temp_bam = 'temp2.bam'
	remap_keep_bam = name + '.' + fa_name + '.int_snps/' + name + fa_name + '.' + bam_suffix + '.bowtie2.remapped.mm10.bowtie2.bam'
	filt_remapped = subprocess.Popen(['python', filter_remapped_reads, to_remap_bam, remap_bam, temp_bam])
	filt_remapped.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', remap_keep_bam, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', remap_keep_bam])
	st_index.wait()

def merge_reads(name, fa_name, fa_to_map_against):
	keep_bam  = name + '.' + fa_name + '.int_snps/' + name + '.' + fa_name + '_mm10lap.bowtie2.keep.bam'
	remap_keep_bam = name + '.' + fa_name + '.int_snps/' + name + '.' + fa_name + '_mm10lap.bowtie2.remapped_keep.bam'
	temp_bam = name + '.' + fa_name + '.temp2.bam'
	merged_bam = name + '.' + fa_name + '.int_snps/' + name + '.' + fa_name + '_mm10lap.bowtie2.merged.bam'
	st_merge = subprocess.Popen([samtools, 'merge', temp_bam, keep_bam, remap_keep_bam])
	st_merge.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', merged_bam, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', merged_bam])
	st_index.wait()

def remove_duplicates_single(name, fa_name, fa_to_map_against):
	merged_bam = name + '.' + fa_name + '.int_snps/' + name + '.' + fa_name + '_mm10lap.bowtie2.merged.bam'
	final_bam = name + '.' + fa_name + '_mm10lap.wasp.bam'
	temp_bam = name + '.' + fa_name + '.temp3.bam'
	rmdup_reads = subprocess.Popen(['python', rmdups_se, merged_bam, temp_bam])
	rmdup_reads.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', name, '-o', final_bam, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', final_bam])
	st_index.wait()

def convert_vcf_to_aaseq_snps(in_vcf_gz, bed, aaseq_txt):
	## vcf to bed  "/opt/bedops-2.4.14/bin/vcf2bed < ${describer}.vcf > ${describer}_pre.bed"
	##vcf to bed
	in_vcf = in_vcf_gz.rsplit('.', 1)[0]
	# print in_vcf
	# print([vcf2bed, '<', in_vcf])
	'''
	with open(bed, 'w') as bed_fh:
		cat_vcf = subprocess.Popen(['cat', in_vcf], stdout=subprocess.PIPE)
		vcf2bed_bo = subprocess.Popen([vcf2bed, '-'], stdin=cat_vcf.stdout, stdout=bed_fh)
		vcf2bed_bo.wait()
	'''
	##bed to asseq file (with added chr?)
	with open(bed, 'r') as bed_fh, open(aaseq_txt, 'w') as txt_fh:
		for line in bed_fh:
			line = line.split(delim)
			line_out = ['chr' + line[0], line[2], line[5], line[6], '\n']
			txt_fh.write(delim.join(line_out))


def run_samtools_fixmate_get_paired(name, fa_name):
	##samtools fixmate ${describer}_C57BL6.bam ${describer}_C57BL6_fixed.bam
	in_bam = name + '.' + fa_name + '_mm10lap.bowtie2.bam'
	print(in_bam)
	temp_bam = name + '.' + fa_name + 'temp.bam'
	# temp2_bam = name + '.' + fa_name + 'temp2.bam'
	out_bam = name + '.' + fa_name + '_mm10lap.bowtie2.fixed.bam'
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-n', '-O', 'bam', '-T', name, '-o', temp_bam, '-@', '10', '-m', '10G', in_bam])
	st_sort_pe.wait()
	st_fix = subprocess.Popen([samtools, 'fixmate', temp_bam, out_bam])
	st_fix.wait()
	# st_view_pe = subprocess.Popen([samtools, 'view', '-h', '-f', '2', '-O', 'bam', '-o', out_bam, '-@', '10', '-m', '10G', temp2_bam])
	# st_view_pe.wait()

def get_allele_specific_reads(name, fa_name, snp_file):
	#Rscript /groups/neuroduo/emi/ChIP-Seq/tools/extractasReads.R $PWD ${describer}_fixed.bam /groups/neuroduo/emi/ChIP-Seq/mouse/external/genomes/PWK_PhJ/mgp.v5.merged.snps.dbSNP142_PASS_PWK_PhJ_foraseq.txt ${describer}_fixed
	in_bam = name + '.' + fa_name + '_mm10lap.bowtie2.fixed.bam'
	out_prefix =  name + '.' + fa_name + '_mm10lap'
	aaseq_asr = subprocess.Popen(['Rscript', extractAsReads_wrapper, working_dir, in_bam, snp_file, out_prefix])
	aaseq_asr.wait()


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
parental_fq_dict =  {'ATAC_Bl6J_m1': ['ATAC-Bl6J_ret_adult_m1.fastq.gz'], 
		'ATAC_Bl6J_m2': ['ATAC-Bl6J_ret_adult_m2.fastq.gz'], 
		'ATAC_Spret_M1': ['ATAC-Spret_ret_adult_M1.fastq.gz'], 
		'ATAC_Spret_M3': ['ATAC-Spret_ret_adult_M3.fastq.gz']}
# parental_fq_dict =  {'ATAC_Bl6J_m1': ['ATAC-Bl6J_ret_adult_m1.fastq.gz']}
'''
for sample in parental_fq_dict:
	##1. map reads
	##using bowtie2
	map_using_bowtie2(sample, parental_fq_dict[sample], bl6_name)
	map_using_bowtie2(sample, parental_fq_dict[sample], spret_name)
	##2. convert to mm10 coordinates 
	align_bams_to_mm10_lapels(sample, bl6_name)
	align_bams_to_mm10_lapels(sample, spret_name)
	##s3. find_intersecting_snps ##only works on spret
	find_intersecting_snps_single_end(sample, spret_name, spret_snp_wasp_dir)
	find_intersecting_snps_single_end(sample, bl6_name, spret_snp_wasp_dir)
	##s4. remap reads that overlap a snp -- repeat against
	remap_fq_single(sample, spret_name, mm10_name_and_dir)
	remap_fq_single(sample, bl6_name, mm10_name_and_dir)
	##s5. filter_remapped_reads
	filter_remapped_reads_sort(sample, spret_name, mm10_fa)
	filter_remapped_reads_sort(sample, bl6_name, mm10_fa)
	##s6. merge all reads
	merge_reads(sample, spret_name, mm10_fa)
	merge_reads(sample, bl6_name, mm10_fa)
	##s7. rmdups
	remove_duplicates_single(sample, spret_name, mm10_fa)
	remove_duplicates_single(sample, bl6_name, mm10_fa)
'''


##for f1 data

atac_f1_fq_dict = {'ATAC-SPREB-1':['ATAC-SPREB-1.R1.fastq.gz', 'ATAC-SPREB-1.R2.fastq.gz'], 
		'ATAC-SPREB-2':['ATAC-SPREB-2.R1.fastq.gz', 'ATAC-SPREB-2.R2.fastq.gz'], 
		'ATAC-SPREB-3':['ATAC-SPREB-3.R1.fastq.gz', 'ATAC-SPREB-3.R2.fastq.gz']}
# atac_f1_fq_dict = {'ATAC-SPREB-1':['ATAC-SPREB-1.R1.fastq.gz', 'ATAC-SPREB-1.R2.fastq.gz']}
# '''
for sample in atac_f1_fq_dict:
	##1. map reads
	'''
	##using bwa
	# map_using_bwa_mem(sample, atac_f1_fq_dict[sample], bl6_fa)
	# map_using_bwa_mem(sample, atac_f1_fq_dict[sample], spret_name)
	'''
	##using bowtie2
	# map_using_bowtie2(sample, atac_f1_fq_dict[sample], bl6_name)
	# map_using_bowtie2(sample, atac_f1_fq_dict[sample], spret_name)
	##2. convert to mm10 coordinates 
	# align_bams_to_mm10_lapels(sample, bl6_name)
	# align_bams_to_mm10_lapels(sample, spret_name)
	##3. fix flags of orphaned reads with no mate
	# run_samtools_fixmate_get_paired(sample, bl6_name)
	# run_samtools_fixmate_get_paired(sample, spret_name)
	##4. get allele specific reads
	# get_allele_specific_reads(sample, bl6_name, spret_passed_aaseq)
	# get_allele_specific_reads(sample, spret_name, spret_passed_aaseq)
	##? split into paired and unpaired (just keeping paired atm) -- ? needed
	##5. find_intersecting_snps ##only works on spret
	# find_intersecting_snps_paired_end(sample, spret_name, '_mm10lap_hap2', spret_snp_wasp_dir)
	# find_intersecting_snps_paired_end(sample, bl6_name, '_mm10lap_hap1', spret_snp_wasp_dir)
	##6. remap reads that overlap a snp to respectuve genome
	remap_fq_paired(sample, spret_name, '_mm10lap_hap2', spret_name)
	remap_fq_paired(sample, bl6_name, '_mm10lap_hap1', bl6_name)
	##7. filter_remapped_reads --- sort out names etc
	# filter_remapped_reads_sort_f1(sample, spret_name, '_mm10lap_hap2', mm10_fa)
	# filter_remapped_reads_sort_f1(sample, bl6_name, '_mm10lap_hap1', mm10_fa)
	##8. remap to mm10
	##8. merge all reads
	##9. rmdups

# '''
