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
analysis of RNASeq for spreb and pwb data

load modules:
module load local_python/2.7.14
module load R/3.3.3



must have environment set up i.e.
add local_python/2.7.14 and rm biobuilds
for lapels to run needed old version of pysam: conda install pysam=0.8.4
for allele specific reads need R/3.3.3 which has asSeq installed
for find_intersected_snps needed pytables version 2 so:conda install pytables=2.4.0
to install/use macs use: conda install -c auto macs 
to use homer: module load homer/4.9.1
run wasp etc -https://github.com/bmvdgeijn/WASP/tree/master/mapping
genomes etc - dl from http://www.csbio.unc.edu/CCstatus/index.py?run=Pseudo
snp vcfs from sanger ftp - ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/
to run tophat2 needed to instal bowtie2: conda install -c bioconda bowtie2
'''

##parameters
delim = '\t'
# working_dir = '/home/atimms/ngs_data/misc/cherry_f1_hybrid_0618'
working_dir = '/home/atimms/ngs_data/misc/cherry_f1_hybrid/spreb_pwb_rnaseq_0918'
os.chdir(working_dir)
threads = '15'


##genomes - dl from http://www.csbio.unc.edu/CCstatus/index.py?run=Pseudo
# star_ref_dir = '/home/atimms/ngs_data/references/star/'
psedogenome_dir = '/home/atimms/ngs_data/references/pseudogenomes/'
bl6_fa = psedogenome_dir + 'C57BL6J_b38.fa'
spret_fa = psedogenome_dir + 'SPRETEiJ_b38_f.fa'
pwk_fa = psedogenome_dir + 'PWKPhJ_b38_f.fa'
bl6_name = bl6_fa.split('/')[-1].split('.')[0]
spret_name = spret_fa.split('/')[-1].split('.')[0]
pwk_name = pwk_fa.split('/')[-1].split('.')[0]

# bl6_index_dir = star_ref_dir + bl6_name
# spret_index_dir = star_ref_dir + spret_name
bl6_hisat_index_prefix = psedogenome_dir + bl6_name
spret_hisat_index_prefix = psedogenome_dir + spret_name
pwk_hisat_index_prefix = psedogenome_dir + pwk_name
# GRCm38_fa = '/home/atimms/ngs_data/references/mm10/GRCm38_68.fa'
# GRCm38_name = 'GRCm38'
# GRCm38_name_and_dir = '/home/atimms/ngs_data/references/mm10/GRCm38_68'
# GRCm38_info = '/home/atimms/ngs_data/references/mm10/GRCm38_68.chromInfo.txt.gz'
# spret_snp_vcf = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz'
# spret_passed_snp_vcf = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.passed.vcf.gz'
# spret_passed_bed = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.passed.bed'
# spret_passed_aaseq = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.passed.txt'
spret_passed_aaseq = psedogenome_dir + 'SPRET_EiJ.mgp.v5.snps.dbSNP142.passed.txt'
# spret_snp_dir = 'spret_vcfs'
spret_snp_wasp_dir = psedogenome_dir + 'spret_wasp_snps'
pwk_snp_vcf = psedogenome_dir + 'PWK_PhJ.mgp.v5.snps.dbSNP142.vcf.gz'
pwk_passed_snp_vcf = psedogenome_dir + 'PWK_PhJ.mgp.v5.snps.dbSNP142.passed.vcf.gz'
pwk_passed_bed = psedogenome_dir + 'PWK_PhJ.mgp.v5.snps.dbSNP142.passed.bed'
pwk_passed_aaseq = psedogenome_dir +'PWK_PhJ.mgp.v5.snps.dbSNP142.passed.txt'
pwk_snp_dir = psedogenome_dir + 'pwk_vcfs'
pwk_snp_wasp_dir = psedogenome_dir + 'pwk_wasp_snps'

mm10_fa = '/home/atimms/ngs_data/references/mm10/mm10.fa'
mm10_name = 'mm10'
mm10_name_and_dir = '/home/atimms/ngs_data/references/mm10/mm10'

##ref files etc
fa_file_path = ['/home/atimms/ngs_data/references/igenomes/', '/genome.fa']
gtf_file_path = ['/home/atimms/ngs_data/references/igenomes/', '/genes.gtf']
ens_hg19_gtf_file = '/home/atimms/ngs_data/references/ensembl/hg19/Homo_sapiens.GRCh37.75.gtf'
# deseq_r_template = ''
# star_bam_suffix = '.Aligned.out.bam'
# sorted_bam_suffix = '.sorted.bam'
genome = 'mm10'
gtf_file = gtf_file_path[0] + genome + gtf_file_path[1]

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
feature_counts = '/home/atimms/programs/subread-1.6.0-Linux-x86_64/bin/featureCounts'

##methods

def make_star_index_files_no_gtf(star_genome_dir, genome_fas, star_threads):
	star_index = subprocess.Popen([star, '--runMode', 'genomeGenerate', '--genomeDir', star_genome_dir, '--genomeFastaFiles', genome_fas, '--runThreadN', star_threads])
	star_index.wait()

def make_hisat_index_files(genome_fasta, index_prefix):
	hisat_index = subprocess.Popen([hisat2_build, genome_fasta, index_prefix])
	hisat_index.wait()

def star_align_paired_end(sample_name, fq_files, star_genome_dir, fa_name):
	r1_fq = fq_files[0]
	r2_fq = fq_files[1]
	print sample_name, r1_fq, r2_fq, star_genome_dir
	star_bam = sample_name + '.' + fa_name + '.Aligned.out.bam'
	sorted_bam = sample_name + '.' + fa_name + '.star.bam'
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.' + fa_name + '.', '--runThreadN', threads])
	star_align.wait()
	##sort bam
	# st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', sample_name, '-o', sorted_bam, '-@', '10', '-m', '10G', star_bam])
	# st_sort_pe.wait()
	# st_index = subprocess.Popen([samtools, 'index', sorted_bam])
	# st_index.wait()

def hisat2_paired_end(sample_name, fq_files, hisat_index_prefix, out_bam):
	r1_fq = fq_files[0]
	r2_fq = fq_files[1]
	print sample_name, r1_fq, r2_fq, hisat_index_prefix, out_bam
	temp_bam = sample_name + '.' + hisat_index_prefix.split('/')[-1] + '.temp.bam'
	hisat2_align = subprocess.Popen([hisat2, '-x', hisat_index_prefix, '-1', r1_fq, '-2', r2_fq, '-p', threads], stdout=subprocess.PIPE)
	st_sam_bam = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-q', '10', '-o', temp_bam, '-'], stdin=hisat2_align.stdout)
	# star_align = subprocess.Popen([hisat2, '-x', hisat_index_prefix, '-1', r1_fq, '-2', r2_fq])
	st_sam_bam.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', sample_name, '-o', out_bam, '-@', '10', '-m', '10G', temp_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', out_bam])
	st_index.wait()

def tophat2_paired_end(sample_name, fq_files, bwt2_index_prefix, out_bam):
	r1_fq = fq_files[0]
	r2_fq = fq_files[1]
	print sample_name, r1_fq, r2_fq, bwt2_index_prefix, out_bam
	temp_dir = sample_name + '.' + bwt2_index_prefix.split('/')[-1] + '.temp.dir'
	temp_bam = temp_dir + '/accepted_hits.bam'
	temp2_bam = temp_dir + '/temp.bam'
	tophat2_align = subprocess.Popen([tophat2,'-p', threads, '--no-mixed', '--no-discordant', '-o', temp_dir, bwt2_index_prefix, r1_fq, r2_fq])
	tophat2_align.wait()
	st_sam_bam = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-q', '10', '-o', temp2_bam, temp_bam])
	st_sam_bam.wait()
	##sort bam
	st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', sample_name, '-o', out_bam, '-@', '10', '-m', '10G', temp2_bam])
	st_sort_pe.wait()
	st_index = subprocess.Popen([samtools, 'index', out_bam])
	st_index.wait()

def star_align_single_end(sample_name, fq_file, star_genome_dir):
	print sample_name, fq_file
	star_align = subprocess.Popen([star,'--genomeDir', star_genome_dir, '--readFilesIn', fq_file, '--outSAMtype', 'BAM', 'Unsorted', '--readFilesCommand', 'zcat', '--outFileNamePrefix', sample_name + '.' + fa_name + '.', '--runThreadN', threads])
	star_align.wait()


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

def convert_vcf_to_aaseq_snps(in_vcf_gz, bed, aaseq_txt):
	## vcf to bed  "/opt/bedops-2.4.14/bin/vcf2bed < ${describer}.vcf > ${describer}_pre.bed"
	##vcf to bed
	bgzip_cmd = subprocess.Popen(['bgzip', '-d', in_vcf_gz])
	bgzip_cmd.wait()
	in_vcf = in_vcf_gz.rsplit('.', 1)[0]
	# print in_vcf
	# print([vcf2bed, '<', in_vcf])
	# '''
	with open(bed, 'w') as bed_fh:
		cat_vcf = subprocess.Popen(['cat', in_vcf], stdout=subprocess.PIPE)
		vcf2bed_bo = subprocess.Popen([vcf2bed, '-'], stdin=cat_vcf.stdout, stdout=bed_fh)
		vcf2bed_bo.wait()
	# '''
	##bed to asseq file (with added chr?)
	with open(bed, 'r') as bed_fh, open(aaseq_txt, 'w') as txt_fh:
		for line in bed_fh:
			line = line.split(delim)
			line_out = ['chr' + line[0], line[2], line[5], line[6], '\n']
			txt_fh.write(delim.join(line_out))


def feature_count_paired_end(bam_suffix, genome_gtf, outfile, samples):
	bam_files = []
	for s in samples:
		bam = s + bam_suffix
		bam_files.append(bam)
	print bam_files, len(bam_files)
	feat_count = subprocess.Popen([feature_counts,'-p', '-t', 'exon', '-g', 'gene_id', '-a', genome_gtf, '-o', outfile, '-T', '10', '-B'] + bam_files)
	feat_count.wait()

def feature_count_paired_end_bam_input(bam_files, genome_gtf, outfile):
	print bam_files, len(bam_files)
	feat_count = subprocess.Popen([feature_counts,'-p', '-t', 'exon', '-g', 'gene_id', '-a', genome_gtf, '-o', outfile, '-T', '10', '-B'] + bam_files)
	feat_count.wait()

def format_feature_counts_file(infile, outfile):
	with open(outfile, "w") as outph, open(infile, "r") as inph:
		line_count = 0
		for line in inph:
			if line[0] != '#':
				line_count += 1
				line = line.strip('\n').split(delim)
				if line_count == 1:
					print line
					samples = line[6:]
					# print samples
					samples = [s.split('.')[0] + '_' +  s.split('.')[1].split('_')[0] for s in samples]
					# print samples
					# samples = [sample_dict_cond[s][0] for s in samples]
					print samples
					header = ['gene'] + samples + ['\n']
					outph.write(delim.join(header))
				else:
					gene = line[0]
					if gene != '':
						lineout = [gene] + line[6:] + ['\n']
						# print lineout
						outph.write(delim.join(lineout))


##make star index files for genomes

##star
# make_star_index_files_no_gtf(bl6_index_dir, bl6_fa, threads)
# make_star_index_files_no_gtf(spret_index_dir, spret_fa, threads)
##hisat
# make_hisat_index_files(bl6_fa, bl6_hisat_index_prefix)
# make_hisat_index_files(spret_fa, spret_hisat_index_prefix)
# make_hisat_index_files(pwk_fa, pwk_hisat_index_prefix)

##make SNP files for PWK
##Convert SNP data to HDF5 format 
# get_passed_snps(pwk_snp_vcf, pwk_passed_snp_vcf)
# split_vcf_by_chr(pwk_passed_snp_vcf, pwk_name, pwk_snp_dir)
#make snp files for wasp (don't use hd5)
# make_txt_files_from_snps(pwk_snp_dir, pwk_snp_wasp_dir)
##get snps for allele specific reads
# convert_vcf_to_aaseq_snps(pwk_passed_snp_vcf, pwk_passed_bed, pwk_passed_aaseq)




##for spreb rnaseq
spreb_rnaseq_fq_dict = {'SPREB1':['SPREB1.R1.fq.gz', 'SPREB1.R2.fq.gz'], 
		'SPREB2':['SPREB2.R1.fq.gz', 'SPREB2.R2.fq.gz'], 'SPREB3':['SPREB3.R1.fq.gz', 'SPREB3.R2.fq.gz']}
##testing
# spreb_rnaseq_fq_dict = {'SPREB3':['SPREB3.R1.fq.gz', 'SPREB3.R2.fq.gz']}

for sample in spreb_rnaseq_fq_dict:
	##file and dir names
	##1 bams
	bam_bl6 = sample + '.' + bl6_name + '.tophat.bam'
	bam_spret = sample + '.' + spret_name + '.tophat.bam'
	##2 bams
	bam_bl6_mm10 = sample + '.' + bl6_name + '.mm10l.tophat.bam'
	bam_spret_mm10 = sample + '.' + spret_name + '.mm10l.tophat.bam'
	##3 bams
	bam_bl6_mm10_fixed = sample + '.' + bl6_name + '.mm10l.tophat.fixed.bam'
	bam_spret_mm10_fixed = sample + '.' + spret_name + '.mm10l.tophat.fixed.bam'
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
	bl6_remap_bam = sample + '.' + bl6_name + '.tophat_remap.bam'
	spret_remap_bam = sample + '.' + spret_name + '.tophat_remap.bam'
	##7 bams converted back to mm10
	bl6_remap_bam_mm10 = sample + '.' + bl6_name + '.tophat_remap.mm10l.bam'
	spret_remap_bam_mm10 = sample + '.' + spret_name + '.tophat_remap.mm10l.bam'
	##8 bam with original reads to map and the the out bams i.e ones to keep
	bl6_to_remap_bam = bl6_intsnp_dir + sample + '.' + bl6_name + '.to.remap.bam'
	spret_to_remap_bam = spret_intsnp_dir + sample + '.' + spret_name + '.to.remap.bam'
	bl6_final_bam = sample + '.' + bl6_name + '.allele_wasp.bam'
	spret_final_bam = sample + '.' + spret_name + '.allele_wasp.bam'



	##1. map reads to both genomes
	'''
	# hisat2_paired_end(sample, spreb_rnaseq_fq_dict[sample], bl6_hisat_index_prefix, bam_bl6)
	# hisat2_paired_end(sample, spreb_rnaseq_fq_dict[sample], spret_hisat_index_prefix, bam_spret)
	tophat2_paired_end(sample, spreb_rnaseq_fq_dict[sample], bl6_hisat_index_prefix, bam_bl6)
	tophat2_paired_end(sample, spreb_rnaseq_fq_dict[sample], spret_hisat_index_prefix, bam_spret)

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
	# hisat2_paired_end(sample, bl6_remap_fqs, bl6_hisat_index_prefix, bl6_remap_bam)
	# hisat2_paired_end(sample, spret_remap_fqs, spret_hisat_index_prefix, spret_remap_bam)

	tophat2_paired_end(sample, bl6_remap_fqs, bl6_hisat_index_prefix, bl6_remap_bam)
	tophat2_paired_end(sample, spret_remap_fqs, spret_hisat_index_prefix, spret_remap_bam)

	##7. convert remaped bam to mm10
	align_bams_to_mm10_lapels(bl6_remap_bam, bl6_remap_bam_mm10, bl6_name)
	align_bams_to_mm10_lapels(spret_remap_bam, spret_remap_bam_mm10, spret_name)
	
	##8. filter_remapped_reads ##no need to merge as no 'kept reads'
	filter_remapped_reads_sort(sample, bl6_name, bl6_fa, bl6_to_remap_bam, bl6_remap_bam_mm10, bl6_final_bam)
	filter_remapped_reads_sort(sample, spret_name, spret_fa, spret_to_remap_bam, spret_remap_bam_mm10, spret_final_bam)
	'''

##rnaseq analysis
##using final post wasp files
project_name = 'spreb_rnaseq_wasp_0918'
# sample_list = ['SPREB1.C57BL6J_b38', 'SPREB1.SPRETEiJ_b38_f', 'SPREB2.C57BL6J_b38', 'SPREB2.SPRETEiJ_b38_f', 
# 		'SPREB3.C57BL6J_b38', 'SPREB3.SPRETEiJ_b38_f' ]
sample_list = ['SPREB1.C57BL6J_b38', 'SPREB1.SPRETEiJ_b38_f', 'SPREB2.C57BL6J_b38', 'SPREB2.SPRETEiJ_b38_f', 'SPREB3.C57BL6J_b38', 'SPREB3.SPRETEiJ_b38_f']
post_wasp_bam_suffix = '.allele_wasp.bam'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'

##get count data then format
# feature_count_paired_end(post_wasp_bam_suffix, gtf_file, feature_count_results_file, sample_list)
# format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)


##using the pre wasp bams
project_name = 'spreb_rnaseq_split_0918'
pre_wasp_bams = ['SPREB1.C57BL6J_b38_hap1.bam', 'SPREB1.SPRETEiJ_b38_f_hap2.bam', 'SPREB2.C57BL6J_b38_hap1.bam', 'SPREB2.SPRETEiJ_b38_f_hap2.bam', 'SPREB3.C57BL6J_b38_hap1.bam', 'SPREB3.SPRETEiJ_b38_f_hap2.bam'] 
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'

##get count data then format
# feature_count_paired_end_bam_input(pre_wasp_bams, gtf_file, feature_count_results_file)
# format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)

##for pwk rnaseq
pwk_rnaseq_fq_dict = {'PWB1':['PWB1.R1.fq.gz', 'PWB1.R2.fq.gz'], 
		'PWB2':['PWB2.R1.fq.gz', 'PWB2.R2.fq.gz'], 'PWB3':['PWB3.R1.fq.gz', 'PWB3.R2.fq.gz']}
# pwk_rnaseq_fq_dict = {'PWB3':['PWB3.R1.fq.gz', 'PWB3.R2.fq.gz']}
##testing
# spreb_rnaseq_fq_dict = {'test':['test.R1.fq.gz', 'test.R2.fq.gz']}


for sample in pwk_rnaseq_fq_dict:
	##file and dir names
	##1 bams
	bam_bl6 = sample + '.' + bl6_name + '.tophat.bam'
	bam_pwk = sample + '.' + pwk_name + '.tophat.bam'
	##2 bams
	bam_bl6_mm10 = sample + '.' + bl6_name + '.mm10l.tophat.bam'
	bam_pwk_mm10 = sample + '.' + pwk_name + '.mm10l.tophat.bam'
	##3 bams
	bam_bl6_mm10_fixed = sample + '.' + bl6_name + '.mm10l.tophat.fixed.bam'
	bam_pwk_mm10_fixed = sample + '.' + pwk_name + '.mm10l.tophat.fixed.bam'
	##4 out_prefix
	bl6_allele = sample + '.' + bl6_name
	pwk_allele = sample + '.' + pwk_name
	##5 in bams and out_dir
	bl6_hap1_bam = bl6_allele + '_hap1.bam'
	pwk_hap2_bam = pwk_allele + '_hap2.bam'
	bl6_intsnp_dir = sample + '.' + bl6_name + '_hap1.int_snps/'
	pwk_intsnp_dir = sample + '.' + pwk_name + '_hap2.int_snps/'
	##6 in fqs and out bam
	bl6_remap_fqs = [bl6_intsnp_dir + sample + '.' + bl6_name + '.remap.fq1.gz', bl6_intsnp_dir + sample + '.' + bl6_name + '.remap.fq2.gz']
	pwk_remap_fqs = [pwk_intsnp_dir + sample + '.' + pwk_name + '.remap.fq1.gz', pwk_intsnp_dir + sample + '.' + pwk_name + '.remap.fq2.gz']
	bl6_remap_bam = sample + '.' + bl6_name + '.hisat2_remap.bam'
	pwk_remap_bam = sample + '.' + pwk_name + '.hisat2_remap.bam'
	##7 bams converted back to mm10
	bl6_remap_bam_mm10 = sample + '.' + bl6_name + '.hisat2_remap.mm10l.bam'
	pwk_remap_bam_mm10 = sample + '.' + pwk_name + '.hisat2_remap.mm10l.bam'
	##8 bam with original reads to map and the the out bams i.e ones to keep
	bl6_to_remap_bam = bl6_intsnp_dir + sample + '.' + bl6_name + '.to.remap.bam'
	pwk_to_remap_bam = pwk_intsnp_dir + sample + '.' + pwk_name + '.to.remap.bam'
	bl6_final_bam = sample + '.' + bl6_name + '.allele_wasp.bam'
	pwk_final_bam = sample + '.' + pwk_name + '.allele_wasp.bam'
	'''
	##1. map reads to both genomes
	tophat2_paired_end(sample, pwk_rnaseq_fq_dict[sample], bl6_hisat_index_prefix, bam_bl6)
	tophat2_paired_end(sample, pwk_rnaseq_fq_dict[sample], pwk_hisat_index_prefix, bam_pwk)

	##2. convert to mm10 coordinates 
	align_bams_to_mm10_lapels(bam_bl6, bam_bl6_mm10, bl6_name)
	align_bams_to_mm10_lapels(bam_pwk, bam_pwk_mm10, pwk_name)

	##3. fix flags of orphaned reads with no mate
	run_samtools_fixmate_get_paired(sample, bl6_name, bam_bl6_mm10, bam_bl6_mm10_fixed)
	run_samtools_fixmate_get_paired(sample, pwk_name, bam_pwk_mm10, bam_pwk_mm10_fixed)

	##4. get allele specific reads
	get_allele_specific_reads(bam_bl6_mm10_fixed, bl6_allele, pwk_passed_aaseq)
	get_allele_specific_reads(bam_pwk_mm10_fixed, pwk_allele, pwk_passed_aaseq)
	
	##5. find_intersecting_snps ##only works on pwk snps
	find_intersecting_snps_paired_end(sample, bl6_name, bl6_hap1_bam, bl6_intsnp_dir, pwk_snp_wasp_dir)
	find_intersecting_snps_paired_end(sample, pwk_name, pwk_hap2_bam, pwk_intsnp_dir, pwk_snp_wasp_dir)
	
	##6. remap reads that overlap a snp to respective genome
	tophat2_paired_end(sample, bl6_remap_fqs, bl6_hisat_index_prefix, bl6_remap_bam)
	tophat2_paired_end(sample, pwk_remap_fqs, pwk_hisat_index_prefix, pwk_remap_bam)

	##7. convert remaped bam to mm10
	align_bams_to_mm10_lapels(bl6_remap_bam, bl6_remap_bam_mm10, bl6_name)
	align_bams_to_mm10_lapels(pwk_remap_bam, pwk_remap_bam_mm10, pwk_name)
	
	##8. filter_remapped_reads ##no need to merge as no 'kept reads'
	filter_remapped_reads_sort(sample, bl6_name, bl6_fa, bl6_to_remap_bam, bl6_remap_bam_mm10, bl6_final_bam)
	filter_remapped_reads_sort(sample, pwk_name, pwk_fa, pwk_to_remap_bam, pwk_remap_bam_mm10, pwk_final_bam)
	'''


##rnaseq analysis
##using final post wasp files
project_name = 'pwb_rnaseq_wasp_1018'
# sample_list = ['SPREB1.C57BL6J_b38', 'SPREB1.SPRETEiJ_b38_f', 'SPREB2.C57BL6J_b38', 'SPREB2.SPRETEiJ_b38_f', 
# 		'SPREB3.C57BL6J_b38', 'SPREB3.SPRETEiJ_b38_f' ]
sample_list = ['PWB1.C57BL6J_b38', 'PWB1.PWKPhJ_b38_f', 'PWB2.C57BL6J_b38', 'PWB2.PWKPhJ_b38_f', 'PWB3.C57BL6J_b38', 'PWB3.PWKPhJ_b38_f']
post_wasp_bam_suffix = '.allele_wasp.bam'
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'

##get count data then format
feature_count_paired_end(post_wasp_bam_suffix, gtf_file, feature_count_results_file, sample_list)
format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)


##using the pre wasp bams
project_name = 'pwb_rnaseq_split_1018'
pre_wasp_bams = ['PWB1.C57BL6J_b38_hap1.bam', 'PWB1.PWKPhJ_b38_f_hap2.bam', 'PWB2.C57BL6J_b38_hap1.bam', 'PWB2.PWKPhJ_b38_f_hap2.bam', 'PWB3.C57BL6J_b38_hap1.bam', 'PWB3.PWKPhJ_b38_f_hap2.bam'] 
feature_count_results_file = project_name + '.feature_counts.all_samples.txt'
deseq_feature_count_file = project_name + '.star_fc.counts.txt'

##get count data then format
feature_count_paired_end_bam_input(pre_wasp_bams, gtf_file, feature_count_results_file)
format_feature_counts_file(feature_count_results_file, deseq_feature_count_file)


