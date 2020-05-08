#!/usr/bin/env python
import sys
import subprocess
import os
import glob


##note
'''
analysis of RNASeq for spreb and pwb data

load modules:
module load biobuilds/2017.11
module load mono/5.10.1.47
module load Pisces/5.1.6.54
module load homer/4.9.1
module load local_python/3.6.4
'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/misc/cherry_atac_analysis_1018'
os.chdir(working_dir)
##files and dirs
fasta = '/home/atimms/ngs_data/references/hg38/hg38.fa'
bwa_bam_suffix = '.bwa_mkdup.bam'
bamlist = 'bams.list'
peak_bed = 'Hu-ret-ATAC_pooled_peaks_conservative_postidr.bed'
pisces_fastq_dir = '/home/atimms/ngs_data/references/illumina_hg38/'

##programs
gatk = '/cm/shared/apps/GATK/3.7/GenomeAnalysisTK.jar'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
pisces = '/cm/shared/apps/Pisces/5.1.6.54/Pisces.exe'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'

##annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,kaviar_20150923,gnomad_genome,avsnp150,rmsk,genomicSuperDups,bed,bed', '-bedfile', 'homer_motifs_hg38_1018.bed,imprint_bed_v2.bed']
av_operation = ['-operation', 'g,f,f,f,f,f,r,r']
# av_options = ['-otherinfo', '-remove', '-nastring', '.']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']


def align_with_bwa(samples, r1_suffix, r2_suffix, final_bam_suffix):
	for sample in samples:
		r1_fq = sample + r1_suffix
		r2_fq = sample + r2_suffix
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
		post_bwa_bam = sample + '.bwa.bam'
		sort_bam = sample + '.bwa_sort.bam'
		mkdup_bam = sample + final_bam_suffix
		##bwa and convert to bam
		bwa_pe = subprocess.Popen(['bwa', 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen(['samtools', 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()
		##sort bam
		st_sort_pe = subprocess.Popen(['samtools', 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam, '-@', '10', '-m', '10G', post_bwa_bam])
		st_sort_pe.wait()
		##mark duplicates
		# picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
		picard_md = subprocess.Popen(['picard', 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
		picard_md.wait()

##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')

def variant_calling_gatk_hc_no_gvcf(bamlist, name_prefix):
	vcf_temp0 = 'temp_gatk0.vcf'
	vcf_raw_snps = 'temp_raw_snps.vcf'
	vcf_filtered_snps = 'temp_filtered_snps.vcf'
	vcf_raw_indels = 'temp_raw_indels.vcf'
	vcf_filtered_indels = 'temp_filtered_indels.vcf'
	vcf_temp1 = 'temp_gatk1.vcf'
	vcf_temp2 = 'temp_gatk2.vcf.gz'
	final_vcf = name_prefix + '.gatkHC.vcf.gz'
	bam_files = []
	##run haplotype caller
	gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', bamlist,'-L', peak_bed, '-o', vcf_temp0])
	gatk_hc.wait()
	##split data in SNPs and indels and apply manual variant filtering (add -L)
	snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_snps, '-selectType', 'SNP'])
	snp_cut.wait()
	snp_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_snps, '-o', vcf_filtered_snps, '--filterExpression ', "QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "snp_filter"])
	snp_vf.wait()
	indel_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_indels, '-selectType', 'INDEL'])
	indel_cut.wait()
	indel_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_indels, '-o', vcf_filtered_indels, '--filterExpression ', "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0", '--filterName', "indel_filter"])
	indel_vf.wait()
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '5', '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', vcf_temp1, '-genotypeMergeOptions', 'UNSORTED'])
	combine_var.wait()
	bgzip_cmd = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip_cmd.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()

def variant_calling_freebayes(bamlist, name_prefix):
	vcf_temp1 = 'temp_fb1.vcf'
	vcf_temp2 = 'temp_fb2.vcf.gz'
	final_vcf = name_prefix + '.freebayes.vcf.gz'
	vcf_fh = open(vcf_temp1, 'w')
	freebayes_run = subprocess.Popen([freebayes,'-f' , fasta, '-L', bamlist, '-t', peak_bed], stdout=vcf_fh)
	freebayes_run.wait()
	vcf_fh.close()
	bgzip_run = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip_run.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()

def variant_calling_pisces(samples, bam_suffix, name_prefix):
	out_dir = name_prefix + '_pisces'
	all_bams = [i + bam_suffix for i in samples]
	print(all_bams)
	##check bai file exists and make if it isn't there
	for bam in all_bams:
		if os.path.isfile(bam + '.bai'):
			print('bam %s alreaded indexed'%bam)
		else:
			print('indexing bam file:', bam)
			st_index = subprocess.Popen(['samtools', 'index', bam])
			st_index.wait()
	##run pisces on all samples
	# '''
	run_pisces = subprocess.Popen(['mono', pisces, '-G', pisces_fastq_dir, '-BamPaths', ','.join(all_bams), '-MinVF', '0.01', '-i', peak_bed, '-t', '16', '-OutFolder', out_dir])
	run_pisces.wait()

def bedtools_intersect_bed_vcf(bed, vcf, out_file):
	with open(out_file, "w") as out_fh: 
		# hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', bed, '-b', vcf, '-wa', '-wb'], stdout=out_fh)
		hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', bed, '-b', vcf, '-wao'], stdout=out_fh)

		hom_bt_intersect.wait()

def count_get_average_allele_freq_per_sample(in_bed, out_bed, filtering_info, var_caller):
	maf_dict = {}
	if var_caller == 'pisces':
		with open(in_bed, "r") as in_fh:
			for line in in_fh:
				line = line.rstrip().split(delim)
				region = '_'.join(line[:3])
				info = line[13]
				if info == '.':
					maf_dict[region] = ['-']
				else:
					alleles = info.split(':')[2].split(',')
					filt = line[10]
					gt = info.split(':')[0]
					if len(alleles) == 2 and gt == '0/1': 
						# print(region, info, filt)
						ref = int(alleles[0])
						alt = float(alleles[1])
						total = ref + alt
						if ref > alt:
							maf =  alt / total
						else:
							maf =  ref / total
						# print(alleles, maf)
						if total >= filtering_info[0] and filt == filtering_info[1]:
							print(region, maf, line)
							if region in maf_dict:
								if maf_dict[region][0] == '-':
									maf_dict[region] = [maf]
								else:
									maf_dict[region].append(maf)
							else:
								maf_dict[region] = [maf]
						else:
							if region not in maf_dict:
								maf_dict[region] = ['-']							
					else:
						if region not in maf_dict:
							maf_dict[region] = ['-']
						print('issue with this line:', line)
	else:
		print('var caller %s not recognized'%variant_caller)
	with open(out_bed, "w") as out_fh:
		for r in maf_dict:
			if maf_dict[r][0] == '-':
				average_maf, var_count = '-', '-'
			else:
				average_maf = str( round(sum(maf_dict[r]) /len(maf_dict[r]),2))
				var_count = str(len(maf_dict[r]))
			pos = r.split('_')
			out_fh.write(delim.join(pos + [average_maf, var_count]) + '\n')

def combine_ind_maf_files(samples, var_caller):
	maf_dict = {}
	if var_caller == 'pisces':
		in_suffix = '.pisces.sum.bed'
		out_file = 'maf_allele_counts.pisces.summary.xls'
		for sample in samples:
			in_bed = sample + in_suffix
			with open(in_bed, "r") as in_fh:
				for line in in_fh:
					line = line.rstrip().split(delim)
					region = '_'.join(line[:3])
					info = line[3:]
					if region in maf_dict:
						maf_dict[region].extend(info)
					else:
						maf_dict[region] = info
	header = ['chr', 'start', 'end']
	for s in samples:
		header.extend([s + ' maf', s + ' var count'])
	header.extend(['average maf', 'total var count'])

	with open(out_file, "w") as out_fh:
		out_fh.write(delim.join(header) + '\n')
		for reg in maf_dict:
			r = reg.split('_')
			counts = maf_dict[reg]
			counts_tidy = list(filter(('-').__ne__, counts))
			if len(counts_tidy) > 0:
				print(counts_tidy)
				mafs = [float(i) for i in counts_tidy[0::2]]
				print(counts_tidy[0::2], mafs)
				a_maf = str(round(sum(mafs)/len(mafs),2))
				var_counts =[int(i) for i in counts_tidy[1::2]]
				t_vc = str(sum(var_counts))
				print(counts_tidy[1::2], var_counts)
			else:
				a_maf = '0'
				t_vc = '0'
			out_fh.write(delim.join(r + counts + [a_maf, t_vc] ) + '\n')


def get_average_maf_per_window(samples, peaks_bed, variant_caller, filter_params):
	for sample in samples:
		if variant_caller == 'pisces':
			snp_vcf = 'ATAC_AI_1018_pisces/' + sample + '.bwa_mkdup.vcf'
			intersected_bed = sample + '.pisces.bt_int.bed'
			counts_bed = sample + '.pisces.sum.bed'
			# snp_bed = sample + '.pisces_snps.bed'
		else:
			print('var caller %s not recognized'%variant_caller)
		##intersect bed file with each vcf
		# bedtools_intersect_bed_vcf(peaks_bed, snp_vcf, intersected_bed)
		# count_get_average_allele_freq_per_sample(intersected_bed, counts_bed, filter_params, variant_caller)
	combine_ind_maf_files(samples, variant_caller)

def homer_find_motifs_hg38(motif_file, out_bed):
	#scanMotifGenomeWide.pl pu1.motif mm9 -bed > pu1.sites.mm9.bed
	with open(out_bed, "w") as out_fh: 
		hom_bt_intersect = subprocess.Popen(['scanMotifGenomeWide.pl', motif_file, 'hg38', '-bed'], stdout=out_fh)
		hom_bt_intersect.wait()


def intersect_vcfs_with_motif_bed(samples, variant_caller, motifs_bed):
	for sample in samples:
		if variant_caller == 'pisces':
			snp_vcf = 'ATAC_AI_1018_pisces/' + sample + '.bwa_mkdup.vcf'
			out_bed = sample + '.' + variant_caller + '.' + motifs_bed.split('.')[0] +  '.bed'

			with open(out_bed, "w") as out_fh: 
				hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', snp_vcf, '-b', motifs_bed, '-wa'], stdout=out_fh)
				hom_bt_intersect.wait()


def annotate_vcf(samples, variant_caller):
	for sample in samples:
		av_prefix = sample + '.' + variant_caller
		if variant_caller == 'pisces':
			snp_vcf = 'ATAC_AI_1018_pisces/' + sample + '.bwa_mkdup.vcf'

		##run_table_annovar(avinput, av_prefix):
		command = [table_annovar] + av_buildver + [snp_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()



def filter_vcf_2_ways(in_vcf, cq_vcf, cqd_vcf):
	print(in_vcf, cq_vcf, cqd_vcf)
	with open(in_vcf, 'r') as in_fh, open(cq_vcf, 'w') as cq_fh, open(cqd_vcf, 'w') as cqd_fh:
		for line in in_fh:
			if line[0] == '#':
				cq_fh.write(line)
				cqd_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				filt = line[6]
				info = line[9]
				ann_info = line[7]
				alleles = info.split(':')[2].split(',')
				coverage = sum([int(a) for a in alleles])
				# print(filt, alleles, coverage)
				if filt == 'PASS' and coverage >= 20:
					cq_fh.write(delim.join(line) + '\n')
					if 'rmsk=.;genomicSuperDups=.' in ann_info and 'avsnp150=rs' in ann_info:
						cqd_fh.write(delim.join(line) + '\n')


def get_maf_per_sample_take2(in_bed, out_bed, var_caller):
	maf_dict = {}
	if var_caller == 'pisces':
		with open(in_bed, "r") as in_fh:
			for line in in_fh:
				line = line.rstrip().split(delim)
				region = '_'.join(line[:3])
				info = line[13]
				ann_info = line[11]
				in_motif = 'no'
				in_imprinted = 'no'
				# print(info)
				if 'bed=N' in ann_info:
					in_motif = 'yes'
				if 'bed2=N' in ann_info:
					in_imprinted = 'yes'
				if info == '.':
					maf_dict[region] = [['-'],[],[]]
				else:
					alleles = info.split(':')[2].split(',')
					gt = info.split(':')[0]
					if len(alleles) == 2 and gt == '0/1': 
						# print(region, info, filt)
						ref = int(alleles[0])
						alt = float(alleles[1])
						total = ref + alt
						if ref > alt:
							maf =  alt / total
						else:
							maf =  ref / total

						# print(alleles, maf)

						if region in maf_dict:
							if maf_dict[region][0][0] == '-':
								maf_dict[region][0] = [maf]
								maf_dict[region][1] = [in_motif]
								maf_dict[region][2] = [in_imprinted]
							else:
								maf_dict[region][0].append(maf)
								maf_dict[region][1].append(in_motif)
								maf_dict[region][2].append(in_imprinted)
						else:
							maf_dict[region] = [[maf],[in_motif],[in_imprinted]]
						
					else:
						if region not in maf_dict:
							maf_dict[region] = [['-'],[],[]]
	# for i in maf_dict:
	# 	print(i, maf_dict[i])

	with open(out_bed, "w") as out_fh:
		for r in maf_dict:
			if maf_dict[r][0][0] == '-':
				average_maf, var_count, motif_var, imprinted = '-', '-', '-', '-'
			else:
				average_maf = str( round(sum(maf_dict[r][0]) /len(maf_dict[r][0]),2))
				var_count = str(len(maf_dict[r][0]))
				if 'yes' in maf_dict[r][1]:
					motif_var = 'yes'
				else:
					motif_var = 'no'
				if 'yes' in maf_dict[r][2]:
					imprinted = 'yes'
					if 'no'in maf_dict[r][2]:
						print('we see yes and no in imprinted', maf_dict[r])
				else:
					imprinted = 'no'				
			pos = r.split('_')
			out_fh.write(delim.join(pos + [average_maf, var_count, motif_var, imprinted]) + '\n')

def combine_ind_maf_files_take2(samples, in_suffix, out_file):
	maf_dict = {}
	for sample in samples:
		in_bed = sample + in_suffix
		with open(in_bed, "r") as in_fh:
			for line in in_fh:
				line = line.rstrip().split(delim)
				region = '_'.join(line[:3])
				info = line[3:]
				if region in maf_dict:
					maf_dict[region].extend(info)
				else:
					maf_dict[region] = info
	header = ['chr', 'start', 'end']
	for s in samples:
		header.extend([s + ' maf', s + ' var count', s + ' var in motif', 'imprinted gene'])
	header.extend(['average maf', 'total var count'])

	with open(out_file, "w") as out_fh:
		out_fh.write(delim.join(header) + '\n')
		for reg in maf_dict:
			r = reg.split('_')
			counts = maf_dict[reg]
			counts_tidy = list(filter(('-').__ne__, counts))
			if len(counts_tidy) > 0:
				print(counts_tidy)
				mafs = [float(i) for i in counts_tidy[0::4]]
				print(counts_tidy[0::4], mafs)
				a_maf = str(round(sum(mafs)/len(mafs),2))
				var_counts =[int(i) for i in counts_tidy[1::4]]
				t_vc = str(sum(var_counts))
				print(counts_tidy[1::2], var_counts)

			else:
				a_maf = '0'
				t_vc = '0'
			out_fh.write(delim.join(r + counts + [a_maf, t_vc] ) + '\n')

def filter_results_file(in_file, out_file):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc == 1:
				header = line + ['average maf var in motif', 'average maf no var in motif', 'average maf difference']
				out_fh.write(delim.join(header) + '\n')
			else:
				total_vars = int(line[36])
				##must be 10 vars over all the samples
				if total_vars >= 10:
					mafs = line[3:32:4]
					in_motif = line[5:34:4]

					##must be one individual with a var in a motif and not all have a var
					if 'yes' in in_motif and 'no' in in_motif:
						yes_indices = [i for i, x in enumerate(in_motif) if x == "yes"]
						yes_mafs = [float(mafs[i]) for i in yes_indices]
						yes_ave_maf = sum(yes_mafs) / len(yes_mafs)
						# print(line)
						# print(mafs, in_motif)
						# print(yes_indices, yes_mafs, yes_ave_maf)
						no_indices = [i for i, x in enumerate(in_motif) if x == "no"]
						no_mafs = [float(mafs[i]) for i in no_indices]
						no_ave_maf = sum(no_mafs) / len(no_mafs)
						# print(no_indices, no_mafs, no_ave_maf)
						maf_diff = no_ave_maf - yes_ave_maf
						line_out = line + [str(yes_ave_maf), str(no_ave_maf), str(maf_diff)]
						out_fh.write(delim.join(line_out) + '\n')

			


def get_average_maf_per_window_with_extras(samples, peaks_bed, variant_caller, motifs_bed):
	for sample in samples:
		# var_in_motif_bed = sample + '.' + variant_caller + '.' + motifs_bed.split('.')[0] +  '.bed'
		annotated_vcf = sample + '.' +  variant_caller + '.hg38_multianno.vcf'
		counts_bed = sample + '.' +  variant_caller + '.sum.bed'
		cov_q_vcf = sample + '.' +  variant_caller + '.hg38.c20_q30.vcf'
		dbsnp_vcf = sample + '.' +  variant_caller + '.hg38.c20_q30_dbsnp.vcf'
		cov_q_int_peaks_bed = cov_q_vcf.rsplit('.', 1)[0] +  '.bt_int_peaks.bed'
		dbsnp_int_peaks_bed = dbsnp_vcf.rsplit('.', 1)[0] +  '.bt_int_peaks.bed'
		cov_q_int_counts = cov_q_vcf.rsplit('.', 1)[0] +  '.sum.bed'
		dbsnp_int_counts = dbsnp_vcf.rsplit('.', 1)[0] +  '.sum.bed'	
		cov_q_int_counts_suffix = '.' +  variant_caller + '.hg38.c20_q30.sum.bed'
		dbsnp_int_counts_suffix = '.' +  variant_caller + '.hg38.c20_q30_dbsn.sum.bed'
		cov_q_int_counts_suffix = '.' +  variant_caller + '.hg38.c20_q30.sum.bed'
		dbsnp_int_counts_suffix = '.' +  variant_caller + '.hg38.c20_q30_dbsnp.sum.bed'
		cov_q_int_counts_final = 'maf_allele_counts.' +  variant_caller + '.hg38.c20_q30.summary.xls'
		dbsnp_int_counts_final = 'maf_allele_counts.' +  variant_caller + '.hg38.c20_q30_dbsnp.summary.xls'	
		cov_q_int_counts_motif = 'maf_allele_counts.' +  variant_caller + '.hg38.c20_q30.motif_filtered.xls'
		dbsnp_int_counts_motif = 'maf_allele_counts.' +  variant_caller + '.hg38.c20_q30_dbsnp.motif_filtered.xls'	
		##filter vars in 3 ways i.e. cov/q, repeat region and in dbsnp
		# filter_vcf_2_ways(annotated_vcf, cov_q_vcf, dbsnp_vcf)
		##intersect bed file with each vcf
		# bedtools_intersect_bed_vcf(peaks_bed, cov_q_vcf, cov_q_int_peaks_bed)
		# bedtools_intersect_bed_vcf(peaks_bed, dbsnp_vcf, dbsnp_int_peaks_bed)

		# get_maf_per_sample_take2(cov_q_int_peaks_bed, cov_q_int_counts, variant_caller)
		# get_maf_per_sample_take2(dbsnp_int_peaks_bed, dbsnp_int_counts, variant_caller)
	# combine_ind_maf_files_take2(samples,cov_q_int_counts_suffix , cov_q_int_counts_final)
	# combine_ind_maf_files_take2(samples,dbsnp_int_counts_suffix , dbsnp_int_counts_final)
	##summarize results
	filter_results_file(cov_q_int_counts_final, cov_q_int_counts_motif)
	filter_results_file(dbsnp_int_counts_final, dbsnp_int_counts_motif)

##run methods

sample_names = ['Hu1_ret_ATAC', 'Hu2_ret_ATAC', 'Hu3_ret_ATAC', 'Hu4_ret_ATAC', 'Hu5_ret_ATAC', 'Hu6_ret_ATAC', 'Hu7_ret_ATAC', 'Hu8_ret_ATAC']
fq_suffixes = ['_R1.fq.gz', '_R2.fq.gz']
project_name = 'ATAC_AI_1018'
homer_db_motifs = 'homer_motifs_1018.txt'
homer_motifs_bed = 'homer_motifs_hg38_1018.bed'
var_filter_params = [20, 'PASS']

##map fqs
# align_with_bwa(sample_names, fq_suffixes[0], fq_suffixes[1], bwa_bam_suffix)

##variant calling
##make file that contains list of all bams
# make_list_of_bams(sample_names, bwa_bam_suffix, bamlist)
# variant_calling_freebayes(bamlist, project_name)
# variant_calling_gatk_hc_no_gvcf(bamlist, project_name)
# variant_calling_pisces(sample_names, bwa_bam_suffix, project_name)

##find motifs in hg38, and then compare to pisces vars
# homer_find_motifs_hg38(homer_db_motifs, homer_motifs_bed)
# intersect_vcfs_with_motif_bed(sample_names, 'pisces', homer_motifs_bed) ##not used added bed to annovar annotation

##count snps in windows etc -- no filtering beyond coverage and qscore
# get_average_maf_per_window(sample_names, peak_bed, 'pisces', var_filter_params)
##annotate with annovar then get maf and compare with vars in motifs
# annotate_vcf(sample_names, 'pisces')
get_average_maf_per_window_with_extras(sample_names, peak_bed, 'pisces', homer_motifs_bed)











