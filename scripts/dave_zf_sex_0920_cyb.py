#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import filtering_annotated
import homozygosity_mapping_cybertron
import math

##parameters
delim = '\t'
thread_number = '20'

##modules 
'''
module load biobuilds/2017.11
module load gcc/6.1.0 ##for homozygosity_mapping_cybertron

'''


##files
fasta = '/home/atimms/ngs_data/references/danRer11/danRer11.fa'
fasta_fai = '/home/atimms/ngs_data/references/danRer11/danRer11.fa.fai'
##first 2 columns of fasta.fai file
bedtools_genome_file  = '/home/atimms/ngs_data/references/danRer11/danRer11.fa.genome'

##programs
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'

##annovar parameters
av_genome = 'danRer11'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,ensGene,rmsk,generic,generic,generic,generic', '-genericdbfile', 'SH_1_A5.avinput,SH_2_B3.avinput,SH_3_m.avinput,SH_4_f.avinput']
av_operation = ['-operation', 'g,g,r,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,', '-vcfinput']
#av_options = ['-otherinfo', '-remove', '-vcfinput']


##methods

def calculate_genome_coverage(samples, bam_suffix, genome_file, prefix):
	##parameters and header
	coverage_required = [1,5,10,20,50,100,200]
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
	output_file = prefix + '.coverage_data.txt'
	print output_file
	with open(output_file, "w") as outfh:
		outfh.write(header)
		##go through bam files and get coverage histogram
		for sample in samples:
			bam = sample + bam_suffix
			print 'calculating coverge for bam file', bam
			##get temp coverage files
			coverage_histogram = sample + '.hist.temp'
			hist_fh = open(coverage_histogram, 'w')
			##just the file
			# bedtools_cov = subprocess.Popen([bedtools, 'genomecov', '-ibam', bam, '-g', genome_file], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen(['bedtools', 'genomecov', '-ibam', bam, '-g', genome_file], stdout=subprocess.PIPE)
			grep_cmd = subprocess.Popen(['grep', '^genome'], stdin=bedtools_cov.stdout, stdout=hist_fh)
			grep_cmd.wait()
			hist_fh.close()
			##calculate numbers from histogram
			with open(coverage_histogram, "r") as cov_temp:
				seq_total = 0
				for line in cov_temp:
					line = line.strip('\n').split(delim)
					target_size = float(line[3])
					if line[1] != '0':
						seq_total += int(line[1]) *  int(line[2])
			cov_stats = []
			for cr in coverage_required:
				coverage_bp = 0
				with open(coverage_histogram, "r") as cov_temp:
					for line in cov_temp:
						line = line.strip('\n').split(delim)
						target_size = float(line[3])
						if int(line[1]) >= cr:
							coverage_bp += int(line[2])
					cov_stats.extend([str(coverage_bp), str((coverage_bp / target_size) * 100)])
			line_out = delim.join([bam, str(target_size), str(seq_total), str(seq_total / target_size)] + cov_stats + ['\n'])
			outfh.write(line_out)



##index fasta file
def index_fasta_with_bwa(fa_file):
	bwa_index = subprocess.Popen(['bwa', 'index', fa_file])
	bwa_index.wait()


##align with bwa
def align_with_bwa(sample_dict, working_dir):
	for sample in sample_dict:
		r1_fq = sample_dict[sample][0]
		r2_fq = sample_dict[sample][1]
		print(sample, r1_fq, r2_fq)
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
		pe_bam = sample + post_bwa_bam
		sort_bam = sample + sorted_bam
		pic_dup_bam = sample + mkdup_bam
		# '''
		bwa_pe = subprocess.Popen(['bwa', 'mem', '-M', '-t', '5', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen(['samtools', 'view', '-q', '20', '-b', '-@', '5', '-o', pe_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()
		st_sort_pe = subprocess.Popen(['samtools', 'sort','-O', 'bam', '-o', sort_bam, '-T', sample,'-@', '10', '-m', '10G', pe_bam])
		st_sort_pe.wait()
		# '''
		##mark duplicates
		picard_md = subprocess.Popen(['picard', 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + pic_dup_bam, 'METRICS_FILE=' + sample + '.metrics', 'TMP_DIR=' + working_dir])
		picard_md.wait()

##make list of all bam files to be analyzed
def make_list_of_bams(sample_names, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample_name in sample_names:
			outfh.write(sample_name + bam_suffix + '\n')


##call samtools on bamfiles
def variant_calling_samtools(project, final_vcf_suffix, bamlist_file):
	vcf_temp1 = project + '.temp_st.vcf.gz'
	vcf_temp2 = project + '.temp_st2.vcf.gz'
	final_vcf = project + final_vcf_suffix
	# '''
	# stmp = subprocess.Popen([samtools,'mpileup', '-ug','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist_file], stdout=subprocess.PIPE)
	# bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '19', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
	stmp = subprocess.Popen(['samtools','mpileup', '-Ou','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist_file], stdout=subprocess.PIPE)
	bcft = subprocess.Popen(['bcftools','call', '-vmO', 'z', '--threads', '9', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
	bcft.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	# '''
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()


##make avinput files
def convert_to_annovar_move_to_annovar_folder(samples, vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
	con_ann.wait()
	for sample in samples:
		temp_av_file = 'temp.' + sample + '.avinput'
		av_file = sample + '.avinput'
		shutil.move(temp_av_file, av_file)
		shutil.copy(av_file, str(av_ref_dir[0]))

def run_table_annovar(vcf, av_prefix):
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def split_info_field(info_list):
	indices = [i for i, s in enumerate(info_list) if 'ANNOVAR_DATE' in s]
	# print indices
	i_count = 0
	final_list = []
	for info in info_list:
		# print info
		if i_count > indices[0] and info != 'ALLELE_END':
			info2 = info.split('=')[1]
			#print info2
			final_list.append(info2)
		i_count += 1
	return final_list

def multianno_to_annotated(out_prefix, samples):
	'''
	##split into individual samples
	post_annovar_vcf = out_prefix + '.' + av_genome + '_multianno.vcf'
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', post_annovar_vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', out_prefix])
	con_ann.wait()
	##format files...
	for sample in samples:
		avinput = out_prefix + '.' + sample + '.avinput'
		outfile = out_prefix + '.' + sample + '.annotated.txt'
		head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Zygosity', 'Qual', 'Coverage', 'Filter', 'Pos', 'Ref2', 'Alt2', 
		'Format', 'Info', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 
		'Func.ensGene', 'Gene.ensGene', 'GeneDetail.ensGene', 'ExonicFunc.ensGene', 'AAChange.ensGene', 'rmsk', 
		'SH_1_A5', 'SH_2_B3', 'SH_3_m', 'SH_4_f']
		head_out = delim.join(head + ['\n'])
		with open(avinput, "r") as av, open(outfile, "w") as final:
			final.write(head_out)
			for line in av:
				line = line.strip('\n').split(delim)
				if '_' not in line[0]:
					stuff = line[0:8] + [line[14]] + [line[9]] + line[11:13]
					info = line[15].split(';')
					info_list = split_info_field(info)
					other_stuff = line[16:]
					line_out = delim.join(stuff + other_stuff + info_list +['\n'])
					final.write(line_out)
	'''
	##for combined avinput
	##make avinput
	# post_annovar_vcf = out_prefix + '.' + av_genome + '_multianno.vcf'
	# con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', post_annovar_vcf, '-includeinfo', '-format', 'vcf4old', '-allallele', '-outfile', out_prefix + '.combined.avinput'])
	# con_ann.wait()
	##annotate
	avinput = out_prefix + '.combined.avinput'
	outfile = out_prefix + '.combined.annotated.txt'
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Qual',
	'Format', 'SH_1_A5.vcf', 'SH_2_B3.vcf', 'SH_3_m.vcf', 'SH_4_f.vcf', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 
	'Func.ensGene', 'Gene.ensGene', 'GeneDetail.ensGene', 'ExonicFunc.ensGene', 'AAChange.ensGene', 'rmsk', 
	'SH_1_A5', 'SH_2_B3', 'SH_3_m', 'SH_4_f']
	head_out = delim.join(head + ['\n'])
	with open(avinput, "r") as av, open(outfile, "w") as final:
		final.write(head_out)
		for line in av:
			line = line.strip('\n').split(delim)
			if '_' not in line[0]:
				stuff = line[0:5] + [line[10]]
				info = line[12].split(';')
				info_list = split_info_field(info)
				other_stuff = line[13:]
				line_out = delim.join(stuff + other_stuff + info_list + ['\n'])
				final.write(line_out)

def filter_parent_snps(in_file, female_snps, male_snps):
	female_snp_outfile = female_snps + '.q50_cov10.txt'
	male_snp_outfile = male_snps + '.q50_cov10.txt'
	female_snp_noalt_outfile = female_snps + '.q50_cov10.no_alt.txt'
	male_snp_noalt_outfile = male_snps + '.q50_cov10.no_alt.txt'
	with open(in_file, 'r') as in_fh, open(female_snp_outfile, 'w') as fout_fh, open(male_snp_outfile, 'w') as mout_fh, open(female_snp_noalt_outfile, 'w') as f2out_fh, open(male_snp_noalt_outfile, 'w') as m2out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc == 1:
				mout_fh.write(delim.join(line) + '\n')
				fout_fh.write(delim.join(line) + '\n')
				m2out_fh.write(delim.join(line) + '\n')
				f2out_fh.write(delim.join(line) + '\n')
			else:
				qual = float(line[5])
				rmsk = line[21]
				male_parent_vcf_info = line[7]
				female_parent_vcf_info = line[8]
				male_pool_vcf_info = line[9]
				female_pool_vcf_info = line[10]
				male_parent_genotype = line[22]
				female_parent_genotype = line[23]
				male_pool_genotype = line[24]
				female_pool_genotype = line[25]
				vcf_infos = line[7:11]
				# coverages = [int(i.split(':')[2]) for i in vcf_infos]
				coverage_info = [i.split(':')[4].split(',') for i in vcf_infos]
				coverages =[]
				for c in coverage_info:
					cov = [int(i) for i in c]
					sum_cov = sum(cov)
					coverages.append(sum_cov)
				print(coverage_info, coverages)
				# print(vcf_infos, coverages)
				##filter vars for coverage, qual and repeats
				if min(coverages) >= 10 and max(coverages) <= 100 and qual >= 50 and rmsk == '.':
					##get male parent is heterozygous (with high confidence) and female is homozygous for ref
					if male_parent_genotype == 'het' and female_parent_genotype == '.':
						mout_fh.write(delim.join(line) + '\n')
						f_gt_coverages = female_parent_vcf_info.split(':')[4].split(',')
						# print(f_gt_coverages)
						if len(f_gt_coverages) == 2 and f_gt_coverages[1] == '0':
							m2out_fh.write(delim.join(line) + '\n')
						else:
							if len(f_gt_coverages) != 2:
								print('multiallele:', line)
					if female_parent_genotype == 'het' and male_parent_genotype == '.':
						fout_fh.write(delim.join(line) + '\n')
						m_gt_coverages = male_parent_vcf_info.split(':')[4].split(',')
						if len(m_gt_coverages) == 2 and m_gt_coverages[1] == '0':
							f2out_fh.write(delim.join(line) + '\n')
						else:
							if len(m_gt_coverages) != 2:
								print('multiallele:', line)


def make_bed_from_filtered_file(infile):
	outbed = infile.rsplit('.', 1)[0] + '.bed'
	with open(outbed, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				genotypes = line[9:11]
				allele_info = [i.split(':')[4].split(',') for i in genotypes]
				nafs = []
				for ai in allele_info:
					alelles = [int(i) for i in ai]
					naf = str(float(alelles[1]) / sum(alelles))
					nafs.append(naf)
				print(line,allele_info, nafs)

				line_out = [line[0], str(int(line[1]) - 1), line[2]] + nafs
				out_fh.write(delim.join(line_out) + '\n')

def make_files_for_graphing(snp_beds, window_size, step_size, genome_fai, working_dir):
	for snp_bed in snp_beds:
		for ws in window_size:
			#make bed file with windows and returns genome name and window size variable
			genome_and_window = homozygosity_mapping_cybertron.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]
			print(genome_and_window)
			window_bed = genome_and_window + '.bed'
			##all shared snps
			out_file_naf = snp_bed.rsplit('.', 1)[0] + '.' + genome_and_window + '.aaf_for_r.txt'
			out_file_combined = snp_bed.rsplit('.', 1)[0] + '.' + genome_and_window + '.combined_for_r.txt'
			##bedtools intersect 
			with open('temp.bed', "w") as naf_fh: 
				# hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', window_bed, '-b', snp_bed, '-c'], stdout=naf_fh)
				# hom_bt_intersect.wait()
				hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', window_bed, '-b', snp_bed, '-wa', '-wb'], stdout=naf_fh)
				hom_bt_intersect.wait()
			with open('temp.bed', "r") as in_fh:
				aaf_dict = {}
				for line in in_fh:
					line = line.rstrip().split(delim)
					chr_start_end = delim.join(line[:3])
					maaf = float(line[6])
					faaf = float(line[7])
					if chr_start_end in aaf_dict:
						aaf_dict[chr_start_end][0].append(maaf)
						aaf_dict[chr_start_end][1].append(faaf)
					else:
						aaf_dict[chr_start_end] = [[maaf],[faaf]]
			##write outfile from dict
			with open(out_file_naf, "w") as outn_fh, open(out_file_combined, "w") as outc_fh:
				outn_fh.write(delim.join(['chr', 'start', 'end', 'test', 'average_aaf']) + '\n')
				outc_fh.write(delim.join(['chr', 'start', 'end', 'test', 'value']) + '\n')

				for window in aaf_dict:
					ave_maaf = sum(aaf_dict[window][0]) / len(aaf_dict[window][0])
					ave_faaf = sum(aaf_dict[window][1]) / len(aaf_dict[window][1])
					aaf_diff = ave_maaf - ave_faaf
					snp_count = len(aaf_dict[window][0])
					log_snp_count = math.log(snp_count, 2)
					outn_fh.write(window + delim + 'male_aaf' + delim + str(ave_maaf) + '\n')
					outn_fh.write(window + delim + 'female_aaf' + delim + str(ave_faaf) + '\n')
					outc_fh.write(window + delim + 'aaf_male' + delim + str(ave_maaf) + '\n')
					outc_fh.write(window + delim + 'aaf_female' + delim + str(ave_faaf) + '\n')
					outc_fh.write(window + delim + 'aaf_difference' + delim + str(aaf_diff) + '\n')
					# outc_fh.write(window + delim + 'snp_count' + delim + str(snp_count) + '\n')
					outc_fh.write(window + delim + 'log2_snp_count' + delim + str(log_snp_count) + '\n')

##run methods
work_dir = '/home/atimms/ngs_data/enu_mapping/dave_zf_sex_0920'
os.chdir(work_dir)

##params
project_name = 'dave_zf_sex_0920'
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
bamlist = 'bams.list'
st_vcf_suffix = '.st.vcf.gz'
result_vcf = project_name + st_vcf_suffix
fq_dict = {'SH_1_A5': ['SH_1_A5_CSFP200005686-1a_HJMKJDSXY_L1_1.fq.gz', 'SH_1_A5_CSFP200005686-1a_HJMKJDSXY_L1_2.fq.gz'],
	'SH_2_B3':['SH_2_B3_CSFP200005687-1a_HJMJNDSXY_L1_1.fq.gz', 'SH_2_B3_CSFP200005687-1a_HJMJNDSXY_L1_2.fq.gz'],
	'SH_3_m':['SH_3_m_CSFP200005688-1a_HJMKJDSXY_L1_1.fq.gz', 'SH_3_m_CSFP200005688-1a_HJMKJDSXY_L1_2.fq.gz'],
	'SH_4_f':['SH_4_f_CSFP200005689-1a_HJLM3DSXY_L1_1.fq.gz', 'SH_4_f_CSFP200005689-1a_HJLM3DSXY_L1_2.fq.gz']}


##index fasta file
# index_fasta_with_bwa(fasta)

##map with bwa and process with samtools etc
# align_with_bwa(fq_dict, work_dir)
# make_list_of_bams(fq_dict, mkdup_bam, bamlist)
# variant_calling_samtools(project_name, st_vcf_suffix, bamlist)

##coverage
# calculate_genome_coverage(fq_dict, mkdup_bam, bedtools_genome_file, project_name)

##annotate with annovar
# convert_to_annovar_move_to_annovar_folder(fq_dict, result_vcf)
# run_table_annovar(result_vcf, project_name)
# multianno_to_annotated(project_name, fq_dict)

##filter vars for graphing - for parent being het and pool being hom reference 
annotated_file = project_name + '.combined.annotated.txt'
mom_het_dad_hom_prefix = project_name + '.mom_het_dad_wt'
dad_het_mom_hom_prefix = project_name + '.dad_het_mom_wt'
##filter 2 ways for both test
# filter_parent_snps(annotated_file, mom_het_dad_hom_prefix, dad_het_mom_hom_prefix)

##make naf bed from filtered files
filtered_files = ['dave_zf_sex_0920.dad_het_mom_wt.q50_cov10.no_alt.txt', 'dave_zf_sex_0920.dad_het_mom_wt.q50_cov10.txt', 
			'dave_zf_sex_0920.mom_het_dad_wt.q50_cov10.no_alt.txt', 'dave_zf_sex_0920.mom_het_dad_wt.q50_cov10.txt']
# filtered_files = ['temp.txt'] ##testing
##make file into bed for garphing
# for filtered_file in filtered_files:
# 	make_bed_from_filtered_file(filtered_file)



##make files for graphing
filtered_beds = ['dave_zf_sex_0920.dad_het_mom_wt.q50_cov10.no_alt.bed', 'dave_zf_sex_0920.mom_het_dad_wt.q50_cov10.no_alt.bed']
window_sizes = [10000000]
step_size = 1000000
make_files_for_graphing(filtered_beds, window_sizes, step_size, fasta_fai, work_dir)







