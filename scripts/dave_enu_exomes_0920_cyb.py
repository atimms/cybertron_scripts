#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import filtering_annotated
import homozygosity_mapping_cybertron

##parameters
delim = '\t'
thread_number = '20'

##modules 
'''
module load biobuilds/2017.11
module load gcc/6.1.0 ##for homozygosity_mapping_cybertron

'''


##files
fasta = '/home/atimms/ngs_data/references/mm10/mm10.fa'
exome_bed_file = '/home/atimms/ngs_data/references/mm10/mm10_refGene_coding_exons.bed'
bedtools_genome_file  = '/home/atimms/ngs_data/references/mm10/mm10.fa.genome'

##prograrm
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'

##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic', '-genericdbfile', 'mgp.v5.merged.snps_all.dbSNP142.avinput,PKD0005.avinput,PKD0006.avinput,PKD0007.avinput,PKD0016.avinput,PKD0017.avinput,PKD0018.avinput,PKD0063.avinput,PKD0064.avinput,PKD0065.avinput,PKD0082.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,,,,,,', '-vcfinput']
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


##get exome coverage for bam files
def calculate_exome_coverage(samples, bam_suffix, exome_bed, prefix):
	##parameters and header
	coverage_required = [1,5,10,20,50,100,200]
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
	output_file = prefix + '.exome.coverage_data.txt'
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
			# bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen(['bedtools', 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=subprocess.PIPE)
			grep_cmd = subprocess.Popen(['grep', '^all'], stdin=bedtools_cov.stdout, stdout=hist_fh)
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
	##for combined avinput
	##make avinput
	# post_annovar_vcf = out_prefix + '.' + av_genome + '_multianno.vcf'
	# con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', post_annovar_vcf, '-includeinfo', '-format', 'vcf4old', '-allallele', '-outfile', out_prefix + '.combined.avinput'])
	# con_ann.wait()
	#re format
	avinput = out_prefix + '.combined.avinput'
	outfile = out_prefix + '.combined.annotated.txt'
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Qual',
	'Format', 'PKD0064_vcf', 'PKD0065_vcf', 'PKD0017_vcf', 'PKD0016_vcf', 'PKD0006_vcf', 'PKD0007_vcf', 'PKD0005_vcf', 
	'PKD0082_vcf', 'PKD0063_vcf', 'PKD0018_vcf', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 
	'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142', 'mgp.v5', 'PKD0005', 'PKD0006', 'PKD0007', 'PKD0016', 'PKD0017', 
	'PKD0018', 'PKD0063', 'PKD0064', 'PKD0065', 'PKD0082']
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


def filter_snps_get_unique(in_file, out_suffix, combined_prefix):
	##filter files q30, cov10 in all samples, not in rmsk, supdups, dbSNP or mgp
	all_filtered_file = combined_prefix + out_suffix
	with open(in_file, 'r') as in_fh, open(all_filtered_file, 'w') as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc == 1:
				out_fh.write(delim.join(line) + '\n')
				samples = line[26:]
			else:
				qual = float(line[5])
				rmsk = line[22]
				supdup = line[23]
				dbsnp = line[24]
				mgp = line[25]
				vcf_infos = line[7:17]
				# coverages = [int(i.split(':')[2]) for i in vcf_infos]
				coverage_info = [i.split(':')[4].split(',') for i in vcf_infos]
				coverages =[]
				for c in coverage_info:
					cov = [int(i) for i in c]
					sum_cov = sum(cov)
					coverages.append(sum_cov)
				# print(coverage_info, coverages)
				# print(vcf_infos, coverages)
				##filter vars for coverage, qual and repeats etc
				if min(coverages) >= 10 and qual >= 50 and rmsk == '.' and supdup == '.' and dbsnp == '.' and mgp == '.':
					out_fh.write(delim.join(line) + '\n')
	for sample in samples:
		outfile = sample + out_suffix
		with open(all_filtered_file, 'r') as aff_fh, open(outfile, 'w') as out_fh:
			lc = 0
			for line in aff_fh:
				lc += 1
				line = line.rstrip().split(delim)
				if lc == 1:
					s_index = line.index(sample)
					out_fh.write(delim.join(line) + '\n')
					# print(sample, s_index)
				else:
					sample_gt = line[s_index]
					gts = line[26:]
					wt_count = gts.count('.')
					print(gts, wt_count, sample_gt)
					if sample_gt != '.' and wt_count == 9:
						out_fh.write(delim.join(line) + '\n')



##run methods
work_dir = '/home/atimms/ngs_data/enu_mapping/dave_enu_exomes_0920'
os.chdir(work_dir)

##params
project_name = 'dave_pkd_exomes_0920'
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
bamlist = 'bams.list'
st_vcf_suffix = '.st.vcf.gz'
result_vcf = project_name + st_vcf_suffix
fq_dict = {'PKD0005': ['ENUxPKD0005_FDHE20H106182-1a-7UDI1353-5UDI1353_1.fq.gz', 'ENUxPKD0005_FDHE20H106182-1a-7UDI1353-5UDI1353_2.fq.gz'],
	'PKD0006': ['ENUxPKD0006_FDHE20H106184-1a-7UDI1337-5UDI1337_1.fq.gz', 'ENUxPKD0006_FDHE20H106184-1a-7UDI1337-5UDI1337_2.fq.gz'],
	'PKD0007': ['ENUxPKD0007_FDHE20H106182-1a-7UDI1356-5UDI1356_1.fq.gz', 'ENUxPKD0007_FDHE20H106182-1a-7UDI1356-5UDI1356_2.fq.gz'],
	'PKD0016': ['ENUxPKD0016_FDHE20H106183-1a-7UDI1338-5UDI1338_1.fq.gz', 'ENUxPKD0016_FDHE20H106183-1a-7UDI1338-5UDI1338_2.fq.gz'],
	'PKD0017': ['ENUxPKD0017_FDHE20H106184-1a-7UDI1341-5UDI1341_1.fq.gz', 'ENUxPKD0017_FDHE20H106184-1a-7UDI1341-5UDI1341_2.fq.gz'],
	'PKD0018': ['ENUxPKD0018_FDHE20H106183-1a-7UDI1343-5UDI1343_1.fq.gz', 'ENUxPKD0018_FDHE20H106183-1a-7UDI1343-5UDI1343_2.fq.gz'],
	'PKD0063': ['ENUxPKD0063_FDHE20H106183-1a-7UDI1344-5UDI1344_1.fq.gz', 'ENUxPKD0063_FDHE20H106183-1a-7UDI1344-5UDI1344_2.fq.gz'],
	'PKD0064': ['ENUxPKD0064_FDHE20H106184-1a-7UDI1346-5UDI1346_1.fq.gz', 'ENUxPKD0064_FDHE20H106184-1a-7UDI1346-5UDI1346_2.fq.gz'],
	'PKD0065': ['ENUxPKD0065_FDHE20H106184-1a-7UDI1357-5UDI1357_1.fq.gz', 'ENUxPKD0065_FDHE20H106184-1a-7UDI1357-5UDI1357_2.fq.gz'],
	'PKD0082': ['ENUxPKD0082_FDHE20H106183-1a-7UDI1348-5UDI1348_1.fq.gz', 'ENUxPKD0082_FDHE20H106183-1a-7UDI1348-5UDI1348_2.fq.gz']}


##index fasta file

##map with bwa and process with samtools etc
# align_with_bwa(fq_dict, work_dir)
# make_list_of_bams(fq_dict, mkdup_bam, bamlist)
# variant_calling_samtools(project_name, st_vcf_suffix, bamlist)

##coverage
# calculate_exome_coverage(fq_dict, mkdup_bam, exome_bed_file, project_name)
calculate_genome_coverage(fq_dict, mkdup_bam, bedtools_genome_file, project_name)


##annotate with annovar
# convert_to_annovar_move_to_annovar_folder(fq_dict, result_vcf)
# run_table_annovar(result_vcf, project_name)
# multianno_to_annotated(project_name, fq_dict)

##filter
annotated = project_name + '.combined.annotated.txt'
filter_suffix = '.q50_cov10_rare.txt'
# filter_snps_get_unique(annotated, filter_suffix, project_name)












