#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import filtering_annotated
import homozygosity_mapping_cybertron


'''
must load 2 modules:
module load java/1.8.0_121 
module load biobuilds/2016.11
module load gcc/6.1.0 ##gcc version used to make bedtools using in homo mapping
'''

##parameters
delim = '\t'
thread_number = '18'

##working dir
# working_dir = '/data/atimms/timon_0317'
working_dir = '/home/atimms/ngs_data/enu_mapping/scott_0718'
os.chdir(working_dir)

##programs and files
gatk = '/cm/shared/apps/GATK/GenomeAnalysisTK.jar'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
fasta = '/home/atimms/ngs_data/references/mm10/mm10.fa'

bwa = 'bwa'
samtools = 'samtools'
bcftools = 'bcftools'
picard = 'picard'
##delly = '/home/atimms/programs/delly_v0.7.6_CentOS5.4_x86_64bit'
bedtools = 'bedtools'
bgzip = 'bgzip'


##files

# st_avinputs = ['K50000022.avinput', 'K50000024.avinput', 'K50000095.avinput', 'combined.avinput']


##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic,generic,generic,generic,generic,generic', '-genericdbfile', 'pleather_sc_st.K416.avinput,daredevil.avinput,bub.avinput,J308.st.avinput,J318.st.avinput,J320.st.avinput,J327.st.avinput,J328.st.avinput,J329.st.avinput,mgp.v5.merged.snps_all.dbSNP142.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,,,,,']
#av_options = ['-otherinfo', '-remove', '-vcfinput']


##filtering and hom mapping parameters
##filtering
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = 'synonymous SNV'
zygosity_col = 24
cov_col = 26
cov_definition = 5
qual_col = 25
qual_definition = 30
##het snp mapping 
genome_fai = '/home/atimms/ngs_data/references/mm10/mm10.fa.fai'
window_size = [1000000,5000000,2000000]
# window_size = [10000000]
step_size = 1000000
info_col = 36
naf_values = [0.8,0.9,0.95]



##methods
def combine_fq_file(r1_to_combine, r2_to_combine, r1_fq, r2_fq):
	print r1_fq, r1_to_combine
	print r2_fq, r2_to_combine
	with open(r1_fq, 'w') as r1_fh:
		cat_files = subprocess.Popen(['cat'] + r1_to_combine, stdout=r1_fh)
		cat_files.wait()
	with open(r2_fq, 'w') as r2_fh:
		cat_files = subprocess.Popen(['cat'] + r2_to_combine, stdout=r2_fh)
		cat_files.wait()


##align with bwa
def align_with_bwa(sample_dict):
	for sample in sample_dict:
		r1_fq = sample_dict[sample][0]
		r2_fq = sample_dict[sample][1]
		print sample, r1_fq, r2_fq
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
		pe_bam = sample + post_bwa_bam
		sort_bam = sample + sorted_bam
		pic_dup_bam = sample + mkdup_bam
		# '''
		bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-q', '20', '-b', '-@', '5', '-o', pe_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()
		st_sort_pe = subprocess.Popen([samtools, 'sort','-O', 'bam', '-o', sort_bam, '-T', sample,'-@', '10', '-m', '10G', pe_bam])
		st_sort_pe.wait()
		# '''
		##mark duplicates
		picard_md = subprocess.Popen([picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + pic_dup_bam, 'METRICS_FILE=' + sample + '.metrics', 'TMP_DIR=' + working_dir])
		picard_md.wait()

##make list of all bam files to be analyzed
def make_list_of_bams(sample_dict, bam_suffix, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in sample_dict:
			bam = sample + bam_suffix
			outfh.write(bam + '\n')




##call samtools on bamfiles
def variant_calling_samtools(samples, bam_suffix, final_vcf_suffix):
	for sample in samples:
		bamfile = sample + bam_suffix
		vcf_temp1 = sample + '.temp_st.vcf.gz'
		vcf_temp2 = sample + '.temp_st2.vcf.gz'
		final_vcf = sample + final_vcf_suffix
		# '''
		stmp = subprocess.Popen([samtools,'mpileup', '-ug','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, bamfile], stdout=subprocess.PIPE)
		bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '19', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
		bcft.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1])
		bcf_index.wait()
		#split multi-allelic variants calls in separate lines
		bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1])
		bcf_norm1.wait()
		# '''
		bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf +'.gz', vcf_temp2])
		bcf_norm2.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
		bcf_index.wait()



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
			# bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=subprocess.PIPE)
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

##make avinput files
def convert_to_annovar(samples, vcf_suffix):
	for sample in samples:
		vcf = sample + vcf_suffix
		con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', 'temp'])
		con_ann.wait()
	temp_files = glob.glob('temp*.avinput')
	for temp_file in temp_files:
		real_file = temp_file[5:]
		os.rename(temp_file, real_file)
		# shutil.copy(real_file, str(av_ref_dir[0]))

def run_table_annovar(samples):
	for sample in samples:
		avinput = sample + '.avinput'
		av_prefix = avinput.rsplit('.',1)[0]
		command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()

def multianno_to_annotated(samples): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142','pleather','daredevil','bub','J308','J318','J320','J327','J328','J329', 'mgp.v5.snps', 'Zygosity', 'Qual', 'Coverage', 'Chr', 'Pos', 'Filter', 'Ref2', 'Alt2','Qual', 'Filter', 'GT_info', 'Format', 'Info']
	head_out = delim.join(head + ['\n'])
	for sample in samples:
		avinput = sample + '.avinput'
		av_prefix = avinput.rsplit('.',1)[0]
		multianno = av_prefix + '.mm10_multianno.txt'
		annotated = av_prefix + '.annotated.txt'
		with open(multianno, "r") as multi, open(annotated, "w") as final:
			final.write(head_out)
			line_count = 0
			for line in multi:
				line_count += 1
				if line_count > 1:
					final.write(line)

def get_shared_snps(samples, file_suffix, out_file):
	snp_dict = {}
	##get all vars and put in a dict to get counts
	for sample in samples:
		in_file = sample + file_suffix
		with open(in_file, "r") as in_fh:
			line_count = 0
			for line in in_fh:
				line = line.split(delim)
				line_count += 1
				if line_count > 1:
					snp = '_'.join(line[:5])
					if snp in snp_dict:
						snp_dict[snp].append(sample)
					else:
						snp_dict[snp] = [sample]
	##get snps in all three sample
	snp_list = []
	for s in snp_dict:
		if len(snp_dict[s]) == 3:
			snp_list.append(s)
			# print s, snp_dict[s]
	print str(len(snp_list)) + ' snps shared in all sample'
	##make final file
	infile = samples[0] + file_suffix
	with open(infile, "r") as infh, open(out_file, "w") as outfh:
		lc = 0
		for line in infh:
			line = line.split(delim)
			lc += 1
			if lc == 1:
				outfh.write(delim.join(line))
			else:
				snp = '_'.join(line[:5])
				if snp in snp_list:
					outfh.write(delim.join(line))

def make_bed_from_ann_txt(infile, outfile):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				out_fh.write(delim.join(line[:3] + ['\n']))





##call methods
##parameters
project_name = 'adelie_0718'
# samples = ['K50000022']
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
st_vcf_suffix = '.st.vcf'
fq_dict = {'ad_02_8': ['ad_02_8_DSW64094-V_HNCT5CCXY_L1_1.fq.gz', 'ad_02_8_DSW64094-V_HNCT5CCXY_L1_2.fq.gz'], 
		'ad_66_6': ['ad_66_6_DSW64095-V_HNCT5CCXY_L1_1.fq.gz', 'ad_66_6_DSW64095-V_HNCT5CCXY_L1_2.fq.gz']}
samples = fq_dict.keys()
bamslist_file = 'bams.list'


##map with bwa and process with samtools etc
# align_with_bwa(fq_dict)
# variant_calling_samtools(fq_dict, mkdup_bam, st_vcf_suffix)
# convert_to_annovar(fq_dict, st_vcf_suffix + '.gz')
# run_table_annovar(fq_dict)
# multianno_to_annotated(fq_dict)




##filter variants for variants,  homozygsity mapping then counts
'''
##filter variants for candidates snps
for sample in samples:
	##exonic_variants
	filtering_annotated.filter(working_dir, "or", sample + '.annotated.txt' , sample + "_1.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", sample + "_1.temp", sample + "_2.temp", [col_function], ['!='], [syn_definition])
	##remove if in dbsnp, sanger, or other mouse line
	filtering_annotated.filter(working_dir, "and", sample + "_2.temp", sample + "_3.temp", [13,14,15,16,17,18,19,20,21,22,23], ['==','==','==','==','==','==','==','==','==','==','=='], ['','','','','','','','','','',''])
	##keep if hom
	filtering_annotated.filter(working_dir, "and", sample + "_3.temp", sample + '.hom_exonic_rare.xls', [zygosity_col], ['=='], ['hom'])
	##filter variants by coverage and quality 
	filtering_annotated.filter(working_dir, "and", sample + '.hom_exonic_rare.xls', sample + '.hom_exonic_rare_qual_filtered.xls', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
'''

	
'''
##homozygosity mapping
# window_size = [100000,500000,1000000,2000000]
window_size = [10000000]
step_size = 100000
for ws in window_size:
	for sample in samples:
		##remove if in dbsnp, sanger, other ped or rmsk
		filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "11.temp", [11,12,13,14,15,16,17,18,19,20,21,22,23], ['==','==','==','==','==','==','==','==','==','==','==','==','=='], ['','','','','','','','','','','','',''])
		# filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "11.temp", [11,12,23], ['==','==','=='], ['','',''])

		##filter variants by coverage and quality 
		filtering_annotated.filter(working_dir, "and", sample + "11.temp", sample + '.hom_temp.txt', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])

		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping_cybertron.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]

		##make bed file from variants
		homozygosity_mapping_cybertron.make_bed_from_ann(working_dir, 'samtools', sample + '.hom_temp.txt', zygosity_col, info_col)
		##hom and het count and hom percentage
		homozygosity_mapping_cybertron.count_and_percentage(working_dir, genome_and_window, sample + '.bed')
		##naf
		homozygosity_mapping_cybertron.naf_in_window(working_dir, genome_and_window, sample + '.bed')
		##total snp number
		homozygosity_mapping_cybertron.total_snp_in_window(working_dir, genome_and_window, sample + '.bed')

		##combine bedgraphs for graphing in r
		homozygosity_mapping_cybertron.combine_bedgraphs_for_r(working_dir, sample, genome_and_window)
'''


##parameters
project_name = 'dartagnan_0718'
# samples = ['K50000022']
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
st_vcf_suffix = '.st.vcf'
fq_dict = {'dar_63_6': ['dar_63_6_DSW64096-V_HNCT5CCXY_L1_1.fq.gz', 'dar_63_6_DSW64096-V_HNCT5CCXY_L1_2.fq.gz']}
samples = fq_dict.keys()
bamslist_file = 'bams.list'


##map with bwa and process with samtools etc
'''
align_with_bwa(fq_dict)
variant_calling_samtools(fq_dict, mkdup_bam, st_vcf_suffix)
convert_to_annovar(fq_dict, st_vcf_suffix + '.gz')
run_table_annovar(fq_dict)
multianno_to_annotated(fq_dict)
'''



##filter variants for variants,  homozygsity mapping then counts
'''
##filter variants for candidates snps
for sample in samples:
	##exonic_variants
	filtering_annotated.filter(working_dir, "or", sample + '.annotated.txt' , sample + "_1.temp", [col_exon, col_exon], ['==','=='], [exon_definition[0],exon_definition[1]])
	##remove synonymous
	filtering_annotated.filter(working_dir, "and", sample + "_1.temp", sample + "_2.temp", [col_function], ['!='], [syn_definition])
	##remove if in dbsnp, sanger, or other mouse line
	filtering_annotated.filter(working_dir, "and", sample + "_2.temp", sample + "_3.temp", [13,14,15,16,17,18,19,20,21,22,23], ['==','==','==','==','==','==','==','==','==','==','=='], ['','','','','','','','','','',''])
	##keep if hom
	filtering_annotated.filter(working_dir, "and", sample + "_3.temp", sample + '.hom_exonic_rare.xls', [zygosity_col], ['=='], ['hom'])
	##filter variants by coverage and quality 
	filtering_annotated.filter(working_dir, "and", sample + '.hom_exonic_rare.xls', sample + '.hom_exonic_rare_qual_filtered.xls', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])
'''

	
# '''
##homozygosity mapping
# window_size = [100000,500000,1000000,2000000]
window_size = [10000000]
step_size = 100000
for ws in window_size:
	for sample in samples:
		##remove if in dbsnp, sanger, other ped or rmsk
		filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "11.temp", [11,12,13,14,15,16,17,18,19,20,21,22,23], ['==','==','==','==','==','==','==','==','==','==','==','==','=='], ['','','','','','','','','','','','',''])
		# filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "11.temp", [11,12,23], ['==','==','=='], ['','',''])

		##filter variants by coverage and quality 
		filtering_annotated.filter(working_dir, "and", sample + "11.temp", sample + '.hom_temp.txt', [cov_col,qual_col], ['>=','>='], [cov_definition,qual_definition])

		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping_cybertron.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]

		##make bed file from variants
		homozygosity_mapping_cybertron.make_bed_from_ann(working_dir, 'samtools', sample + '.hom_temp.txt', zygosity_col, info_col)
		##hom and het count and hom percentage
		homozygosity_mapping_cybertron.count_and_percentage(working_dir, genome_and_window, sample + '.bed')
		##naf
		homozygosity_mapping_cybertron.naf_in_window(working_dir, genome_and_window, sample + '.bed')
		##total snp number
		homozygosity_mapping_cybertron.total_snp_in_window(working_dir, genome_and_window, sample + '.bed')

		##combine bedgraphs for graphing in r
		homozygosity_mapping_cybertron.combine_bedgraphs_for_r(working_dir, sample, genome_and_window)
# '''






