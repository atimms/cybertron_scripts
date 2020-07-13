#!/usr/bin/env python
import subprocess
import os
import glob
import shutil
import filtering_annotated
import homozygosity_mapping_cybertron
from collections import OrderedDict

##parameters
delim = '\t'
thread_number = '20'

##modules 
'''
module load biobuilds/2017.11
module load gcc/6.1.0 ##gcc version used to make bedtools using in homo mapping
'''

##programs and files
gatk = '/cm/shared/apps/GATK/GenomeAnalysisTK.jar'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
fasta = '/home/atimms/ngs_data/references/mm10/mm10.fa'
vcf_bed = '/home/atimms/programs/svtools/vcfToBedpe'

##first 2 columns of fasta.fai file
bedtools_genome_file  = '/home/atimms/ngs_data/references/mm10/mm10.fa.genome'
bwa = 'bwa'
samtools = 'samtools'
bcftools = 'bcftools'
picard = 'picard'
bedtools = 'bedtools'
bgzip = 'bgzip'
delly = '/home/atimms/programs/delly/src/delly'
manta_config = '/home/atimms/programs/manta-1.6.0.centos6_x86_64/bin/configManta.py'


##files
bamlist = 'bams.list'
# st_avinputs = ['K50000022.avinput', 'K50000024.avinput', 'K50000095.avinput', 'combined.avinput']
exome_bed_file = '/home/atimms/ngs_data/references/mm10/mm10_refGene_coding_exons.bed'
delly_exclude_regions = '/home/atimms/programs/delly/excludeTemplates/mouse.mm10.excl.tsv'

##file suffix
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
st_vcf_suffix = '.st.vcf.gz'

##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic,generic,generic,generic', '-genericdbfile', 'mgp.v5.merged.snps_all.dbSNP142.avinput,0620.kc1_wt.avinput,0620.kc2_ko.avinput,129P2_OlaHsd.mgp.v5.snps.dbSNP142.avinput,129S1_SvImJ.mgp.v5.snps.dbSNP142.avinput,129S5SvEvBrd.mgp.v5.snps.dbSNP142.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,,,,', '-vcfinput']
#av_options = ['-otherinfo', '-remove', '-vcfinput']

##filtering and hom mapping parameters
##filtering
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = 'synonymous SNV'
zygosity_col = 28
cov_col = 30
cov_definition = 10
qual_col = 29
qual_definition = 30
##het snp mapping 
genome_fai = '/home/atimms/ngs_data/references/mm10/mm10.fa.fai'
window_size = [1000000,5000000,2000000]
# window_size = [10000000]
step_size = 1000000
info_col = 40

##methods

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

def align_with_bwa_no_filter(sample_dict):
	for sample in sample_dict:
		r1_fq = sample_dict[sample][0]
		r2_fq = sample_dict[sample][1]
		print sample, r1_fq, r2_fq
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
		pe_bam = sample + '.temp1.bam'
		sort_bam = sample + '.temp2.bam'
		pic_dup_bam = sample + '.bwa_mkdup.no_filter.bam'
		# '''
		bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-o', pe_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()
		st_sort_pe = subprocess.Popen([samtools, 'sort','-O', 'bam', '-o', sort_bam, '-T', sample,'-@', '10', '-m', '10G', pe_bam])
		st_sort_pe.wait()
		# '''
		##mark duplicates
		picard_md = subprocess.Popen([picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + pic_dup_bam, 'METRICS_FILE=' + sample + '.metrics', 'TMP_DIR=' + working_dir])
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
	stmp = subprocess.Popen([samtools,'mpileup', '-Ou','-t', 'DP,DV,DPR', '-q', '20', '-C', '50', '-f', fasta, '-b', bamlist_file], stdout=subprocess.PIPE)
	bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '10', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
	bcft.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen([bcftools, 'norm', '-m-both', '-o', vcf_temp2, vcf_temp1])
	bcf_norm1.wait()
	# '''
	bcf_norm2 = subprocess.Popen([bcftools, 'norm', '-f', fasta, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()


##make avinput files
def convert_to_annovar_move_to_annovar_folder(samples, vcf):
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', '0620'])
	con_ann.wait()
	for sample in samples:
		av_file = '0620.' + sample + '.avinput'
		shutil.copy(av_file, str(av_ref_dir[0]))

def run_table_annovar(vcf, av_prefix):
	command = [table_annovar] + av_buildver + [vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def multianno_to_annotated(av_prefix): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'snp142', 'mgp.v5.snps', 'kc1_wt', 'kc2_ko', '129P2_OlaHsd', '129S1_SvImJ', '129S5SvEvBrd', 'vcf_info', 'vcf_format', 'vcf_kc1_wt', 'vcf_kc2_ko']
	head_out = delim.join(head + ['\n'])

	multianno = av_prefix + '.mm10_multianno.txt'
	annotated = av_prefix + '.annotated.txt'
	with open(multianno, "r") as multi, open(annotated, "w") as final:
		final.write(head_out)
		line_count = 0
		for line in multi:
			line = line.split(delim)
			line_count += 1
			if line_count > 1:
				line_out = line[:18] + line[28:]
				final.write(delim.join(line_out))


def filter_ann_file(file_prefix):
	##remove if in rmsk, segdup
	# filtering_annotated.filter(working_dir, "and", file_prefix + '.annotated.txt', file_prefix + "11.temp", [11,12], ['==','=='], ['.','.'])
	##in chr17 region
	# filtering_annotated.filter(working_dir, "and", file_prefix + '11.temp', file_prefix + '.chr17_20-30mb.xls', [1,3,3], ['==','>=','<='], ['chr17',20000000,30000000])
	# filtering_annotated.filter(working_dir, "and", file_prefix + '11.temp', file_prefix + '.chr17_20-40mb.xls', [1,3,3], ['==','>=','<='], ['chr17',20000000,40000000])
	# filtering_annotated.filter(working_dir, "and", file_prefix + '11.temp', file_prefix + '.chr17_20-60mb.xls', [1,3,3], ['==','>=','<='], ['chr17',20000000,60000000])
	# filtering_annotated.filter(working_dir, "and", file_prefix + '11.temp', file_prefix + '.chr17.xls', [1], ['=='], ['chr17'])

	# filtering_annotated.filter(working_dir, "and", file_prefix + '11.temp', file_prefix + '.chr17_20-30mb.xls', [1,3,3], ['==','>=','<='], ['chr17','20000000','30000000'])

	##unique to ko
	# filtering_annotated.filter(working_dir, "and", file_prefix + '.chr17_20-30mb.xls', file_prefix + '.chr17_20-30mb.ko_only.xls', [15], ['=='], ['.'])
	# filtering_annotated.filter(working_dir, "and", file_prefix + '.chr17_20-40mb.xls', file_prefix + '.chr17_20-40mb.ko_only.xls', [15], ['=='], ['.'])
	# filtering_annotated.filter(working_dir, "and", file_prefix + '.chr17_20-60mb.xls', file_prefix + '.chr17_20-60mb.ko_only.xls', [15], ['=='], ['.'])
	# filtering_annotated.filter(working_dir, "and", file_prefix + '.chr17.xls', file_prefix + '.chr17.ko_only.xls', [15], ['=='], ['.'])
	##in 129
	filtering_annotated.filter(working_dir, "or", file_prefix + '.chr17.ko_only.xls', file_prefix + '.chr17.ko_only.129.xls', [17,18,19], ['!=','!=','!='], ['.','.','.'])


def run_manta_svs(samples):
	outdir = 'manta_run'
	bams = []
	for sample in samples:
		in_bam = sample + '.bwa_mkdup.bam'
		bams.append(in_bam)
	##call svs
	manta_config_setup = subprocess.Popen([manta_config, '--bam', bams[0], '--bam', bams[1], '--referenceFasta', fasta, '--runDir', outdir])
	manta_config_setup.wait()
	manta_run = subprocess.Popen([outdir + '/runWorkflow.py'])
	manta_run.wait()

def filter_manta_results(in_vcf, out_file, bedfile):
	temp_bed = 'temp.bed'
	temp2_bed = 'temp2.bed'
	final_header = ['CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B', 'END_B', 'ID', 'QUAL', 'STRAND_A', 'STRAND_B', 'TYPE', 'FILTER', 'INFO', 'FORMAT', 'kc1_wt', 'kc2_ko', 'chr_gene', 'tx_start', 'tx_end', 'gene_name', 'overlaps']
	##make vcf into bed
	manta_config = subprocess.Popen([vcf_bed, '-i', in_vcf, '-o', temp_bed])
	manta_config.wait()	
	##use bedtools to add genename
	with open(temp2_bed, "w") as out_fh:
		manta_config = subprocess.Popen(['bedtools', 'intersect', '-a', temp_bed, '-b', bedfile, '-wao'], stdout=out_fh)
		manta_config.wait()		
	with open(temp2_bed, "r") as in_fh, open(out_file, "w") as out_fh:
		out_fh.write(delim.join(final_header) + '\n')
		for line in in_fh:
			line = line.split(delim)
			c1 = line[0]
			c2 = line[3]
			start = int(line[1])
			if c1 == c2:
				end = int(line[5])
			else:
				end = int(line[2])
			# if c1 == 'chr17' and end > 25000000 and start < 27000000:
			if c1 == 'chr17':
				out_fh.write(delim.join(line))



def make_bed_from_ann_txt(infile, outfile):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				out_fh.write(delim.join(line[:3]) + '\n')


def make_file_for_r(all_file, ko_only_file, out_file):
	out_dict = OrderedDict()
	# out_dict = {}
	with open(all_file, "r") as all_fh:
		lc = 0
		for line in all_fh:
			lc += 1
			if lc >1:
				line = line.rstrip().split(delim)
				line = [line[0].replace('chr','')] + line[1:] #removes chr from start of line
				start = int(line[1])
				end = int(line[2])
				count = int(line[3])
				midpoint = start + ((end - start) / 2)
				info = '_'.join(line[:3])
				out_dict[info] = [midpoint, count]
	with open(ko_only_file, "r") as ko_fh:
		lc = 0
		for line in ko_fh:
			lc += 1
			if lc >1:
				line = line.rstrip().split(delim)
				line = [line[0].replace('chr','')] + line[1:] #removes chr from start of line
				count = float(line[3])
				info = '_'.join(line[:3])
				out_dict[info].append(count)
	with open(out_file, "w") as out_fh:
		out_fh.write(delim.join(['chr', 'start', 'end', 'position', 'type', 'count', '\n']))
		for r in out_dict:
			region = r.split('_')
			print(out_dict[r][2], out_dict[r][1])
			if out_dict[r][1] == 0:
				ko_perc = 0
			else:
				ko_perc = out_dict[r][2] / out_dict[r][1]
			out_fh.write(delim.join(region + [str(out_dict[r][0]), 'total_snps', str(out_dict[r][1])]) + '\n')
			out_fh.write(delim.join(region + [str(out_dict[r][0]), 'ko_snps', str(int(out_dict[r][2]))]) + '\n')
			out_fh.write(delim.join(region + [str(out_dict[r][0]), 'ko_percentage', str(ko_perc)]) + '\n')


##run_methods
working_dir = '/home/atimms/ngs_data/misc/kenny_mouse_0620/wgs'
os.chdir(working_dir)

project_name = 'kenny_wgs_0620'
fq_dict = {'kc1_wt':['kc1_wt.r1.fastq.gz', 'kc1_wt.r2.fastq.gz'], 'kc2_ko':['kc2_ko.r1.fastq.gz', 'kc2_ko.r2.fastq.gz']}
samples = fq_dict.keys()
result_vcf = project_name + st_vcf_suffix

##align all samples
# align_with_bwa(fq_dict)

##align ko without the q score filter
fq_dict_ko = {'kc2_ko':['kc2_ko.r1.fastq.gz', 'kc2_ko.r2.fastq.gz']}
# align_with_bwa_no_filter(fq_dict_ko)

##call vars
# make_list_of_bams(samples, mkdup_bam, bamlist)
# variant_calling_samtools(project_name, st_vcf_suffix, bamlist)

##annotate and filter
# convert_to_annovar_move_to_annovar_folder(samples, result_vcf)
# run_table_annovar(result_vcf, project_name)
# multianno_to_annotated(project_name)
# filter_ann_file(project_name)

##get SVs
# run_manta_svs(fq_dict)
##filter and add genes to manta results, if called in affected in regions
manta_vcf = 'manta_run/results/variants/diploidSV.vcf'
manta_filter_file = 'kenny_wgs_0620.chr17_25-27mb.manta_svs.0620.xls'
refgene_bed = '/home/atimms/ngs_data/references/mm10/mm10.refGene_genes.bed'
# filter_manta_results(manta_vcf, manta_filter_file,refgene_bed)
manta_filter_file = 'kenny_wgs_0620.chr17.manta_svs.0620.xls'
# filter_manta_results(manta_vcf, manta_filter_file,refgene_bed)


##make graphs from ann.txt files
##mapping parameters
genome_fai = '/home/atimms/ngs_data/references/mm10/mm10.fa.fai'
window_size = [200000, 500000, 1000000, 2000000]
# window_size = [10000000]
step_size = 200000
# info_col = 40
files_to_graph = ['kenny_wgs_0620.chr17_20-30mb.xls', 'kenny_wgs_0620.chr17_20-30mb.ko_only.xls',
		'kenny_wgs_0620.chr17_20-40mb.xls', 'kenny_wgs_0620.chr17_20-40mb.ko_only.xls',
		'kenny_wgs_0620.chr17_20-60mb.xls', 'kenny_wgs_0620.chr17_20-60mb.ko_only.xls',
		'kenny_wgs_0620.chr17.xls', 'kenny_wgs_0620.chr17.ko_only.xls', 'kenny_wgs_0620.chr17.ko_only.129.xls']
beds_to_graph = []
# '''
for file_to_graph in files_to_graph:
	bed_to_graph = file_to_graph.rsplit('.', 1)[0] + '.bed'
	beds_to_graph.append(bed_to_graph)
	make_bed_from_ann_txt(file_to_graph, bed_to_graph)

for bed_to_graph in beds_to_graph:
	for ws in window_size:
		#make bed file with windows and returns genome name and window size variable
		genome_and_window = homozygosity_mapping_cybertron.make_windows(working_dir, genome_fai, ws, step_size).split('.')[0]
		print genome_and_window
		window_bed = genome_and_window + '.bed'
		##all shared snps
		out_bed = bed_to_graph.rsplit('.', 1)[0] + '.' + genome_and_window + '.bed'
		##bedtools intersect 
		with open('temp.bed', "w") as naf_fh: 
			hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', window_bed, '-b', bed_to_graph, '-c'], stdout=naf_fh)
			hom_bt_intersect.wait()
		##filter for region of interest
		if '30mb' in bed_to_graph:
			filtering_annotated.filter(working_dir, "and", 'temp.bed', 'temp2.bed', [1,3,3], ['==','>','<='], ['chr17',20000000,30000000])
		elif '40mb' in bed_to_graph:
			filtering_annotated.filter(working_dir, "and", 'temp.bed', 'temp2.bed', [1,3,3], ['==','>','<='], ['chr17',20000000,40000000])
		elif '60mb' in bed_to_graph:
			filtering_annotated.filter(working_dir, "and", 'temp.bed', 'temp2.bed', [1,3,3], ['==','>','<='], ['chr17',20000000,60000000])

		else:
			filtering_annotated.filter(working_dir, "and", 'temp.bed', 'temp2.bed', [1], ['=='], ['chr17'])

		##add header
		with open(out_bed, "w") as out_fh, open('temp2.bed', "r") as in_fh:
			out_fh.write(delim.join(['chr', 'start', 'end', 'snp_number', '\n']))
			for line in in_fh:
				if line.startswith('chr17'):
					out_fh.write(line)


for ws in window_size:
	ws = int(ws/1000)
	ss = int(step_size/1000)
	all_snps = 'kenny_wgs_0620.chr17_20-30mb.mm10_' + str(ws) + 'kb_200kb.bed'
	ko_only_snps = 'kenny_wgs_0620.chr17_20-30mb.ko_only.mm10_' + str(ws) + 'kb_200kb.bed'
	combined_file = 'kenny_wgs_0620.chr17_20-30mb.' + str(ws) + 'kb.file_for_graphing.txt'
	make_file_for_r(all_snps, ko_only_snps, combined_file)
	all_snps = 'kenny_wgs_0620.chr17_20-40mb.mm10_' + str(ws) + 'kb_200kb.bed'
	ko_only_snps = 'kenny_wgs_0620.chr17_20-40mb.ko_only.mm10_' + str(ws) + 'kb_200kb.bed'
	combined_file = 'kenny_wgs_0620.chr17_20-40mb.' + str(ws) + 'kb.file_for_graphing.txt'
	make_file_for_r(all_snps, ko_only_snps, combined_file)
	all_snps = 'kenny_wgs_0620.chr17_20-60mb.mm10_' + str(ws) + 'kb_200kb.bed'
	ko_only_snps = 'kenny_wgs_0620.chr17_20-60mb.ko_only.mm10_' + str(ws) + 'kb_200kb.bed'
	combined_file = 'kenny_wgs_0620.chr17_20-60mb.' + str(ws) + 'kb.file_for_graphing.txt'
	make_file_for_r(all_snps, ko_only_snps, combined_file)
	all_snps = 'kenny_wgs_0620.chr17.mm10_' + str(ws) + 'kb_200kb.bed'
	ko_only_snps = 'kenny_wgs_0620.chr17.ko_only.mm10_' + str(ws) + 'kb_200kb.bed'
	combined_file = 'kenny_wgs_0620.chr17.' + str(ws) + 'kb.file_for_graphing.txt'
	make_file_for_r(all_snps, ko_only_snps, combined_file)
	all_snps = 'kenny_wgs_0620.chr17.mm10_' + str(ws) + 'kb_200kb.bed'
	ko_only_snps = 'kenny_wgs_0620.chr17.ko_only.129.mm10_' + str(ws) + 'kb_200kb.bed'
	combined_file = 'kenny_wgs_0620.chr17.129_ko.' + str(ws) + 'kb.file_for_graphing.txt'
	make_file_for_r(all_snps, ko_only_snps, combined_file)
# '''




