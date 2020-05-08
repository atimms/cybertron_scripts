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
for standad analysis:
module load biobuilds/2017.11
module load gcc/6.1.0 ##for homozygosity_mapping_cybertron

for sv analsis with sve:
start interactive session
load modules:
module load local_python/2.7.14
module load R/3.3.3

install packages, just once:
conda install htseq
conda install -c bioconda bx-python
conda install -c bioconda crossmap
conda install -c bioconda mygene
install root:
dl from https://root.cern.ch/downloading-root
mv to cluster
tar xvzf root_v5.34.20.Linux-slc5_amd64-gcc4.3.tar.gz
export ROOTSYS=/home/atimms/programs/root/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
install SVE:
git clone --recursive https://github.com/TheJacksonLaboratory/SVE.git
cd SVE
make

'''


##working dir
working_dir = '/home/atimms/ngs_data/enu_mapping/warman_enu_0919'
os.chdir(working_dir)

##programs and files
gatk = '/cm/shared/apps/GATK/GenomeAnalysisTK.jar'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
fasta = '/home/atimms/ngs_data/references/mm10/mm10.fa'
##first 2 columns of fasta.fai file
bedtools_genome_file  = '/home/atimms/ngs_data/references/mm10/mm10.fa.genome'
bwa = 'bwa'
samtools = 'samtools'
bcftools = 'bcftools'
picard = 'picard'
bedtools = 'bedtools'
bgzip = 'bgzip'
freebayes = '/home/atimms/programs/freebayes/bin/freebayes'
manta_config = '/home/atimms/programs/manta-1.6.0.centos6_x86_64/bin/configManta.py'

##files
bamslist_file = 'bams.list'
# st_avinputs = ['K50000022.avinput', 'K50000024.avinput', 'K50000095.avinput', 'combined.avinput']


##annovar parameters
av_genome = 'mm10'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
##filter on daredevil, bub, pleather and some of jabiers mice
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,snp142,generic,generic,generic', '-genericdbfile', 'F1446.avinput,F1477.avinput,mgp.v5.merged.snps_all.dbSNP142.avinput']
av_operation = ['-operation', 'g,r,r,f,f,f,f']
av_options = ['-otherinfo', '-remove','-arg', '-splicing 10 ,,,,,,']
#av_options = ['-otherinfo', '-remove', '-vcfinput']


##filtering and hom mapping parameters
##filtering
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = 'synonymous SNV'
zygosity_col = 24
cov_col = 26
cov_definition = 10
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
	print r1_fq, r1_to_combine, len(r1_to_combine)
	print r2_fq, r2_to_combine, len(r2_to_combine)
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
		shutil.copy(real_file, str(av_ref_dir[0]))

def run_table_annovar(samples):
	for sample in samples:
		avinput = sample + '.avinput'
		av_prefix = avinput.rsplit('.',1)[0]
		command = [table_annovar] + av_buildver + [avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()

def multianno_to_annotated(samples): 
	head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 
	'genomicSuperDups', 'snp142','F1446.avinput','F1477.avinput', 'mgp.v5.snps', 'Zygosity', 'Qual', 'Coverage', 'Chr', 'Pos', 'Filter', 
	'Ref2', 'Alt2','Qual', 'Filter', 'GT_info', 'Format', 'Info']
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

def make_bed_from_ann_txt(infile, outfile):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count >1:
				line = line.rstrip().split(delim)
				out_fh.write(delim.join(line[:3]) + '\n')

def make_list_of_bams(bam_files, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_files:
			outfh.write(bam + '\n')

def genotype_vars_with_freebayes(bedfile, final_vcf, bams):
	bamlist = 'bams.list'
	make_list_of_bams(bams, bamlist)
	with open(final_vcf, 'w') as vcf_fh:
		run_freebayes = subprocess.Popen([freebayes, '-f', fasta, '-L', bamlist, '-t', bedfile, '--haplotype-length', '0', '--min-alternate-count', '1', '--min-alternate-fraction', '0', '--pooled-continuous', '--report-monomorphic'], stdout=vcf_fh)
		run_freebayes.wait()


def run_manta_svs(bam_files, project_name, working_dir):
	bams = []
	for bam_file in bam_files:
		in_bam = ['--bam', bam_file]
		bams.extend(in_bam)
	##call svs
	run_manta_config = subprocess.Popen([manta_config] + bams + ['--referenceFasta', fasta, '--runDir', working_dir + '/' + project_name ])
	run_manta_config.wait()
	manta_run = subprocess.Popen([working_dir + '/' + project_name + '/runWorkflow.py'])
	manta_run.wait()

##call methods
##parameters
project_name = 'warman_mapping_0919'
# samples = ['K50000022']
post_bwa_bam = '.bwa.bam'
sorted_bam = '.bwa_sorted.bam'
mkdup_bam = '.bwa_mkdup.bam'
st_vcf_suffix = '.st.vcf'

##dicts... note F1446 in the mutant and F1477 the wt
fq_combine_dict = {'F1446': 'DKWGS190823001','F1477': 'DKWGS190823002' }
fq_dict = {}
##combine fq files and make fq dict
for sample in fq_combine_dict:
	r1_fq = sample + '_1.fq.gz'
	r2_fq = sample + '_2.fq.gz'
	r1s = sorted(glob.glob(fq_combine_dict[sample] + '*_1.fq.gz'))
	r2s = sorted(glob.glob(fq_combine_dict[sample] + '*_2.fq.gz'))
	fq_dict[sample] = [r1_fq, r2_fq]
	# combine_fq_file(r1s, r2s, r1_fq, r2_fq)

##map with bwa and process with samtools etc
'''
print(fq_dict)
align_with_bwa(fq_dict)
variant_calling_samtools(fq_dict, mkdup_bam, st_vcf_suffix)
convert_to_annovar(fq_dict, st_vcf_suffix + '.gz')
run_table_annovar(fq_dict)
multianno_to_annotated(fq_dict)
'''

##coverage
'''
calculate_genome_coverage(samples, mkdup_bam, bedtools_genome_file, project_name)
calculate_exome_coverage(samples, mkdup_bam, exome_bed_file, project_name)
'''

##genotype markers
bams_to_genotype = ['F1446.bwa_mkdup.bam', 'F1477.bwa_mkdup.bam']
genotype_vcf = project_name + '.genotype_markers.vcf'
marker_bed = 'test.bed'
# genotype_vars_with_freebayes(marker_bed, genotype_vcf, bams_to_genotype)
genotype_vcf = project_name + '.genotype_markers_0120.vcf'
marker_bed = '0120.bed'
# genotype_vars_with_freebayes(marker_bed, genotype_vcf, bams_to_genotype)

##filter and graphing data
##parameters
##filtering -- column nubers start at 1
col_exon = 6
exon_definition = ['exonic', 'splicing']
col_function = 9
syn_definition = 'synonymous SNV'
zygosity_col = 17
cov_col = 19
cov_definition = [20,200]
qual_col = 25
qual_definition = 30
##het snp mapping 
genome_fai = '/home/atimms/ngs_data/references/mm10/mm10.fa.fai'
# window_size = [1000000,5000000,2000000]
# window_size = [5000000]
step_size = 1000000
info_col = 40
naf_values = [0.8,0.9,0.95]

##filter vars
'''
for sample in ['F1446']:
	##filter variants by coverage and quality 
	filtering_annotated.filter(working_dir, "and", sample + '.annotated.txt', sample + "2.temp", [cov_col,cov_col,qual_col], ['>=','<=','>='], [cov_definition[0], cov_definition[1],qual_definition])
	##remove if in rpt region
	filtering_annotated.filter(working_dir, "and", sample + '2.temp', sample + '.all_vars.xls', [11,12], ['==','=='], ['',''])
	##filter variants by all affected being hom
	filtering_annotated.filter(working_dir, "and", sample + '.all_vars.xls', sample + '.hom_vars.xls', [14], ['=='], ['hom'])
	##all vars not in any unaffected
	filtering_annotated.filter(working_dir, "and", sample + '.all_vars.xls', sample +'.not_unaffected.all_vars.xls', [15], ['=='], [''])
	#hom vars not hom in unaffected
	filtering_annotated.filter(working_dir, "and", sample + '.hom_vars.xls', sample +'.not_hom_unaffected.hom_vars.xls', [15], ['!='], ['hom'])
	filtering_annotated.filter(working_dir, "and", sample + '.hom_vars.xls', sample +'.not_unaffected.hom_vars.xls', [15], ['=='], [''])
'''

##make files to graph the data
files_to_graph = ['F1446.all_vars.xls', 'F1446.hom_vars.xls', 'F1446.not_unaffected.all_vars.xls', 'F1446.not_hom_unaffected.hom_vars.xls', 'F1446.not_unaffected.hom_vars.xls']
# files_to_graph = ['Smo90.in_all_aff_not_unaffected.all_vars.xls', 'Smo90.in_all_aff_not_hom_unaffected.hom_vars.xls']

##make bed from enu vars
'''
beds_to_graph = []
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
		##add header
		with open(out_bed, "w") as out_fh, open('temp.bed', "r") as in_fh:
			out_fh.write(delim.join(['chr', 'start', 'end', 'snp_number', '\n']))
			for line in in_fh:
				out_fh.write(line)
'''

##run SVs on manta
run_manta_svs(bams_to_genotype, project_name, working_dir)

