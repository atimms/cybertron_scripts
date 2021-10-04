#!/usr/bin/env python
import subprocess
import os

##info

'''
qsub -Iq cdbrmq -l mem=120gb,ncpus=10 -P 0f241fcd-0dd3-491c-acac-15c6cff1c880
module load java/1.8.0_202
module load R #for vardict
'''



##parameters
delim = '\t'
##programs
picard = '/home/atimms/programs/picard_2.25.6/picard.jar'
bwa = '/home/atimms/programs/bwa-0.7.17/bwa'
fgbio = '/home/atimms/programs/fgbio_1.3/fgbio-1.3.0.jar'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
samtools = '/home/atimms/programs/samtools-1.11/bin/samtools'
# vardict_dir = '/home/atimms/programs/VarDictJava/VarDict/'
vardict_dir = '/home/atimms/programs/VarDictJava/dist/VarDict-1.8.0/bin/'
vardict = vardict_dir + 'VarDict'
teststrandbias = vardict_dir + 'teststrandbias.R'
var2vcf_valid = vardict_dir + 'var2vcf_valid.pl'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'

##files 
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,gnomad211_exome,gnomad211_genome,dbnsfp41a,clinvar_20210123']
av_operation = ['-operation', 'g,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput' ]

##params
vardict_maf_req = '0.0025'


##methods
def get_sample_dict(in_file):
	s_dict = {}
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			sample = line[0]
			fqs = line[1:]
			if sample in s_dict:
				print(sample, ' seen multiple times')
			else:
				s_dict[sample] = fqs
	return(s_dict)

def make_unmapped_bam(s_dict):
	for sample in s_dict:
		unmapped_bam = sample + '_unmapped_temp.bam'
		fq1 = s_dict[sample][0]
		fq2 = s_dict[sample][1]
		#java -jar picard.jar FastqToSam FASTQ=unmapped_R1.fastq.gz FASTQ2=unmapped_R2.fastq.gz O=unmapped.bam SM=sample
		picard_ms = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'FastqToSam', 'FASTQ=' + fq1, 'FASTQ2=' + fq2, 'OUTPUT=' + unmapped_bam, 'SM=' + sample])
		picard_ms.wait()

def add_umi_unmapped_bam(s_dict):
	for sample in s_dict:
		unmapped_bam = sample + '_unmapped_temp.bam'
		unmapped_bam_with_umi = sample + '_unmapped_with_umi_temp.bam'
		#java -jar fgbio.jar ExtractUmisFromBam --input=unmapped.bam --output=unmapped.withUMI.bam --read-structure=8M143T 8M143T --molecular-index-tags=ZA ZB --single-tag=RX
		fgbio_eu = subprocess.Popen(['java', '-jar', fgbio, 'ExtractUmisFromBam', '--input=' + unmapped_bam, '--output=' + unmapped_bam_with_umi, '--read-structure=8M143T', '8M143T', '--molecular-index-tags=ZA', 'ZB', '--single-tag=RX'])
		fgbio_eu.wait()

def align_reads(s_dict):
	for sample in s_dict:
		unmapped_bam_with_umi = sample + '_unmapped_with_umi_temp.bam'
		temp_mapped_bam = sample + '_bwa_temp.bam'
		mapped_bam = sample + '_bwa.bam'
		
		#java -jar picard.jar SamToFastq I=unmapped.withUMI.bam F=/dev/stdout
		picard_sf = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'SamToFastq', 'I=' + unmapped_bam_with_umi, 'F=/dev/stdout', 'INTERLEAVE=true'], stdout=subprocess.PIPE)
		#bwa mem -p -t 4 hg38.fa /dev/stdin
		with open(temp_mapped_bam, "w") as tmb_fh:
			bwa_int = subprocess.Popen([bwa, 'mem', '-p', '-t', '4', fasta, '/dev/stdin'], stdin=picard_sf.stdout, stdout=tmb_fh)
			bwa_int.wait()
		#java -jar picard.jar MergeBamAlignment UNMAPPED=unmapped.withUMI.bam ALIGNED=/dev/stdin O=mapped.bam R=hg38.fa \ SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
		picard_mba = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'MergeBamAlignment', 'UNMAPPED=' + unmapped_bam_with_umi, 'ALIGNED=' + temp_mapped_bam, 'O=' + mapped_bam, 'R=' + fasta, 'SO=coordinate', 'ALIGNER_PROPER_PAIR_FLAGS=true', 'MAX_GAPS=-1', 'ORIENTATIONS=FR', 'VALIDATION_STRINGENCY=SILENT', 'CREATE_INDEX=true'])
		picard_mba.wait()



def error_correct_collapse(s_dict):
	for sample in s_dict:
		mapped_bam = sample + '_bwa.bam'
		mapped_bam_fixed_umi = sample + '_bwa_fixedumi_temp.bam'
		rejected_bam = sample + '_rejected_temp.bam'
		metrics = sample + '.correct_umi_metrics_temp.txt'
		grouped_bam = sample + '_grouped_temp.bam'
		consensus_unmapped_bam = sample + '_consensus_unmapped_temp.bam'
		temp_mapped_bam = sample + '_consensus_mapped_temp.bam'
		consensus_mapped_bam = sample + '_consensus_mapped_temp.bam'
		consensus_mapped_filtered_bam = sample + '_consensus_mapped_filtered_temp.bam'
		##correct umis
		#java -jar fgbio.jar CorrectUmis -i mapped.bam -o mapped.fixedumi.bam --max-mismatches=3 --min-distance=1 -M metrics.txt -r rejected.bam -t RX
		fgbio_cu = subprocess.Popen(['java', '-jar', fgbio, 'CorrectUmis', '-i', mapped_bam, '-o', mapped_bam_fixed_umi, '--max-mismatches=3', '--min-distance=1', '-M', metrics , '-r', rejected_bam, '-t', 'RX', '-u', 'GAGACGAT', 'TTCCAAGG', 'CGCATGAT', 'ACGGAACA', 'CGGCTAAT', 'GCTATCCT', 'TGGACTCT', 'ATCCAGAG', 'CTTAGGAC', 'GTGCCATA', 'TCGCTGTT', 'TTCGTTGG', 'AAGCACTG', 'GTCGAAGA', 'ACCACGAT', 'GATTACCG', 'GCACAACT', 'GCGTCATT', 'GAAGGAAG', 'ACTGAGGT', 'TGAAGACG', 'GTTACGCA', 'AGCGTGTT', 'GATCGAGT', 'TTGCGAAG', 'CTGTTGAC', 'GATGTGTG', 'ACGTTCAG', 'TTGCAGAC', 'CAATGTGG', 'ACGACTTG', 'ACTAGGAG'])
		fgbio_cu.wait()
		##Group read families 
		#java -jar fgbio.jar GroupReadsByUmi  --input=mapped.fixedumi.bam --output=grouped.bam  --strategy=paired --edits=0 --min-map-q=20
		fgbio_gr = subprocess.Popen(['java', '-jar', fgbio, 'GroupReadsByUmi', '-i', mapped_bam_fixed_umi, '-o', grouped_bam, '--strategy=paired', '--edits=0', '--min-map-q=20'])
		fgbio_gr.wait()
		##Collapse combined read families
		#java -jar fgbio.jar CallDuplexConsensusReads  --input=grouped.bam --output=consensus.unmapped.bam  --error-rate-pre-umi=45 --error-rate-post-umi=30  --min-input-base-quality=30
		fgbio_cdcr = subprocess.Popen(['java', '-jar', fgbio, 'CallDuplexConsensusReads', '-i', grouped_bam, '-o', consensus_unmapped_bam, '--error-rate-pre-umi=45', '--error-rate-post-umi=30', '--min-input-base-quality=30'])
		fgbio_cdcr.wait()
		##Align collapsed combined families
		#java -jar picard.jar SamToFastq I=consensus.unmapped.bam  F=/dev/stdout
		picard_sf = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'SamToFastq', 'I=' + consensus_unmapped_bam, 'F=/dev/stdout', 'INTERLEAVE=true'], stdout=subprocess.PIPE)
		#bwa mem -p -t 4 hg38.fa /dev/stdin
		with open(temp_mapped_bam, "w") as tmb_fh:
			bwa_int = subprocess.Popen([bwa, 'mem', '-p', '-t', '4', fasta, '/dev/stdin'], stdin=picard_sf.stdout, stdout=tmb_fh)
			bwa_int.wait()
		#java -jar picard.jar MergeBamAlignment UNMAPPED=unmapped.withUMI.bam ALIGNED=/dev/stdin O=consensus.mapped.bam R=hg38.fa \ SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
		picard_mba = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'MergeBamAlignment', 'UNMAPPED=' + consensus_unmapped_bam, 'ALIGNED=' + temp_mapped_bam, 'O=' + consensus_mapped_bam, 'R=' + fasta, 'SO=coordinate', 'ALIGNER_PROPER_PAIR_FLAGS=true', 'MAX_GAPS=-1', 'ORIENTATIONS=FR', 'VALIDATION_STRINGENCY=SILENT', 'CREATE_INDEX=true'])
		picard_mba.wait()
		##Filter collapsed combined families
		#java -jar fgbio.jar FilterConsensusReads --input=consensus.mapped.bam --output=consensus.mapped.filtered.bam --ref=hg38.fa --min-reads=2 1 1 --max-read-error-rate=0.05 --max-base-error-rate=0.1 --min-base-quality=50 --max-no-call-fraction=0.05 --require-single-strand-agreement=true
		fgbio_cdcr = subprocess.Popen(['java', '-jar', fgbio, 'FilterConsensusReads', '-i', consensus_mapped_bam, '-o', consensus_mapped_filtered_bam, '--ref=' + fasta, '--min-reads=2', '1', '1', '--max-read-error-rate=0.05', '--max-base-error-rate=0.1', '--min-base-quality=50', '--max-no-call-fraction=0.05', '--require-single-strand-agreement=true'])
		fgbio_cdcr.wait()		

def clip_overlap(s_dict, inbam_prefix, outbam_prefix):
	for sample in s_dict:
		consensus_mapped_filtered_bam = sample + inbam_prefix
		consensus_mapped_filtered_clipped_bam = sample + outbam_prefix
		#java -jar fgbio.jar ClipBam --input=consensus.mapped.filtered.bam --output=consensus.mapped.filtered.clipped.bam --ref=hg38.fa --soft-clip=false --clip-overlapping-reads=true
		fgbio_cb = subprocess.Popen(['java', '-jar', fgbio, 'ClipBam', '-i', consensus_mapped_filtered_bam, '-o', consensus_mapped_filtered_clipped_bam, '--ref=' + fasta, '--clipping-mode=Hard', '--clip-overlapping-reads=true'])
		fgbio_cb.wait()			

def calcluate_coverage(samples, bam_suffix, exome_bed, prefix):
	##parameters and header
	coverage_required = [10,50,100,200,500,1000]
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
	output_file = prefix + '.coverage_data.txt'
	print(output_file)
	with open(output_file, "w") as outfh:
		outfh.write(header)
		##go through bam files and get coverage histogram
		for sample in samples:
			bam = sample + bam_suffix
			print('calculating coverge for bam file', bam)
			##get temp coverage files
			coverage_histogram = sample + '.hist.temp'
			hist_fh = open(coverage_histogram, 'w')
			bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist'], stdout=subprocess.PIPE)
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

def calcluate_coverage_per_region(samples, bam_suffix, exome_bed, prefix):
	##parameters and header
	header = delim.join(['chr', 'start', 'end', 'target size', 'total bp', 'coverage']) + '\n'
	##open results file and write header
	
	##go through bam files and get coverage histogram
	for sample in samples:
		bam = sample + bam_suffix
		output_file = prefix + '.' + sample + '.coverage_data_per_region.txt'
		print('calculating coverge for bam file', bam)
		##get temp file for counts/bed region
		coverage_histogram = sample + '.hist.temp'
		with open(coverage_histogram, "w") as out_hist_fh:
			bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist'], stdout=out_hist_fh)
			bedtools_cov.wait()		
		##calculate numbers from histogram
		cov_dict = {}
		with open(coverage_histogram, "r") as cov_temp:
			for line in cov_temp:
				line = line.strip('\n').split(delim)
				if line[0] != 'all':
					region = '_'.join(line[:3])
					target_size = float(line[5])
					if line[3] != '0':
						seq_total = int(line[3]) *  int(line[4])
					else:
						seq_total = 0
					if region in cov_dict:
						cov_dict[region][0] += seq_total
					else:
						cov_dict[region] = [seq_total, target_size]
					# if line[1] == '115251155':
					# 	print(region, seq_total, cov_dict[region])
		with open(output_file, "w") as outfh:
			outfh.write(header)
			for r in cov_dict:
				t_size = cov_dict[r][1]
				s_size = cov_dict[r][0]
				coverage = s_size / float(t_size)
				line_out = r.split('_') + [str(t_size), str(s_size), str(coverage)]
				outfh.write(delim.join(line_out) + '\n')

def var_calling_vardict(samples, bed, bam_suffix, vcf_suffix):
	for sample in samples:
		bam = sample + bam_suffix
		out_vcf = sample + vcf_suffix
		with open(out_vcf, 'w') as out_fh:
			run_vardict = subprocess.Popen([vardict, '-G', fasta, '-f', vardict_maf_req, '-N', sample, '-b', bam, '-c', '1', '-S', '2', '-E', '3', '-th', '1', '--nosv', bed], stdout=subprocess.PIPE)
			run_teststrandbias = subprocess.Popen([teststrandbias], stdin=run_vardict.stdout, stdout=subprocess.PIPE)
			run_var2vcf = subprocess.Popen([var2vcf_valid, '-N', sample, '-E', '-f', vardict_maf_req], stdin=run_teststrandbias.stdout, stdout=out_fh)
			run_var2vcf.wait()

def annotate_vcf(samples, vcf_suffix, out_file):
	with open(out_file, "w") as out_fh:
		head = ['Sample', 'Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'DamagePredCount', 'SIFT_pred', 'SIFT4G_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'PROVEAN_pred', 'VEST4_score', 'MetaSVM_pred', 'MetaLR_pred', 'M-CAP_pred', 'REVEL_score', 'MutPred_score', 'MVP_score', 'MPC_score', 'PrimateAI_pred', 'DEOGEN2_pred', 'BayesDel_addAF_pred', 'BayesDel_noAF_pred', 'ClinPred_pred', 'LIST-S2_pred', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_pred', 'fathmm-XF_coding_pred', 'Eigen-raw_coding', 'Eigen-phred_coding', 'Eigen-PC-raw_coding', 'Eigen-PC-phred_coding', 'GenoCanyon_score', 'integrated_fitCons_score', 'GM12878_fitCons_score', 'H1-hESC_fitCons_score', 'HUVEC_fitCons_score', 'LINSIGHT', 'GERP++_NR', 'GERP++_RS', 'phyloP100way_vertebrate', 'phyloP30way_mammalian', 'phyloP17way_primate', 'phastCons100way_vertebrate', 'phastCons30way_mammalian', 'phastCons17way_primate', 'bStatistic', 'Interpro_domain', 'GTEx_V8_gene', 'GTEx_V8_tissue', 'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG',  'zygosity', 'quality', 'coverage', 'vcf_format', 'vcf_info']
		out_fh.write(delim.join(head) + '\n')
		for sample in samples:
			input_vcf = sample + vcf_suffix
			out_prefix = sample
			# '''
			##annotate vcf with annovar i.e. run_table_annovar
			command = [table_annovar] + av_buildver + [input_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
			annovar = subprocess.Popen(command)
			annovar.wait()
			# '''
			##combine all variants
			multianno = sample +'.hg19_multianno.txt'
			lc = 0
			with open(multianno, "r") as in_fh:
				for line in in_fh:
					line = line.split(delim)
					lc +=1
					if lc > 1:
						out_fh.write(delim.join([sample] + line[:103] + line[112:]))

def dedup_by_position(s_dict):
	for sample in s_dict:
		bwa_bam = sample + '_bwa.bam'
		mark_dup_bam = sample + '.method_a_mkdup_temp.bam'
		rm_dup_bam = sample + '.method_a_rmdup_temp.bam'
		metrics = sample + '.markduplicates_a.metrics.txt'
		#java -jar picard.jar MarkDuplicates I=mapped.bam O=markduplicates.bam M=markduplicates.metrics.txt
		##markdups
		picard_ms = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'MarkDuplicates', 'I=' + bwa_bam, 'O=' + mark_dup_bam, 'M=' + metrics])
		picard_ms.wait()
		##rm dups
		picard_ms = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'MarkDuplicates', 'I=' + bwa_bam, 'O=' + rm_dup_bam, 'M=' + metrics, 'REMOVE_DUPLICATES=true'])
		picard_ms.wait()

def dedup_by_umi(s_dict):
	for sample in s_dict:
		bwa_bam = sample + '_bwa.bam'
		mark_dup_bam = sample + '.method_b_mkdup_temp.bam'
		rm_dup_bam = sample + '.method_b_rmdup_temp.bam'
		metrics = sample + '.markduplicates_b.metrics.txt'
		#java -jar picard.jar MarkDuplicates I=mapped.bam O=markduplicates.bam M=markduplicates.metrics.txt
		##markdups
		picard_ms = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'MarkDuplicates', 'I=' + bwa_bam, 'O=' + mark_dup_bam, 'M=' + metrics, 'BARCODE_TAG=RX'])
		picard_ms.wait()
		##rm dups
		picard_ms = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'MarkDuplicates', 'I=' + bwa_bam, 'O=' + rm_dup_bam, 'M=' + metrics, 'REMOVE_DUPLICATES=true', 'BARCODE_TAG=RX'])
		picard_ms.wait()


def error_correct_collapse_srf(s_dict):
	for sample in s_dict:
		mapped_bam = sample + '_bwa.bam'
		mapped_bam_fixed_umi = sample + '_bwa_fixedumi_temp_c.bam'
		rejected_bam = sample + '_rejected_temp_c.bam'
		metrics = sample + '.correct_umi_metrics_c.txt'
		grouped_bam = sample + '_grouped_temp_c.bam'
		consensus_unmapped_bam = sample + '_consensus_unmapped_temp_c.bam'
		temp_mapped_bam = sample + '_consensus_mapped_temp_c.bam'
		consensus_mapped_bam = sample + '_consensus_mapped_temp_c.bam'
		consensus_mapped_filtered_bam = sample + '_consensus_mapped_filtered_temp_c.bam'
		##correct umis
		#java -jar fgbio.jar CorrectUmis -i mapped.bam -o mapped.fixedumi.bam --max-mismatches=3 --min-distance=1 -M metrics.txt -r rejected.bam -t RX
		fgbio_cu = subprocess.Popen(['java', '-jar', fgbio, 'CorrectUmis', '-i', mapped_bam, '-o', mapped_bam_fixed_umi, '--max-mismatches=3', '--min-distance=1', '-M', metrics , '-r', rejected_bam, '-t', 'RX', '-u', 'GAGACGAT', 'TTCCAAGG', 'CGCATGAT', 'ACGGAACA', 'CGGCTAAT', 'GCTATCCT', 'TGGACTCT', 'ATCCAGAG', 'CTTAGGAC', 'GTGCCATA', 'TCGCTGTT', 'TTCGTTGG', 'AAGCACTG', 'GTCGAAGA', 'ACCACGAT', 'GATTACCG', 'GCACAACT', 'GCGTCATT', 'GAAGGAAG', 'ACTGAGGT', 'TGAAGACG', 'GTTACGCA', 'AGCGTGTT', 'GATCGAGT', 'TTGCGAAG', 'CTGTTGAC', 'GATGTGTG', 'ACGTTCAG', 'TTGCAGAC', 'CAATGTGG', 'ACGACTTG', 'ACTAGGAG'])
		fgbio_cu.wait()
		##Group read families 
		#java -jar fgbio.jar GroupReadsByUmi  --input=mapped.fixedumi.bam --output=grouped.bam  --strategy=paired --edits=0 --min-map-q=20
		fgbio_gr = subprocess.Popen(['java', '-jar', fgbio, 'GroupReadsByUmi', '-i', mapped_bam_fixed_umi, '-o', grouped_bam, '--strategy=paired', '--edits=0', '--min-map-q=20'])
		fgbio_gr.wait()
		##Collapse single read families
		#java -jar fgbio.jar CallDuplexConsensusReads  --input=grouped.bam --output=consensus.unmapped.bam  --error-rate-pre-umi=45 --error-rate-post-umi=30  --min-input-base-quality=30
		fgbio_cdcr = subprocess.Popen(['java', '-jar', fgbio, 'CallMolecularConsensusReads', '-i', grouped_bam, '-o', consensus_unmapped_bam, '--error-rate-pre-umi=45', '--error-rate-post-umi=30', '--min-input-base-quality=30', '--min-consensus-base-quality=40', '--min-reads=1'])
		fgbio_cdcr.wait()
		##Align collapsed combined families
		#java -jar picard.jar SamToFastq I=consensus.unmapped.bam  F=/dev/stdout
		picard_sf = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'SamToFastq', 'I=' + consensus_unmapped_bam, 'F=/dev/stdout', 'INTERLEAVE=true'], stdout=subprocess.PIPE)
		#bwa mem -p -t 4 hg38.fa /dev/stdin
		with open(temp_mapped_bam, "w") as tmb_fh:
			bwa_int = subprocess.Popen([bwa, 'mem', '-p', '-t', '4', fasta, '/dev/stdin'], stdin=picard_sf.stdout, stdout=tmb_fh)
			bwa_int.wait()
		#java -jar picard.jar MergeBamAlignment UNMAPPED=unmapped.withUMI.bam ALIGNED=/dev/stdin O=consensus.mapped.bam R=hg38.fa \ SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
		picard_mba = subprocess.Popen(['java', '-Xmx4g', '-jar', picard, 'MergeBamAlignment', 'UNMAPPED=' + consensus_unmapped_bam, 'ALIGNED=' + temp_mapped_bam, 'O=' + consensus_mapped_bam, 'R=' + fasta, 'SO=coordinate', 'ALIGNER_PROPER_PAIR_FLAGS=true', 'MAX_GAPS=-1', 'ORIENTATIONS=FR', 'VALIDATION_STRINGENCY=SILENT', 'CREATE_INDEX=true'])
		picard_mba.wait()
		##Filter collapsed combined families
		fgbio_cdcr = subprocess.Popen(['java', '-jar', fgbio, 'FilterConsensusReads', '-i', consensus_mapped_bam, '-o', consensus_mapped_filtered_bam, '--ref=' + fasta, '--min-reads=3', '--max-read-error-rate=0.05', '--max-base-error-rate=0.1', '--min-base-quality=40', '--max-no-call-fraction=0.1'])
		fgbio_cdcr.wait()		


def mini_van_analysis_master(in_file, project_name, bed, hotspot_bed):
	sample_dict = get_sample_dict(in_file)
	# '''
	##make unmapped bam from fqs
	make_unmapped_bam(sample_dict)
	##add umi
	add_umi_unmapped_bam(sample_dict)
	##align reads
	align_reads(sample_dict)
	##Error correct by collapsing combined read families
	error_correct_collapse(sample_dict)
	##Clip overlap between read pairs
	clip_overlap(sample_dict)
	##calculate raw and collapsed coverage
	calcluate_coverage(sample_dict, '_bwa.bam', bed, project_name + '_all')
	calcluate_coverage(sample_dict, '_bwa_collapsed.bam', bed, project_name + '_collapsed')
	calcluate_coverage_per_region(sample_dict, '_bwa_collapsed.bam', bed, project_name + '_collapsed')
	calcluate_coverage(sample_dict, '_bwa_collapsed.bam', hotspot_bed, project_name + '_collapsed_hotspot')
	calcluate_coverage_per_region(sample_dict, '_bwa_collapsed.bam', hotspot_bed, project_name + '_collapsed_hotspot')
	# '''
	##variant calling
	var_calling_vardict(sample_dict, bed)

def minivan_analysis_master_v2(in_file, project_name, bed, hotspot_bed, error_correction_methods):
	sample_dict = get_sample_dict(in_file)
	print(sample_dict)
	# '''
	##for all samples
	##make unmapped bam from fqs
	make_unmapped_bam(sample_dict)
	##add umi
	add_umi_unmapped_bam(sample_dict)
	##align reads
	align_reads(sample_dict)
	##calculate raw coverage
	calcluate_coverage(sample_dict, '_bwa.bam', bed, project_name + '_all')
	# '''
	##error correction A
	if 'a' in error_correction_methods:
		# '''
		##Deduplicate by start-stop position
		dedup_by_position(sample_dict)
		##Clip overlap between read pairs
		clip_overlap(sample_dict, '.method_a_mkdup_temp.bam', '.method_a_mkdup.bam')
		clip_overlap(sample_dict, '.method_a_rmdup_temp.bam', '.method_a_rmdup_temp2.bam')
		##calculate dedup coverage
		calcluate_coverage(sample_dict, '.method_a_rmdup_temp2.bam', bed, project_name + '_method_a')
		calcluate_coverage_per_region(sample_dict, '.method_a_rmdup_temp2.bam', bed, project_name + '_method_a')
		calcluate_coverage(sample_dict, '.method_a_rmdup_temp2.bam', hotspot_bed, project_name + '_method_a_hotspot')
		calcluate_coverage_per_region(sample_dict, '.method_a_rmdup_temp2.bam', hotspot_bed, project_name + '_method_a_hotspot')
		# '''
		##variant calling
		var_calling_vardict(sample_dict, bed, '.method_a_mkdup.bam', '.method_a_mkdup.vardict.vcf')
		##annovar
		annotate_vcf(sample_dict, '.method_a_mkdup.vardict.vcf', project_name + '.method_a_mkdup.vardict.xls')
	##error correction B
	if 'b' in error_correction_methods:
		# '''
		##Deduplicate by umi
		dedup_by_umi(sample_dict)
		##Clip overlap between read pairs
		clip_overlap(sample_dict, '.method_b_mkdup_temp.bam', '.method_b_mkdup.bam')
		clip_overlap(sample_dict, '.method_b_rmdup_temp.bam', '.method_b_rmdup_temp2.bam')
		##calculate dedup coverage
		calcluate_coverage(sample_dict, '.method_b_rmdup_temp2.bam', bed, project_name + '_method_b')
		calcluate_coverage_per_region(sample_dict, '.method_b_rmdup_temp2.bam', bed, project_name + '_method_b')
		calcluate_coverage(sample_dict, '.method_b_rmdup_temp2.bam', hotspot_bed, project_name + '_method_b_hotspot')
		calcluate_coverage_per_region(sample_dict, '.method_b_rmdup_temp2.bam', hotspot_bed, project_name + '_method_b_hotspot')
		# '''		
		##variant calling
		var_calling_vardict(sample_dict, bed, '.method_b_mkdup.bam', '.method_b_mkdup.vardict.vcf')
		##annovar
		annotate_vcf(sample_dict, '.method_b_mkdup.vardict.vcf', project_name + '.method_b_mkdup.vardict.xls')

	##error correction C
	if 'c' in error_correction_methods:
		# '''
		##Error correct by collapsing single read families
		error_correct_collapse_srf(sample_dict)
		##Clip overlap between read pairs
		clip_overlap(sample_dict, '_consensus_mapped_filtered_temp_c.bam', '.method_c_collapsed.bam')
		##calculate collapsed coverage
		calcluate_coverage(sample_dict, '.method_c_collapsed.bam', bed, project_name + '_method_c')
		calcluate_coverage_per_region(sample_dict, '.method_c_collapsed.bam', bed, project_name + '_method_c')
		calcluate_coverage(sample_dict, '.method_c_collapsed.bam', hotspot_bed, project_name + '_method_c_hotspot')
		calcluate_coverage_per_region(sample_dict, '.method_c_collapsed.bam', hotspot_bed, project_name + '_method_c_hotspot')
		# '''
		##variant calling
		var_calling_vardict(sample_dict, bed, '.method_c_collapsed.bam', '.method_c_collapsed.vardict.vcf')
		##annovar
		annotate_vcf(sample_dict, '.method_c_collapsed.vardict.vcf', project_name + '.method_c_collapsed.vardict.xls')

	##error correction D
	if 'd' in error_correction_methods:
		# '''
		##Error correct by collapsing combined read families
		error_correct_collapse(sample_dict)
		##Clip overlap between read pairs
		clip_overlap(sample_dict, '_consensus_mapped_filtered_temp.bam', '.method_d_collapsed.bam')
		##calculate collapsed coverage
		calcluate_coverage(sample_dict, '.method_d_collapsed.bam', bed, project_name + '_method_d')
		calcluate_coverage_per_region(sample_dict, '.method_d_collapsed.bam', bed, project_name + '_method_d')
		calcluate_coverage(sample_dict, '.method_d_collapsed.bam', hotspot_bed, project_name + '_method_d_hotspot')
		calcluate_coverage_per_region(sample_dict, '.method_d_collapsed.bam', hotspot_bed, project_name + '_method_d_hotspot')
		# '''		
		##variant calling
		var_calling_vardict(sample_dict, bed, '.method_d_collapsed.bam', '.method_d_collapsed.vardict.vcf')
		##annovar
		annotate_vcf(sample_dict, '.method_d_collapsed.vardict.vcf', project_name + '.method_d_collapsed.vardict.xls')

##run methods
# working_dir = '/home/atimms/ngs_data/targetted/jimmy_minivan_0621'
# os.chdir(working_dir)

##params
# minvan_bed = 'minivan_hg19_nochr_0621.bed'
# hotspot_bed = 'minivan_hotspots_hg19_0621.bed'
minvan_bed = '/home/atimms/ngs_data/references/hg19/minivan_hg19_nochr_0621.bed'
hotspot_bed = '/home/atimms/ngs_data/references/hg19/minivan_hotspots_hg19_0621.bed'

##minivan test 062821
##3 columns, sample name, fq1, fq2
info_file = 'minivan_test_062821.txt'
project = info_file.split('.')[0]
# mini_van_analysis_master(info_file, project, minvan_bed, hotspot_bed)

##minivan test 070121
##3 columns, sample name, fq1, fq2
info_file = 'minivan_test_070121.txt'
project = info_file.split('.')[0]
# mini_van_analysis_master(info_file, project, minvan_bed, hotspot_bed)

##combine the two minivan tests 0721
##3 columns, sample name, fq1, fq2
info_file = 'minivan_test_0721.txt'
project = info_file.split('.')[0]
error_correction_wanted = ['a', 'b', 'c', 'd']
error_correction_wanted = ['c']
# minivan_analysis_master_v2(info_file, project, minvan_bed, hotspot_bed, error_correction_wanted)

##new run 0821
# working_dir = '/home/atimms/ngs_data/targetted/jimmy_minivan_0821'
# os.chdir(working_dir)
##3 columns, sample name, fq1, fq2
info_file = 'minivan_test_0821.txt'
project = info_file.split('.')[0]
error_correction_wanted = ['a', 'b', 'c', 'd']
# minivan_analysis_master_v2(info_file, project, minvan_bed, hotspot_bed, error_correction_wanted)

##new run 0821
# working_dir = '/home/atimms/ngs_data/targetted/jimmy_minivan_nextseq_0821'
# os.chdir(working_dir)
##3 columns, sample name, fq1, fq2
info_file = 'minivan_nextseq_0821.txt'
project = info_file.split('.')[0]
error_correction_wanted = ['a', 'b', 'c', 'd']
# minivan_analysis_master_v2(info_file, project, minvan_bed, hotspot_bed, error_correction_wanted)


##new run 0921
working_dir = '/home/atimms/ngs_data/targetted/jimmy_minivan_0921'
os.chdir(working_dir)
##3 columns, sample name, fq1, fq2
info_file = 'minivan_0921.txt'
project = info_file.split('.')[0]
error_correction_wanted = ['a', 'b', 'c', 'd']
minivan_analysis_master_v2(info_file, project, minvan_bed, hotspot_bed, error_correction_wanted)
