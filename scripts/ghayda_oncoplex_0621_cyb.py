#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
qsub -Iq cdbrmq -l mem=120gb,ncpus=10 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load java/1.8.0_202 
module load mono/5.10.1.47
'''

##set input variables and parameters
delim = '\t'

##programs
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
freebayes = '/home/atimms/programs/freebayes-1.3.4-linux-static-AMD64'
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
bgzip = '/home/atimms/programs/samtools-1.11/bin/bgzip'
samtools = '/home/atimms/programs/samtools-1.11/bin/samtools'
pisces = '/cm/shared/apps/Pisces/5.1.6.54/Pisces.exe'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,gnomad211_exome,gnomad211_genome,dbnsfp41a,clinvar_20210123,generic', '-genericdbfile', 'Geno2MP.avinput']
av_operation = ['-operation', 'g,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput' ]

##files
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
pisces_fastq_dir = '/home/atimms/ngs_data/references/illumina_human_g1k_v37/'
onco_bed = '/home/atimms/ngs_data/references/hg19/OncoPlex_v6_112018_no_chr_MSI.bed'

##methods
def variant_calling_pisces(bam_dict, bed):
	for sample in bam_dict:
		out_dir = sample + '_pisces'
		bam = bam_dict[sample]
		pisces_vcf = out_dir + '/' + bam.rsplit('.', 1)[0] + '.vcf'
		vcf_temp1 = sample + 'temp_pis1.vcf.gz'
		final_vcf = sample + '.pisces.cov20_q20.vcf.gz'
		##check bai file exists and make if it isn't there
		if os.path.isfile(bam + '.bai'):
			print('bam %s alreaded indexed'%bam)
		else:
			print('indexing bam file:', bam)
			st_index = subprocess.Popen([samtools, 'index', bam])
			st_index.wait()
		##run pisces on all samples
		run_pisces = subprocess.Popen(['mono', pisces, '-G', pisces_fastq_dir, '-BamPaths', bam, '-i', bed, '-t', '10', '-OutFolder', out_dir])
		run_pisces.wait()
		bgzip_cmd = subprocess.Popen([bgzip, pisces_vcf])
		bgzip_cmd.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', pisces_vcf + '.gz'])
		bcf_index.wait()	
		##filter
		bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>20 & MIN(FORMAT/DP)>20", '-V', 'indels','--threads', '10', '-o', vcf_temp1, '-O', 'z', pisces_vcf + '.gz'])
		bcftools_filter.wait()
		##decomposing and subsetting vcf
		zless_vcf = subprocess.Popen(['zcat', vcf_temp1], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
		bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta,'--threads', '20', '-o', final_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
		bcftools_norm.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
		bcf_index.wait()

def annotate_vcf(samples, combined_exonic):
	with open(combined_exonic, "w") as out_fh:
		head = ['Sample', 'Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'DamagePredCount', 'SIFT_pred', 'SIFT4G_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'PROVEAN_pred', 'VEST4_score', 'MetaSVM_pred', 'MetaLR_pred', 'M-CAP_pred', 'REVEL_score', 'MutPred_score', 'MVP_score', 'MPC_score', 'PrimateAI_pred', 'DEOGEN2_pred', 'BayesDel_addAF_pred', 'BayesDel_noAF_pred', 'ClinPred_pred', 'LIST-S2_pred', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_pred', 'fathmm-XF_coding_pred', 'Eigen-raw_coding', 'Eigen-phred_coding', 'Eigen-PC-raw_coding', 'Eigen-PC-phred_coding', 'GenoCanyon_score', 'integrated_fitCons_score', 'GM12878_fitCons_score', 'H1-hESC_fitCons_score', 'HUVEC_fitCons_score', 'LINSIGHT', 'GERP++_NR', 'GERP++_RS', 'phyloP100way_vertebrate', 'phyloP30way_mammalian', 'phyloP17way_primate', 'phastCons100way_vertebrate', 'phastCons30way_mammalian', 'phastCons17way_primate', 'bStatistic', 'Interpro_domain', 'GTEx_V8_gene', 'GTEx_V8_tissue', 'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG', 'Geno2MP',  'zygosity', 'quality', 'coverage', 'vcf_format', 'vcf_info']
		out_fh.write(delim.join(head) + '\n')
		for sample in samples:
			input_vcf = sample + '.pisces.cov20_q20.vcf.gz'
			out_prefix = sample
			'''
			##annotate vcf with annovar i.e. run_table_annovar
			command = [table_annovar] + av_buildver + [input_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
			annovar = subprocess.Popen(command)
			annovar.wait()
			# '''
			##combine all exonic variants
			multianno = sample +'.hg19_multianno.txt'
			lc = 0
			with open(multianno, "r") as in_fh:
				for line in in_fh:
					line = line.split(delim)
					lc +=1
					if lc >= 1:
						func_refgne = line[5]
						ex_func = line[8]
						af_max = line[11]
						if af_max == '.':
							af_max = '0'
						if func_refgne == 'exonic' or func_refgne == 'splicing':
							if ex_func != 'synonymous SNV':
								if float(af_max) < 0.01:
									out_fh.write(delim.join([sample] + line[:104] + line[112:]))


##run methods
working_dir = '/home/atimms/ngs_data/targetted/ghayda_oncoplex_0621'
os.chdir(working_dir)
##first set
sample_dict = {'LR19-241_1906010': '51531_G04_OPXv6_ND0461.final.bam', 'LR11-366_1903019': '51532_H04_OPXv6_ND0461.final.bam', 
	'LR11-366_1903024': '51534_B05_OPXv6_ND0461.final.bam', 'LR11-366_1903025': '51535_C05_OPXv6_ND0461.final.bam', 
	'LR15-251_1509030': '51744_D06_OPXv6_NA0551.final.bam', 'LR16-242_1606101': '51745_E06_OPXv6_NA0551.final.bam', 
	'LR16-413_1705105': '51746_F06_OPXv6_NA0551.final.bam', 'LR16-470_1708086': '51747_G06_OPXv6_NA0551.final.bam', 
	'LR18-024_1801071': '51748_H06_OPXv6_NA0551.final.bam', 'LR18-024_1802009': '51749_A07_OPXv6_NA0551.final.bam', 
	'LR16-313_s1912049': '51750_B07_OPXv6_NA0551.final.bam', 'LR20-019_s2002001': '51751_C07_OPXv6_NA0551.final.bam', 
	'LR20-019_s2002002': '51752_D07_OPXv6_NA0551.final.bam', 'LR20-033_S2002014': '51753_E07_OPXv6_NA0551.final.bam', 
	'LR20-010a1_S2002028': '51754_F07_OPXv6_NA0551.final.bam', 'LR13-129_S2010021': '51755_G07_OPXv6_NA0551.final.bam', 
	'LR13-129_S2010022': '51756_H07_OPXv6_NA0551.final.bam', 'LR13-129_S2010024': '51757_A08_OPXv6_NA0551.final.bam', 
	'LR21-012_S2101015': '51758_B08_OPXv6_NA0551.final.bam', 'LR17-514_S2104017': '51759_C08_OPXv6_NA0551.final.bam'}
##second set
sample_dict_2 = {'LR12-317_1209034': '51737_E05_OPXv6_NA0551.final.bam', 'LR13-129_1305080': '51738_F05_OPXv6_NA0551.final.bam', 
	'LR12-246_1401036': '51739_G05_OPXv6_NA0551.final.bam', 'LR14-046_1402058': '51740_H05_OPXv6_NA0551.final.bam', 
	'LR14-155_1405044': '51741_A06_OPXv6_NA0551.final.bam', 'LR14-155_1405045': '51742_B06_OPXv6_NA0551.final.bam', 
	'LR15-251_1509029': '51743_C06_OPXv6_NA0551.final.bam'}


##variant calling
# variant_calling_pisces(sample_dict, onco_bed)
# variant_calling_pisces(sample_dict_2, onco_bed)

##combine the 2 dictionaries
sample_dict_2.update(sample_dict)

##annotation and combine into single fils
final_results = 'oncoplex_061521.exonic_rare_cov20_q20.xls'
annotate_vcf(sample_dict_2, final_results)



