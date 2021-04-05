#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

'''
set up:
interactive session i.e.
##for mapping/annotation
qsub -Iq cdbrmq -l mem=200gb,ncpus=20 -P 19833a08-f6fb-4bea-8526-8a79069da878
qsub -Iq cdbrmq -l mem=20gb,ncpus=1 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load java/1.8.0_202 
module load mono/5.10.1.47

##for graphing
conda activate pd_np_plt_etc

'''
##parameters
delim = '\t'

##files and dirs
fasta = '/home/atimms/ngs_data/references/hg38/hg38.fa'
pisces_fastq_dir = '/home/atimms/ngs_data/references/illumina_hg38/'

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
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'rmsk,genomicSuperDups,gnomad211_genome,avsnp150,generic', '-genericdbfile', 'dsvm.random.all_scores.avinput']
av_operation = ['-operation', 'r,r,f,f,f']
# av_options = ['-otherinfo', '-remove', '-nastring', '.']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']



##methods

##call freebayes in all sample listed in list of bamfiles
def variant_calling_freebayes(samples, bam_suffix, dvsm_bed):
	for sample in samples:
		vcf_temp1 = sample + 'temp_fb1.vcf'
		vcf_temp2 = sample + 'temp_fb2.vcf.gz'
		final_vcf = sample + '.fb.cov20_q20.vcf.gz'
		bam = sample + bam_suffix
		##call vars
		freebayes_run = subprocess.Popen([freebayes,'-f', fasta, '-t', dvsm_bed, bam, '-v', vcf_temp1])
		freebayes_run.wait()
		bgzip_cmd = subprocess.Popen([bgzip, vcf_temp1])
		bgzip_cmd.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
		bcf_index.wait()	
		##filter
		bcftools_filter = subprocess.Popen([bcftools, 'view', '-i', "QUAL>20 & MIN(FORMAT/DP)>20", '-V', 'indels','--threads', '20', '-o', vcf_temp2, '-O', 'z', vcf_temp1 + '.gz'])
		bcftools_filter.wait()
		##decomposing and subsetting vcf
		zless_vcf = subprocess.Popen(['zcat', vcf_temp2], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
		bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '--threads', '20', '-o', final_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
		bcftools_norm.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
		bcf_index.wait()	

def variant_calling_gatk(samples, bam_suffix, dvsm_bed):
	for sample in samples:
		vcf_temp0 = sample + 'temp_gatk0.vcf'
		vcf_temp1 = sample + 'temp_gatk1.vcf'
		vcf_temp2 = sample + 'temp_gatk2.vcf'
		vcf_temp3 = sample + 'temp_gatk3.vcf.gz'
		final_vcf = sample + '.gatk.cov20_q20.vcf.gz'
		bam = sample + bam_suffix
		##call vars
		##run haplotype caller
		gatk_hc = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', bam,'-L', dvsm_bed, '-o', vcf_temp0])
		gatk_hc.wait()
		##split data in SNPs and indels and apply manual variant filtering (add -L)
		snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-nt', '5', '-V', vcf_temp0, '-o', vcf_temp1, '-selectType', 'SNP'])
		snp_cut.wait()
		snp_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_temp1, '-o', vcf_temp2, '--filterExpression ', "QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "snp_filter"])
		snp_vf.wait()
		bgzip_cmd = subprocess.Popen([bgzip, vcf_temp2])
		bgzip_cmd.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp2 + '.gz'])
		bcf_index.wait()	
		##filter
		bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>20 & MIN(FORMAT/DP)>20", '-V', 'indels', '-o', vcf_temp3, '-O', 'z', vcf_temp2 + '.gz'])
		bcftools_filter.wait()
		##decomposing and subsetting vcf
		zless_vcf = subprocess.Popen(['zcat', vcf_temp3], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
		bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta,'--threads', '20', '-o', final_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
		bcftools_norm.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
		bcf_index.wait()	

def variant_calling_samtools(samples, bam_suffix, dvsm_bed):
	for sample in samples:
		vcf_temp1 = sample + 'temp_st1.vcf.gz'
		vcf_temp2 = sample + 'temp_st2.vcf.gz'
		final_vcf = sample + '.st.cov20_q20.vcf.gz'
		bam = sample + bam_suffix
		##call vars
		stmp = subprocess.Popen([bcftools,'mpileup', '-Ou', '-q', '20', '-C', '50', '-a', 'FORMAT/AD,FORMAT/DP', '--threads', '10', '-T', dvsm_bed, '-f', fasta, bam], stdout=subprocess.PIPE)
		bcft = subprocess.Popen([bcftools,'call', '-vmO', 'z', '--threads', '10', '-V', 'indels', '-o', vcf_temp1], stdin=stmp.stdout)
		bcft.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1])
		bcf_index.wait()	
		##filter
		bcftools_filter = subprocess.Popen([bcftools, 'view', '-i', "QUAL>20 & MIN(FORMAT/DP)>20", '-V', 'indels','--threads', '20', '-o', vcf_temp2, '-O', 'z', vcf_temp1])
		bcftools_filter.wait()
		##decomposing and subsetting vcf
		zless_vcf = subprocess.Popen(['zcat', vcf_temp2], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
		bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta,'--threads', '20', '-o', final_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
		bcftools_norm.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
		bcf_index.wait()	

def variant_calling_pisces(samples, bam_suffix, dvsm_bed):
	for sample in samples:
		out_dir = sample + '_pisces'
		bam = sample + bam_suffix
		pisces_vcf = out_dir + '/' + sample + '.bwa_mkdup.vcf'
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
		# '''
		run_pisces = subprocess.Popen(['mono', pisces, '-G', pisces_fastq_dir, '-BamPaths', bam, '-i', dvsm_bed, '-t', '18', '-OutFolder', out_dir])
		run_pisces.wait()
		bgzip_cmd = subprocess.Popen([bgzip, pisces_vcf])
		bgzip_cmd.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', pisces_vcf + '.gz'])
		bcf_index.wait()	
		##filter
		bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>20 & MIN(FORMAT/DP)>20", '-V', 'indels','--threads', '20', '-o', vcf_temp1, '-O', 'z', pisces_vcf + '.gz'])
		bcftools_filter.wait()
		##decomposing and subsetting vcf
		zless_vcf = subprocess.Popen(['zcat', vcf_temp1], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
		bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta,'--threads', '20', '-o', final_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
		bcftools_norm.wait()
		bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
		bcf_index.wait()	

def annotate_vcfs(samples, callers):
	for sample in samples:
		for caller in callers:
			in_vcf = sample + '.' + caller + '.cov20_q20.vcf.gz'
			out_prefix = sample + '.' + caller
			##annotate vcfs with annovar i.e. run_table_annovar
			command = [table_annovar] + av_buildver + [in_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
			annovar = subprocess.Popen(command)
			annovar.wait()

def filter_make_files_for_graphing_and_graph(samples, callers, out_prefix):
	for caller in callers:
		out_file = out_prefix + caller + '.for_scatter.txt'
		ab50_pdf = 'ab50_vs_dsvm.' + caller + '.for_scatter.pdf'
		with open(out_file, 'w') as out_fh:
			header = ['var', 'dSVM', 'AB', 'AB50', 'sample']
			out_fh.write(delim.join(header) + '\n')
			for sample in samples:
				multianno = sample + '.' + caller + '.hg38_multianno.txt'
				with open(multianno, 'r') as in_fh:
					lc = 0
					for line in in_fh:
						line = line.rstrip().split(delim)
						lc += 1
						if lc > 1:
							var = '_'.join([sample] + line[:5])
							rmsk = line[5]
							sedup = line[6]
							dsvm = line[25]
							gt = line[26]
							dbsnp = line[24]
							cov = int(line[28])
							##remove if in repeat no dvsm score, not in dbsnp or homozygote
							if rmsk == '.' and sedup == '.' and dsvm != '.' and gt == '0.5' and dbsnp != '.':
							# if rmsk == '.' and sedup == '.' and dsvm != '.' and gt == '0.5' and dbsnp != '.' and cov >= 40:
								if caller == 'fb' or caller == 'pisces':
									alleles = line[38].split(':')[2].split(',')
								elif caller == 'gatk':
									alleles = line[38].split(':')[1].split(',')
								elif caller == 'st':
									alleles = line[38].split(':')[3].split(',')
								alleles = [float(i) for i in alleles]
								ab = alleles[1] / sum(alleles)
								##remove the real homozoygous vars
								if ab < 0.95 and ab > 0.05:
									if ab <= 0.5:
										ab50 = ab
									else:
										ab50 = 1 - ab
									if len(alleles) != 2:
										print('multiallelic: ', line)
									line_out = [var, dsvm, str(ab), str(ab50), sample]
									out_fh.write(delim.join(line_out) + '\n')

		##make dict into a panda's dataframe and get 
		int_data = pd.read_table(out_file, index_col=0 )
		sns.set_style('whitegrid')
		correlation = round(int_data['AB'].corr(int_data['dSVM']),3)
		sns.regplot(x = int_data['AB'], y = int_data['dSVM'],  ci = None, line_kws = {'color': 'red'})
		plt.title("Pearson's r = " + str(correlation))
		# print(round(int_data['AB'].corr(int_data['dSVM']),3))
		pdf_name = out_file.rsplit('.', 1)[0] + '.pdf' 
		plt.savefig(pdf_name)
		plt.close()
		sns.set_style('whitegrid')
		correlation = round(int_data['AB50'].corr(int_data['dSVM']),3)
		sns.regplot(x = int_data['AB50'], y = int_data['dSVM'],  ci = None, line_kws = {'color': 'red'})
		plt.title("Pearson's r = " + str(correlation))
		# print(round(int_data['AB'].corr(int_data['dSVM']),3))
		plt.savefig(ab50_pdf)
		plt.close()


def make_box_plots(callers, in_prefix):
	for caller in callers:
		in_file = in_prefix + caller + '.for_scatter.txt'
		boxplot_file = in_prefix + caller + '.for_boxplot.txt'
		boxplot_pdf = in_prefix + caller + '.boxplot.pdf'
		violinplot_pdf = in_prefix + caller + '.violinplot.pdf'
		with open(in_file, 'r') as in_fh, open(boxplot_file, 'w') as out_fh:
			header = ['var', 'dSVM', 'AB', 'AB50', 'sample', 'dSVM_type']
			out_fh.write(delim.join(header) + '\n')
			lc = 0
			for line in in_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc > 1:
					dsvm = float(line[1])
					ab = float(line[2])
					ab50 = float(line[3])
					if dsvm < -4:
						line_out = line + ['negative_dsvm']
						out_fh.write(delim.join(line_out) + '\n')
					elif dsvm > 4:
						line_out = line + ['positive_dsvm']
						out_fh.write(delim.join(line_out) + '\n')
					elif dsvm > -0.05 and dsvm < 0.05:
						line_out = line + ['neutral_dsvm']
						out_fh.write(delim.join(line_out) + '\n')
		##make boxplot and violin plot
		int_data = pd.read_table(boxplot_file, index_col=0 )
		sns.boxplot( x=int_data["dSVM_type"], y=int_data["AB50"], order = ['negative_dsvm', 'neutral_dsvm', 'positive_dsvm'] )
		plt.savefig(boxplot_pdf)
		plt.close()
		sns.violinplot( x=int_data["dSVM_type"], y=int_data["AB50"], order = ['negative_dsvm', 'neutral_dsvm', 'positive_dsvm'] )
		plt.savefig(violinplot_pdf)
		plt.close()
		##calculate ttest
		cat1 = int_data[int_data['dSVM_type']=='negative_dsvm']
		cat2 = int_data[int_data['dSVM_type']=='neutral_dsvm']
		cat3 = int_data[int_data['dSVM_type']=='positive_dsvm']
		tt_neg_neu = ttest_ind(cat1['AB50'], cat2['AB50'], equal_var=False)
		print(boxplot_file, tt_neg_neu)






##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_vaf_dSVM_comparison_0221'
os.chdir(working_dir)

sample_names = ['Hu1_ret_ATAC', 'Hu2_ret_ATAC', 'Hu3_ret_ATAC', 'Hu4_ret_ATAC', 'Hu5_ret_ATAC', 'Hu6_ret_ATAC', 'Hu7_ret_ATAC', 'Hu8_ret_ATAC']
# sample_names = ['Hu1_ret_ATAC']
bwa_bam_suffix = '.bwa_mkdup.bam'
##dsvm scores
dsvm_leah = 'dsvm_simplified.txt'
dsvm_random_all_bed = 'dsvm.random.all_scores.chr_added_merged.bed'
##var caller names
var_callers = ['fb', 'gatk', 'pisces', 'st']
project_name = 'ab_vs_dsvm.'


##variant calling
# variant_calling_freebayes(sample_names, bwa_bam_suffix, dsvm_random_all_bed)
# variant_calling_gatk(sample_names, bwa_bam_suffix, dsvm_random_all_bed)
# variant_calling_pisces(sample_names, bwa_bam_suffix, dsvm_random_all_bed)
# variant_calling_samtools(sample_names, bwa_bam_suffix, dsvm_random_all_bed)

##annotate with annovar: dbsnp, gnomad, rpts, dvsm
##make avinput: awk '{print $1"\t"$2+1"\t"$2+1"\t"$4"\t"$5"\t"$8;}' dsvm_simplified.txt > dsvm.random.all_scores.avinput
# annotate_vcfs(sample_names, var_callers)

##filter vars and make files for graphing, then graph
# sample_names = ['Hu1_ret_ATAC', 'Hu2_ret_ATAC', 'Hu3_ret_ATAC', 'Hu4_ret_ATAC']
# var_callers = ['fb']
filter_make_files_for_graphing_and_graph(sample_names, var_callers, project_name)

##make box plot for high low and neutral dsvm
make_box_plots(var_callers, project_name)

