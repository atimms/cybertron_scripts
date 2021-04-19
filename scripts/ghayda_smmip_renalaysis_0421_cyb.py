#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
qsub -Iq cdbrmq -l mem=30gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
'''

##set input variables and parameters
delim = '\t'

##programs
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
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

##methods

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


def proccess_vcfs_etc(input_vcfs):
	for input_vcf in input_vcfs:
		out_prefix = input_vcf.split('.')[0]
		temp1_vcf = out_prefix + 'temp1.vcf.gz'
		# '''
		##decomposing and subsetting vcf
		# zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
		zless_vcf = subprocess.Popen(['cat', input_vcf], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
		bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', temp1_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
		bcftools_norm.wait()
		##annotate vcf with annovar i.e. run_table_annovar
		command = [table_annovar] + av_buildver + [temp1_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
		annovar = subprocess.Popen(command)
		annovar.wait()
		##split into individual samples
		post_annovar_vcf = out_prefix + '.hg19_multianno.vcf'
		con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', post_annovar_vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', out_prefix])
		con_ann.wait()
		# '''
		avinputs = glob.glob(out_prefix + '.L*put')
		print(avinputs)
		# '''
		##format files into annotated
		ann_files = []
		for avinput in avinputs:
			outfile = avinput.rsplit('.',1)[0] + '.annotated.xls'
			ann_files.append(outfile)
			head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Zygosity', 'Qual', 'Coverage', 'Format', 'Info', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'DamagePredCount', 'SIFT_pred', 'SIFT4G_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'PROVEAN_pred', 'VEST4_score', 'MetaSVM_pred', 'MetaLR_pred', 'M-CAP_pred', 'REVEL_score', 'MutPred_score', 'MVP_score', 'MPC_score', 'PrimateAI_pred', 'DEOGEN2_pred', 'BayesDel_addAF_pred', 'BayesDel_noAF_pred', 'ClinPred_pred', 'LIST-S2_pred', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_pred', 'fathmm-XF_coding_pred', 'Eigen-raw_coding', 'Eigen-phred_coding', 'Eigen-PC-raw_coding', 'Eigen-PC-phred_coding', 'GenoCanyon_score', 'integrated_fitCons_score', 'GM12878_fitCons_score', 'H1-hESC_fitCons_score', 'HUVEC_fitCons_score', 'LINSIGHT', 'GERP++_NR', 'GERP++_RS', 'phyloP100way_vertebrate', 'phyloP30way_mammalian', 'phyloP17way_primate', 'phastCons100way_vertebrate', 'phastCons30way_mammalian', 'phastCons17way_primate', 'bStatistic', 'Interpro_domain', 'GTEx_V8_gene', 'GTEx_V8_tissue', 'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG', 'Geno2MP']
			head_out = delim.join(head + ['\n'])
			with open(avinput, "r") as av, open(outfile, "w") as final:
				final.write(head_out)
				for line in av:
					line = line.strip('\n').split(delim)
					stuff = line[0:8] + line[16:]
					info = line[15].split(';')
					# print(info)
					info_list = split_info_field(info)
					line_out = delim.join(stuff + info_list +['\n'])
					final.write(line_out)
		##combine all pathogenic variants
		combined_pathogenic = 'combined.' + out_prefix + '.pathogenic.xls'
		with open(combined_pathogenic, "w") as out_fh:
			fc = 0
			for ann_file in ann_files:
				sample = ann_file.split('.')[1]
				fc += 1
				lc = 0
				with open(ann_file, "r") as in_fh:
					for line in in_fh:
						line = line.split(delim)
						lc +=1
						if lc == 1:
							if fc == 1:
								out_fh.write(delim.join(['sample'] + line))
						else:
							clinsig = line[104]
							if 'Pathogenic' in clinsig or 'Likely_pathogenic' in clinsig:
								out_fh.write(delim.join([sample] + line))

		# '''

##run methods
working_dir = '/home/atimms/ngs_data/targetted/ghayda_smmip_renalaysis_0421'
os.chdir(working_dir)
vcfs = ['mcd16_0517_collapsed_freebayes.vcf', 'mcd16_0517_collapsed_gatkUG.vcf', 'mcd16_0517_total_freebayes.vcf', 
	'mcd16_p3_0517_collapsed_freebayes.vcf', 'mcd16_p3_0517_collapsed_gatkUG.vcf', 'mcd16_p3_0517_total_freebayes.vcf', 
	'mcd17_0518.combined_pisces.vcf', 'meg_16_0217.combined_pisces.vcf', 'meg16_p1r_collapsed_freebayes.vcf', 
	'meg16_p1r_collapsed_gatkUG.vcf', 'meg16_p1r_total_freebayes.vcf', 'meg16_p2_collapsed_freebayes.vcf', 
	'meg16_p2_collapsed_gatkUG.vcf', 'meg16_p2_total_freebayes.vcf', 'meg16_p3_collapsed_freebayes.vcf', 
	'meg16_p3_collapsed_gatkUG.vcf', 'meg16_p3_total_freebayes.vcf', 'meg16_p4_collapsed_freebayes.vcf', 
	'meg16_p4_collapsed_gatkUG.vcf', 'meg16_p4_total_freebayes.vcf']

# vcfs = ['mcd16_0517_collapsed_freebayes.vcf']
proccess_vcfs_etc(vcfs)





