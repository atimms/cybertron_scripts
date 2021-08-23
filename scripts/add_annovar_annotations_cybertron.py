#!/usr/bin/python
import os
import subprocess
import glob
import shutil

##set input variables and parameters
delim = '\t'

##programs and ref files
vt = 'vt'
# gemini = '/home/atimms/programs/gemini/data/anaconda/bin/gemini'
gemini = '/home/atimms/scripts/gemini'
bcftools_12 = 'bcftools'
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
# dbdb_data = ref_dir + 'gemini/dbdb.gene_association.txt'
# mgi_data = ref_dir + 'gemini/mgi.abnormal_brain.txt'
# omim_data = ref_dir + 'gemini/omim_genemap2_0816.txt'
# previous_exomes_data = ref_dir + 'gemini/all_exome_data_std_pipeline.1217.xls'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
ann_var = '/home/atimms/programs/annovar/annotate_variation.pl'
prep_ann = '/home/atimms/programs/annovar/prepare_annovar_user.pl'

##annovar parameters
# av_genome = 'hg19'
av_genome = 'hg38'
# av_genome = 'danRer11'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,kaviar_20150923,gnomad_genome,popfreq_all_20150413']
av_operation = ['-operation', 'g,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput']



##methods
def download_annovar_db_from_annovar(file_to_download):
	annovar_dl = subprocess.Popen([ann_var, '-buildver', av_genome, '-downdb', file_to_download, av_ref_dir[0], '-webfrom', 'annovar'])
	annovar_dl.wait()

def download_annovar_db_not_from_annovar(file_to_download):
	annovar_dl = subprocess.Popen([ann_var, '-buildver', av_genome, '-downdb', file_to_download, av_ref_dir[0]])
	annovar_dl.wait()

def filter_vcf_by_region(invcf, outvcf, region_bed):
	# bgzip_vcf = subprocess.Popen(['bgzip', invcf])
	# bgzip_vcf.wait()
	# bcf_index = subprocess.Popen(['bcftools', 'index', invcf + '.gz'])
	# bcf_index.wait()
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-R', region_bed, '-o', outvcf, '-O', 'z', invcf])
	bcftools_filter.wait()

def annotate_vcf(vcf, av_prefix):
	##convert vcf file to individual avinput file
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', vcf, '-allsample', '--withfreq', '-outfile', av_prefix + '.avinput'])
	con_ann.wait()
	##run_table_annovar(avinput, av_prefix):
	command = [table_annovar] + av_buildver + [av_prefix + '.avinput'] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', av_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()


##run methods
##parameters
# working_dir = '/home/atimms/ngs_data/exomes/ghayda_gw_1117/add_annovar_LR07-201_0118'
# os.chdir(working_dir)
# region_file = 'LR07-201_vars.bed'
# std_vcf = 'LR07-201.int.norm.vcf.gz'
# filtered_vcf = 'LR07-201.region.vcf.gz'
# file_prefix = 'LR07-201'

##get vcf with vars of interest
# filter_vcf_by_region(std_vcf, filtered_vcf, region_file)
# annotate_vcf(filtered_vcf, file_prefix)


##download annovar databases tp be used for hg38
# ftds = ['refGene', 'avsnp150', 'gnomad_genome', 'kaviar_20150923']
# for ftd in ftds:
# 	download_annovar_db_from_annovar(ftd)
# ftds = ['rmsk', 'genomicSuperDups']
# for ftd in ftds:
# 	download_annovar_db_not_from_annovar(ftd)
##more dbs for hg38 genomes anlysis 0219
# ftds = ['dbnsfp35a', 'exac03', 'gnomad_exome', 'clinvar_20180603', 'popfreq_max_20150413']
# for ftd in ftds:
# 	download_annovar_db_from_annovar(ftd)
# ftds = ['rmsk', 'genomicSuperDups']
# for ftd in ftds:
# 	download_annovar_db_not_from_annovar(ftd)

##new hg19 annotation
# ftds = ['clinvar_20190305']
# for ftd in ftds:
# 	download_annovar_db_from_annovar(ftd)
# ftds = ['gnomad211_exome', 'gnomad211_genome']
# ftds = ['dbnsfp35c']
# ftds = ['dbnsfp31a_interpro']
# ftds = ['dbnsfp35a', 'gnomad211_exome', 'gnomad211_genome', 'mcap13', 'revel', 'mcap14']
# for ftd in ftds:
# 	download_annovar_db_from_annovar(ftd)
# ftds = ['gerp++gt2']
# for ftd in ftds:
# 	download_annovar_db_from_annovar(ftd)
##0121
# ftds = ['dbnsfp41a', 'clinvar_20200316']
# ftds = ['refGene']
# ftds = ['clinvar_20210123']
# for ftd in ftds:
# 	download_annovar_db_from_annovar(ftd)


##new hg38 files
# ftds = ['clinvar_20190305', 'dbnsfp31a_interpro']
# ftds = ['gnomad211_genome', 'gnomad211_exome', 'gnomad30_genome']
# ftds = ['clinvar_20200316']
# ftds = ['clinvar_20210123', 'dbnsfp41a']
ftds = ['knownGene', 'ensGene']
for ftd in ftds:
	download_annovar_db_from_annovar(ftd)

##make cosmic dbs - 0320
##download latest verion (90) i.e. 4 files from cosmic and de compress (make sure you have correct version)
##run commands - hg38
#/home/atimms/programs/annovar/prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.vcf > hg38_cosmic90_coding.txt
#/home/atimms/programs/annovar/prepare_annovar_user.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingVariants.vcf > hg38_cosmic90_noncoding.txt
##run commands - hg19
#/home/atimms/programs/annovar/prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.vcf > hg19_cosmic90_coding.txt
#/home/atimms/programs/annovar/prepare_annovar_user.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingVariants.vcf > hg19_cosmic90_noncoding.txt


##zebrafish annotations
# ftds = ['refGene', 'rmsk', 'ensGene']
# # ftds = ['rmsk', 'ensGene']
# ftds = ['ensGene']


# for ftd in ftds:
# 	download_annovar_db_not_from_annovar(ftd)






