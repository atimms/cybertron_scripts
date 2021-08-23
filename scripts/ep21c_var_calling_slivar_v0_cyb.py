#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
requires:
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878

#java 1.8 for gatk
module load java/1.8.0_202 


'''

##set input variables and parameters
delim = '\t'
threads = '5'

##programs
bgzip = '/home/atimms/programs/samtools-1.11/bin/bgzip'
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
slivar = '/home/atimms/programs/slivar_0121/slivar'
# slivar = '/home/atimms/programs/slivar_0121/slivar_dev'
pslivar = '/home/atimms/programs/slivar_0121/pslivar'
slivar_functions = '/home/atimms/programs/slivar_0121/slivar-0.2.1/js/slivar-functions_at.js'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
freebayes = '/home/atimms/programs/freebayes-1.3.4-linux-static-AMD64'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
tabix = '/home/atimms/programs/samtools-1.11/bin/tabix'
# plink = '/home/atimms/programs/plink'
somalier = '/home/atimms/programs/somalier_0721/somalier'

##files
fasta = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens_assembly38.fasta'
ens_gff = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens.GRCh38.103.gff3.gz'
gnomad_gnotate = '/home/atimms/ngs_data/references/slivar/gnomad.hg38.genomes.v3.fix.zip'
#wget -qO - https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | zcat | cut -f 1,21,24 | tail -n+2 | awk '{ printf("%s\tpLI=%.3g;oe_lof=%.5g\n", $1, $2, $3)}' > pli.lookup
pli_lookup = '/home/atimms/ngs_data/references/slivar/pli.lookup'
#grep -v '^#' omim_genemap2_0321.txt | cut -f9,13 | grep -v '^\s' > /home/atimms/ngs_data/references/slivar/omim.lookup
omim_lookup = '/home/atimms/ngs_data/references/slivar/omim.lookup'
dbsnp_common_vcf = '/home/atimms/ngs_data/references/hg38_gatk/common_all_20180418.chr_added.vcf.gz'
somalier_sites_vcf = '/home/atimms/ngs_data/references/somalier/sites.hg38.vcf.gz'

##ped types
trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
duo_types = ['duo', 'duo*', 'parent_sibship']
single_types = ['singleton', 'singleton*', 'sibship']
duo_single_types = duo_types + single_types
multiplex_types = ['multiplex']

##annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,knownGene,ensGene,dbnsfp41a,clinvar_20210123,hg38_gnomad30_genome']
av_operation = ['-operation', 'g,g,g,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput', '-arg', '-splicing 10 ,,,,,' ]

##methods




















##run methods
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_genedx_0821/'
exome_info_file = 'ghayda_genedx_0821.txt'
run_sarek_by_cohort(work_dir, exome_info_file)




