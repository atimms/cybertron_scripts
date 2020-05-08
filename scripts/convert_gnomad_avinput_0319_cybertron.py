#!/usr/bin/env python
import subprocess
import os

##note
'''
module load local_python/3.6.4
module load biobuilds/2017.11
'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/references/gnomad_data'
os.chdir(working_dir)



def convert_vcf_annovar_with_freq_info(in_vcf, out_file):
	##bcftools query -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%pab_max\n' gnomad.exomes.r2.1.sites.chr22.vcf.gz
	bt_query = subprocess.Popen(['bcftools','query', '-f', '%CHROM\t%POS\t%END\t%REF\t%ALT\t%AC\t%AN\t%AF\t%nhomalt\t%non_neuro_AC\t%non_neuro_AN\t%non_neuro_AF\t%non_neuro_nhomalt\n', '-o', out_file, in_vcf])
	bt_query.wait()


##get files using gsutil

gnomad_genome_vcf = 'gnomad.genomes.r2.1.sites.vcf.gz'
gnomad_genome_avinput = 'gnomad.genomes.r2.1.avinput'
gnomad_exome_vcf = 'gnomad.exomes.r2.1.1.sites.vcf.gz'
gnomad_exome_avinput = 'gnomad.exomes.r2.1.1.avinput'
genomad_test_vcf = 'gnomad.exomes.r2.1.sites.chr22.vcf.gz'
genomad_test_avinput = 'test.avinput'
gnomad_genome_exome_int_vcf = 'gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.gz'
gnomad_genome_exome_int_avinput = 'gnomad.genomes.r2.1.1.exome_calling_intervals.sites.avinput'

# convert_vcf_annovar_with_freq_info(gnomad_exome_vcf, gnomad_exome_avinput)
# convert_vcf_annovar_with_freq_info(gnomad_genome_vcf, gnomad_genome_avinput)
convert_vcf_annovar_with_freq_info(gnomad_genome_exome_int_vcf, gnomad_genome_exome_int_avinput)
# convert_vcf_annovar_with_freq_info(genomad_test_vcf, genomad_test_avinput)