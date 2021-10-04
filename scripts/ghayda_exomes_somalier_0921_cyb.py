#!/usr/bin/env python
import os
import subprocess
import glob

'''
requires:
qsub -Iq cdbrmq -l mem=30gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878

'''

##set input variables and parameters
delim = '\t'

somalier = '/home/atimms/programs/somalier_0721/somalier'
somalier_sites_vcf = '/home/atimms/ngs_data/references/somalier/sites.GRCh37.vcf.gz'
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'


def run_somalier(vcfs, ped_file):
	for vcf in vcfs:
		#somalier extract -d extracted/ --sites sites.vcf.gz -f /data/human/g1k_v37_decoy.fa $cohort.vcf.gz
		som_extract = subprocess.Popen([somalier, 'extract', '-d', 'somalier_extract/', '--sites', somalier_sites_vcf, '-f', fasta, vcf])
		som_extract.wait()
	#somalier relate --ped $pedigree extracted/*.somalier
	# som_relate = subprocess.Popen([somalier, 'relate', '--ped', ped, 'somalier_extract/*.somalier'])
	som_relate = subprocess.Popen([somalier, 'relate', '--infer', '--ped', ped_file, 'somalier_extract/*.somalier'])
	som_relate.wait()


##somalier
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalyze_exomes_0121/ghayda_exomes_somalier_0921/'
os.chdir(working_dir)
combined_ped = 'combined_peds.ped'
vcf_files = glob.glob('*intersect.bcftools.GRCh37_87.vcf.gz')
run_somalier(vcf_files, combined_ped)







