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
# ens_gff = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens.GRCh38.103.gff3.gz'
ens_gff = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens.GRCh38.103.chr_added.gff3.gz'
gnomad_gnotate = '/home/atimms/ngs_data/references/slivar/gnomad.hg38.genomes.v3.fix.zip'
#wget -qO - https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | zcat | cut -f 1,21,24 | tail -n+2 | awk '{ printf("%s\tpLI=%.3g;oe_lof=%.5g\n", $1, $2, $3)}' > pli.lookup
pli_lookup = '/home/atimms/ngs_data/references/slivar/pli.lookup'
#grep -v '^#' omim_genemap2_0321.txt | cut -f9,13 | grep -v '^\s' > /home/atimms/ngs_data/references/slivar/omim.lookup
omim_lookup = '/home/atimms/ngs_data/references/slivar/omim.lookup'
#cut -f2,3,4,6,8 dbdb.gene_association.txt | awk 'BEGIN{FS="\t"}{ printf("%s\tinheritance=%s;phenotype=%s;syndrome=%s;loe=%s\n", $1, $2, $3, $4, $5)}' > ~/ngs_data/references/slivar/dbdb.lookup
dbdb_lookup = '/home/atimms/ngs_data/references/slivar/dbdb.lookup'
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
av_protocol = ['-protocol', 'refGene,dbnsfp41a,clinvar_20210123,gnomad30_genome']
av_operation = ['-operation', 'g,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput' ]

##methods

def combine_snp_indel_vcfs(vcf_filtered_snps, vcf_filtered_indels, comb_vcf):
	##index in vcfs
	tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', vcf_filtered_snps])
	tabix_vcf.wait()
	tabix2_vcf = subprocess.Popen([tabix, '-p', 'vcf', vcf_filtered_indels])
	tabix2_vcf.wait()
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '-nt', '5', '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', comb_vcf, '-genotypeMergeOptions', 'UNSORTED'])
	combine_var.wait()


def run_somalier(vcf, ped):
	#somalier extract -d extracted/ --sites sites.vcf.gz -f /data/human/g1k_v37_decoy.fa $cohort.vcf.gz
	# som_extract = subprocess.Popen([somalier, 'extract', '-d', 'somalier_extract/', '--sites', somalier_sites_vcf, '-f', fasta, vcf])
	# som_extract.wait()
	#somalier relate --ped $pedigree extracted/*.somalier
	# som_relate = subprocess.Popen([somalier, 'relate', '--ped', ped, 'somalier_extract/*.somalier'])
	som_relate = subprocess.Popen([somalier, 'relate', '--infer', '--ped', ped, 'somalier_extract/*.somalier'])
	som_relate.wait()

def process_annotate_vcf(input_vcf, out_vcf):
	temp_vcf =input_vcf.rsplit('.',2)[0] + '.temp.vcf.gz'
	##dcompress, change a gatk header thing,  (should be redundant)
	##so just and decompse and normalize
	#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
	# '''
	bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', temp_vcf, '-O', 'z', input_vcf])
	bcftools_norm.wait()
	tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', temp_vcf])
	tabix_vcf.wait()
	# '''
	##annotate with bcftools, have to use unphased version
	#bcftools csq -f hs37d5.fa -g Homo_sapiens.GRCh37.82.gff3.gz in.vcf -Ob -o out.bcf
	bcftools_csq = subprocess.Popen([bcftools, 'csq', '-g', ens_gff, '-f', fasta, '-o', out_vcf, '-O', 'z', '--local-csq', '-s', '-', temp_vcf])
	bcftools_csq.wait()
	tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', out_vcf])
	tabix_vcf.wait()


def slivar_std_analysis_on_trio(in_vcf, ped_file):
	slivar_std_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.vcf'
	slivar_cp_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.vcf'
	slivar_std_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.xls'
	slivar_cp_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.xls'
	slivar_std_annovar = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp'
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.hg38_multianno.txt'
	slivar_std_final = ped_file.rsplit('.', 1)[0] + '.std_analysis.xls'
	##get std file
	# '''
	slivar_trio = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
		'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
		'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
		'--family-expr', 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001',
		'--family-expr', 'recessive:fam.every(segregating_recessive)',
		'--family-expr', 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001',
		'--family-expr', 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)',
		'--trio', 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < 10',
		'-o', slivar_std_vcf])
	slivar_trio.wait()	
	##comp het from std vcf
	#slivar compound-hets -v vcfs/$cohort.vcf --sample-field comphet_side --sample-field denovo -p $ped > vcfs/$cohort.ch.vcf
	slivar_cp = subprocess.Popen([slivar, 'compound-hets', '--vcf', slivar_std_vcf, '--ped', ped_file,
				'--sample-field', 'comphet_side', '--sample-field', 'denovo', '-o', slivar_cp_vcf])
	slivar_cp.wait()
	##different params for diff annotations
	c_param = 'BCSQ'
	##convert to tsv (add extra annotation?)
	slivar_trio_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param,
		'-o', slivar_std_tsv, '-s', 'denovo', '-s', 'x_denovo', '-s', 'recessive', '-s', 'x_recessive',
		'-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt', '-g', pli_lookup,
		'-g', dbdb_lookup,'-g', omim_lookup, slivar_std_vcf])
	slivar_trio_tsv.wait()
	slivar_ch_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param, '-o', slivar_cp_tsv, 
		'-s', 'slivar_comphet', '-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt',
		'-g', pli_lookup, '-g', dbdb_lookup,'-g', omim_lookup, slivar_cp_vcf])
	slivar_ch_tsv.wait()
	# '''
	##if no variants don't run annovar
	with open(slivar_std_tsv, "r") as sst_fh:
		std_lc = 0
		for line in sst_fh:
			std_lc += 1

	if std_lc > 1:
		##annotate vcfs with annovar i.e. run_table_annovar
		command = [table_annovar] + av_buildver + [slivar_std_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_std_annovar]
		annovar = subprocess.Popen(command)
		annovar.wait()
		##combine slivar files with annovar annotation
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final, slivar_std_vcf)

def slivar_std_analysis_on_duo(in_vcf, ped_file):
	slivar_std_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.vcf'
	slivar_cp_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.vcf'
	slivar_std_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.xls'
	slivar_cp_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.xls'
	slivar_std_annovar = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp'
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.hg38_multianno.txt'
	slivar_std_final = ped_file.rsplit('.', 1)[0] + '.std_analysis.xls'
	##get std file
	slivar_duo = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
		'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
		'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
		'--family-expr', 'potential_denovo:fam.every(function(s) {return s.het == s.affected && s.hom_ref == !s.affected && s.GQ >=15 && s.DP >= 10 && s.AB > 0.20 == s.affected && s.AB < 0.1 == !s.affected}) && INFO.gnomad_popmax_af < 0.001',
		'--family-expr', 'recessive:fam.every(segregating_recessive)',
		'--family-expr', 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)',
		##gets all het and hom_ref vars in all inds, filtered by comp het
		'--family-expr', 'comphet_side:fam.every(function(s) {return (s.het || s.hom_ref) && s.GQ >= 15 && s.DP >= 10})',
		##are more specific alternative would be
		#'--family-expr', 'comphet_side:fam.every(function(s) {return (s.het || s.hom_ref) && s.GQ >= 20 }) && fam.some(function(s) { return s.het && s.affected })',
		'-o', slivar_std_vcf])
	slivar_duo.wait()
	##comp het from std vcf
	#slivar compound-hets -v vcfs/$cohort.vcf --sample-field comphet_side --sample-field denovo -p $ped > vcfs/$cohort.ch.vcf
	slivar_cp = subprocess.Popen([slivar, 'compound-hets', '--vcf', slivar_std_vcf, '--ped', ped_file,
				'-s', 'comphet_side', '-s', 'potential_denovo', '--allow-non-trios', '-o', slivar_cp_vcf])
	slivar_cp.wait()
	# '''
	##different params for diff annotations
	c_param = 'BCSQ'
	##convert to tsv (add extra annotation?)
	slivar_trio_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param,
		'-o', slivar_std_tsv, '-s', 'potential_denovo', '-s', 'recessive', '-s', 'x_recessive',
		'-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt', '-g', pli_lookup,
		'-g', dbdb_lookup,'-g', omim_lookup, slivar_std_vcf])
	slivar_trio_tsv.wait()
	slivar_ch_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param, '-o', slivar_cp_tsv, 
		'-s', 'slivar_comphet', '-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt',
		'-g', pli_lookup, '-g', dbdb_lookup,'-g', omim_lookup, slivar_cp_vcf])
	slivar_ch_tsv.wait()
	##annotate vcfs with annovar i.e. run_table_annovar
	command = [table_annovar] + av_buildver + [slivar_std_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_std_annovar]
	annovar = subprocess.Popen(command)
	annovar.wait()
	##combine slivar files with annovar annotation
	merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final, slivar_std_vcf)

	# '''


def slivar_std_analysis_on_multiplex(in_vcf, ped_file):
	slivar_std_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.vcf'
	slivar_cp_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.vcf'
	slivar_std_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.xls'
	slivar_cp_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.xls'
	slivar_std_annovar = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp'
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.hg38_multianno.txt'
	slivar_std_final = ped_file.rsplit('.', 1)[0] + '.std_analysis.xls'
	##get std file
	slivar_multi = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
		'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
		'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
		##all affexted family member are het
		'--family-expr', 'aff_het:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.het && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.2)})',
		'--family-expr', 'aff_hom:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.hom_alt && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.8)})',
		'--family-expr', 'comphet_side:fam.every(function(s) {return (s.het || s.hom_ref) && s.GQ >= 15 && s.DP >= 10 })',

		'-o', slivar_std_vcf])
	slivar_multi.wait()
	#slivar compound-hets -v vcfs/$cohort.vcf --sample-field comphet_side --sample-field denovo -p $ped > vcfs/$cohort.ch.vcf
	slivar_cp = subprocess.Popen([slivar, 'compound-hets', '--vcf', slivar_std_vcf, '--ped', ped_file,
				'-s', 'comphet_side', '-s', 'aff_het', '--allow-non-trios', '-o', slivar_cp_vcf])
	slivar_cp.wait()
	##different params for diff annotations
	c_param = 'BCSQ'
	##convert to tsv (add extra annotation?)
	slivar_trio_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param,
		'-o', slivar_std_tsv, '-s', 'aff_het', '-s', 'aff_hom',
		'-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt', '-g', pli_lookup,
		'-g', dbdb_lookup,'-g', omim_lookup, slivar_std_vcf])
	slivar_trio_tsv.wait()
	slivar_ch_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param, '-o', slivar_cp_tsv, 
		'-s', 'slivar_comphet', '-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt',
		'-g', pli_lookup, '-g', dbdb_lookup,'-g', omim_lookup, slivar_cp_vcf])
	slivar_ch_tsv.wait()
	##annotate vcfs with annovar i.e. run_table_annovar
	command = [table_annovar] + av_buildver + [slivar_std_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_std_annovar]
	annovar = subprocess.Popen(command)
	annovar.wait()
	##combine slivar files with annovar annotation
	merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final, slivar_std_vcf)


def slivar_std_analysis_on_trio_all_vars(in_vcf, ped_file):
	slivar_std_vcf = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_std.temp.vcf'
	slivar_cp_vcf = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_comphet.temp.vcf'
	slivar_std_tsv = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_std.temp.xls'
	slivar_cp_tsv = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_comphet.temp.xls'
	slivar_std_annovar = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_std.temp'
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_std.temp.hg38_multianno.txt'
	slivar_std_final = ped_file.rsplit('.', 1)[0] + '.all_vars.std_analysis.xls'
	##get std file
	# '''
	slivar_trio = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
		'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
		'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
		'--family-expr', 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001',
		'--family-expr', 'recessive:fam.every(segregating_recessive)',
		'--family-expr', 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001',
		'--family-expr', 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)',
		'--trio', 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < 10',
		'-o', slivar_std_vcf])
	slivar_trio.wait()	
	##comp het from std vcf
	#slivar compound-hets -v vcfs/$cohort.vcf --sample-field comphet_side --sample-field denovo -p $ped > vcfs/$cohort.ch.vcf
	slivar_cp = subprocess.Popen([slivar, 'compound-hets', '--vcf', slivar_std_vcf, '--ped', ped_file,
				'--sample-field', 'comphet_side', '--sample-field', 'denovo', '-o', slivar_cp_vcf])
	slivar_cp.wait()
	##different params for diff annotations
	c_param = 'BCSQ'
	##convert to tsv (add extra annotation?)
	slivar_trio_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param,
		'-o', slivar_std_tsv, '-s', 'denovo', '-s', 'x_denovo', '-s', 'recessive', '-s', 'x_recessive',
		'-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt', '-g', pli_lookup,
		'-g', dbdb_lookup,'-g', omim_lookup, slivar_std_vcf])
	slivar_trio_tsv.wait()
	slivar_ch_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param, '-o', slivar_cp_tsv, 
		'-s', 'slivar_comphet', '-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt',
		'-g', pli_lookup, '-g', dbdb_lookup,'-g', omim_lookup, slivar_cp_vcf])
	slivar_ch_tsv.wait()
	# '''
	##if no variants don't run annovar
	with open(slivar_std_tsv, "r") as sst_fh:
		std_lc = 0
		for line in sst_fh:
			std_lc += 1

	if std_lc > 1:
		##annotate vcfs with annovar i.e. run_table_annovar
		command = [table_annovar] + av_buildver + [slivar_std_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_std_annovar]
		annovar = subprocess.Popen(command)
		annovar.wait()
		##combine slivar files with annovar annotation
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final, slivar_std_vcf)


def slivar_std_analysis_on_duo_all_vars(in_vcf, ped_file):
	slivar_std_vcf = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_std.temp.vcf'
	slivar_cp_vcf = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_comphet.temp.vcf'
	slivar_std_tsv = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_std.temp.xls'
	slivar_cp_tsv = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_comphet.temp.xls'
	slivar_std_annovar = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_std.temp'
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.all_vars.slivar_std.temp.hg38_multianno.txt'
	slivar_std_final = ped_file.rsplit('.', 1)[0] + '.all_vars.std_analysis.xls'
	##get std file
	slivar_duo = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
		'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
		'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
		'--family-expr', 'potential_denovo:fam.every(function(s) {return s.het == s.affected && s.hom_ref == !s.affected && s.GQ >=15 && s.DP >= 10 && s.AB > 0.20 == s.affected && s.AB < 0.1 == !s.affected}) && INFO.gnomad_popmax_af < 0.001',
		'--family-expr', 'recessive:fam.every(segregating_recessive)',
		'--family-expr', 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)',
		##gets all het and hom_ref vars in all inds, filtered by comp het
		'--family-expr', 'comphet_side:fam.every(function(s) {return (s.het || s.hom_ref) && s.GQ >= 15 && s.DP >= 10})',
		##are more specific alternative would be
		#'--family-expr', 'comphet_side:fam.every(function(s) {return (s.het || s.hom_ref) && s.GQ >= 20 }) && fam.some(function(s) { return s.het && s.affected })',
		'-o', slivar_std_vcf])
	slivar_duo.wait()
	##comp het from std vcf
	#slivar compound-hets -v vcfs/$cohort.vcf --sample-field comphet_side --sample-field denovo -p $ped > vcfs/$cohort.ch.vcf
	slivar_cp = subprocess.Popen([slivar, 'compound-hets', '--vcf', slivar_std_vcf, '--ped', ped_file,
				'-s', 'comphet_side', '-s', 'potential_denovo', '--allow-non-trios', '-o', slivar_cp_vcf])
	slivar_cp.wait()
	# '''
	##different params for diff annotations
	c_param = 'BCSQ'
	##convert to tsv (add extra annotation?)
	slivar_trio_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param,
		'-o', slivar_std_tsv, '-s', 'potential_denovo', '-s', 'recessive', '-s', 'x_recessive',
		'-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt', '-g', pli_lookup,
		'-g', dbdb_lookup,'-g', omim_lookup, slivar_std_vcf])
	slivar_trio_tsv.wait()
	slivar_ch_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param, '-o', slivar_cp_tsv, 
		'-s', 'slivar_comphet', '-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt',
		'-g', pli_lookup, '-g', dbdb_lookup,'-g', omim_lookup, slivar_cp_vcf])
	slivar_ch_tsv.wait()
	##annotate vcfs with annovar i.e. run_table_annovar
	command = [table_annovar] + av_buildver + [slivar_std_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_std_annovar]
	annovar = subprocess.Popen(command)
	annovar.wait()
	##combine slivar files with annovar annotation
	merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final, slivar_std_vcf)

	# '''

def merge_std_silvar_annovar(std_tsv, cp_tsv, multianno, outfile, std_vcf):
	##get vcf names fro std vcf
	with open(std_vcf, "r") as sv_fh:
		for line in sv_fh:
			if line[:2] != '##' and line[0] == '#':
				line = line.rstrip().split(delim)
				vcf_samples = line[9:]
	##make dict from annovar multianno
	multianno_dict = {}
	with open(multianno, "r") as ma_fh:
		lc = 0
		for line in ma_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc == 1:
				##add header info
				extra_header = line[8:10] + line[13:15] + line[34:36] + line[49:51] + line[61:79] + ['vcf_info', 'vcf_format'] + vcf_samples
			else:
				chr_pos_ref_alt = ':'.join(line[82:84] + line[85:87])
				##keep refgene, clivar, cadd, gerp, polyphe, extra gnomad, vcf info
				info_to_keep = line[8:10] + line[13:15] + line[34:36] + line[49:51] + line[61:79] + line[89:]				
				multianno_dict[chr_pos_ref_alt] = info_to_keep
	with open(outfile, "w") as out_fh:
		with open(std_tsv, "r") as std_fh:
			lc = 0
			for line_1 in std_fh:	
				line_1 = line_1.rstrip('\n').split(delim)
				lc += 1
				if lc == 1:
					header = line_1[:13] + ['pLI', 'dbdb', 'omim'] + extra_header
					out_fh.write(delim.join(header) + '\n')
				else:
					cpra = line_1[3]
					line1_out = line_1[:16] + multianno_dict[cpra]
					out_fh.write(delim.join(line1_out) + '\n')
		with open(cp_tsv, "r") as cp_fh:
			lc = 0
			for line_2 in cp_fh:
				line_2 = line_2.rstrip('\n').split(delim)
				lc += 1
				if lc > 1:
					cpra = line_2[3]
					line2_out = line_2[:16] + multianno_dict[cpra]
					out_fh.write(delim.join(line2_out) + '\n')


def silvar_by_ped_type(ped_info_file, project_prefix, annotated_vcf):
	##get sample counts for the three ped types and make ped files
	ped_type_dict = {'trios': 0, 'single_duos': 0, 'multiplex': 0}
	trio_ped_file = project_prefix + '.trios.ped' 
	single_duos_ped_file = project_prefix + '.single_duos.ped' 
	multiplex_ped_file = project_prefix + '.multiplex.ped' 
	with open(ped_info_file, 'r') as pi_fh, open(trio_ped_file, 'w') as tpf_fh, open(single_duos_ped_file, 'w') as sdpf_fh, open(multiplex_ped_file, 'w') as mpf_fh:
		lc = 0
		for line in pi_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc == 1:
				tpf_fh.write(delim.join(line[:6]) + '\n')
				sdpf_fh.write(delim.join(line[:6]) + '\n')
				mpf_fh.write(delim.join(line[:6]) + '\n')
			else:		
				ped_type = line[6]
				if ped_type in trio_types:
					ped_type_dict['trios'] += 1
					tpf_fh.write(delim.join(line[:6]) + '\n')
				elif ped_type in duo_single_types:
					ped_type_dict['single_duos'] += 1
					sdpf_fh.write(delim.join(line[:6]) + '\n')
				elif ped_type in multiplex_types:
					ped_type_dict['multiplex'] += 1
					mpf_fh.write(delim.join(line[:6]) + '\n')

				else:
					print('ped type: ' + ped_type + ' not recognized')
	print(ped_type_dict)
	##filter with slivar
	for ped_type in ped_type_dict:
		if ped_type_dict[ped_type] > 0:
			ped_file = project_prefix + '.' + ped_type + '.ped' 
			# print(annotated_vcf, ped_file, formatted_ped_type)
			if ped_type == 'trios':
				# slivar_std_analysis_on_trio(annotated_vcf, ped_file)
				slivar_std_analysis_on_trio_all_vars(annotated_vcf, ped_file)
			elif ped_type == 'single_duos':
				# slivar_std_analysis_on_duo(annotated_vcf, ped_file)
				slivar_std_analysis_on_duo_all_vars(annotated_vcf, ped_file)
			elif ped_type == 'multiplex':
				# slivar_std_analysis_on_multiplex(annotated_vcf, ped_file)
				pass
			else:
				print(ptype, ' ped type not recognized')


def slivar_analysis_master(ped_info_file, project_prefix, in_vcf, annotated_vcf):
	##normalize and annotate the combined vcf
	# process_annotate_vcf(in_vcf, annotated_vcf)
	##analyze in ped type groups
	silvar_by_ped_type(ped_info_file, project_prefix, annotated_vcf)





##run methods
working_dir = '/home/atimms/ngs_data/genomes/kim_genomes_0721'
os.chdir(working_dir)

##params
project_name = 'kim_genomes_0721'
snp_vcf = 'Joint.HC.VQSR.snp.g.vcf.gz'
indel_vcf = 'Joint.HC.VQSR.indel.g.vcf.gz'
combined_vcf = 'Joint.HC.VQSR.all.vcf.gz'
combined_ped = 'combined.ped'
annotated_comb_vcf = project_name + '.bcftools.GRCh38.103.vcf.gz'
ped_info_file = 'kim_genomes_0721.ped_info.txt'
# ped_info_file = 'test.ped_info.txt'

##combine the snp and indel vcfs
# combine_snp_indel_vcfs(snp_vcf, indel_vcf, combined_vcf)

##run somalier on combined vcf --- ped just 1st 6 cols of pedinfo
# run_somalier(combined_vcf, combined_ped)

##run silvar
slivar_analysis_master(ped_info_file, project_name, combined_vcf, annotated_comb_vcf)


