#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
requires:
qsub -Iq cdbrmq -l mem=30gb,ncpus=2 -P 9aa67182-09d6-405c-b7f7-93b0320828c1
qsub -Iq cdbrmq -l mem=100gb,ncpus=10 -P 9aa67182-09d6-405c-b7f7-93b0320828c1
#java 1.8 for snpsift
module load java/1.8.0_202 
'''

##set input variables and parameters
delim = '\t'
threads = '5'

##programs
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
snpeff_jar = '/home/atimms/programs/snpEff_0121/snpEff.jar'
snpsift_jar  = '/home/atimms/programs/snpEff_0121/SnpSift.jar'
# slivar = '/home/atimms/programs/slivar_0121/slivar'
slivar = '/home/atimms/programs/slivar_0121/slivar_dev'
pslivar = '/home/atimms/programs/slivar_0121/pslivar'
slivar_functions = '/home/atimms/programs/slivar_0121/slivar-0.2.1/js/slivar-functions_at.js'
slivar_functions_nodp = '/home/atimms/programs/slivar_0121/slivar-0.2.1/js/slivar-functions_nodp.js'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
freebayes = '/home/atimms/programs/freebayes-1.3.4-linux-static-AMD64'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
plink = '/home/atimms/programs/plink'
dbsnp_common_vcf = '/home/atimms/ngs_data/references/hg19/common_all_20180423.vcf.gz'
tabix = '/home/atimms/programs/samtools-1.11/bin/tabix'


##files
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
ens_gff = ref_dir + 'Homo_sapiens.GRCh37.87.gff3.gz'
gnomad_gnotate = '/home/atimms/ngs_data/references/slivar/gnomad.hg37.zip'
#wget -qO - https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | zcat | cut -f 1,21,24 | tail -n+2 | awk '{ printf("%s\tpLI=%.3g;oe_lof=%.5g\n", $1, $2, $3)}' > pli.lookup
pli_lookup = '/home/atimms/ngs_data/references/slivar/pli.lookup'
#cut -f2,3,4,6,8 dbdb.gene_association.txt | awk 'BEGIN{FS="\t"}{ printf("%s\tinheritance=%s;phenotype=%s;syndrome=%s;loe=%s\n", $1, $2, $3, $4, $5)}' > ~/ngs_data/references/slivar/dbdb.lookup
dbdb_lookup = '/home/atimms/ngs_data/references/slivar/dbdb.lookup'
#grep -v '^#' omim_genemap2_0816.txt | cut -f9,13 | grep -v '^\s' > /home/atimms/ngs_data/references/slivar/omim.lookup
omim_lookup = '/home/atimms/ngs_data/references/slivar/omim.lookup'
hg19_selfchain = '/home/atimms/ngs_data/references/slivar/selfchain-LCR.hg19.bed'
hg19_refgene_exons = '/home/atimms/ngs_data/references/slivar/hg19.refflat_exons.nochr.bed'
decipher_delevopmental_file = '/home/atimms/ngs_data/references/slivar/DDG2P_11_2_2021.csv'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,gnomad211_exome,gnomad211_genome,dbnsfp41a,clinvar_20210123,generic', '-genericdbfile', 'Geno2MP.avinput']
av_operation = ['-operation', 'g,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput' ]


##run methods
def plink_relatadness_check(vcf, file_prefix, seq_type):
	## add dbsnp to vcf and filter for known snps and pop freq
	# java -jar SnpSift.jar annotate dbSnp132.vcf input.vcf > output_annotated.vcf
	# '''
	snpsift1_vcf = file_prefix + '.snpsift_temp1.vcf'
	with open (snpsift1_vcf, 'w') as ss1_fh:
		snpeff_ann = subprocess.Popen(['java', '-Xmx20G', '-jar', snpsift_jar, 'annotate', dbsnp_common_vcf, vcf], stdout=ss1_fh)
		snpeff_ann.wait()
	##make file (id_dot.txt) with a dot in it to filter against
	with open('id_dot.txt', 'w') as id_fh:
		id_fh.write('.')
	##correct filtering?? - use id_dot.txt to remove non-dbsnp snps
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(FORMAT/DP)>20 & ID!=@id_dot.txt", '-V', 'indels', '-o', 'temp_plink_genome.vcf.gz', '-O', 'z', snpsift1_vcf])
	bcftools_filter.wait()
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(INFO/DP)>50", '-V', 'indels', '-o', 'temp_plink_sex.vcf.gz', '-O', 'z', vcf])
	bcftools_filter.wait()
	##generate plink file from vcf
	make_plink = subprocess.Popen([plink, '--vcf', 'temp_plink_genome.vcf.gz', '--vcf-filter', '--const-fid', '--out', 'temp_genome.pass_q50_dp50'])
	make_plink.wait()
	make_plink2 = subprocess.Popen([plink, '--vcf', 'temp_plink_sex.vcf.gz', '--vcf-filter', '--const-fid', '--out', 'temp_sex.pass_q50_dp50'])
	make_plink2.wait()
	##check sex -results in .sexcheck file
	plink_sex = subprocess.Popen([plink, '--bfile', 'temp_sex.pass_q50_dp50', '--make-bed', '--split-x', 'hg19', '--out', 'temp_sex.pass_q50_dp50_splitx'])
	plink_sex.wait()	
	plink_sex = subprocess.Popen([plink, '--bfile', 'temp_sex.pass_q50_dp50_splitx', '--check-sex', '--out', file_prefix + '.pass_q50_dp50'])
	plink_sex.wait()
	# '''
	##ibd check
	if seq_type == 'exome':
		plink_indep = subprocess.Popen([plink, '--bfile', 'temp_genome.pass_q50_dp50', '--indep', '50', '5', '2'])
		plink_indep.wait()
	elif seq_type == 'genome':
		plink_indep = subprocess.Popen([plink, '--bfile', 'temp_genome.pass_q50_dp50', '--indep', '50', '5', '2', '--maf', '0.10'])
		plink_indep.wait()
	plink_extract = subprocess.Popen([plink, '--bfile', 'temp_genome.pass_q50_dp50', '--extract', 'plink.prune.in', '--make-bed', '--out', 'temp_genome.pass_q50_dp50_ldprune'])
	plink_extract.wait()	
	plink_ibd = subprocess.Popen([plink,  '--bfile', 'temp_genome.pass_q50_dp50_ldprune', '--genome', '--out', file_prefix + '.pass_q50_dp50'])
	plink_ibd.wait()



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
				extra_header = line[8:10] + line[96:100] + ['Geno2MP'] + [line[69]] + line[83:85] + line[47:49] + line[50:52] + [line[58]] + line[10:44] + ['vcf_format'] + vcf_samples
			else:
				chr_pos_ref_alt = ':'.join(line[104:106] + line[107:109])
				##keep refgene, clivar/Geno2MP, cadd, gerp, polyphe, mutataster/mutass,REVEL, extra gnomad, vcf info
				info_to_keep = line[8:10] + line[96:101] + [line[69]] + line[83:85] + line[47:49]  + line[50:52] + [line[58]] + line[10:44] + line[112:]
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

def process_annotate_vcf(bcftools_vcf, input_vcf, file_prefix):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	temp_vcf = file_prefix + '.temp1.vcf.gz'
	##decomposing and subsetting vcf
	# zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
	# sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
	# bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', temp_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
	# bcftools_norm.wait()
	##using single or 10 threads with no sed command as ad good in vcf
	bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', temp_vcf, '-O', 'z', input_vcf])
	# bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '--threads', '10', '-w', '10000', '-f', fasta, '-o', temp_vcf, '-O', 'z', input_vcf])
	bcftools_norm.wait()
	##annotate with bcftools, have to use unphased version
	#bcftools csq -f hs37d5.fa -g Homo_sapiens.GRCh37.82.gff3.gz in.vcf -Ob -o out.bcf
	##single thread
	bcftools_csq = subprocess.Popen([bcftools, 'csq', '-g', ens_gff, '-f', fasta, '-o', bcftools_vcf, '-O', 'z', '--local-csq', '-s', '-', temp_vcf])
	##using 10 threads
	# bcftools_csq = subprocess.Popen([bcftools, 'csq', '-g', ens_gff, '-f', fasta, '-o', bcftools_vcf, '--threads', '10', '-O', 'z', '--local-csq', '-s', '-', temp_vcf])
	bcftools_csq.wait()
	# tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', bcftools_vcf])
	# tabix_vcf.wait()

def slivar_std_analysis_on_trio(in_vcf, ped_file, analysis_type):
	slivar_std_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.vcf'
	slivar_cp_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.vcf'
	slivar_std_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.xls'
	slivar_cp_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.xls'
	slivar_std_annovar = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp'
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.hg19_multianno.txt'
	slivar_std_final = ped_file.rsplit('.', 1)[0] + '.std_analysis.xls'
	##get std file
	# '''
	if analysis_type == 'std':
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
	elif analysis_type == 'nodp':
		slivar_trio = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions_nodp, '--info',
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
		# '''
		##combine slivar files with annovar annotation
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final, slivar_std_vcf)
	

def slivar_std_analysis_on_duo(in_vcf, ped_file, analysis_type):
	slivar_std_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.vcf'
	slivar_cp_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.vcf'
	slivar_std_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.xls'
	slivar_cp_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.xls'
	slivar_std_annovar = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp'
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.hg19_multianno.txt'
	slivar_std_final = ped_file.rsplit('.', 1)[0] + '.std_analysis.xls'
	##get std file
	# '''
	if analysis_type == 'std':
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
	elif analysis_type == 'nodp':
		slivar_duo = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
			'--family-expr', 'potential_denovo:fam.every(function(s) {return s.het == s.affected && s.hom_ref == !s.affected && s.GQ >=15 && s.AB > 0.20 == s.affected && s.AB < 0.1 == !s.affected}) && INFO.gnomad_popmax_af < 0.001',
			'--family-expr', 'recessive:fam.every(segregating_recessive)',
			'--family-expr', 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)',
			##gets all het and hom_ref vars in all inds, filtered by comp het
			'--family-expr', 'comphet_side:fam.every(function(s) {return (s.het || s.hom_ref) && s.GQ >= 15})',
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
	
def slivar_std_analysis_on_multiplex(in_vcf, ped_file, analysis_type):
	slivar_std_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.vcf'
	slivar_cp_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.vcf'
	slivar_std_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.xls'
	slivar_cp_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.xls'
	slivar_std_annovar = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp'
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.hg19_multianno.txt'
	slivar_std_final = ped_file.rsplit('.', 1)[0] + '.std_analysis.xls'
	##get std file
	# '''
	if analysis_type == 'std':
		slivar_multi = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
			##all affexted family member are het
			'--family-expr', 'aff_het:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.het && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.2)})',
			'--family-expr', 'aff_hom:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.hom_alt && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.8)})',
			'--family-expr', 'comphet_side:fam.every(function(s) {return (s.het || s.hom_ref) && s.GQ >= 15 && s.DP >= 10 })',

			'-o', slivar_std_vcf])
		slivar_multi.wait()
	elif analysis_type == 'nodp':
		slivar_multi = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
			##all affexted family member are het
			'--family-expr', 'aff_het:fam.every(function(s) {return (!s.affected && s.GQ >= 15) || (s.affected && s.het && s.GQ >= 15 && s.AB > 0.2)})',
			'--family-expr', 'aff_hom:fam.every(function(s) {return (!s.affected && s.GQ >= 15) || (s.affected && s.hom_alt && s.GQ >= 15 && s.AB > 0.8)})',
			'--family-expr', 'comphet_side:fam.every(function(s) {return (s.het || s.hom_ref) && s.GQ >= 15})',

			'-o', slivar_std_vcf])
		slivar_multi.wait()
	# '''
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
	# '''

##not used now
def get_specific_ped_file_get_samples(in_ped, ped_name, out_ped):
	sample_list = []
	with open(in_ped, 'r') as in_fh, open(out_ped, 'w') as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc == 1:
				out_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				fam_id = line[0]
				ind_id = line[1]
				if fam_id == ped_name:
					out_fh.write(delim.join(line) + '\n')
					sample_list.append(ind_id)
	return(sample_list)	

##not used now
def get_list_of_samples(infile):
	samples = []
	with open(infile, 'r') as in_fh:
		for line in in_fh:
			sample = line.rstrip()
			samples.append(sample)
	return(samples)
##not used now
def make_sample_file(samples, outfile, samples_to_include):
	with open(outfile, "w") as out_fh:
		for sample in samples:
			if sample in samples_to_include:
				out_fh.write(sample + '\n')

##not used now
def get_vcf_per_ped(in_vcf, sample_list, vcf_samples, ped_vcf, ped):
	sample_file = ped + '.samples_temp.txt'
	make_sample_file(sample_list, sample_file, vcf_samples)
	##get the samples we want, and remove when we don't see a call, using one or multiple threads
	# bcftools_view = subprocess.Popen([bcftools, 'view', '-a', '-Ou', '-S', sample_file, in_vcf], stdout=subprocess.PIPE)
	# bcftools_view2 = subprocess.Popen([bcftools, 'view', '-m', '2', '-O', 'z', '-o', ped_vcf, '-'], stdin=bcftools_view.stdout)
	bcftools_view = subprocess.Popen([bcftools, 'view', '-a', '--threads', '4', '-Ou', '-S', sample_file, in_vcf], stdout=subprocess.PIPE)
	bcftools_view2 = subprocess.Popen([bcftools, 'view', '-m', '2', '--threads', '4', '-O', 'z', '-o', ped_vcf, '-'], stdin=bcftools_view.stdout)
	bcftools_view2.wait()

def split_ped_file_into_ped_types(infile, ped_dict, out_prefix):
	for ptype in ped_dict:
		out_ped = out_prefix + '.' + ptype + '.ped'
		with open(infile, 'r') as in_fh, open(out_ped, 'w') as out_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				if lc == 1:
					out_fh.write(line)
				else:
					line = line.rstrip().split(delim)
					fam_id = line[0]
					if fam_id in ped_dict[ptype]:
						out_fh.write(delim.join(line) + '\n')



def standard_slivar_protocol_v2(ped_info_file, project_name, combined_vcf, combined_ped, anal_type):
	trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
	duo_types = ['duo', 'duo*', 'parent_sibship']
	single_types = ['singleton', 'singleton*', 'sibship']
	duo_single_types = duo_types + single_types
	multiplex_types = ['multiplex']
	##annoate the combined vcf
	vcf_prefix = combined_vcf.rsplit('.',2)[0]
	annotated_vcf = vcf_prefix + '.bcftools.GRCh37_87.vcf.gz'
	process_annotate_vcf(annotated_vcf, combined_vcf, vcf_prefix)
	##get ped ids for the three groups
	ped_type_dict = {'trios': [], 'single_duos': [], 'multiplex':[]}
	with open(ped_info_file, 'r') as pi_fh:
		for line in pi_fh:
			line = line.rstrip().split(delim)
			ped = line[0]
			ped_type = line[1]
			if ped_type in trio_types:
				ped_type_dict['trios'].append(ped)
			elif ped_type in duo_single_types:
				ped_type_dict['single_duos'].append(ped)
			elif ped_type in multiplex_types:
				ped_type_dict['multiplex'].append(ped)
			else:
				print('ped type: ' + ped_type + ' not recognized')
	##make ped file for the three groups
	split_ped_file_into_ped_types(combined_ped, ped_type_dict, project_name)
	##filter with slivar
	for ped_type in ped_type_dict:
		ped_file = project_name + '.' + ped_type + '.ped' 
		# print(annotated_vcf, ped_file, formatted_ped_type)
		# '''
		if ped_type == 'trios':
			slivar_std_analysis_on_trio(annotated_vcf, ped_file, anal_type)
		elif ped_type == 'single_duos':
			slivar_std_analysis_on_duo(annotated_vcf, ped_file, anal_type)
		elif ped_type == 'multiplex':
			slivar_std_analysis_on_multiplex(annotated_vcf, ped_file, anal_type)
		else:
			print(ptype, ' ped type not recognized')
		# '''
		# if ped_type.lower() in trio_types or ped_type.lower() in multiplex_types or ped_type.lower() in duo_types:
		# 	slivar_duo_del(annotated_vcfs, ped_file, formatted_ped_type)




##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/daniela_slivar_analysis_0221'
os.chdir(working_dir)

##ped types for slivar analysis
trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
duo_single_types = ['duo', 'duo*', 'parent_sibship', 'singleton', 'singleton*', 'sibship']
multiplex_types = ['multiplex']

##exomes
exome_vcf = 'microtia_exomes_recall_0221.gatkHC.vcf.gz'
exome_ped = 'exomes_0221.ped'
exome_prefix = 'exomes_0221'
exome_ped_info = 'exomes_0221.ped_info.txt'

##check sex and ibd sharing
# plink_relatadness_check(exome_vcf, exome_prefix, 'exome')

##run slivar after splitting ped and vcf file
# standard_slivar_protocol_v2(exome_ped_info, exome_prefix, exome_vcf, exome_ped, 'std')

##genomes -- combined all genomes (these are odd with depth missing in some cases) 
genome_vcf = 'luquetti_grc_wgs_combined.sample_ids.vcf.gz'
genome_ped = 'genome_0221.ped'
genome_prefix = 'genome_0221'
genome_ped_info = 'genome_0221.ped_info.txt'
genome_sample_in_vcf = 'genome_original_vcf_ids.txt'

##check sex and ibd sharing
# plink_relatadness_check(genome_vcf, genome_prefix, 'genome')

##run slivar after splitting ped and vcf file
# standard_slivar_protocol_v2(genome_ped_info, genome_prefix, genome_vcf, genome_ped, 'nodp')


##genomes 0219
##maunally switched ids to our using bcftools
genome_vcf = 'luquetti_external_uwcmg_cfm_1.HF.sample_ids.vcf.gz'
genome_ped = 'genome_0219.ped'
genome_prefix = 'genome_0219'
genome_ped_info = 'genome_0219.ped_info.txt'

##run slivar after splitting ped and vcf file
standard_slivar_protocol_v2(genome_ped_info, genome_prefix, genome_vcf, genome_ped, 'std')



##genomes 0620
##maunally switched ids to our using bcftools
genome_vcf = 'luquetti_grc_wgs_2.HF.sample_ids.vcf.gz'
genome_ped = 'genome_0620.ped'
genome_prefix = 'genome_0620'
genome_ped_info = 'genome_0620.ped_info.txt'

##run slivar after splitting ped and vcf file
# standard_slivar_protocol_v2(genome_ped_info, genome_prefix, genome_vcf, genome_ped, 'std')




