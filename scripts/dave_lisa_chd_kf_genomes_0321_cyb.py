#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
requires:
qsub -Iq cdbrmq -l mem=20gb,ncpus=1 -P 2cf4ba3b-f9ef-4cfc-a2c1-be85f5db6730
#java 11 for snpeff
module load java/11.0.10
#biobuilds for tabix etc
module load biobuilds/2017.11

using peddy:
conda activate peddy


'''

##set input variables and parameters
delim = '\t'
threads = '5'

##programs
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
snpeff_jar = '/home/atimms/programs/snpEff_0121/snpEff.jar'
slivar = '/home/atimms/programs/slivar_0121/slivar'
# slivar = '/home/atimms/programs/slivar_0121/slivar_dev'
pslivar = '/home/atimms/programs/slivar_0121/pslivar'
slivar_functions = '/home/atimms/programs/slivar_0121/slivar-0.2.1/js/slivar-functions_at.js'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
freebayes = '/home/atimms/programs/freebayes-1.3.4-linux-static-AMD64'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
tabix = '/home/atimms/programs/samtools-1.11/bin/tabix'
snpsift_jar  = '/home/atimms/programs/snpEff_0121/SnpSift.jar'
plink = '/home/atimms/programs/plink'

##files
fasta = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens_assembly38.fasta'
ens_gff = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens.GRCh38.103.gff3.gz'
gnomad_gnotate = '/home/atimms/ngs_data/references/slivar/gnomad.hg38.genomes.v3.fix.zip'
#wget -qO - https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | zcat | cut -f 1,21,24 | tail -n+2 | awk '{ printf("%s\tpLI=%.3g;oe_lof=%.5g\n", $1, $2, $3)}' > pli.lookup
pli_lookup = '/home/atimms/ngs_data/references/slivar/pli.lookup'
#grep -v '^#' omim_genemap2_0321.txt | cut -f9,13 | grep -v '^\s' > /home/atimms/ngs_data/references/slivar/omim.lookup
omim_lookup = '/home/atimms/ngs_data/references/slivar/omim.lookup'
dbsnp_common_vcf = '/home/atimms/ngs_data/references/hg38_gatk/common_all_20180418.chr_added.vcf.gz'

##annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,dbnsfp41a,clinvar_20210123']
av_operation = ['-operation', 'g,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput' ]

##methods
def make_ped_files(input_file):
	ped_dict = {}
	with open(input_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc == 1:
				header = line
			else:
				ped_name = line[0]
				if ped_name in ped_dict:
					ped_dict[ped_name].append(line)
				else:
					ped_dict[ped_name] = [line]
	for ped in ped_dict:
		ped_file = ped + '.ped'
		ad_ped_file = ped + '.ad.ped'
		with open(ped_file, "w") as pf_fh, open(ad_ped_file, "w") as apf_fh:
			pf_fh.write(delim.join(header) + '\n')
			apf_fh.write(delim.join(header) + '\n')
			for outline in ped_dict[ped]:
				pf_fh.write(delim.join(outline) + '\n')
				aff_status = outline[5]
				if aff_status == '1':
					dom_outline = outline[:5] + ['2']
					apf_fh.write(delim.join(dom_outline) + '\n')
				else:
					apf_fh.write(delim.join(outline) + '\n')


def plink_relatadness_check(input_file):
	with open(input_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			file_prefix = line[0]
			vcf = line[1]
			tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', vcf])
			tabix_vcf.wait()
			## add dbsnp to vcf and filter for known snps and pop freq
			##not needed as have rs numbers
			# java -jar SnpSift.jar annotate dbSnp132.vcf input.vcf > output_annotated.vcf
			'''
			snpsift1_vcf = file_prefix + '.snpsift_temp1.vcf'
			
			with open (snpsift1_vcf, 'w') as ss1_fh:
				snpeff_ann = subprocess.Popen(['java', '-Xmx20G', '-jar', snpsift_jar, 'annotate', dbsnp_common_vcf, vcf], stdout=ss1_fh)
				snpeff_ann.wait()
			##make file (id_dot.txt) with a dot in it to filter against
			with open('id_dot.txt', 'w') as id_fh:
				id_fh.write('.')
			##correct filtering?? - use id_dot.txt to remove non-dbsnp snps
			bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(FORMAT/DP)>20 & ID!=@id_dot.txt", '-V', 'indels', '-o', file_prefix + 'temp_plink_genome.vcf.gz', '-O', 'z', vcf])
			bcftools_filter.wait()
			bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(INFO/DP)>50", '-V', 'indels', '-o', file_prefix + 'temp_plink_sex.vcf.gz', '-O', 'z', vcf])
			bcftools_filter.wait()
			##generate plink file from vcf
			make_plink = subprocess.Popen([plink, '--vcf', file_prefix + 'temp_plink_genome.vcf.gz', '--vcf-filter', '--const-fid', '--out', file_prefix + 'temp_genome.pass_q50_dp50'])
			make_plink.wait()
			make_plink2 = subprocess.Popen([plink, '--vcf', file_prefix + 'temp_plink_sex.vcf.gz', '--vcf-filter', '--const-fid', '--out', file_prefix + 'temp_sex.pass_q50_dp50'])
			make_plink2.wait()
			##check sex -results in .sexcheck file
			plink_sex = subprocess.Popen([plink, '--bfile', file_prefix + 'temp_sex.pass_q50_dp50', '--make-bed', '--split-x', 'hg19', '--out', file_prefix + 'temp_sex.pass_q50_dp50_splitx'])
			plink_sex.wait()	
			plink_sex = subprocess.Popen([plink, '--bfile', file_prefix + 'temp_sex.pass_q50_dp50_splitx', '--check-sex', '--out', file_prefix + '.pass_q50_dp50'])
			plink_sex.wait()
			##ibd check
			plink_indep = subprocess.Popen([plink, '--bfile', file_prefix + 'temp_genome.pass_q50_dp50', '--indep', '50', '5', '2', '--maf', '0.10'])
			plink_indep.wait()
			plink_extract = subprocess.Popen([plink, '--bfile', file_prefix + 'temp_genome.pass_q50_dp50', '--extract', 'plink.prune.in', '--make-bed', '--out', file_prefix + 'temp_genome.pass_q50_dp50_ldprune'])
			plink_extract.wait()	
			plink_ibd = subprocess.Popen([plink,  '--bfile', file_prefix + 'temp_genome.pass_q50_dp50_ldprune', '--genome', '--out', file_prefix + '.pass_q50_dp50'])
			plink_ibd.wait()
			'''
			# '''
			##correct filtering?? - use id_dot.txt to remove non-dbsnp snps
			# bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(INFO/DP)>60 & ID!=@id_dot.txt", '-V', 'indels', '-o', file_prefix + 'temp_plink.vcf.gz', '-O', 'z', vcf])
			# bcftools_filter.wait()
			##generate plink file from vcf
			# make_plink = subprocess.Popen([plink, '--vcf', file_prefix + 'temp_plink.vcf.gz', '--vcf-filter', '--const-fid', '--out', file_prefix + 'temp.pass_q50_dp50'])
			# make_plink.wait()
			##check sex -results in .sexcheck file
			# plink_sex1 = subprocess.Popen([plink, '--bfile', file_prefix + 'temp.pass_q50_dp50', '--make-bed', '--split-x', 'hg38', '--out', file_prefix + 'temp.pass_q50_dp50_splitx'])
			# plink_sex1.wait()	
			# plink_sex2 = subprocess.Popen([plink, '--bfile', file_prefix + 'temp.pass_q50_dp50', '--check-sex', '--out', file_prefix + '.pass_q50_dp50'])
			# plink_sex2.wait()
			##ibd check
			plink_indep = subprocess.Popen([plink, '--bfile', file_prefix + 'temp.pass_q50_dp50', '--indep', '50', '5', '2', '--maf', '0.10'])
			plink_indep.wait()
			plink_extract = subprocess.Popen([plink, '--bfile', file_prefix + 'temp.pass_q50_dp50', '--extract', 'plink.prune.in', '--make-bed', '--out', file_prefix + 'temp.pass_q50_dp50_ldprune'])
			plink_extract.wait()	
			plink_ibd = subprocess.Popen([plink,  '--bfile', file_prefix + 'temp.pass_q50_dp50_ldprune', '--genome', '--out', file_prefix + '.pass_q50_dp50'])
			plink_ibd.wait()
			# '''

def peddy_relatadness_check(input_file):
	with open(input_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			ped_name = line[0]
			ped_file = ped_name + '.ped'
			vcf = line[1]
			tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', vcf])
			tabix_vcf.wait()
			#python -m peddy --plot --sites hg38 --prefix ceph-1463 data/ceph1463.peddy.vcf.gz data/ceph1463.ped
			run_peddy = subprocess.Popen(['python', '-m', 'peddy',  '--plot', '--sites', 'hg38', '--prefix', ped_name + '.peddy', vcf, ped_file])
			run_peddy.wait()

def process_annotate_vcf(ped, input_vcf, out_vcf):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	##so just normalize now and use vep annotation that's already there
	#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
	bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', out_vcf, '-O', 'z', input_vcf])
	bcftools_norm.wait()
	tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', out_vcf])
	tabix_vcf.wait()
	##annotate with bcftools, have to use unphased version
	'''
	#bcftools csq -f hs37d5.fa -g Homo_sapiens.GRCh37.82.gff3.gz in.vcf -Ob -o out.bcf
	bcftools_csq = subprocess.Popen([bcftools, 'csq', '-g', ens_gff, '-f', fasta, '-o', bcftools_vcf, '-O', 'z', '--local-csq', '-s', '-', temp_vcf])
	bcftools_csq.wait()
	tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', bcftools_vcf])
	tabix_vcf.wait()
	'''

def merge_std_silvar_annovar(std_tsv, cp_tsv, multianno, outfile):
	##make dict from annovar multianno
	multianno_dict = {}
	with open(multianno, "r") as ma_fh:
		lc = 0
		for line in ma_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc == 1:
				##add header info
				extra_header = line[8:10] + line[13:15] + line[34:36] + line[49:51] + ['vcf_info']
			else:
				chr_pos_ref_alt = ':'.join(line[69:71] + line[72:74])
				##keep refgene, clivar/Geno2MP, cadd, gerp, polyphe, mutataster/mutass,REVEL, extra gnomad, vcf info
				info_to_keep = line[8:10] + line[13:15] + line[34:36] + line[49:51] + line[77:]
				multianno_dict[chr_pos_ref_alt] = info_to_keep
	with open(outfile, "w") as out_fh:
		with open(std_tsv, "r") as std_fh:
			lc = 0
			for line_1 in std_fh:	
				line_1 = line_1.rstrip('\n').split(delim)
				lc += 1
				if lc == 1:
					header = line_1[:13] + ['pLI', 'omim'] + extra_header
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

def slivar_std_analysis_on_trio(in_vcf, ped_file):
	slivar_std_vcf = in_vcf.split('.')[0] + '.slivar_std.temp.vcf'
	slivar_cp_vcf = in_vcf.split('.')[0] + '.slivar_comphet.temp.vcf'
	slivar_std_tsv = in_vcf.split('.')[0] + '.slivar_std.temp.xls'
	slivar_cp_tsv = in_vcf.split('.')[0] + '.slivar_comphet.temp.xls'
	slivar_std_annovar = in_vcf.split('.')[0] + '.slivar_std.temp'
	slivar_std_multi = in_vcf.split('.')[0] + '.slivar_std.temp.hg38_multianno.txt'
	slivar_std_final = in_vcf.split('.')[0] + '_slivar.std_analysis.xls'
	##get std file
	# '''
	slivar_trio = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
		'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
		'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
		'--family-expr', 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001',
		'--family-expr', 'recessive:fam.every(segregating_recessive)',
		'--family-expr', 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001',
		'--family-expr', 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)',
		##added these
		'--family-expr', 'aff_het:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.het && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.2)})',
		'--family-expr', 'aff_hom:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.hom_alt && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.8)})',
		'--trio', 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < 10',

		'-o', slivar_std_vcf])
	slivar_trio.wait()
	##comp het from std vcf
	#slivar compound-hets -v vcfs/$cohort.vcf --sample-field comphet_side --sample-field denovo -p $ped > vcfs/$cohort.ch.vcf
	slivar_cp = subprocess.Popen([slivar, 'compound-hets', '--vcf', slivar_std_vcf, '--ped', ped_file,
				'--sample-field', 'comphet_side', '--sample-field', 'denovo', '-o', slivar_cp_vcf])
	slivar_cp.wait()
	##convert to tsv (add extra annotation?) -- using ANN for vep annotation
	slivar_trio_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', 'ANN',
		'-o', slivar_std_tsv, '-s', 'denovo', '-s', 'x_denovo', '-s', 'recessive', '-s', 'x_recessive',
		 '-s', 'aff_het', '-s', 'aff_hom',
		'-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt', '-g', pli_lookup,
		'-g', omim_lookup, slivar_std_vcf])
	slivar_trio_tsv.wait()
	slivar_ch_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', 'ANN', '-o', slivar_cp_tsv, 
		'-s', 'slivar_comphet', '-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt',
		'-g', pli_lookup,'-g', omim_lookup, slivar_cp_vcf])
	slivar_ch_tsv.wait()
	##if no variants don't run annovar
	var_lc = 0
	with open(slivar_std_tsv, "r") as sst_fh, open(slivar_cp_tsv, "r") as sct_fh:
		for line in sst_fh:
			var_lc += 1
		for line in sct_fh:
			var_lc += 1			
	if var_lc > 2:
		##annotate vcfs with annovar i.e. run_table_annovar
		command = [table_annovar] + av_buildver + [slivar_std_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_std_annovar]
		annovar = subprocess.Popen(command)
		annovar.wait()
	# '''
	##combine slivar files with annovar annotation
	merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final)

def slivar_ad_analysis_on_trio(in_vcf, ped_file):
	slivar_dom_vcf = in_vcf.split('.')[0] + '.slivar_dominant.vcf'
	slivar_dom_tsv = in_vcf.split('.')[0] + '.slivar_dominant.temp.xls'
	slivar_dom_annovar = in_vcf.split('.')[0] + '.slivar_dominant.temp'
	slivar_dom_multi = in_vcf.split('.')[0] + '.slivar_dominant.temp.hg38_multianno.txt'
	slivar_dom_final = in_vcf.split('.')[0] + '_slivar.dominant_analysis.xls'
	##dominant -- not working
	# slivar_dom = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
	# 	'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
	# 	'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
	# 	'--family-expr', 'dominant:fam.every(segregating_dominant)',
	# 	'-o', slivar_dom_vcf])
	# slivar_dom.wait()
	slivar_multi = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
		'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
		'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
		##all affected family member are het or hom
		'--family-expr', 'aff_het:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.het && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.2)})',
		'--family-expr', 'aff_hom:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.hom_alt && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.8)})',
		'-o', slivar_dom_vcf])
	slivar_multi.wait()
	slivar_dom_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', 'ANN', '-o', slivar_dom_tsv, 
		'-s', 'aff_het','-s', 'aff_hom', '-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt',
		'-g', pli_lookup,'-g', omim_lookup, slivar_dom_vcf])
	slivar_dom_tsv.wait()
	##if no variants don't run annovar
	var_lc = 0
	with open(slivar_std_tsv, "r") as sst_fh, open(slivar_cp_tsv, "r") as sct_fh:
		for line in sst_fh:
			var_lc += 1
		for line in sct_fh:
			var_lc += 1			
	if var_lc > 2:
		##annotate vcfs with annovar i.e. run_table_annovar
		command = [table_annovar] + av_buildver + [slivar_dom_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_dom_annovar]
		annovar = subprocess.Popen(command)
		annovar.wait()
	# '''
	##combine slivar files with annovar annotation
	# merge_std_silvar_annovar_dominant(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final)

def slivar_analysis_master(infile):
	with open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			ped = line[0]
			vcf = line[1]
			ped_file = ped + '.ped'
			ad_ped_file = ped + '.ad.ped'
			annotated_norm_vcf = ped + '.norm.vcf.gz'
			##process vcf files, so have normalized and annotated vcf 
			process_annotate_vcf(ped, vcf, annotated_norm_vcf)
			##slivar std analysis
			slivar_std_analysis_on_trio(annotated_norm_vcf, ped_file)
			##dom analysis --not used, incorperated in std analysis
			# slivar_ad_analysis_on_trio(annotated_norm_vcf, ped_file)


def combine_std_files_get_gene_specific(info_file, genelist_1, filename_1, genelist_2, filename_2):
	first_outfile = 'combined.' + filename_1 + '.xls'
	second_outfile = 'combined.' + filename_2 + '.xls'
	all_outfile = 'combined.all_genes.xls'
	with open(info_file, "r") as in_fh, open(first_outfile, "w") as f_fh, open(second_outfile, "w") as s_fh, open(all_outfile, "w") as a_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			ped = line[0]
			std_file = ped + '_slivar.std_analysis.xls'
			if os.path.isfile(std_file):
				with open(std_file, "r") as std_fh:
					lc2 = 0
					for line2 in std_fh:
						line2 = line2.rstrip().split(delim)
						lc2 += 1
						if lc2 == 1:
							if lc == 1:
								f_fh.write(delim.join(line2) + '\n')
								s_fh.write(delim.join(line2) + '\n')
								a_fh.write(delim.join(line2) + '\n')
						else:
							genes = line2[8].split(';')
							a_fh.write(delim.join(line2) + '\n')
							if bool(set(genes) & set(genelist_1)):
							# if genes in genelist_1:
								f_fh.write(delim.join(line2) + '\n')
							if bool(set(genes) & set(genelist_2)):
								s_fh.write(delim.join(line2) + '\n')

def filter_vars_refgene_type(in_file, out_file):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc == 1:
				out_fh.write(delim.join(line) + '\n')
			else:
				var_type = line[0]
				ex_func = line[15]
				if var_type == 'aff_het' and ex_func != '.' and ex_func != 'synonymous SNV' and ex_func != 'unknown':
					out_fh.write(delim.join(line) + '\n')

def count_vars_per_ped(in_file, gene_list, lisa_gene_peds, out_file):
	ped_dict = {}
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc > 1:
				ped_id = line[1]
				genes = line[8].split(';')
				chd_gene, lisa_ped = 'no', 'no'
				##if a chd gene
				if bool(set(genes) & set(gene_list)):
					chd_gene = 'yes'
					chd_gene_names = list(set(genes).intersection(gene_list))
				##if ped has a variant in one of lisa's genes
				if ped_id in lisa_gene_peds:
					lisa_ped = 'yes'
				##make dict
				if ped_id in ped_dict:
					ped_dict[ped_id][1] += 1
					ped_dict[ped_id][3].extend(genes)
					if chd_gene == 'yes':
						ped_dict[ped_id][2] += 1
						ped_dict[ped_id][4].extend(chd_gene_names)

				else:
					if chd_gene == 'yes':
						ped_dict[ped_id] = [lisa_ped, 1, 1, genes, chd_gene_names]
					else:
						ped_dict[ped_id] = [lisa_ped, 1, 0, genes, []]
	# for p in ped_dict:
	# 	print(p, ped_dict[p])
	with open(out_file, "w") as out_fh:
		# header = ['ped', 'has var in gene', 'var count', 'chd var count', 'gene names']
		header = ['ped', 'has var in gene', 'var count', 'chd var count', 'CHD gene names']
		out_fh.write(delim.join(header) + '\n')
		for p in ped_dict:
			# line_out = [p, ped_dict[p][0], str(ped_dict[p][1]), str(ped_dict[p][2]), '_'.join(ped_dict[p][3])]
			line_out = [p, ped_dict[p][0], str(ped_dict[p][1]), str(ped_dict[p][2]), ', '.join(ped_dict[p][4])]
			out_fh.write(delim.join(line_out) + '\n')





##run methods
working_dir = '/home/atimms/ngs_data/genomes/dave_lisa_chd_kf_genomes_0321'
os.chdir(working_dir)
project_name = 'chd_genomes_0321'
pedigree_file_info = project_name + '.ped_file_info.txt'
ped_info_file = project_name + '.ped_info.txt' ##all peds

##make 2x ped files per family (ad peds not used)
# make_ped_files(pedigree_file_info)

##check relationships -- plink not working properly, so use peddy
# ped_info_file = 'genome_info.2.txt'
# plink_relatadness_check(ped_info_file)
# peddy_relatadness_check(ped_info_file)

##run analysis
# ped_info_file = 'genomes_info_aa'
# ped_info_file = 'genomes_info_ab'
# ped_info_file = 'genomes_info_ac'
# ped_info_file = 'genomes_info_ad'
# ped_info_file = 'genomes_info_ae'
# slivar_analysis_master(ped_info_file)

##combined files and get gene specific vars
lisa_genes = ['POMP', 'PSMD6', 'PSMA6', 'PSMD3', 'PSMA7']
lisa_filename = 'lisa_genes'
chd_genes = ['ABL1', 'ACTC1', 'ACVR1', 'ACVR2B', 'ADAMTS10', 'AFF4', 'ANKRD11', 'ARID1A', 'ARID1B', 'B3GAT3', 'BCOR', 
		'BMPR2', 'BRAF', 'CDK13', 'CFC1', 'CHD4', 'CHD7', 'CHST14', 'CITED2', 'CREBBP', 'CRELD1', 'DLL4', 'DNAH11', 
		'DOCK6', 'EFTUD2', 'EHMT1', 'ELN', 'EP300', 'ESCO2', 'EVC', 'EVC2', 'FBN1', 'FGFR2', 'FLNA', 'FLT4', 'FOXC1', 
		'FOXC2', 'FOXH1', 'FOXP1', 'GATA4', 'GATA5', 'GATA6', 'GDF1', 'GJA1', 'GLI3', 'GPC3', 'HAND1', 'HAND2', 'HDAC8', 
		'HNRNPK', 'HRAS', 'INVS', 'JAG1', 'KANSL1', 'KAT6A', 'KAT6B', 'KDM6A', 'KMT2A', 'KMT2D', 'KRAS', 'KYNU', 'MAP2K1', 
		'MAP2K2', 'MAP3K7', 'MED12', 'MED13L', 'MEIS2', 'MESP1', 'MYBPC3', 'MYH11', 'MYH6', 'MYH7', 'NF1', 'NIPBL', 
		'NKX2-5', 'NKX2-6', 'NODAL', 'NONO', 'NOTCH1', 'NOTCH2', 'NPHP3', 'NPHP4', 'NR2F2', 'NRAS', 'NSD1', 'NUP188', 
		'PBX1', 'PIGL', 'PIGV', 'PITX2', 'PKD1L1', 'PRDM6', 'PRKD1', 'PTPN11', 'RAB23', 'RAD21', 'RAF1', 'RBFOX2', 'RERE', 'RIT1', 'SALL1', 'SALL4', 'SF3B4', 'SHOC2', 'SMAD2', 'SMAD3', 'SMAD4', 'SMAD6', 'SMARCA4', 'SMARCB1', 'SMARCE1', 'SMC1A', 'SMC3', 'SMG9', 'SON', 'SOS1', 'STRA6', 'TAB2', 'TBX1', 'TBX20', 'TBX5', 'TFAP2B', 'TGFBR1', 'TGFBR2', 'TLL1', 'TRAF7', 'TXNL4A', 'UBR1', 'WASHC5', 'ZEB2', 'ZFPM2', 'ZIC3']
chd_filename = 'chd_genes'
# combine_std_files_get_gene_specific(ped_info_file, lisa_genes, lisa_filename, chd_genes, chd_filename)

##look for burden/other genes in relation to lisa's gene
##filter vars i.e. only keep aff_hom and in ref gene
all_vars_file = 'combined.all_genes.xls'
all_vars_file_filtered = 'combined.all.refGene_affHet.xls'
lisa_vars_file = 'combined.lisa_genes.xls'
lisa_vars_file_filtered = 'combined.lisa_genes.refGene_affHet.xls'
chd_vars_file = 'combined.chd_genes.xls'
chd_vars_file_filtered = 'combined.chd_genes.refGene_affHet.xls'
# filter_vars_refgene_type(chd_vars_file, chd_vars_file_filtered) 
# filter_vars_refgene_type(lisa_vars_file, lisa_vars_file_filtered) 
# filter_vars_refgene_type(all_vars_file, all_vars_file_filtered) 

##get counts/gene info
lisa_genes_peds = ['FM_3H7YG3NN', 'FM_5B9CHB1E', 'FM_5PV0WQTJ', 'FM_6X9F8KJM', 'FM_BCG60FZS', 'FM_CNJVXHYD', 'FM_EE5MNZ82', 
	'FM_EHMVXHZ5', 'FM_G2H13CNH', 'FM_KHACJ1N0', 'FM_MQTN1D84', 'FM_PB8EYGK1', 'FM_PVP8RFNR', 'FM_RPY6F5F8', 'FM_RY1NW5R1', 
	'FM_TCQPAHMZ']
counts_per_ped_file = 'counts_per_ped.refGene_affHet.xls'
count_vars_per_ped(all_vars_file_filtered, chd_genes, lisa_genes_peds, counts_per_ped_file)


