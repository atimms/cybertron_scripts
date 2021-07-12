#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
requires:
qsub -Iq cdbrmq -l mem=30gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
#java 1.8 for snpeff
module load java/1.8.0_202 
#biobuilds for tabix etc
module load biobuilds/2017.11
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

##files
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
ens_gff = ref_dir + 'Homo_sapiens.GRCh37.87.gff3.gz'
exome_capture_bed = ref_dir + 'dobyns_exome.in_any_padded_target.1015.bed'
gnomad_gnotate = '/home/atimms/ngs_data/references/slivar/gnomad.hg37.zip'
#wget -qO - https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | zcat | cut -f 1,21,24 | tail -n+2 | awk '{ printf("%s\tpLI=%.3g;oe_lof=%.5g\n", $1, $2, $3)}' > pli.lookup
pli_lookup = '/home/atimms/ngs_data/references/slivar/pli.lookup'
#cut -f2,3,4,6,8 dbdb.gene_association.txt | awk 'BEGIN{FS="\t"}{ printf("%s\tinheritance=%s;phenotype=%s;syndrome=%s;loe=%s\n", $1, $2, $3, $4, $5)}' > ~/ngs_data/references/slivar/dbdb.lookup
dbdb_lookup = '/home/atimms/ngs_data/references/slivar/dbdb.lookup'
#grep -v '^#' omim_genemap2_0321.txt | cut -f9,13 | grep -v '^\s' > /home/atimms/ngs_data/references/slivar/omim.lookup
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

def process_annotate_vcf(ped, input_vcf):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	temp_vcf = ped + '.temp1.vcf.gz'
	# normalized_vcf = ped + '.int.norm.vcf'
	bcftools_vcf = ped + '.bcftools.GRCh37_87.vcf.gz'
	snpeff1_vcf = ped + '.snpeff.GRCh37_87.vcf'
	snpeff2_vcf = ped + '.snpeff.hg19.vcf'
	snpeff3_vcf = ped + '.snpeff.hg19kg.vcf'
	##decomposing and subsetting vcf
	zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
	sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
	#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
	bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', temp_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
	bcftools_norm.wait()
	##annotate with bcftools, have to use unphased version
	#bcftools csq -f hs37d5.fa -g Homo_sapiens.GRCh37.82.gff3.gz in.vcf -Ob -o out.bcf
	bcftools_csq = subprocess.Popen([bcftools, 'csq', '-g', ens_gff, '-f', fasta, '-o', bcftools_vcf, '-O', 'z', '--local-csq', '-s', '-', temp_vcf])
	bcftools_csq.wait()
	tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', bcftools_vcf])
	tabix_vcf.wait()
	##annotate with snpeff and compress and index (try 2x references) -- just use bcftools
	'''
	with open (snpeff1_vcf, 'w') as se_fh:
		snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'GRCh37.75', '-noStats', temp_vcf], stdout=se_fh)
		snpeff_vcf.wait()
	with open (snpeff2_vcf, 'w') as se_fh:
		snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'hg19', '-noStats', temp_vcf], stdout=se_fh)
		snpeff_vcf.wait()
	##weird reference, so don't use
	# with open (snpeff3_vcf, 'w') as se_fh:
	# 	snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'hg19kg', '-noStats', temp_vcf], stdout=se_fh)
	# 	snpeff_vcf.wait()
	for normalized_vcf in [snpeff1_vcf, snpeff2_vcf]:
		bgzip_vcf = subprocess.Popen(['bgzip', normalized_vcf])
		bgzip_vcf.wait()
		tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', normalized_vcf + '.gz'])
		tabix_vcf.wait()
	'''
	##remove intermediate files
	# os.remove(temp_vcf)


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
				extra_header = line[8:10] + line[96:100] + ['Geno2MP'] + [line[69]] + line[83:85] + line[47:49] + line[50:52] + [line[58]] + line[10:44] + ['vcf_info']
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

def merge_dn_silvar_annovar(std_tsv, multianno, outfile):
	##make dict from annovar multianno
	multianno_dict = {}
	with open(multianno, "r") as ma_fh:
		lc = 0
		for line in ma_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc == 1:
				##add header info
				extra_header = line[8:10] + line[96:100] + ['Geno2MP'] + [line[69]] + line[83:85] + line[47:49] + line[50:52] + [line[58]] + line[10:44] + ['vcf_info']
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



def slivar_std_analysis_on_trio(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg19_multianno.txt'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.xls'
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
		# '''
		##different params for diff annotations
		ann_program = in_vcf.split('.')[2]
		if ann_program == 'bcftools':
			c_param = 'BCSQ'
		elif ann_program == 'snpeff':
			c_param = 'ANN'
		else:
			print(ann_program, 'not recognized as annotation program')
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
		else:
			with open(slivar_std_final, "w") as out_fh:
				with open(slivar_std_tsv, "r") as std_fh:
					lc = 0
					for line_1 in std_fh:	
						line_1 = line_1.rstrip('\n').split(delim)
						lc += 1
						if lc == 1:
							header = line_1[:13] + ['pLI', 'dbdb', 'omim']
							out_fh.write(delim.join(header) + '\n')		

def slivar_std_analysis_on_trio_gatk_all(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg19_multianno.txt'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.not_filtered_analysis.xls'
		##get std file
		# '''
		slivar_trio = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.ALT[0] != "*"',
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
		# '''
		##different params for diff annotations
		ann_program = in_vcf.split('.')[2]
		if ann_program == 'bcftools':
			c_param = 'BCSQ'
		elif ann_program == 'snpeff':
			c_param = 'ANN'
		else:
			print(ann_program, 'not recognized as annotation program')
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
		else:
			with open(slivar_std_final, "w") as out_fh:
				with open(slivar_std_tsv, "r") as std_fh:
					lc = 0
					for line_1 in std_fh:	
						line_1 = line_1.rstrip('\n').split(delim)
						lc += 1
						if lc == 1:
							header = line_1[:13] + ['pLI', 'dbdb', 'omim']
							out_fh.write(delim.join(header) + '\n')		

def slivar_genic_dn_on_trio(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_dn_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_exonic_dn.vcf'
		slivar_dn_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_dn.temp.xls'
		slivar_dn_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_dn.temp'
		slivar_dn_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_dn.temp.hg19_multianno.txt'
		slivar_dn_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.genic_dn_analysis.xls'
		##get std file
		# '''
		##all exonic de novo, i.e. use genic instead of impactful
		slivar_dn = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.genic && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
			'--family-expr', 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001',
			'-o', slivar_dn_vcf])
		slivar_dn.wait()
		# '''
		##different params for diff annotations
		ann_program = in_vcf.split('.')[2]
		if ann_program == 'bcftools':
			c_param = 'BCSQ'
		elif ann_program == 'snpeff':
			c_param = 'ANN'
		else:
			print(ann_program, 'not recognized as annotation program')
		##convert to tsv (add extra annotation?)
		slivar_trio_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param,
			'-o', slivar_dn_tsv, '-s', 'denovo',
			'-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt', '-g', pli_lookup,
			'-g', dbdb_lookup,'-g', omim_lookup, slivar_dn_vcf])
		slivar_trio_tsv.wait()
		##if no variants don't run annovar
		var_lc = 0
		with open(slivar_dn_tsv, "r") as sst_fh:
			for line in sst_fh:
				var_lc += 1		
		if var_lc > 1:
			##annotate vcfs with annovar i.e. run_table_annovar
			command = [table_annovar] + av_buildver + [slivar_dn_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_dn_annovar]
			annovar = subprocess.Popen(command)
			annovar.wait()
			##combine slivar files with annovar annotation
			merge_dn_silvar_annovar(slivar_dn_tsv, slivar_dn_multi, slivar_dn_final)
		else:
			##make dummy file with just header is no var
			with open(slivar_dn_final, "w") as out_fh:
				with open(slivar_dn_tsv, "r") as std_fh:
					lc = 0
					for line_1 in std_fh:	
						line_1 = line_1.rstrip('\n').split(delim)
						lc += 1
						if lc == 1:
							header = line_1[:13] + ['pLI', 'dbdb', 'omim']
							out_fh.write(delim.join(header) + '\n')		

def slivar_std_analysis_on_duo(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg19_multianno.txt'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.xls'
		##get std file
		# '''
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
		ann_program = in_vcf.split('.')[2]
		##test, not used as standard
		# if ann_program == 'bcftools':
		# 	c_param = 'BCSQ'
		# elif ann_program == 'snpeff':
		# 	c_param = 'ANN'
		# else:
		# 	print(ann_program, 'not recognized as annotation program')
		# ##get all vars from std vcf, for checking comphet
		# tmp_xls = in_vcf.split('.')[0] + '.temp,xls'
		# slivar_trio_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param,
		# 	'-o', tmp_xls, '-s', 'potential_denovo', '-s', 'comphet_side',
		# 	'-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt', '-g', pli_lookup,
		# 	'-g', dbdb_lookup,'-g', omim_lookup, slivar_std_vcf])
		# slivar_trio_tsv.wait()

		# '''
		##different params for diff annotations
		ann_program = in_vcf.split('.')[2]
		if ann_program == 'bcftools':
			c_param = 'BCSQ'
		elif ann_program == 'snpeff':
			c_param = 'ANN'
		else:
			print(ann_program, 'not recognized as annotation program')
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
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final)

def slivar_std_analysis_on_duo_gatk_all(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg19_multianno.txt'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.not_filtered_analysis.xls'
		##get std file
		# '''
		slivar_duo = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.ALT[0] != "*"',
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
		ann_program = in_vcf.split('.')[2]
		# '''
		##different params for diff annotations
		ann_program = in_vcf.split('.')[2]
		if ann_program == 'bcftools':
			c_param = 'BCSQ'
		elif ann_program == 'snpeff':
			c_param = 'ANN'
		else:
			print(ann_program, 'not recognized as annotation program')
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
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final)

def slivar_std_analysis_on_multiplex(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg19_multianno.txt'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.xls'
		##get std file
		'''
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
		'''
		ann_program = in_vcf.split('.')[2]
		##different params for diff annotations
		ann_program = in_vcf.split('.')[2]
		if ann_program == 'bcftools':
			c_param = 'BCSQ'
		elif ann_program == 'snpeff':
			c_param = 'ANN'
		else:
			print(ann_program, 'not recognized as annotation program')
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
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final)
		# '''

def slivar_std_analysis_on_multiplex_gatk_all(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg19_multianno.txt'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.not_filtered_analysis.xls'
		##get std file
		# '''
		slivar_multi = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.ALT[0] != "*"',
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
		# '''
		ann_program = in_vcf.split('.')[2]
		##different params for diff annotations
		ann_program = in_vcf.split('.')[2]
		if ann_program == 'bcftools':
			c_param = 'BCSQ'
		elif ann_program == 'snpeff':
			c_param = 'ANN'
		else:
			print(ann_program, 'not recognized as annotation program')
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
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final)
		# '''

def slivar_duo_del(in_vcfs, ped_file, ped_type):
	##make dict from decipher data
	decipher_dict = {}
	with open(decipher_delevopmental_file, "r") as decipher_fh:
		for line in decipher_fh:
			line = line.rstrip().split(',') 
			genename = line[0]
			data = line
			decipher_dict[genename] = ','.join(data)
	# for g in decipher_dict:
	# 	print(g, decipher_dict[g])
	##go through 
	for in_vcf in in_vcfs:
		out_file = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.duo_del.xls'
		temp1_file = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.duo_del_temp1.xls'
		temp1_bed = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.duo_del_temp1.bed'
		temp2_bed = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.duo_del_temp2.bed'
		with open(temp1_file, "w") as t1_fh:
			run_slivar_duo_del = subprocess.Popen([slivar, 'duo-del', '--ped', ped_file, '--exclude', hg19_selfchain,
			'-a', in_vcf], stdout=t1_fh)
			run_slivar_duo_del.wait()
		##make bed from duo_del
		with open(temp1_file, "r") as t1_fh, open(temp1_bed, "w") as tb_fh:
			lc = 0
			for line in t1_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc > 1:
					tb_fh.write(delim.join(line[:3]) + '\n')
		##use bedtools to add genename 
		with open(temp2_bed, "w") as out_fh:
			manta_config = subprocess.Popen([bedtools, 'intersect', '-a', temp1_bed, '-b', hg19_refgene_exons, '-wa', '-wb'], stdout=out_fh)
			manta_config.wait()	
		##make dict from bed
		bed_dict = {}
		with open(temp2_bed, "r") as t2b_fh:
			for line in t2b_fh:
				line = line.rstrip().split(delim)
				cse = '_'.join(line[:3])
				gene = line[6].split('_')[0]
				if cse in bed_dict:
					if gene not in bed_dict[cse][0]:
						bed_dict[cse][0].append(gene)
						if gene in decipher_dict:
							bed_dict[cse][1].append(decipher_dict[gene])

				else:
					if gene in decipher_dict:
						bed_dict[cse] = [[gene],[decipher_dict[gene]]]
					else:
						bed_dict[cse] = [[gene],[]]

		##now combine all together
		with open(temp1_file, "r") as t1_fh, open(out_file, "w") as out_fh:
			lc = 0
			for line in t1_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc == 1:
					out_fh.write(delim.join(line + ['genes', 'decipher_development_phenotypes']) + '\n')
				else:
					cse = '_'.join(line[:3])
					if cse in bed_dict:
						line_out = line + ['; '.join(bed_dict[cse][0]), '; '.join(bed_dict[cse][1])]
						out_fh.write(delim.join(line_out) + '\n')
					else:
						out_fh.write(delim.join(line + ['.', '.']) + '\n')





def standard_slivar_protocol(working_directory, pedigree, ped_type):
	os.chdir(working_directory)
	in_vcf = pedigree + '.intersected_vcfs/0002.vcf.gz'
	gatk_vcf = pedigree + '.gatkHC.vcf.gz'
	ped_file = pedigree + '.ped'
	# annotated_vcfs = [pedigree + '.intersected.bcftools.GRCh37_87.vcf.gz', pedigree + '.intersected.snpeff.GRCh37_87.vcf.gz',
	# 		pedigree + '.intersected.snpeff.hg19.vcf.gz', pedigree + '.gatk.bcftools.GRCh37_87.vcf.gz', 
	# 		pedigree + '.gatk.snpeff.GRCh37_87.vcf.gz', pedigree + '.gatk.snpeff.hg19.vcf.gz']
	annotated_vcfs = [pedigree + '.intersected.bcftools.GRCh37_87.vcf.gz', pedigree + '.gatk.bcftools.GRCh37_87.vcf.gz']
	##process vcf files, so have normalized and annotated vcf 
	'''
	process_annotate_vcf(pedigree + '.intersected', in_vcf)
	process_annotate_vcf(pedigree + '.gatk', gatk_vcf)
	# '''
	##filter with slivar
	trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
	duo_types = ['duo', 'duo*', 'parent_sibship']
	single_types = ['singleton', 'singleton*', 'sibship']
	multiplex_types = ['multiplex']
	##testing other ped types with parents
	if '*' in ped_type:
		formatted_ped_type = ped_type.replace('*', 's')
	else:
		formatted_ped_type = ped_type
	# print(pedigree, ped_type, formatted_ped_type)
	'''
	if ped_type.lower() in trio_types:
		slivar_std_analysis_on_trio(annotated_vcfs, ped_file, formatted_ped_type)
	elif ped_type.lower() in duo_types or ped_type.lower() in single_types:
		slivar_std_analysis_on_duo(annotated_vcfs, ped_file, formatted_ped_type)
	elif ped_type.lower() in multiplex_types:
		slivar_std_analysis_on_multiplex(annotated_vcfs, ped_file, formatted_ped_type)
	else:
		print(ped_type, ' ped type not recognized')
	# '''
	if ped_type.lower() in trio_types or ped_type.lower() in multiplex_types or ped_type.lower() in duo_types:
		slivar_duo_del(annotated_vcfs, ped_file, formatted_ped_type)



def slivar_exome_inputs(work_directory, infile):
	os.chdir(work_directory)
	with open(infile, "r") as in_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc > 1:
				ped = line[0]
				w_dir = line[1]
				ped_type = line[2]
				print(w_dir, ped, ped_type)
				standard_slivar_protocol(w_dir, ped, ped_type)


def get_var_counts_per_ped(work_directory, info_file, in_file_suffices):
	os.chdir(work_directory)
	for in_file_suffix in in_file_suffices:
		outfile = 'var_count_summary' + in_file_suffix
		var_count_dict = {}
		with open(info_file, "r") as in_fh:
			lc = 0
			for line in in_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc > 1:
					ped = line[0]
					ped_type = line[2]
					std_file = ped + in_file_suffix
					with open(std_file, "r") as std_fh:
						lc2 = 0
						for line2 in std_fh:
							line2 = line2.rstrip().split(delim)
							lc2 += 1
							if lc2 > 1:
								if 'slivar' in in_file_suffix:
									var_type = line2[0]
								else:
									var_type = line2[-1]
								if var_type.startswith('slivar_comphet'):
									var_type = 'slivar_comphet'
								if ped in var_count_dict:
									var_count_dict[ped][0].append(var_type)
								else:
									var_count_dict[ped] = [[var_type], ped_type]
		##get list of all var types we see with these std files
		all_var_types = []
		for p in var_count_dict:
			all_var_types.extend(var_count_dict[p][0])
		all_var_types = list(set(all_var_types))
		# print(in_file_suffix, all_var_types)
		##count var types and print out
		with open(outfile, "w") as out_fh:
			out_fh.write(delim.join(['ped', 'ped_type'] + all_var_types) + '\n')

			for ped in var_count_dict:
				p_type =  var_count_dict[ped][1]
				counts = []
				for var in all_var_types:
					count = var_count_dict[ped][0].count(var)
					counts.append(count)
				counts = [str(i) for i in counts]
				print(ped, p_type, all_var_types, counts)
				# print(var_count_dict[ped][0])
				out_fh.write(delim.join([ped, p_type] + counts) + '\n')

def get_var_counts(work_directory, info_file, analysis_types):
	os.chdir(work_directory)
	for analysis_type in analysis_types:
		count_outfile = 'var_counts' + analysis_type + 'xls'
		var_count_dict = {}
		with open(info_file, "r") as in_fh:
			lc = 0
			for line in in_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc > 1:
					ped = line[0]
					ped_type = line[2]
					if '*' in ped_type:
						formatted_ped_type = ped_type.replace('*', 's')
					else:
						formatted_ped_type = ped_type
					std_file = ped + analysis_type + formatted_ped_type + '.std_analysis.xls'
					if os.path.isfile(std_file):
						with open(std_file, "r") as std_fh:
							lc2 = 0
							for line2 in std_fh:
								line2 = line2.rstrip().split(delim)
								lc2 += 1
								if lc2 > 1:
									if 'slivar' in std_file:
										var_type = line2[0]
									else:
										var_type = line2[-1]
									if var_type.startswith('slivar_comphet'):
										var_type = 'slivar_comphet'
									if ped in var_count_dict:
										var_count_dict[ped][0].append(var_type)
									else:
										var_count_dict[ped] = [[var_type], ped_type]
		##get list of all var types we see with these std files
		all_var_types = []
		for p in var_count_dict:
			all_var_types.extend(var_count_dict[p][0])
		all_var_types = list(set(all_var_types))
		# print(in_file_suffix, all_var_types)
		##count var types and print out
		with open(count_outfile, "w") as cout_fh:
			cout_fh.write(delim.join(['ped', 'ped_type'] + all_var_types) + '\n')

			for ped in var_count_dict:
				p_type =  var_count_dict[ped][1]
				counts = []
				for var in all_var_types:
					count = var_count_dict[ped][0].count(var)
					counts.append(count)
				counts = [str(i) for i in counts]
				print(ped, p_type, all_var_types, counts)
				# print(var_count_dict[ped][0])
				cout_fh.write(delim.join([ped, p_type] + counts) + '\n')


def combine_var_counts(work_directory, info_file, analysis_types):
	os.chdir(work_directory)
	for analysis_type in analysis_types:
		single_duo_outfile = 'combined' + analysis_type + 'single_duo.all.xls'
		trio_quad_outfile = 'combined' + analysis_type + 'trio_quad.all.xls'
		multiplex_outfile = 'combined' + analysis_type + 'multiplex.all.xls'
		trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
		single_duo_types = ['duo', 'duo*', 'parent_sibship','singleton', 'singleton*', 'sibship']
		multiplex_types = ['multiplex']
		with open(info_file, "r") as in_fh, open(single_duo_outfile, "w") as sd_fh, open(trio_quad_outfile, "w") as tq_fh, open(multiplex_outfile, "w") as m_fh:
			lc = 0
			for line in in_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc > 1:
					ped = line[0]
					ped_type = line[2]
					if '*' in ped_type:
						formatted_ped_type = ped_type.replace('*', 's')
					else:
						formatted_ped_type = ped_type
					std_file = ped + analysis_type + formatted_ped_type + '.std_analysis.xls'
					if os.path.isfile(std_file):
						with open(std_file, "r") as std_fh:
							lc2 = 0
							for line2 in std_fh:
								line2 = line2.rstrip().split(delim)
								lc2 += 1
								if lc2 == 1:
									if lc == 2:
										sd_fh.write(delim.join(line2) + '\n')
										tq_fh.write(delim.join(line2) + '\n')
										m_fh.write(delim.join(line2) + '\n')
								else:
									if ped_type in trio_types:
										tq_fh.write(delim.join(line2) + '\n')
									elif ped_type in single_duo_types:
										sd_fh.write(delim.join(line2) + '\n')
									elif ped_type in multiplex_types:
										m_fh.write(delim.join(line2) + '\n')
									else:
										print(ped_type, ' not recognized as a ped type')

def combine_unfiltered_gatk_slivar_files(work_directory, info_file, file_suffix):
	os.chdir(work_directory)
	single_duo_outfile = 'combined.gatk38_not_filtered.single_duo.all.xls'
	trio_quad_outfile = 'combined.gatk38_not_filtered.trio_quad.all.xls'
	multiplex_outfile = 'combined.gatk38_not_filtered.multiplex.all.xls'
	trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
	single_duo_types = ['duo', 'duo*', 'parent_sibship','singleton', 'singleton*', 'sibship']
	multiplex_types = ['multiplex']
	with open(info_file, "r") as in_fh, open(single_duo_outfile, "w") as sd_fh, open(trio_quad_outfile, "w") as tq_fh, open(multiplex_outfile, "w") as m_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc > 1:
				ped = line[0]
				ped_type = line[2]
				if '*' in ped_type:
					formatted_ped_type = ped_type.replace('*', 's')
				else:
					formatted_ped_type = ped_type
				std_file = ped + '.gatk38_slivar.' + formatted_ped_type + file_suffix
				if os.path.isfile(std_file):
					with open(std_file, "r") as std_fh:
						lc2 = 0
						for line2 in std_fh:
							line2 = line2.rstrip().split(delim)
							lc2 += 1
							if lc2 == 1:
								if lc == 2:
									sd_fh.write(delim.join(line2) + '\n')
									tq_fh.write(delim.join(line2) + '\n')
									m_fh.write(delim.join(line2) + '\n')
							else:
								if ped_type in trio_types:
									tq_fh.write(delim.join(line2) + '\n')
								elif ped_type in single_duo_types:
									sd_fh.write(delim.join(line2) + '\n')
								elif ped_type in multiplex_types:
									m_fh.write(delim.join(line2) + '\n')
								else:
									print(ped_type, ' not recognized as a ped type')


def filter_clinsig_vars(work_directory, in_files, out_suffix, path_definitions):
	os.chdir(work_directory)
	for in_file in in_files:
		out_file = in_file.rsplit('.', 2)[0] + out_suffix
		with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
			lc = 0
			for line in in_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc == 1:
					out_fh.write(delim.join(line) + '\n')
				else:
					clinsig = line[21]
					if clinsig in path_definitions:
						out_fh.write(delim.join(line) + '\n')


def filter_parental_mosiac(work_directory, in_files, out_suffix):
	os.chdir(work_directory)
	for in_file in in_files:
		out_file = in_file.rsplit('.', 2)[0] + out_suffix
		with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
			lc = 0
			for line in in_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc == 1:
					out_fh.write(delim.join(line) + '\n')
				else:
					mode = line[0]
					if mode == 'denovo' or mode == 'x_denovo':
						print(line[11])
						dad_ab = float(line[11].split(',')[1])
						mom_ab = float(line[11].split(',')[2])
						if dad_ab > 0 or mom_ab > 0:
							out_fh.write(delim.join(line) + '\n')						

##make list of all bam files to be analyzed
def make_list_of_bams(bam_list, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for bam in bam_list:
			outfh.write(bam + '\n')

def variant_calling_gatk_hc_no_gvcf(bamlist, name_prefix):
	vcf_temp0 = name_prefix + 'temp_gatk0.vcf'
	vcf_raw_snps = name_prefix + 'temp_raw_snps.vcf'
	vcf_filtered_snps = name_prefix + 'temp_filtered_snps.vcf'
	vcf_raw_indels = name_prefix + 'temp_raw_indels.vcf'
	vcf_filtered_indels = name_prefix + 'temp_filtered_indels.vcf'
	vcf_temp1 = name_prefix + 'temp_gatk1.vcf'
	final_vcf = name_prefix + '.gatk38hc.vcf.gz'
	##run haplotype caller
	gatk_hc = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fasta, '-I', bamlist,'-L', exome_capture_bed, '-o', vcf_temp0])
	gatk_hc.wait()
	##split data in SNPs and indels and apply manual variant filtering (add -L)
	snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-V', vcf_temp0, '-o', vcf_raw_snps, '-selectType', 'SNP'])
	snp_cut.wait()
	snp_vf = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_snps, '-o', vcf_filtered_snps, '--filterExpression ', "QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "snp_filter"])
	snp_vf.wait()
	indel_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fasta, '-V', vcf_temp0, '-o', vcf_raw_indels, '-selectType', 'INDEL'])
	indel_cut.wait()
	indel_vf = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fasta, '-V', vcf_raw_indels, '-o', vcf_filtered_indels, '--filterExpression ', "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0", '--filterName', "indel_filter"])
	indel_vf.wait()
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fasta, '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', vcf_temp1, '-genotypeMergeOptions', 'UNSORTED'])
	combine_var.wait()
	bgzip_cmd = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip_cmd.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	##decomposing and subsetting vcf - AD looks good with gatk38, so sed not needed
	# zless_vcf = subprocess.Popen(['zless', vcf_temp1 + '.gz'], stdout=subprocess.PIPE)
	# sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
	#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
	# bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', final_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
	bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', final_vcf, '-O', 'z', vcf_temp1 + '.gz'])
	bcftools_norm.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()	

##call freebayes in all sample listed in list of bamfiles
def variant_calling_freebayes(bamlist, name_prefix):
	vcf_temp1 = name_prefix + 'temp_fb1.vcf'
	final_vcf = name_prefix + '.freebayes13.vcf.gz'
	freebayes_run = subprocess.Popen([freebayes,'-f' , fasta, '-L', bamlist, '-t', exome_capture_bed, '-v', vcf_temp1])
	freebayes_run.wait()
	##decomposing and subsetting vcf
	zless_vcf = subprocess.Popen(['cat', vcf_temp1], stdout=subprocess.PIPE)
	sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
	#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
	bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', final_vcf, '-O', 'z'], stdin=sed_vcf.stdout)
	bcftools_norm.wait()
	bcf_index = subprocess.Popen([bcftools, 'index', final_vcf])
	bcf_index.wait()	

def intersect_two_vcf_files(vcf1, vcf2, output_dir):
	bcf_isec = subprocess.Popen([bcftools, 'isec', vcf1, vcf2, '-p', output_dir ])
	bcf_isec.wait()

def standard_slivar_protocol_v2(working_directory, pedigree, ped_type, gatk_vcf, int_vcf):
	os.chdir(working_directory)
	ped_file = pedigree + '.ped'
	annotated_vcfs = [pedigree + '.intersect.bcftools.GRCh37_87.vcf.gz', pedigree + '.gatk38.bcftools.GRCh37_87.vcf.gz']
	##process vcf files, so have normalized and annotated vcf 
	'''
	process_annotate_vcf(pedigree + '.intersect', int_vcf)
	process_annotate_vcf(pedigree + '.gatk38', gatk_vcf)
	'''
	##filter with slivar
	trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
	duo_types = ['duo', 'duo*', 'parent_sibship']
	single_types = ['singleton', 'singleton*', 'sibship']
	multiplex_types = ['multiplex']
	##change ped type for final file name
	if '*' in ped_type:
		formatted_ped_type = ped_type.replace('*', 's')
	else:
		formatted_ped_type = ped_type
	# print(pedigree, ped_type, formatted_ped_type)
	'''
	##standard analyses
	if ped_type.lower() in trio_types:
		slivar_std_analysis_on_trio(annotated_vcfs, ped_file, formatted_ped_type)
	elif ped_type.lower() in duo_types or ped_type.lower() in single_types:
		slivar_std_analysis_on_duo(annotated_vcfs, ped_file, formatted_ped_type)
	elif ped_type.lower() in multiplex_types:
		slivar_std_analysis_on_multiplex(annotated_vcfs, ped_file, formatted_ped_type)
	else:
		print(ped_type, ' ped type not recognized')
	
	##run duo_del
	if ped_type.lower() in trio_types or ped_type.lower() in multiplex_types or ped_type.lower() in duo_types:
		slivar_duo_del(annotated_vcfs, ped_file, formatted_ped_type)
	
	##get genic denovo vars
	if ped_type.lower() in trio_types:
		slivar_genic_dn_on_trio(annotated_vcfs, ped_file, formatted_ped_type)
	
	##slivar on all gatk vars analyses
	if ped_type.lower() in trio_types:
		slivar_std_analysis_on_trio_gatk_all([annotated_vcfs[1]], ped_file, formatted_ped_type)
	elif ped_type.lower() in duo_types or ped_type.lower() in single_types:
		slivar_std_analysis_on_duo_gatk_all([annotated_vcfs[1]], ped_file, formatted_ped_type)
	elif ped_type.lower() in multiplex_types:
		slivar_std_analysis_on_multiplex_gatk_all([annotated_vcfs[1]], ped_file, formatted_ped_type)
	else:
		print(ped_type, ' ped type not recognized')
	# '''



def combine_files(work_directory, var_types, file_suffix):
	os.chdir(work_directory)
	for var_type in var_types:
		infiles = glob.glob('*' + var_type + '*' + file_suffix)
		outfile = 'combined' + var_type + file_suffix
		infiles.remove(outfile)
		print(var_type, len(infiles))
		# print(infiles)
		fc = 0
		with open(outfile, "w") as out_fh:
			for infile in infiles:
				with open(infile, "r") as in_fh:
					fc += 1
					lc = 0
					for line in in_fh:
						lc += 1
						if lc == 1:
							if fc == 1:
								out_fh.write(line)
						else:
							out_fh.write(line)



def repeat_var_calling_0221(work_directory, infile):
	rss_dir = '/archive/.snapshot/202005281800-Dobyns_W-Archive/dobyns_w/exome_data/all_exome_files/'
	rss_dir2 = '/archive/.snapshot/202005281800-Dobyns_W-Archive/dobyns_w/exome_data/exomes_after_0616/jimmy_trios_0916/'
	rss_dir3 = '/archive/.snapshot/202005281800-Dobyns_W-Archive/dobyns_w/exome_data/batch1_and_2_exomes_2016/acc_fp1/'
	rss_dir4 = '/archive/.snapshot/202005281800-Dobyns_W-Archive/dobyns_w/exome_data/exomes_after_1217/ghayda_gw_0218/'
	rss_dir5 = '/archive/.snapshot/202005281800-Dobyns_W-Archive/dobyns_w/exome_data/exomes_after_1217/ghayda_gw_1117/'
	new_dir = '/archive/mirzaa_g/exomes/exomes_after_0920/ghayda_broad_exomes_0920/'
	os.chdir(work_directory)
	with open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			ped = line[0]
			alt_ped = line[1]
			ped_type = line[2]
			bam_location = line[3]		
			samples = line[4:]
			if bam_location == 'rss':
				bams = [rss_dir + alt_ped + '/' + s + '.bwa_gatk.bam' for s in samples]
			elif bam_location == 'rss2':
				bams = [rss_dir2 + s + '.bwa_gatk.bam' for s in samples]
			elif bam_location == 'rss3':
				bams = [rss_dir3 + s + '.bwa_gatk.bam' for s in samples]
			elif bam_location == 'rss4':
				bams = [rss_dir4 + s + '.bwa_gatk.bam' for s in samples]
			elif bam_location == 'rss5':
				bams = [rss_dir5 + s + '.bwa_gatk.bam' for s in samples]
			elif bam_location == '920':
				bams = [new_dir + s + '.bwa_gatk.bam' for s in samples]
			bams_list = ped + '.bams.list'
			make_list_of_bams(bams, bams_list)
			gatk_vcf = ped + '.gatk38hc.vcf.gz'
			fb_vcf = ped + '.freebayes13.vcf.gz'
			int_vcf = ped + '.intersect_vcfs/0002.vcf'
			##call vars and intersect
			variant_calling_gatk_hc_no_gvcf(bams_list, ped)
			variant_calling_freebayes(bams_list, ped)
			intersect_two_vcf_files(gatk_vcf, fb_vcf, ped + '.intersect_vcfs')
			##run slivar
			standard_slivar_protocol_v2(work_directory, ped, ped_type, gatk_vcf, int_vcf)								

def slivar_analysis_master_0221(work_directory, infile):
	os.chdir(work_directory)
	with open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			ped = line[0]
			ped_type = line[2]
			gatk_vcf = ped + '.gatk38hc.vcf.gz'
			int_vcf = ped + '.intersect_vcfs/0002.vcf'
			##run slivar
			standard_slivar_protocol_v2(work_directory, ped, ped_type, gatk_vcf, int_vcf)


def slivar_std_analysis_on_trio_split(in_vcfs, ped_file, ped_type, ped_name):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg19_multianno.txt'
		slivar_std_final = ped_name + '.' + in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.xls'
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
		# '''
		##different params for diff annotations
		ann_program = in_vcf.split('.')[2]
		if ann_program == 'bcftools':
			c_param = 'BCSQ'
		elif ann_program == 'snpeff':
			c_param = 'ANN'
		else:
			print(ann_program, 'not recognized as annotation program')
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
		else:
			with open(slivar_std_final, "w") as out_fh:
				with open(slivar_std_tsv, "r") as std_fh:
					lc = 0
					for line_1 in std_fh:	
						line_1 = line_1.rstrip('\n').split(delim)
						lc += 1
						if lc == 1:
							header = line_1[:13] + ['pLI', 'dbdb', 'omim']
							out_fh.write(delim.join(header) + '\n')		

def slivar_std_analysis_on_duo_split(in_vcfs, ped_file, ped_type, ped_name):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg19_multianno.txt'
		slivar_std_final = ped_name + '.' + in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.xls'
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
		ann_program = in_vcf.split('.')[2]

		##different params for diff annotations
		ann_program = in_vcf.split('.')[2]
		if ann_program == 'bcftools':
			c_param = 'BCSQ'
		elif ann_program == 'snpeff':
			c_param = 'ANN'
		else:
			print(ann_program, 'not recognized as annotation program')
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
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final)

def standard_split_slivar_protocol(working_directory, ped_file, ped_type, gatk_vcf, int_vcf, new_ped_name):
	os.chdir(working_directory)
	annotated_vcfs = [int_vcf, gatk_vcf]
	##filter with slivar
	trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
	duo_types = ['duo', 'duo*', 'parent_sibship']
	single_types = ['singleton', 'singleton*', 'sibship']
	multiplex_types = ['multiplex']
	##change ped type for final file name
	if '*' in ped_type:
		formatted_ped_type = ped_type.replace('*', 's')
	else:
		formatted_ped_type = ped_type
	# print(pedigree, ped_type, formatted_ped_type)
	# '''
	##standard analyses
	if ped_type.lower() in trio_types:
		slivar_std_analysis_on_trio_split(annotated_vcfs, ped_file, formatted_ped_type, new_ped_name)
	elif ped_type.lower() in duo_types or ped_type.lower() in single_types:
		slivar_std_analysis_on_duo_split(annotated_vcfs, ped_file, formatted_ped_type, new_ped_name)
	else:
		print(ped_type, ' ped type not recognized')



def slivar_analysis_split_ped(work_directory, infile):
	os.chdir(work_directory)
	with open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			ped = line[0]
			ped_split = line[1]
			ped_type = line[2]
			ped_file = ped_split + '.ped'
			gatk_vcf = ped + '.gatk38.bcftools.GRCh37_87.vcf.gz'
			int_vcf = ped + '.intersect.bcftools.GRCh37_87.vcf.gz'
			##run slivar
			standard_split_slivar_protocol(work_directory, ped_file, ped_type, gatk_vcf, int_vcf, ped_split)		

##run on original gatk and intersected files
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalyze_exomes_0121/all_exome_data'

##run as a batch
# input_exomes_file = 'exomes_1.txt'
# input_exomes_file = 'exomes_2.txt'
# input_exomes_file = 'exomes_3.txt'
# input_exomes_file = 'exomes_4.txt'
# input_exomes_file = 'exomes_5.txt'
# input_exomes_file = 'exomes_6.txt'
# input_exomes_file = 'exomes_all.txt'
# input_exomes_file = 'exomes_s2.txt'

##input files
# slivar_exome_inputs(working_dir, input_exomes_file)

##ind peds
# standard_slivar_protocol(working_dir, "LR02-027", 'trio')
# standard_slivar_protocol(working_dir, "LR05-203", 'multiplex')
# standard_slivar_protocol(working_dir, "LP99-100", 'sibship')


#rm *temp* *snpeff* *bcftools*

##get var counts per ped for each version of slivar data
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalyze_exomes_0121/analysis_files_0121'
std_file_suffices = ['.std_analysis.xls', '.gatk.bcftools.GRCh37_87.slivar.std_analysis.xls', 
		'.gatk.snpeff.GRCh37_87.slivar.std_analysis.xls', '.gatk.snpeff.hg19.slivar.std_analysis.xls',
		'.intersected.bcftools.GRCh37_87.slivar.std_analysis.xls', '.intersected.snpeff.GRCh37_87.slivar.std_analysis.xls',
		'.intersected.snpeff.hg19.slivar.std_analysis.xls']
ped_info_file = 'all_exome_info.txt'
# ped_info_file = 'test_exome_info.txt'
# get_var_counts_per_ped(working_dir, ped_info_file, std_file_suffices)

##latest versions of analysis files - 0221
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalyze_exomes_0121/analysis_files_0221'
var_types = ['.gatk_slivar.', '.intersected_slivar.']
##get counts
# get_var_counts(working_dir, ped_info_file, var_types)
##combine all single/duo trio and mutiplex files
# combine_var_counts(working_dir, ped_info_file, var_types)
##filter the combined files
combined_var_files = ['combined.gatk_slivar.multiplex.all.xls', 'combined.gatk_slivar.single_duo.all.xls', 
		'combined.gatk_slivar.trio_quad.all.xls', 'combined.intersected_slivar.multiplex.all.xls', 
		'combined.intersected_slivar.single_duo.all.xls', 'combined.intersected_slivar.trio_quad.all.xls']
pathogenic_definitions = ['Pathogenic', 'Pathogenic/Likely_pathogenic', 'Likely_pathogenic']
# filter_clinsig_vars(working_dir, combined_var_files, '.clinvar.xls', pathogenic_definitions)

##parental mosaic from denovo
trio_files = ['combined.gatk_slivar.trio_quad.all.xls', 'combined.intersected_slivar.trio_quad.all.xls']
# filter_parental_mosiac(working_dir, trio_files, '.parental_mosaic.xls')

##cnv analysis



##repeat var calling/slivar analysis
# exome_ped_file = 'exome_test' #test
exome_ped_file = 'exome_info.txt' ##all the peds
##split 5 ways
# exome_ped_file = 'exome_infoaa'
# exome_ped_file = 'exome_infoab'
# exome_ped_file = 'exome_infoac'
# exome_ped_file = 'exome_infoad'
# exome_ped_file = 'exome_infoae'

##for combining files
ped_info_file = 'all_exome_info.txt'

working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalyze_exomes_0121/recall_vars_0221'
##repeat var calling and slivar
# repeat_var_calling_0221(working_dir, exome_ped_file)

##redo slivar with new omim and repeats analyses etc
# slivar_analysis_master_0221(working_dir, exome_ped_file)

##latest versions of analysis files - 0221
var_types = ['.gatk38_slivar.', '.intersect_slivar.']
##get counts
# get_var_counts(working_dir, ped_info_file, var_types)

##combine all single/duo trio and mutiplex files for std analysis files
# combine_var_counts(working_dir, ped_info_file, var_types)
##also combine the non filtered gatk slivar files
# combine_unfiltered_gatk_slivar_files(working_dir, exome_ped_file, '.not_filtered_analysis.xls')


##filter the combined files
##clivar
combined_var_files = ['combined.gatk38_slivar.multiplex.all.xls', 'combined.gatk38_slivar.single_duo.all.xls', 
		'combined.gatk38_slivar.trio_quad.all.xls', 'combined.intersect_slivar.multiplex.all.xls', 
		'combined.intersect_slivar.single_duo.all.xls', 'combined.intersect_slivar.trio_quad.all.xls', 
		'combined.gatk38_not_filtered.multiplex.all.xls', 'combined.gatk38_not_filtered.single_duo.all.xls', 
		'combined.gatk38_not_filtered.trio_quad.all.xls']
pathogenic_definitions = ['Pathogenic', 'Pathogenic/Likely_pathogenic', 'Likely_pathogenic']
# filter_clinsig_vars(working_dir, combined_var_files, '.clinvar.xls', pathogenic_definitions)

##parental mosaic from denovo
trio_files = ['combined.gatk38_slivar.trio_quad.all.xls', 'combined.intersect_slivar.trio_quad.all.xls']
# filter_parental_mosiac(working_dir, trio_files, '.parental_mosaic.xls')

##combine duo_del files, need to run on all samples first
# combine_files(working_dir, var_types, 'duo_del.xls')

##combine genic de novos
# combine_files(working_dir, var_types, 'genic_dn_analysis.xls')

##split quads/quints etc
##need to split ped and make ped analysis file  i.e. ped, ped split, ped split type
# exome_ped_file = 'exome_split_1.txt'
# exome_ped_file = 'exome_split_2.txt'
# exome_ped_file = 'exome_split_3.txt'
# exome_ped_file = 'exome_split_4.txt'
# exome_ped_file = 'exome_split_5.txt'
# exome_ped_file = 'exome_split_6.txt'
##redo slivar with
# slivar_analysis_split_ped(working_dir, exome_ped_file)


def slivar_std_analysis_on_trio_dom(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_dn_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_dom.vcf'
		slivar_dn_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_dom.temp.xls'
		slivar_dn_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_dom.temp'
		slivar_dn_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_dom.temp.hg19_multianno.txt'
		slivar_dn_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.dominant_analysis.xls'
		##get std file
		'''
		##all dom
		slivar_dn = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
			'--family-expr', 'aff_het:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.het && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.2)})',
			'-o', slivar_dn_vcf])
		slivar_dn.wait()

		'''
		##different params for diff annotations
		ann_program = in_vcf.split('.')[2]
		if ann_program == 'bcftools':
			c_param = 'BCSQ'
		elif ann_program == 'snpeff':
			c_param = 'ANN'
		else:
			print(ann_program, 'not recognized as annotation program')
		##convert to tsv (add extra annotation?)
		slivar_trio_tsv = subprocess.Popen([slivar, 'tsv', '--ped', ped_file, '-c', c_param,
			'-o', slivar_dn_tsv, '-s', 'aff_het',
			'-i', 'gnomad_popmax_af', '-i', 'gnomad_popmax_af_filter', '-i', 'gnomad_nhomalt', '-g', pli_lookup,
			'-g', dbdb_lookup,'-g', omim_lookup, slivar_dn_vcf])
		slivar_trio_tsv.wait()
		##if no variants don't run annovar
		var_lc = 0
		with open(slivar_dn_tsv, "r") as sst_fh:
			for line in sst_fh:
				var_lc += 1		
		if var_lc > 1:
			##annotate vcfs with annovar i.e. run_table_annovar
			command = [table_annovar] + av_buildver + [slivar_dn_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_dn_annovar]
			annovar = subprocess.Popen(command)
			annovar.wait()
			##combine slivar files with annovar annotation
			merge_dn_silvar_annovar(slivar_dn_tsv, slivar_dn_multi, slivar_dn_final)
		else:
			##make dummy file with just header is no var
			with open(slivar_dn_final, "w") as out_fh:
				with open(slivar_dn_tsv, "r") as std_fh:
					lc = 0
					for line_1 in std_fh:	
						line_1 = line_1.rstrip('\n').split(delim)
						lc += 1
						if lc == 1:
							header = line_1[:13] + ['pLI', 'dbdb', 'omim']
							out_fh.write(delim.join(header) + '\n')	

def standard_dom_slivar_protocol(working_directory, ped_file, ped_type, gatk_vcf, int_vcf):
	os.chdir(working_directory)
	annotated_vcfs = [int_vcf, gatk_vcf]
	##filter with slivar
	trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
	duo_types = ['duo', 'duo*', 'parent_sibship']
	single_types = ['singleton', 'singleton*', 'sibship']
	multiplex_types = ['multiplex']
	##change ped type for final file name
	if '*' in ped_type:
		formatted_ped_type = ped_type.replace('*', 's')
	else:
		formatted_ped_type = ped_type
	# print(pedigree, ped_type, formatted_ped_type)
	# '''
	##standard analyses
	if ped_type.lower() in trio_types:
		slivar_std_analysis_on_trio_dom(annotated_vcfs, ped_file, formatted_ped_type)
	else:
		print(ped_type, ' ped type not recognized')


def slivar_dominant_analysis(work_directory, infile):
	os.chdir(work_directory)
	with open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			ped = line[0]
			ped_type = line[1]
			ped_file = ped + '.ped'
			gatk_vcf = ped + '.gatk38.bcftools.GRCh37_87.vcf.gz'
			int_vcf = ped + '.intersect.bcftools.GRCh37_87.vcf.gz'
			##run slivar
			standard_dom_slivar_protocol(work_directory, ped_file, ped_type, gatk_vcf, int_vcf)

##dominant analysis on specific ped/peds
# exome_ped_file = 'dom_analysis1.txt'
exome_ped_file = 'dom_analysis2.txt'
slivar_dominant_analysis(working_dir, exome_ped_file)


