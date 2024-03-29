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
or use daniela_genomes.pbs i.e.
qsub -q cdbrmq daniela_genomes.pbs
'''

##set input variables and parameters
delim = '\t'
threads = '5'

##programs
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
slivar = '/home/atimms/programs/slivar_1221/slivar'
# slivar_functions = '/home/atimms/programs/slivar_1221/slivar-0.2.7/js/slivar-functions_adj.js'
slivar_functions = '/home/atimms/programs/slivar_1221/slivar-0.2.7/js/slivar-functions.js'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
freebayes = '/home/atimms/programs/freebayes-1.3.4-linux-static-AMD64'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
plink = '/home/atimms/programs/plink'
dbsnp_common_vcf = '/home/atimms/ngs_data/references/hg19/common_all_20180423.vcf.gz'
tabix = '/home/atimms/programs/samtools-1.11/bin/tabix'
somalier = '/home/atimms/programs/somalier_0721/somalier'

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
somalier_sites_vcf = '/home/atimms/ngs_data/references/somalier/sites.GRCh37.vcf.gz'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,gnomad211_exome,gnomad211_genome,dbnsfp41a,clinvar_20210123,generic', '-genericdbfile', 'Geno2MP.avinput']
av_operation = ['-operation', 'g,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput' ]


##run methods
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
				extra_header = line[8:10] + line[96:100] + ['Geno2MP'] + [line[69]] + line[83:85] + line[47:49] + line[50:52] + [line[58]] + line[10:44] + ['vcf_info', 'vcf_format'] + vcf_samples
			else:
				chr_pos_ref_alt = ':'.join(line[104:106] + line[107:109])
				##keep refgene, clivar/Geno2MP, cadd, gerp, polyphen, mutataster/mutass,REVEL, extra gnomad, vcf info
				info_to_keep = line[8:10] + line[96:101] + [line[69]] + line[83:85] + line[47:49]  + line[50:52] + [line[58]] + line[10:44] + line[111:]
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
	##decomposing and subsetting vcf... not needed
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

def slivar_std_analysis_on_trio(in_vcf, ped_file):
	slivar_std_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.vcf'
	slivar_cp_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.vcf'
	slivar_std_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.xls'
	slivar_cp_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.xls'
	slivar_std_annovar = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp'
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.hg19_multianno.txt'
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
	

def slivar_std_analysis_on_duo(in_vcf, ped_file):
	slivar_std_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.vcf'
	slivar_cp_vcf = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.vcf'
	slivar_std_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.xls'
	slivar_cp_tsv = ped_file.rsplit('.', 1)[0] + '.slivar_comphet.temp.xls'
	slivar_std_annovar = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp'
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.hg19_multianno.txt'
	slivar_std_final = ped_file.rsplit('.', 1)[0] + '.std_analysis.xls'
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
	slivar_std_multi = ped_file.rsplit('.', 1)[0] + '.slivar_std.temp.hg19_multianno.txt'
	slivar_std_final = ped_file.rsplit('.', 1)[0] + '.std_analysis.xls'
	##get std file
	# '''
	slivar_multi = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
		'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
		'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
		##all affexted family member are het
		'--family-expr', 'aff_het:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.het && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.2)})',
		'--family-expr', 'aff_hom:fam.every(function(s) {return (!s.affected && s.GQ >= 15 && s.DP >= 10) || (s.affected && s.hom_alt && s.GQ >= 15 && s.DP >= 10 && s.AB > 0.8)})',
		'--family-expr', 'comphet_side:fam.every(function(s) {return (s.het || s.hom_ref) && s.GQ >= 15 && s.DP >= 10 })',

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



def standard_slivar_protocol(ped_info_file, project_name, combined_vcf, combined_ped):
	##ped types for slivar analysis
	trio_types = ['trio', 'quad', 'trio_with_sib', 'trio_ms', 'quint']
	duo_single_types = ['duo', 'duo_ms', 'parent_sibship', 'singleton', 'singleton_ms', 'sibship']
	multiplex_types = ['multiplex']
	##annoate the combined vcf
	vcf_prefix = combined_vcf.rsplit('.',2)[0]
	annotated_vcf = vcf_prefix + '.bcftools.GRCh37_87.vcf.gz'
	'''
	process_annotate_vcf(annotated_vcf, combined_vcf, vcf_prefix)
	'''
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
		'''
		if ped_type == 'trios':
			if len(ped_type_dict['trios']) > 0:
				slivar_std_analysis_on_trio(annotated_vcf, ped_file)
		elif ped_type == 'single_duos':
			if len(ped_type_dict['single_duos']) > 0:
				slivar_std_analysis_on_duo(annotated_vcf, ped_file)
		elif ped_type == 'multiplex':
			if len(ped_type_dict['multiplex']) > 0:
				slivar_std_analysis_on_multiplex(annotated_vcf, ped_file)
		else:
			print(ptype, ' ped type not recognized')
		'''
		##temp for just analyzing multiplex ped
		if ped_type == 'multiplex':
			if len(ped_type_dict['multiplex']) > 0:
				slivar_std_analysis_on_multiplex(annotated_vcf, ped_file)
		# if ped_type.lower() in trio_types or ped_type.lower() in multiplex_types or ped_type.lower() in duo_types:
		# 	slivar_duo_del(annotated_vcfs, ped_file, formatted_ped_type)

def reheader_vcf(in_vcf, out_vcf, rehead_sample_file):
	bcftools_rehead = subprocess.Popen([bcftools, 'reheader', '-s', rehead_sample_file, '-o', out_vcf, in_vcf])
	bcftools_rehead.wait()
	bcftools_index = subprocess.Popen([bcftools, 'index', out_vcf])
	bcftools_index.wait()

def somalier_check(vcf, ped):
	#somalier extract -d extracted/ --sites sites.vcf.gz -f /data/human/g1k_v37_decoy.fa $cohort.vcf.gz
	som_extract = subprocess.Popen([somalier, 'extract', '-d', 'somalier_extract/', '--sites', somalier_sites_vcf, '-f', fasta, vcf])
	som_extract.wait()
	#somalier relate --ped $pedigree extracted/*.somalier
	# som_relate = subprocess.Popen([somalier, 'relate', '--ped', ped, 'somalier_extract/*.somalier'])
	som_relate = subprocess.Popen([somalier, 'relate', '--infer', '--ped', ped, 'somalier_extract/*.somalier'])
	som_relate.wait()



def get_var_counts_per_ped(info_file, in_files):
	##make dict from info file
	info_dict = {}
	with open(info_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			ped = line[0]
			ped_info = line[1:]
			if ped in info_dict:
				print('ped seen multiple times')
			else:
				info_dict[ped] = ped_info

	for in_file in in_files:
		outfile = in_file.rsplit(".", 1)[0] + '.counts.xls'
		##make counts dict
		var_count_dict = {}
		with open(in_file, "r") as std_fh:
			lc2 = 0
			for line2 in std_fh:
				line2 = line2.rstrip().split(delim)
				lc2 += 1
				if lc2 > 1:
					var_type = line2[0]
					ped = line2[1]
					if var_type.startswith('slivar_comphet'):
						var_type = 'slivar_comphet'
					if ped in var_count_dict:
						var_count_dict[ped].append(var_type)
					else:
						var_count_dict[ped] = [var_type]
		##get list of all var types we see with these std files
		all_var_types = []
		for p in var_count_dict:
			all_var_types.extend(var_count_dict[p])
		all_var_types = list(set(all_var_types))
		# print(in_file_suffix, all_var_types)
		##count var types and print out
		with open(outfile, "w") as out_fh:
			out_fh.write(delim.join(['ped', 'ped_type', 'seq_batch'] + all_var_types) + '\n')
			for ped in var_count_dict:
				ped_tb = info_dict[ped]
				counts = []
				for var in all_var_types:
					count = var_count_dict[ped].count(var)
					counts.append(count)
				counts = [str(i) for i in counts]
				# print(ped, ped_tb, all_var_types, counts)
				# print(var_count_dict[ped][0])
				out_fh.write(delim.join([ped] + ped_tb + counts) + '\n')


##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/daniela_combined_analysis_0222'
os.chdir(working_dir)



##names etc
original_vcf = 'luquetti_grc_combined_2.HF.final.vcf.gz'
vcf_sample_names = 'luquetti_grc_combined_2.sample_names.txt'
combined_vcf = 'luquetti_grc_combined_2.HF.sample_ids.vcf.gz'
combined_prefix = 'cfm_combined_0222'
combined_ped = combined_prefix + '.ped'
## 3 columns: ped id, ped type, batch
combined_ped_info = combined_prefix + '.ped_info.txt'
##result files 
combined_var_files = 'cfm_combined_0222.multiplex.std_analysis.xls', 'cfm_combined_0222.single_duos.std_analysis.xls', 'cfm_combined_0222.trios.std_analysis.xls'


##change sample names in original vcf
# reheader_vcf(original_vcf, combined_vcf, vcf_sample_names)

##check sex and ibd sharing using somalier
# somalier_check(combined_vcf, combined_ped)

##run slivar after splitting ped and vcf file
standard_slivar_protocol(combined_ped_info, combined_prefix, combined_vcf, combined_ped)

##counts per ped
get_var_counts_per_ped(combined_ped_info, combined_var_files)

##additional filtering etc






