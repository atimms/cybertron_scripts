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
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'
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
# av_protocol = ['-protocol', 'refGene,gnomad211_exome,gnomad211_genome,dbnsfp41a,clinvar_20210123,generic', '-genericdbfile', 'Geno2MP.avinput']
# av_operation = ['-operation', 'g,f,f,f,f,f']
# av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput' ]
av_protocol = ['-protocol', 'refGene']
av_operation = ['-operation', 'g']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput' ]


def merge_std_silvar_annovar(std_tsv, cp_tsv, avinput, outfile):
	##make dict from slivar vars
	var_list = []
	with open(std_tsv, "r") as std_fh:
		lc = 0
		for line_1 in std_fh:	
			line_1 = line_1.rstrip('\n').split(delim)
			lc += 1
			if lc > 1:
				cpra = line_1[3]
				var_list.append(cpra)

	with open(cp_tsv, "r") as cp_fh:
		lc = 0
		for line_2 in cp_fh:	
			line_2 = line_2.rstrip('\n').split(delim)
			lc += 1
			if lc > 1:
				cpra = line_2[3]
				var_list.append(cpra)


	with open(avinput, "r") as av_fh, open(outfile, "w") as out_fh:
		lc = 0
		for line in av_fh:
			line = line.split(delim)
			lc += 1
			if lc > 1:
				chr_pos_ref_alt = ':'.join(line[8:10] + line[11:13])
				if chr_pos_ref_alt in var_list:
					out_fh.write(delim.join(line))




def slivar_std_analysis_on_trio(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_avinput = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.avinput'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.avinput'
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
		# '''
		##if no variants don't run annovar
		var_lc = 0
		with open(slivar_std_tsv, "r") as sst_fh, open(slivar_cp_tsv, "r") as sct_fh:
			for line in sst_fh:
				var_lc += 1
			for line in sct_fh:
				var_lc += 1			
		if var_lc > 2:
			##annotate vcfs with annovar i.e. run_table_annovar just to get avinput file
			command = [table_annovar] + av_buildver + [slivar_std_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_std_annovar]
			annovar = subprocess.Popen(command)
			annovar.wait()
			##get avinput vars if in slivar files
			merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_avinput, slivar_std_final)
		else:
			with open(slivar_std_final, "w") as out_fh:
				pass


def slivar_std_analysis_on_duo(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_avinput = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.avinput'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.avinput'
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
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_avinput, slivar_std_final)

def slivar_std_analysis_on_multiplex(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_avinput = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.avinput'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.avinput'
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
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_avinput, slivar_std_final)
		# '''


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
	print(pedigree, ped_type, formatted_ped_type)
	##standard analyses
	if ped_type.lower() in trio_types:
		slivar_std_analysis_on_trio(annotated_vcfs, ped_file, formatted_ped_type)
	elif ped_type.lower() in duo_types or ped_type.lower() in single_types:
		slivar_std_analysis_on_duo(annotated_vcfs, ped_file, formatted_ped_type)
	elif ped_type.lower() in multiplex_types:
		slivar_std_analysis_on_multiplex(annotated_vcfs, ped_file, formatted_ped_type)
	else:
		print(ped_type, ' ped type not recognized')
	

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




def combine_var_files(work_directory, info_file, analysis_types):
	os.chdir(work_directory)
	for analysis_type in analysis_types:
		all_outfile = 'combined' + analysis_type + 'all.avinput'
		trio_quad_outfile = 'combined' + analysis_type + 'trio_quad.avinput'
		trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
		# single_duo_types = ['duo', 'duo*', 'parent_sibship','singleton', 'singleton*', 'sibship']
		# multiplex_types = ['multiplex']
		with open(info_file, "r") as in_fh, open(all_outfile, "w") as all_fh, open(trio_quad_outfile, "w") as tq_fh:
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
					std_file = ped + analysis_type + formatted_ped_type + '.std_analysis.avinput'
					# if os.path.isfile(std_file):
					with open(std_file, "r") as std_fh:
						for line2 in std_fh:
							line2 = line2.split(delim)
							line_out = delim.join(line2[:5] + [ped] + line2[5:])
							all_fh.write(line_out)
							if ped_type in trio_types:
								tq_fh.write(line_out)



##for combining files


working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalyze_exomes_0121/recall_vars_0221'

##repeat var calling/slivar analysis
# exome_ped_file = 'exome_test1' #test
# exome_ped_file = 'exome_info.txt' ##all the peds
##split 5 ways
# exome_ped_file = 'exome_infoaa'
# exome_ped_file = 'exome_infoab'
# exome_ped_file = 'exome_infoac'
# exome_ped_file = 'exome_infoad'
exome_ped_file = 'exome_infoae'

##redo slivar with new omim and repeats analyses etc
# slivar_analysis_master_0221(working_dir, exome_ped_file)

##latest versions of analysis files - 0221
var_types = ['.gatk38_slivar.', '.intersect_slivar.']


##combine avinput files, maybe in trios and all
ped_info_file = 'all_exome_info.txt'
combine_var_files(working_dir, ped_info_file, var_types)

##then manually make the 'solved' set of candidates



