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
pslivar = '/home/atimms/programs/slivar_0121/pslivar'
slivar_functions = '/home/atimms/programs/slivar_0121/slivar-0.2.1/js/slivar-functions_at.js'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'

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

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,gnomad211_exome,gnomad211_genome,dbnsfp41a,clinvar_20200316,generic', '-genericdbfile', 'Geno2MP.avinput']

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
	##annotate with snpeff and compress and index (try 3x references)
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




def slivar_std_analysis_on_trio(in_vcfs, ped_file):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg19_multianno.txt'
		slivar_std_final = in_vcf.rsplit('.', 2)[0] + '.slivar.std_analysis.xls'
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
		##annotate vcfs with annovar i.e. run_table_annovar
		command = [table_annovar] + av_buildver + [slivar_std_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', slivar_std_annovar]
		annovar = subprocess.Popen(command)
		annovar.wait()
		# '''
		##combine slivar files with annovar annotation
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final)

		##extra analyses
		'''
		slivar_dn_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_exonic_dn.vcf'
		slivar_dom_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_dominant.vcf'
		
		##all exonic de novo, i.e. use genic instead of impactful
		slivar_dn = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.genic && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
			'--family-expr', 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001',
			'-o', slivar_dn_vcf])
		slivar_dn.wait()
		##dominant -- not working
		slivar_dom = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
			'--family-expr', 'dominant:fam.every(segregating_dominant)',
			'-o', slivar_dom_vcf])
		slivar_dom = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'-g', gnomad_gnotate, '--js', slivar_functions, '--info',
			'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"',
			'-o', slivar_dom_vcf])
		slivar_dom.wait()
		'''		


def standard_slivar_protocol(working_directory, pedigree, ped_type):
	os.chdir(working_directory)
	in_vcf = pedigree + '.intersected_vcfs/0002.vcf.gz'
	gatk_vcf = pedigree + '.gatkHC.vcf.gz'
	ped_file = pedigree + '.ped'
	annotated_vcfs = [pedigree + '.intersected.bcftools.GRCh37_87.vcf.gz', pedigree + '.intersected.snpeff.GRCh37_87.vcf.gz',
			pedigree + '.intersected.snpeff.hg19.vcf.gz', pedigree + '.gatk.bcftools.GRCh37_87.vcf.gz', 
			pedigree + '.gatk.snpeff.GRCh37_87.vcf.gz', pedigree + '.gatk.snpeff.hg19.vcf.gz']
	# annotated_vcfs = [pedigree + '.intersected.bcftools.GRCh37_87.vcf.gz']
	##process vcf files, so have normalized and annotated vcf 
	process_annotate_vcf(pedigree + '.intersected', in_vcf)
	process_annotate_vcf(pedigree + '.gatk', gatk_vcf)
	##filter with slivar
	# if ped_type.lower() == 'trio':
	# 	slivar_std_analysis_on_trio(annotated_vcfs, ped_file)


	##testing other ped types
	# else:
	# 	slivar_std_analysis_on_trio(annotated_vcfs, ped_file)


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

##tests on 0920 plus exomes
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalyze_exomes_0121/exomes_after_0920'
input_exomes_file = 'exomes_after_0920.txt'
# input_exomes_file = 'exomes_after_0920_test.txt'

##input files
# slivar_exome_inputs(working_dir, input_exomes_file)

##individual peds
# standard_slivar_protocol(working_dir, "LR08-081", 'trio')
# standard_slivar_protocol(working_dir, "LR05-298", 'trio')
# standard_slivar_protocol(working_dir, "LR14-331", 'multiplex')
# standard_slivar_protocol(working_dir, "LR18-427", 'duo')
# standard_slivar_protocol(working_dir, "LR18-545", 'trio_plus')

##exomes from rss
working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalyze_exomes_0121/exomes_from_rss_0420'
# input_exomes_file = 'exomes_rss_0420_a.txt'
# input_exomes_file = 'exomes_rss_0420_b.txt'
input_exomes_file = 'exomes_rss_0420_c.txt'

##input files
# slivar_exome_inputs(working_dir, input_exomes_file)

##ind samples
standard_slivar_protocol(working_dir, "LP97-070", 'trio')
#rm *temp* *snpeff* *bcftools*


