#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys
import ep21_pisces_mosaic_calling_v0


'''
requires:
##without pisces
qsub -Iq cdbrmq -l mem=30gb,ncpus=2 -P 19833a08-f6fb-4bea-8526-8a79069da878
##with pisces
qsub -Iq cdbrmq -l mem=200gb,ncpus=20 -P 19833a08-f6fb-4bea-8526-8a79069da878
#java 1.8 for gatk
module load java/1.8.0_202 


'''

##set input variables and parameters
delim = '\t'
threads = '5'

##programs
bgzip = '/home/atimms/programs/samtools-1.11/bin/bgzip'
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
slivar = '/home/atimms/programs/slivar_0921/slivar'
# slivar = '/home/atimms/programs/slivar_0121/slivar_dev'
# pslivar = '/home/atimms/programs/slivar_0121/pslivar'
slivar_functions = '/home/atimms/programs/slivar_0921/slivar-0.2.5/js/slivar-functions_at.js'
# slivar_functions = '/home/atimms/programs/slivar_0921/slivar-0.2.5/js/slivar-functions.js'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
freebayes = '/home/atimms/programs/freebayes-1.3.4-linux-static-AMD64'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
tabix = '/home/atimms/programs/samtools-1.11/bin/tabix'
# plink = '/home/atimms/programs/plink'
somalier = '/home/atimms/programs/somalier_0721/somalier'

##files
fasta = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens_assembly38.fasta'
exome_capture_bed = '/home/atimms/ngs_data/references/exome_beds_hg38/hg38_targets_combined_padded_0721.bed'
##sorted specifically for this reference
hg38_refseq_exons = '/home/atimms/ngs_data/references/hg38/hg38_RefSeq_exons.sorted.bed'
ens_gff = '/home/atimms/ngs_data/references/hg38_gatk/Homo_sapiens.GRCh38.103.chr_added.gff3.gz'
gnomad_gnotate = '/home/atimms/ngs_data/references/slivar/gnomad.hg38.genomes.v3.fix.zip'
topmed_gnotate = '/home/atimms/ngs_data/references/slivar/topmed.hg38.dbsnp.151.zip'
#wget -qO - https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | zcat | cut -f 1,21,24 | tail -n+2 | awk '{ printf("%s\tpLI=%.3g;oe_lof=%.5g\n", $1, $2, $3)}' > pli.lookup
pli_lookup = '/home/atimms/ngs_data/references/slivar/pli.lookup'
#grep -v '^#' omim_genemap2_0321.txt | cut -f9,13 | grep -v '^\s' > /home/atimms/ngs_data/references/slivar/omim.lookup
omim_lookup = '/home/atimms/ngs_data/references/slivar/omim.lookup'
#cut -f2,3,4,6,8 dbdb.gene_association.txt | awk 'BEGIN{FS="\t"}{ printf("%s\tinheritance=%s;phenotype=%s;syndrome=%s;loe=%s\n", $1, $2, $3, $4, $5)}' > ~/ngs_data/references/slivar/dbdb.lookup
dbdb_lookup = '/home/atimms/ngs_data/references/slivar/dbdb.lookup'
hg38_selfchain = '/home/atimms/ngs_data/references/slivar/selfchain-LCR.hg38.bed'
hg38_refseq_bed = '/home/atimms/ngs_data/references/slivar/hg38.refseq.bed'

decipher_delevopmental_file = '/home/atimms/ngs_data/references/slivar/DDG2P_11_2_2021.csv'
dbsnp_common_vcf = '/home/atimms/ngs_data/references/hg38_gatk/common_all_20180418.chr_added.vcf.gz'
somalier_sites_vcf = '/home/atimms/ngs_data/references/somalier/sites.hg38.vcf.gz'

##ped types
trio_types = ['trio', 'quad', 'trio_with_sib', 'trio*', 'quint']
duo_types = ['duo', 'duo*', 'parent_sibship']
single_types = ['singleton', 'singleton*', 'sibship']
duo_single_types = duo_types + single_types
multiplex_types = ['multiplex']

##file sufffixes
gatk_vcf_suffix = '.gatk38hc.vcf.gz'
fb_vcf_suffix = '.freebayes13.vcf.gz'
int_vcf_dir = '.intersect_vcfs'
int_vcf_suffix = '.intersect_vcfs/0002.vcf'
gatk_bcftools_vcf_suffix = '.gatk38hc.bcftools.GRCh38_103.vcf.gz'
int_bcftools_vcf_suffix = '.intersect.bcftools.GRCh38_103.vcf.gz'

##annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,knownGene,ensGene,dbnsfp41a,clinvar_20210123,gnomad30_genome']
av_operation = ['-operation', 'g,g,g,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput', '-arg', '-splicing 10 ,,,,,' ]


##methods
def make_pedigree_dict(input_file):
	ped_dict = {}
	with open(input_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				ped_name = line[0]
				sample_name = line[1]
				ped_type = line[6]
				sequenced = line[7]
				mosaic = line[10]
				# sex = line[4]
				##add data to dict
				if sequenced == 'yes':
					if ped_name in ped_dict:
						ped_dict[ped_name][0].append(sample_name)
						ped_dict[ped_name][2].append(mosaic)
					else:
						ped_dict[ped_name] = [[sample_name],ped_type, [mosaic]]
	##retrun analysis info
	return ped_dict


##make list of all bam files to be analyzed
def make_list_of_bams(samples, bamlist_file):
	with open(bamlist_file, "w") as outfh:
		for sample in samples:
			bam = 'Preprocessing/' + sample + '/Recalibrated/' + sample + '.recal.bam'
			outfh.write(bam + '\n')

def variant_calling_gatk_hc_no_gvcf(bamlist, name_prefix):
	vcf_temp0 = name_prefix + 'temp_gatk0.vcf'
	vcf_raw_snps = name_prefix + 'temp_raw_snps.vcf'
	vcf_filtered_snps = name_prefix + 'temp_filtered_snps.vcf'
	vcf_raw_indels = name_prefix + 'temp_raw_indels.vcf'
	vcf_filtered_indels = name_prefix + 'temp_filtered_indels.vcf'
	vcf_temp1 = name_prefix + 'temp_gatk1.vcf'
	final_vcf = name_prefix + gatk_vcf_suffix
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
	bgzip_cmd = subprocess.Popen([bgzip, vcf_temp1])
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
	final_vcf = name_prefix + fb_vcf_suffix
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

def var_calling(ped_dict):
	for ped in ped_dict:
		bams_list = ped + '.bams.list'
		samples = ped_dict[ped][0]
		make_list_of_bams(samples, bams_list)
		gatk_vcf = ped + gatk_vcf_suffix
		fb_vcf = ped + fb_vcf_suffix
		##call vars and intersect
		variant_calling_gatk_hc_no_gvcf(bams_list, ped)
		variant_calling_freebayes(bams_list, ped)
		intersect_two_vcf_files(gatk_vcf, fb_vcf, ped + int_vcf_dir)


def process_annotate_vcf(input_vcf, out_vcf):
	temp_vcf = input_vcf.rsplit('.',2)[0] + '.temp.vcf.gz'
	##dcompress, change a gatk header thing,  (should be redundant)
	##so just and decompse and normalize
	#bcftools norm -m - -w 10000 -f $fasta -O b -o $clean_bcf
	bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', temp_vcf, '-O', 'z', input_vcf])
	bcftools_norm.wait()
	tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', temp_vcf])
	tabix_vcf.wait()
	##annotate with bcftools, have to use unphased version
	#bcftools csq -f hs37d5.fa -g Homo_sapiens.GRCh37.82.gff3.gz in.vcf -Ob -o out.bcf
	bcftools_csq = subprocess.Popen([bcftools, 'csq', '-g', ens_gff, '-f', fasta, '-o', out_vcf, '-O', 'z', '--local-csq', '-s', '-', temp_vcf])
	bcftools_csq.wait()
	tabix_vcf = subprocess.Popen([tabix, '-p', 'vcf', out_vcf])
	tabix_vcf.wait()

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
			run_slivar_duo_del = subprocess.Popen([slivar, 'duo-del', '--ped', ped_file, '--exclude', hg38_selfchain, '-a', in_vcf], stdout=t1_fh)
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
			manta_config = subprocess.Popen([bedtools, 'intersect', '-a', temp1_bed, '-b', hg38_refseq_bed, '-wa', '-wb'], stdout=out_fh)
			manta_config.wait()	
		##make dict from bed
		bed_dict = {}
		with open(temp2_bed, "r") as t2b_fh:
			for line in t2b_fh:
				line = line.rstrip().split(delim)
				cse = '_'.join(line[:3])
				# gene = line[6].split('_')[0]
				gene = line[6]
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

def slivar_std_analysis_on_trio(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg38_multianno.txt'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.xls'
		##get std file
		# '''
		slivar_trio = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '-g', topmed_gnotate, '--js', slivar_functions, '--info',
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
			merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final, slivar_std_vcf)
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

def slivar_std_analysis_on_duo(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg38_multianno.txt'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.xls'
		##get std file
		# '''
		slivar_duo = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '-g', topmed_gnotate, '--js', slivar_functions, '--info',
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
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final, slivar_std_vcf)

def slivar_std_analysis_on_multiplex(in_vcfs, ped_file, ped_type):
	for in_vcf in in_vcfs:
		slivar_std_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.vcf'
		slivar_cp_vcf = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.vcf'
		slivar_std_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.xls'
		slivar_cp_tsv = in_vcf.rsplit('.', 2)[0] + '.slivar_comphet.temp.xls'
		slivar_std_annovar = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp'
		slivar_std_multi = in_vcf.rsplit('.', 2)[0] + '.slivar_std.temp.hg38_multianno.txt'
		slivar_std_final = in_vcf.rsplit('.', 4)[0] + '_slivar.' + ped_type + '.std_analysis.xls'
		##get std file
		# '''
		slivar_multi = subprocess.Popen([slivar, 'expr', '--vcf', in_vcf, '--ped', ped_file,
			'--pass-only', '-g', gnomad_gnotate, '-g', topmed_gnotate, '--js', slivar_functions, '--info',
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
		merge_std_silvar_annovar(slivar_std_tsv, slivar_cp_tsv, slivar_std_multi, slivar_std_final, slivar_std_vcf)


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
				extra_header = line[8:18] + line[71:75] + [line[45]] + line[59:61] + line[76:88] + ['vcf_info', 'vcf_format'] + vcf_samples
			else:
				chr_pos_ref_alt = ':'.join(line[92:94] + line[95:97])
				##keep refgene, clivar, cadd, gerp, polyphe, extra gnomad, vcf info
				info_to_keep = line[8:18] + line[71:75] + [line[45]] + line[59:61] + line[76:88] + line[99:]				
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


def standard_slivar_protocol(ped_dict):
	for ped in ped_dict:
		ped_type = ped_dict[ped][1]
		ped_file = ped + '.ped'
		gatk_vcf = ped + gatk_vcf_suffix
		int_vcf = ped + int_vcf_suffix
		gatk_bcftools_vcf = ped + gatk_bcftools_vcf_suffix
		int_bcftools_vcf = ped + int_bcftools_vcf_suffix	
		annotated_vcfs = [gatk_bcftools_vcf, int_bcftools_vcf]
		##process vcf files, so have annotated vcf 
		# '''
		process_annotate_vcf(int_vcf, int_bcftools_vcf)
		process_annotate_vcf(gatk_vcf, gatk_bcftools_vcf)
		# '''
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
		print(ped, ped_type, formatted_ped_type)
		##standard analyses
		# '''
		if ped_type.lower() in trio_types:
			slivar_std_analysis_on_trio(annotated_vcfs, ped_file, formatted_ped_type)
		elif ped_type.lower() in duo_types or ped_type.lower() in single_types:
			slivar_std_analysis_on_duo(annotated_vcfs, ped_file, formatted_ped_type)
		elif ped_type.lower() in multiplex_types:
			slivar_std_analysis_on_multiplex(annotated_vcfs, ped_file, formatted_ped_type)
		else:
			print(ped_type, ' ped type not recognized')
		# '''
		##run duo_del
		# '''
		if ped_type.lower() in trio_types or ped_type.lower() in multiplex_types or ped_type.lower() in duo_types:
			slivar_duo_del(annotated_vcfs, ped_file, formatted_ped_type)
		# '''


def run_somalier(ped_dict, ped_file):
	for ped in ped_dict:
		gatk_vcf = ped + gatk_vcf_suffix
		#somalier extract -d extracted/ --sites sites.vcf.gz -f /data/human/g1k_v37_decoy.fa $cohort.vcf.gz
		som_extract = subprocess.Popen([somalier, 'extract', '-d', 'somalier_extract/', '--sites', somalier_sites_vcf, '-f', fasta, gatk_vcf])
		som_extract.wait()
	#somalier relate --ped $pedigree extracted/*.somalier
	# som_relate = subprocess.Popen([somalier, 'relate', '--ped', ped, 'somalier_extract/*.somalier'])
	som_relate = subprocess.Popen([somalier, 'relate', '--infer', '--ped', ped_file, 'somalier_extract/*.somalier'])
	som_relate.wait()


##get exome coverage for bam files
def calculate_exome_coverage(ped_dict):
	for ped in ped_dict:
		# print(ped, ped_dict[ped])
		##parameters and header
		coverage_required = [1,5,10,20,50,100,200]
		cov_head = []
		for cr in coverage_required:
			cr1 = 'total bp >=' + str(cr)
			cr2 = 'percentage covered >=' + str(cr)
			cov_head.extend([cr1, cr2])
		header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
		##open results file and write header
		output_file = ped + '.coverage_data.txt'
		print(output_file)
		samples = ped_dict[ped][0]
		with open(output_file, "w") as outfh:
			outfh.write(header)
			##go through bam files and get coverage histogram
			for sample in samples:
				bam = 'Preprocessing/' + sample + '/Recalibrated/' + sample + '.recal.bam'
				print('calculating coverge for bam file', bam)
				##get temp coverage files
				coverage_histogram = sample + '.hist.temp'
				hist_fh = open(coverage_histogram, 'w')
				# bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=hist_fh)
				# bedtools_cov.wait()
				bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', hg38_refseq_exons, '-b', bam, '-hist', '-sorted'], stdout=subprocess.PIPE)
				grep_cmd = subprocess.Popen(['grep', '^all'], stdin=bedtools_cov.stdout, stdout=hist_fh)
				grep_cmd.wait()
				hist_fh.close()
				##calculate numbers from histogram
				with open(coverage_histogram, "r") as cov_temp:
					seq_total = 0
					for line in cov_temp:
						line = line.strip('\n').split(delim)
						target_size = float(line[3])
						if line[1] != '0':
							seq_total += int(line[1]) *  int(line[2])
				cov_stats = []
				for cr in coverage_required:
					coverage_bp = 0
					with open(coverage_histogram, "r") as cov_temp:
						for line in cov_temp:
							line = line.strip('\n').split(delim)
							target_size = float(line[3])
							if int(line[1]) >= cr:
								coverage_bp += int(line[2])
						cov_stats.extend([str(coverage_bp), str((coverage_bp / target_size) * 100)])
				line_out = delim.join([bam, str(target_size), str(seq_total), str(seq_total / target_size)] + cov_stats + ['\n'])
				outfh.write(line_out)

def combine_results(project):
	results_dir = project + '_results'
	os.mkdir(results_dir)
	files_to_go = glob.glob('*coverage_data.txt') + glob.glob('*.std_analysis.xls') + glob.glob('*.pisces.xls') + ['somalier.html']
	for file_to_go in files_to_go:
		shutil.copy(file_to_go, results_dir)



def post_bam_master(working_dir, info_file):
	os.chdir(working_dir)
	project_name = info_file.rsplit('.',1)[0]
	##get ped info, samples/pedtype/mosaic per ped
	pedigree_dict = make_pedigree_dict(info_file)
	# print(pedigree_dict)
	'''
	##call vars on gatk and freeabyes and intersect
	var_calling(pedigree_dict)
	##slivar
	standard_slivar_protocol(pedigree_dict)
	##somalier
	run_somalier(pedigree_dict, project_name + '.ped')
	##coverage
	calculate_exome_coverage(pedigree_dict)
	##mosaic
	ep21_pisces_mosaic_calling_v0.run_mosaic_variant_calling(working_dir, pedigree_dict)
	'''
	##combine results files into new directory
	combine_results(project_name)















##run methods
##ghayda 0821 2x singles plus extras
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_genedx_0821/'
exome_info_file = 'ghayda_genedx_0821.txt'
# exome_info_file = 'LR10-064_0821.txt' ##test on just LR10-064 (issue with rg)
# post_bam_master(work_dir, exome_info_file)

##kim 0921 34 singles
work_dir = '/home/atimms/ngs_data/exomes/working/kim_exomes_0621/'
exome_info_file = 'kim_exomes_0621.txt'
post_bam_master(work_dir, exome_info_file)
 

