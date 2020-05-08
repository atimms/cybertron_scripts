#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'
star_threads = '10'

##programs
picard = '/home/atimms/programs/picard_2.19/picard.jar'
qorts = '/home/atimms/programs/hartleys-QoRTs-39cd1fc/QoRTs.jar'
gatk = '/cm/shared/apps/GATK/3.7/GenomeAnalysisTK.jar'
table_annovar = '/home/atimms/programs/annovar_0618/table_annovar.pl'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'

##ref files etc
star_genome_dir = '/home/atimms/ngs_data/references/star/hg19'
fa_file = '/home/atimms/ngs_data/references/igenomes/hg19/genome.fa'
gtf_file = '/home/atimms/ngs_data/references/igenomes/hg19/genes.gtf'
ref_dir = '/home/atimms/ngs_data/references/hg19/'
indels_mills = ref_dir + 'Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
dbsnp = ref_dir + 'dbsnp_138.hg19.vcf'
exome_capture_bed = ref_dir + 'agilent_v7_hg19_regions_padded.bed'
star_bam_suffix = '.Aligned.out.bam'

#annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp35a,dbnsfp31a_interpro,avsnp147,gnomad211_genome,gnomad211_exome,exac03']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.','-vcfinput']

##sample info dict

def get_info_from_file(sample_file):
	sample_info_dict = {}
	with open(sample_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				sample_name = line[0]
				##add data to ped file dict
				sample_info = line[1:]
				if sample_name in sample_info_dict:
					print('sample %s seen multiple times'%sample_name)
				else:
					sample_info_dict[sample_name] = sample_info
	return sample_info_dict		

def convert_bam_fastq_bedtools(bamfile, fq_names):
	read_sorted_bam = bamfile[:-4] + 'n_sorted.bam'
	r1_fastq = fq_names[0]
	r2_fastq = fq_names[1]
	st_n_sort = subprocess.Popen(['samtools', 'sort', '-nO', 'bam', '-o', read_sorted_bam, '-@', '10', '-m', '10G', '-T', bamfile[:-4] + 'tempy', bamfile])
	st_n_sort.wait()
	bam_fq = subprocess.Popen(['bedtools', 'bamtofastq', '-i', read_sorted_bam, '-fq', r1_fastq, '-fq2', r2_fastq])
	bam_fq.wait()
	os.remove(read_sorted_bam)

def star_align_paired_end_2_pass(sample_name, fq_names):
	r1_fq = fq_names[0]
	r2_fq = fq_names[1]
	# star_align = subprocess.Popen(['STAR','--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--outSAMmapqUnique', '60', '--readFilesCommand', 'zcat', '--twopassMode', 'Basic', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align = subprocess.Popen(['STAR','--genomeDir', star_genome_dir, '--readFilesIn', r1_fq, r2_fq, '--outSAMtype', 'BAM', 'Unsorted', '--outSAMmapqUnique', '60', '--twopassMode', 'Basic', '--outFileNamePrefix', sample_name + '.', '--runThreadN', star_threads])
	star_align.wait()
	os.remove(r1_fq)
	os.remove(r2_fq)



def proccess_bam_files(sample_name, proccessed_bam):
	in_bam = sample_name + star_bam_suffix
	arg_bam = sample_name + '.with_rg.bam'
	mkdup_bam = sample_name + '.mkdup.bam'
	reorder_bam = sample_name + '.reorder.bam'
	split_bam = sample_name + '.split.bam'
	final_bam = proccessed_bam
	picard_arg = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'AddOrReplaceReadGroups', 'INPUT=' + in_bam, 'OUTPUT=' + arg_bam, 'SO=coordinate', 'RGID=' + sample_name, 'RGLB=' + sample_name, 'RGPL=Illumina', 'RGPU=machine1', 'RGSM=' + sample_name ])
	picard_arg.wait()	
	picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=SILENT', 'INPUT=' + arg_bam, 'OUTPUT=' + mkdup_bam,'CREATE_INDEX=true', 'METRICS_FILE=' + sample_name + '.metrics'])
	picard_md.wait()
	picard_rs = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'ReorderSam', 'INPUT=' + mkdup_bam, 'OUTPUT=' + reorder_bam, 'CREATE_INDEX=true', 'REFERENCE=' + fa_file])
	picard_rs.wait()
	# gatk_snt = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'SplitNCigarReads',  '-R', fa_file, '-I', reorder_bam, '-o', split_bam, '-rf', 'ReassignOneMappingQuality', '-RMQF', '255', '-RMQT', '60', '-U', 'ALLOW_N_CIGAR_READS'])
	gatk_snt = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'SplitNCigarReads',  '-R', fa_file, '-I', reorder_bam, '-o', split_bam, '-U', 'ALLOW_N_CIGAR_READS'])
	gatk_snt.wait()
	gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '6', '-R', fa_file, '-I', split_bam, '-knownSites', indels_mills, '-knownSites', dbsnp, '-o', sample_name + '.recal_data.table'])
	gatk_br.wait()
	gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '6', '-R', fa_file, '-I', split_bam, '-BQSR', sample_name + '.recal_data.table', '-o', final_bam])
	gatk_pr.wait()
	##remove intermeidate files
	files_to_go = [arg_bam, mkdup_bam,reorder_bam,split_bam]
	print 'removing files:'
	for f in files_to_go:
		os.remove(f)
		print(f)

def haplotype_caller_gvcf(sample, in_bam):
	##using exome only regions
	# gatk_hc = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fa_file, '-I', in_bam,'-L', exome_capture_bed, '-dontUseSoftClippedBases', '--emitRefConfidence', 'GVCF', '--variant_index_type', 'LINEAR', '--variant_index_parameter', '128000', '--max_alternate_alleles', '3', '-o', sample + '.g.vcf'])
	gatk_hc = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'HaplotypeCaller', '-R', fa_file, '-I', in_bam, '-dontUseSoftClippedBases', '-stand_call_conf', '20', '--emitRefConfidence', 'GVCF', '--variant_index_type', 'LINEAR', '--variant_index_parameter', '128000', '--max_alternate_alleles', '3', '-o', sample + '.g.vcf'])
	gatk_hc.wait()

def combine_gvcfs(samples, out_vcf):
	gvcf_files = []
	for sample in samples:
		gvcf = ['-V', sample + '.g.vcf']
		gvcf_files.extend(gvcf)
	print gvcf_files
	gatk_hc = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineGVCFs', '-R', fa_file] + gvcf_files + ['-o', out_vcf])
	gatk_hc.wait()

def genotype_gvcfs_and_filter(batch_gvcfs, outfile_prefix):
	gvcf_cmds = []
	vcf_temp0 = outfile_prefix + 'temp0.vcf'
	final_vcf = outfile_prefix + '.gatk_vars.vcf'
	for batch_gvcf in batch_gvcfs:
		gvcf = ['-V', batch_gvcf]
		gvcf_cmds.extend(gvcf)
	print gvcf_cmds
	##genotype g.vcfs
	##with nt and exomes
	# command = ['java', '-Xmx100g', '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', fa_file, '-nt', '10','-L', exome_capture_bed] + gvcf_cmds + ['-o', vcf_temp0]
	# command = ['java', '-Xmx100g', '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', fa_file, '-L', exome_capture_bed] + gvcf_cmds + ['-o', vcf_temp0]
	##straight up
	# '''
	command = ['java', '-Xmx100g', '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', fa_file] + gvcf_cmds + ['-o', vcf_temp0]
	gatk_gg = subprocess.Popen(command)
	gatk_gg.wait()
	# '''
	# gatk_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fa_file, '-V', vcf_temp0, '-o', final_vcf, '-window', '35', '-cluster', '3', '-filterName', 'FS','-filter', 'FS > 30.0', '-filterName', 'QD', '-filter', 'QD < 2.0'])
	gatk_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fa_file, '-V', vcf_temp0, '-o', final_vcf, '-filterName', 'FS','-filter', 'FS > 30.0', '-filterName', 'QD', '-filter', 'QD < 2.0'])

	gatk_vf.wait()

def split_info_field(info_list):
	indices = [i for i, s in enumerate(info_list) if 'ANNOVAR_DATE' in s]
	# print indices
	i_count = 0
	final_list = []
	for info in info_list:
		# print info
		if i_count > indices[0] and info != 'ALLELE_END':
			info2 = info.split('=')[1]
			#print info2
			final_list.append(info2)
		i_count += 1
	return final_list

def annotate_vars_from_vcf(out_prefix, in_vcf, norm_vcf, samples):
	'''
	bgzip_cmd = subprocess.Popen(['bgzip', in_vcf])
	bgzip_cmd.wait()
	vcf_temp1 = 'temp111.vcf'
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-o', vcf_temp1, in_vcf + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fa_file, '-O', 'z', '-o', norm_vcf, vcf_temp1])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', norm_vcf])
	bcf_index.wait()
	
	##annotate vcf file
	command = [table_annovar] + av_buildver + [norm_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()
	##split into individual samples
	post_annovar_vcf = out_prefix + '.' + av_genome + '_multianno.vcf'
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', post_annovar_vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', out_prefix])
	con_ann.wait()
	'''
	##format files...
	for sample in samples:
		avinput = out_prefix + '.' + sample + '.avinput'
		outfile = out_prefix + '.' + sample + '.annotated.txt'
		head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Zygosity', 'Qual', 'Coverage', 'Filter', 'Pos', 'Ref2', 'Alt2', 'Format', 'Info', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'SIFT_score', 'SIFT_converted_rankscore', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_rankscore', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_converted_rankscore', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_converted_rankscore', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_score_rankscore', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_converted_rankscore', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_converted_rankscore', 'PROVEAN_pred', 'VEST3_score', 'VEST3_rankscore', 'MetaSVM_score', 'MetaSVM_rankscore', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_rankscore', 'MetaLR_pred', 'M-CAP_score', 'M-CAP_rankscore', 'M-CAP_pred', 'REVEL_score', 'REVEL_rankscore', 'MutPred_score', 'MutPred_rankscore', 'CADD_raw', 'CADD_raw_rankscore', 'CADD_phred', 'DANN_score', 'DANN_rankscore', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_rankscore', 'fathmm-MKL_coding_pred', 'Eigen_coding_or_noncoding', 'Eigen-raw', 'Eigen-PC-raw', 'GenoCanyon_score', 'GenoCanyon_score_rankscore', 'integrated_fitCons_score', 'integrated_fitCons_score_rankscore', 'integrated_confidence_value', 'GERP++_RS', 'GERP++_RS_rankscore', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP20way_mammalian', 'phyloP20way_mammalian_rankscore', 'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons20way_mammalian', 'phastCons20way_mammalian_rankscore', 'SiPhy_29way_logOdds', 'SiPhy_29way_logOdds_rankscore', 'Interpro_domain', 'GTEx_V6p_gene', 'GTEx_V6p_tissue', 'Interpro_domain', 'avsnp147', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax','ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS']
		head_out = delim.join(head + ['\n'])
		with open(avinput, "r") as av, open(outfile, "w") as final:
			final.write(head_out)
			for line in av:
				line = line.strip('\n').split(delim)
				stuff = line[0:8] + [line[14]] + [line[9]] + line[11:13]
				info = line[15].split(';')
				info_list = split_info_field(info)
				other_stuff = line[16:]
				line_out = delim.join(stuff + other_stuff + info_list +['\n'])
				final.write(line_out)

def combine_filter_vars(out_prefix, s_dict):
	var_dict = {}
	case_count, ctl_count = 0,0
	cases_temp_file = out_prefix + '.cases_temp1.xls'
	cases_outfile = out_prefix + '.cases_filtered.xls'
	controls_temp_file = out_prefix + '.controls_temp1.xls'
	controls_outfile = out_prefix + '.controls_filtered.xls'
	with open(cases_temp_file, 'w') as case_tf, open(controls_temp_file, 'w') as ctl_tf:
		fc = 0
		for sample in s_dict:
			fc += 1
			ann_text = out_prefix + '.' + sample + '.annotated.txt'
			sample_info = s_dict[sample][:4]
			dx = sample_info[2]
			if dx == 'Control':
				ctl_count += 1
			else:
				case_count += 1
			print sample_info, sample
			with open(ann_text, 'r') as in_fh:
				lc = 0
				for line in in_fh:
					line = line.strip().split(delim)
					lc += 1
					if lc ==1:
						if fc == 1:
							header = delim.join(['sample', 'ID', 'sex', 'dx', 'cohort'] + line) + '\n'
							case_tf.write(header)
							ctl_tf.write(header)
					else:
						##only keep filtered vars:
						# print line[:8]
						cov = line[7]
						if cov == '.':
							cov = 0
						gatk_filter = line[8]
						rmsk = line[19]
						supdup = line[20]
						if int(cov) >= 10 and gatk_filter == 'PASS' and rmsk == '.' and supdup == '.':
							line_out = delim.join([sample] + sample_info + line) + '\n'
							var = '_'.join(line[:5])
							if dx == 'Control':
								ctl_tf.write(line_out)
								
								if var in var_dict:
									var_dict[var][1] += 1
								else:
									var_dict[var] = [0,1]
							else:
								case_tf.write(line_out)
								
								if var in var_dict:
									var_dict[var][0] += 1
								else:
									var_dict[var] = [1,0]
	# print var_dict
	with open(cases_temp_file, 'r') as case_tf, open(cases_outfile, 'w') as case_out:
		lc = 0
		for line in case_tf:
			line = line.strip().split(delim)
			lc += 1
			if lc == 1:
				header = delim.join(line[:5]+ ['case_count', 'control_count', 'genename'] + line[5:]) + '\n'
				case_out.write(header)
			else:				
				v = '_'.join(line[5:10])
				v_info = [str(i) for i in var_dict[v]]
				genename = line[20]
				line_out = delim.join(line[:5] + v_info + [genename] + line[5:]) + '\n'
				case_out.write(line_out)				
	with open(controls_temp_file, 'r') as ctl_tf, open(controls_outfile, 'w') as ctl_out:
		lc = 0
		for line in ctl_tf:
			line = line.strip().split(delim)
			lc += 1
			if lc == 1:
				header = delim.join(line[:5]+ ['case_count', 'control_count', 'genename'] + line[5:]) + '\n'
				ctl_out.write(header)
			else:				
				v = '_'.join(line[5:10])
				v_info = [str(i) for i in var_dict[v]]
				genename = line[20]
				line_out = delim.join(line[:5] + v_info + [genename] + line[5:]) + '\n'
				ctl_out.write(line_out)	
	print case_count, ctl_count

def call_all_rnaseq_methods(working_dir, sample_file):
	os.chdir(working_dir)
	sample_dict = get_info_from_file(sample_file)
	project = sample_file.split('.')[0]
	combined_batch_gvcf = project + '.g.vcf'
	# print(sample_dict)
	for sample in sample_dict:
		in_bam = sample_dict[sample][4]
		read1_fq = sample + '.r1.fastq'
		read2_fq = sample + '.r2.fastq'
		out_bam = sample + '.star_gatk.bam'
		# convert_bam_fastq_bedtools(in_bam, [read1_fq, read2_fq])
		# star_align_paired_end_2_pass(sample, [read1_fq, read2_fq])
		# proccess_bam_files(sample, out_bam)
		##can be run on few cores and lower memory --- call on batches
		# haplotype_caller_gvcf(sample, out_bam)
	##call on batches
	# combine_gvcfs(sample_dict, combined_batch_gvcf)
	##do seperatley and just once
	b_gvcfs = glob.glob('*b*g.vcf')
	print(b_gvcfs)
	# genotype_gvcfs_and_filter(b_gvcfs, project)
	##annotatate variants
	variant_vcf = project + '.gatk_vars.vcf'
	normalized_vcf = project + '.gatk_vars_normalized.vcf.gz'
	# sample_dict = ['163499', '163500']
	annotate_vars_from_vcf(project, variant_vcf, normalized_vcf, sample_dict)
	combine_filter_vars(project, sample_dict)


def filter_var_files_for_exonic(infile, ex_out, exu_out, exf_out):
	with open(infile, "r") as in_fh, open(ex_out, "w") as ex_fh, open(exu_out, "w") as exu_fh, open(exf_out, "w") as exf_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc == 1:
				ex_fh.write(line)
				exu_fh.write(line)
				exf_fh.write(line)
			else:
				line = line.split(delim)
				func_rg = line[22].split("\\")[0]
				func_exonic = line[25]
				gnomad_exome_AF = line[103]
				gnomad_genome_AF = line[86]
				if gnomad_exome_AF == '.':
					gnomad_exome_AF = 0
				else:
					gnomad_exome_AF = float(gnomad_exome_AF)
				if gnomad_genome_AF == '.':
					gnomad_genome_AF = 0
				else:
					gnomad_genome_AF = float(gnomad_genome_AF)
				if func_rg == 'UTR3' or func_rg == 'UTR5' or func_rg == 'exonic' or func_rg == 'splicing': 
					if func_exonic != 'synonymous_SNV':
						exu_fh.write(delim.join(line))
				if func_rg == 'exonic' or func_rg == 'splicing':
					if func_exonic != 'synonymous_SNV':
						ex_fh.write(delim.join(line))
				if func_rg == 'exonic' or func_rg == 'splicing':
					if func_exonic != 'synonymous_SNV':
						if gnomad_genome_AF <= 0.05 and gnomad_exome_AF <= 0.05:
							exf_fh.write(delim.join(line))

def filter_var_files_for_gene(infile, outfile, genelist):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc == 1:
				out_fh.write(line)
			else:
				line = line.split(delim)
				genenames = line[7].split("x3b")
				for genename in genenames:
					if genename in genelist:
						out_fh.write(delim.join(line))

def filter_var_files_for_aaf_exonic(infile, aaf_out, aaf_exonic_out, aaf_req):
	with open(infile, "r") as in_fh, open(aaf_out, "w") as aaf_fh, open(aaf_exonic_out, "w") as aafex_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc == 1:
				aaf_fh.write(line)
				aafex_fh.write(line)
			else:
				line = line.split(delim)
				func_rg = line[22].split("\\")[0]
				func_exonic = line[25]
				gnomad_exome_AF = line[118]
				gnomad_genome_AF = line[101]
				if gnomad_exome_AF == '.':
					gnomad_exome_AF = 0
				else:
					gnomad_exome_AF = float(gnomad_exome_AF)
				if gnomad_genome_AF == '.':
					gnomad_genome_AF = 0
				else:
					gnomad_genome_AF = float(gnomad_genome_AF)
				if gnomad_genome_AF <= aaf_req and gnomad_exome_AF <= aaf_req:
					aaf_fh.write(delim.join(line))
					if func_rg == 'exonic' or func_rg == 'splicing':
						if func_exonic != 'synonymous_SNV':
							aafex_fh.write(delim.join(line))




##sample files
info_file = 'cunn_rnaseq_0419.txt'
b1_file = 'cunn_rnaseq_0419_b1.txt'
b2_file = 'cunn_rnaseq_0419_b2.txt'
b3_file = 'cunn_rnaseq_0419_b3.txt'
b4_file = 'cunn_rnaseq_0419_b4.txt'
b5_file = 'cunn_rnaseq_0419_b5.txt'
b6_file = 'cunn_rnaseq_0419_b6.txt'
b7_file = 'cunn_rnaseq_0419_b7.txt'
b8_file = 'cunn_rnaseq_0419_b8.txt'
b9_file = 'cunn_rnaseq_0419_b9.txt'
b10_file = 'cunn_rnaseq_0419_b10.txt'
work_dir = '/home/atimms/ngs_data/rnaseq/cunningham_rnaseq_0419'
os.chdir(work_dir)
##batches
# call_all_rnaseq_methods(work_dir, b1_file)
# call_all_rnaseq_methods(work_dir, b2_file)
# call_all_rnaseq_methods(work_dir, b3_file)
# call_all_rnaseq_methods(work_dir, b4_file)
# call_all_rnaseq_methods(work_dir, b5_file)
# call_all_rnaseq_methods(work_dir, b6_file)
# call_all_rnaseq_methods(work_dir, b7_file)
# call_all_rnaseq_methods(work_dir, b8_file)
# call_all_rnaseq_methods(work_dir, b9_file)
# call_all_rnaseq_methods(work_dir, b10_file)

##all samples -- used to genotype the gvcfs and then annotate etc
# call_all_rnaseq_methods(work_dir, info_file)

##test set
# call_all_rnaseq_methods(work_dir, 'cunn_rnaseq_0419.temp.txt')

##filter variant files
all_case_vars = 'cunn_rnaseq_0419.cases_filtered.xls'
all_control_vars = 'cunn_rnaseq_0419.controls_filtered.xls'
# exonic_cases = 'cunn_rnaseq_0419.cases_filtered.exonic.xls'
# exonic_controls = 'cunn_rnaseq_0419.controls_filtered.exonic.xls'
# exonic_utr_cases = 'cunn_rnaseq_0419.cases_filtered.exonic_utr.xls'
# exonic_utr_controls = 'cunn_rnaseq_0419.controls_filtered.exonic_utr.xls'
# exonic_5perc_cases = 'cunn_rnaseq_0419.cases_filtered.exonic.05af.xls'
# exonic_5perc_controls = 'cunn_rnaseq_0419.controls_filtered.exonic.05af.xls'

##get exonic vars
# os.chdir(work_dir)
# filter_var_files_for_exonic(all_case_vars, exonic_cases, exonic_utr_cases, exonic_5perc_cases)
# filter_var_files_for_exonic(all_control_vars, exonic_controls, exonic_utr_controls, exonic_5perc_controls)

##get gene specific variants
# genes_to_get = ['PIEZO1']
# cases_for_genes = all_case_vars.rsplit('.',1)[0] + '.PIEZO1.xls'
# controls_for_genes = all_control_vars.rsplit('.',1)[0] + '.PIEZO1.xls'
# filter_var_files_for_gene(all_case_vars, cases_for_genes, genes_to_get)
# filter_var_files_for_gene(all_control_vars, controls_for_genes, genes_to_get)
# genes_to_get = ['ADAMTSL4']
# cases_for_genes = all_case_vars.rsplit('.',1)[0] + '.ADAMTSL4.xls'
# controls_for_genes = all_control_vars.rsplit('.',1)[0] + '.ADAMTSL4.xls'
# filter_var_files_for_gene(all_case_vars, cases_for_genes, genes_to_get)
# filter_var_files_for_gene(all_control_vars, controls_for_genes, genes_to_get)

##get one file per gene, must just be one gene per list
# genelists_to_get = [['PIEZO1'], ['FLNA'], ['FLNB'], ['FLNC'], ['AXL']]
# for genelist in genelists_to_get:
# 	cases_for_genes = all_case_vars.rsplit('.',1)[0] + '.' + genelist[0] + '.xls'
# 	controls_for_genes = all_control_vars.rsplit('.',1)[0] + '.' + genelist[0] + '.xls'
# 	filter_var_files_for_gene(all_case_vars, cases_for_genes, genelist)
# 	filter_var_files_for_gene(all_control_vars, controls_for_genes, genelist)

##filter by allele feq and if exonic
allele_freqs_to_get = [0.02,0.01,0.001]
for allele_freq in allele_freqs_to_get:
	aaf_cases = all_case_vars.rsplit('.',1)[0] + '.' + str(allele_freq) + '_aaf.xls'
	aaf_exonic_cases = all_case_vars.rsplit('.',1)[0] + '.' + str(allele_freq) + '_aaf.exonic.xls'
	aaf_controls = all_control_vars.rsplit('.',1)[0] + '.' + str(allele_freq) + '_aaf.xls'
	aaf_exonic_controls = all_control_vars.rsplit('.',1)[0] + '.' + str(allele_freq) + '_aaf.exonic.xls'
	filter_var_files_for_aaf_exonic(all_case_vars, aaf_cases, aaf_exonic_cases, allele_freq)
	filter_var_files_for_aaf_exonic(all_control_vars, aaf_controls, aaf_exonic_controls, allele_freq)



