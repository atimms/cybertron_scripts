#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import gzip

'''
##load modules required for analysis
module load biobuilds/2016.11
module load vt/0.5772-60f436c3 
'''


##parameters
delim = '\t'

##programs and ref files
vt = 'vt'
gemini = '/home/atimms/scripts/gemini'
snpeff_jar = '/home/atimms/programs/snpEff/snpEff.jar'
bcftools_12 = 'bcftools'
ref_dir = '/home/atimms/ngs_data/references/hg19/'
plink = '/home/atimms/programs/plink'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'
ann_var = '/home/atimms/programs/annovar_1019/annotate_variation.pl'
convert_2_annovar = '/home/atimms/programs/annovar_1019/convert2annovar.pl'
manta_config = '/home/atimms/programs/manta-1.6.0.centos6_x86_64/bin/configManta.py'

##unique for this project
fasta = ref_dir + 'hs37d5.fa'
gemini_ref_dir = '/home/atimms/ngs_data/references/gemini/'
# rvis_data = gemini_ref_dir + 'RVIS_ExAC_4KW.txt'
# gdi_data = gemini_ref_dir + 'GDI_full_10282015.txt'
# mgi_data = gemini_ref_dir + 'mgi.abnormal_outer_ear_morhology.txt'
# hpo_data = gemini_ref_dir + 'genes_for_HP_0000356.csv'
# sbse_rnaseq_data = gemini_ref_dir + 'esra.gene_exp.diff'
# hmx1_rnaseq_data = gemini_ref_dir + 'jess.gene_exp.diff'
# human_ear_rnaseq_data = gemini_ref_dir + 'human_ear.gene_exp.diff'
# human_syndrome_data = gemini_ref_dir + 'microtia_human_syndromes.txt'
# hoxa2_data = gemini_ref_dir + 'hoxa2.expression_data.txt'
# ba_exp_data = gemini_ref_dir + 'ba_expression_papers_0816.txt'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,gerp++gt2,rmsk,genomicSuperDups,gnomad211_genome,gnomad211_exome']
av_operation = ['-operation', 'g,f,r,r,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.']


##lists of peds
duo_list = ['CL3003', 'CL3023', 'CL3051', 'F040', 'PM1010037']
single_list = ['CFM-MOS-11', 'F168']


def make_analysis_dicts(input_file):
	ped_file_dict, analysis_dict = {}, {}
	with open(input_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				ped_name = line[0]
				sample_name = line[1]
				seq = line[6]
				anal_type = line[7]
				# print ped_name, line
				##add data to ped file dict
				ped_file_info = line[:6]
				# print ped_file_info
				if ped_name in ped_file_dict:
					ped_file_dict[ped_name].append(ped_file_info)
				else:
					ped_file_dict[ped_name] = [ped_file_info]
				##add data to analysis_dict
				if seq == 'yes':
					if ped_name in analysis_dict:
						analysis_dict[ped_name][0].append(sample_name)
					else:
						analysis_dict[ped_name] = [[sample_name], anal_type]


	return ped_file_dict, analysis_dict		

def make_ped_files(input_dict):
	for ped in input_dict:
		# print ped, input_dict[ped]
		outfile = ped + '.ped'
		header = ['#family_id', 'name', 'paternal_id', 'maternal_id', 'sex', 'phenotype', '\n']
		with open(outfile, "w") as out_fh:
			out_fh.write(delim.join(header))
			for outline in input_dict[ped]:
				out_fh.write(delim.join(outline) + '\n')

def make_sample_file(samples, outfile):
	with open(outfile, "w") as out_fh:
		for sample in samples:
			out_fh.write(sample + '\n')

def get_vcf_per_ped(pedigree_dict, in_vcf):
	for pedigree in pedigree_dict:
		sample_list = pedigree_dict[pedigree][0]
		sample_file = pedigree + '.samples_temp.txt'
		ped_vcf = pedigree + '.vcf.gz'
		make_sample_file(sample_list, sample_file)
		##get the samples we want, and remove when we don't see a call
		bcftools_view = subprocess.Popen(['bcftools', 'view', '-a', '--threads', '10', '-Ou', '-S', sample_file, in_vcf], stdout=subprocess.PIPE)
		bcftools_view2 = subprocess.Popen(['bcftools', 'view', '-m', '2', '--threads', '10', '-O', 'z', '-o', ped_vcf, '-'], stdin=bcftools_view.stdout)
		bcftools_view2.wait()
		

def cp_gemini_db(in_db, out_db):
	cp_command = subprocess.Popen(['cp', in_db, out_db])
	cp_command.wait()

def add_new_ped_file_to_gemini_db(pedfile, geminidb):
	gemini_amend = subprocess.Popen([gemini, 'amend', '--sample', pedfile, geminidb])
	gemini_amend.wait()

def make_ad_files(in_ped, out_ped, in_db, out_db):
	##make ad ped, so convert 1's to 0's in ped
	with open(out_ped, "w") as out_fh, open(in_ped, "r") as in_fh:
		for line in in_fh:
			line = line.strip('\n').split(delim)
			pheno = line[5]
			the_rest = line[:5]
			if pheno == '1':
				pheno = '0'
			out_fh.write(delim.join(the_rest + [pheno, '\n']))
	##cp gemini database then add ad ped file
	cp_gemini_db(in_db, out_db)
	add_new_ped_file_to_gemini_db(out_ped, out_db)

def load_single_vcf_into_gemini(ped, db_name, input_vcf, ped_file):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	temp_vcf = ped + 'temp1.vcf'
	normalized_vcf = ped + '.int.norm.vcf'
	# ped_file = ped + '.ped'
	with open (temp_vcf, 'w') as tvcf_fh:
		zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		vt_decompose = subprocess.Popen([vt, 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
		vt_normalize = subprocess.Popen([vt, 'normalize', '-r', fasta, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
		vt_normalize.wait()
	##annotate with snpeff and compress and index
	with open (normalized_vcf, 'w') as nvcf_fh:
		snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'GRCh37.75', '-v', '-formatEff', '-classic', temp_vcf], stdout=nvcf_fh)
		snpeff_vcf.wait()
		bgzip_vcf = subprocess.Popen(['bgzip', normalized_vcf])
		bgzip_vcf.wait()
		tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', normalized_vcf + '.gz'])
		tabix_vcf.wait()
	##add file to gemini, requires ped file
	gemini_temp_dir = './' + ped + '.gemini_load_dir'
	gemini_load = subprocess.Popen([gemini, 'load', '--cores', '20', '-t', 'snpEff', '-p', ped_file, '--tempdir', gemini_temp_dir,  '-v', normalized_vcf + '.gz', db_name])
	# gemini_load = subprocess.Popen([gemini, 'load', '--cores', '5', '-t', 'snpEff', '-p', ped_file, '--tempdir', gemini_temp_dir,  '-v', normalized_vcf + '.gz', db_name])
	gemini_load.wait()
	##remove intermediate files
	os.remove(temp_vcf)

def gemini_comp_het(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.comp_hets.xls'
	with open (outfile, 'w') as out_fh:
		##compound het: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		# gemini_query = subprocess.Popen([gemini, 'comp_hets', '--max-priority', '2', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter',  "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf) , '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'comp_hets', '--max-priority', '2', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter',  "impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf) , '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_recessive(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.autosomal_recessive.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		# gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile

def gemini_denovo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		# gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile

##not used atm
def gemini_dominant(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.autosomal_dominant.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		# gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile

def gemini_inherited(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.inher.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		# gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile


def gemini_xlinked(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.x_linked_recessive.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		# gemini_query = subprocess.Popen([gemini, 'x_linked_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'x_linked_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile

def gemini_xlinked_de_novo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.x_linked_de_novo.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		# gemini_query = subprocess.Popen([gemini, 'x_linked_de_novo', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'x_linked_de_novo', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile

def gemini_de_novo_syn(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo_syn.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		# gemini_query = subprocess.Popen([gemini, 'de_novo', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity == 'LOW' and is_coding == 1 and gerp_bp_score >= 2 AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "impact_severity == 'LOW' and is_coding == 1 and gerp_bp_score >= 2 AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile

def gemini_recessive_syn(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.auto_rec_syn.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		# gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity == 'LOW' and is_coding == 1 and gerp_bp_score >= 2 AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "impact_severity == 'LOW' and is_coding == 1 and gerp_bp_score >= 2 AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile


def combine_gemini_results(ped, file_list, add_annovar_pos, out_suffix, position_to_finish):
	##file names
	outfile = ped + out_suffix
	temp_file = 'combined.temp.xls'
	# temp2_file = 'combined.temp2.xls'
	print 'combining files:', file_list
	##combine all files, making sure col lengths are all the same
	with open(temp_file, "w") as temp_fh:
		file_count, total_line_count = 0, 0
		for filename in file_list:
			# print filename
			analysis_type = filename.split('.')[-2]
			# file_count += 1
			# print analysis_type, file_count
			##add header from first file
			with open(filename, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line_count +=1
					total_line_count += 1
					line = line.rstrip().split(delim)
					if line_count == 1:
						file_count += 1
						if file_count == 1:
							header = delim.join(line[:position_to_finish] + ['analysis', '\n'])
							temp_fh.write(header)
					else:
						line_out = delim.join(line[:position_to_finish] + [analysis_type, '\n'])
						temp_fh.write(line_out)
	if total_line_count >=2:
		add_annovar_ts_info(ped, temp_file, outfile, add_annovar_pos)
		##remove intermediate files
		os.remove(temp_file)

def add_annovar_ts_info(pedigree, in_file, out_file, pos_to_insert):
	##files
	normalized_vcf = pedigree + '.int.norm.vcf.gz'
	region_file = pedigree + '_temp.regions'
	temp_vcf = pedigree + 'temp0.vcf.gz'
	av_file = pedigree + '.avinput'
	multianno = pedigree + '.hg19_multianno.txt'
	## make_list_of_regions(var_file, region_file):
	with open(in_file, 'r') as var_fh, open(region_file, 'w') as reg_fh:
		line_count = 0
		for line in var_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count >1:
				chrom = line[0][3:]
				# start = int(line[1]) - 2
				# end = int(line[2]) + 2
				start = line[1]
				end = line[2]
				reg_fh.write(delim.join([chrom, str(start), str(end), '\n']))
				# print chrom, start_end
	## get_var_from_regions_in_vcf(in_vcf, regions, out_vcf):
	bcftools_filter = subprocess.Popen([bcftools_12, 'view', '-R', region_file, '-o', temp_vcf, '-O', 'z', normalized_vcf])
	bcftools_filter.wait()
	##convert vcf file to individual avinput file
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', temp_vcf, '-allsample', '--withfreq', '-outfile', av_file])
	con_ann.wait()
	##run_table_annovar(avinput, av_prefix):
	command = [table_annovar] + av_buildver + [av_file] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', pedigree]
	annovar = subprocess.Popen(command)
	annovar.wait()
	## add_annovar_info(var_file, multianno, outfile, pos_to_insert):
	##make dict with annovar data
	ann_dict = {}
	with open(multianno, 'r') as multi_fh:
		for line in multi_fh:
			line = line.rstrip().split(delim)
			chrom = line[0]
			pos = line[2]
			ch_pos = '_'.join([chrom,pos])
			ref = line[3]
			alt = line[4]
			ts_info = line[9]
			if ts_info == '.':
				ts_info = line[7]
			# print chrom, pos, ts_info
			ann_dict[ch_pos] = ts_info

	with open(in_file, 'r') as in_fh, open(out_file, 'w') as out_fh:
		line_count = 0
		for line2 in in_fh:
			line2 = line2.split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join(line2[:pos_to_insert] + ['annovar_ts_info'] + line2[pos_to_insert:]))
			else:
				chrom2 = line2[0][3:]
				pos2 = line2[2]
				ch_pos2 = '_'.join([chrom2,pos2])
				if ch_pos2 in ann_dict:
					# print line2[:5],ch_pos2, ann_dict[ch_pos2]
					out_fh.write(delim.join(line2[:pos_to_insert] + [ann_dict[ch_pos2]] + line2[pos_to_insert:]))
				else:
					print 'not found:', line2[:5], ch_pos2

def format_autosomal_dominant_analysis(ped, file_list, add_annovar_pos, position_to_finish):
	##file names
	outfile = ped + '.auto_dominant_analysis.xls'
	temp_file = 'combined.temp.xls'
	# temp2_file = 'combined.temp2.xls'
	# temp3_file = 'combined.temp3.xls'
	print 'combining files:', file_list
	##combine all files, making sure col lengths are all the same
	with open(temp_file, "w") as temp_fh:
		file_count = 0
		for filename in file_list:
			print filename
			##add header from first file
			with open(filename, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line_count +=1
					line = line.rstrip().split(delim)
					if line_count == 1:
						file_count += 1
						if file_count == 1:
							header = delim.join(line[:position_to_finish] + ['\n'])
							temp_fh.write(header)
					else:
						line_out = delim.join(line[:position_to_finish] + ['\n'])
						temp_fh.write(line_out)

	add_annovar_ts_info(ped, temp_file, outfile, add_annovar_pos)
	##remove intermediate files
	os.remove(temp_file)

##run all other methods
def standard_gemini_protocol(pedigree, ped_type):
	gemini_db = pedigree + '.gemini.db'
	ad_db = pedigree + '.gemini_ad.db'
	in_vcf = pedigree + '.vcf.gz'
	std_ped_file = pedigree + '.ped'
	ad_ped_file = pedigree + '_ad.ped'
	# freq_req = 0.01
	freq_req = 0.1
	freq_req_recessive = 0.1
	# coverage_req = 10
	coverage_req = 5
	##load vcf file into db
	# load_single_vcf_into_gemini(pedigree, gemini_db, in_vcf, std_ped_file)
	##make auto dom ped and gemini db for duo and trio analysis
	if ped_type != 'singleton':
		make_ad_files(std_ped_file, ad_ped_file, gemini_db, ad_db)
	# '''
	##run queries and get files to combine
	files_to_combine = []
	if ped_type == 'singleton':
		files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_dominant(pedigree, gemini_db, freq_req, coverage_req))
	elif ped_type == 'duo':
		files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_dominant(pedigree, ad_db, freq_req, coverage_req))
	elif ped_type == 'trio':
		files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_denovo(pedigree, gemini_db, freq_req, coverage_req))
		files_to_combine.append(gemini_xlinked(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_xlinked_de_novo(pedigree, gemini_db, freq_req, coverage_req))
		auto_dom_file = gemini_dominant(pedigree, ad_db, freq_req, coverage_req)
		format_autosomal_dominant_analysis(pedigree, [auto_dom_file], 13, 22)
	else:
		print 'ped type %s not recognized'%ped_type
	##combine results
	combine_gemini_results(pedigree, files_to_combine, 6, '.std_analysis.xls', 22)
	# '''



def get_vars_in_candidate_regions(analysis_dict, reg_bed, out_prefix):
	combined_outfile = out_prefix + '.vars_in_regions.xls'
	combined_avinput = out_prefix + '.vars_in_regions.avinput'
	ann_txt = out_prefix + '.hg19_multianno.txt'
	filtered_vars = out_prefix + '.vars.rare_gerp_rpts.xls'
	'''
	region_list = []
	with open(reg_bed, "r") as bed_fh:
		for line in bed_fh:
			line = line.rstrip().split(delim)
			region = line[0] + ':' + line[1] + '-' + line[2]
			region_list.append(region)
	for ped in analysis_dict:
		outfile = ped + '.vars_in_regions.xls'
		with open(outfile, "w") as out_fh:
			rc = 0
			for reg in region_list:
				rc += 1
				# gemini_query = subprocess.Popen([gemini, 'query', '--region', reg, '-q', '"select *,(gts).(*) from variants"', ped + '.gemini.db'], stdout=out_fh)
				if rc == 1:
					gemini_query = subprocess.Popen([gemini, 'query', '--region', reg, '-q', "select *,(gts).(*) from variants", '--header', ped + '.gemini.db'], stdout=out_fh)
					gemini_query.wait()
				else:
					gemini_query = subprocess.Popen([gemini, 'query', '--region', reg, '-q', "select *,(gts).(*) from variants", ped + '.gemini.db'], stdout=out_fh)
					gemini_query.wait()

	##combined the file
	with open(combined_outfile, "w") as cout_fh:
		pc = 0
		for ped in analysis_dict:
			infile = ped + '.vars_in_regions.xls'
			pc += 1
			lc = 0
			with open(infile, "r") as in_fh:
				for line in in_fh:
					lc += 1
					if lc == 1:
						if pc == 1:
							cout_fh.write('ped' + '\t' + line)
					else:
						cout_fh.write(ped + '\t' + line)
	'''
	##make avinput
	with open(combined_outfile, "r") as cin_fh, open(combined_avinput, "w") as cout_fh:
		lc = 0
		for line in cin_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc >1:
				cout_fh.write(delim.join([line[1][3:], str(int(line[2]) + 1), line[3]] + line[7:9] + [line[0], ','.join(line[152:])]) + '\n')
	##run_table_annovar(avinput, av_prefix):
	command = [table_annovar] + av_buildver + [combined_avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()
	
	with open(ann_txt, "r") as in_fh, open(filtered_vars, "w") as out_fh:
		lc, fc = 0, 0
		for line in in_fh:
			lc += 1
			line = line.split(delim)
			rmsk = line[11]
			segdup = line[12]
			gerp = line[10]
			af = line[13]
			if lc == 1:
				out_fh.write(delim.join(line))
			else:
				if af == '.':
					af = 0
				else:
					af = float(af)
				# if rmsk == '.' and segdup == '.' and gerp != '.' and af <= 0.01:
				if rmsk == '.' and segdup == '.' and gerp != '.' and af <= 0.05:
					out_fh.write(delim.join(line))
					fc += 1
	print(lc,fc)

def get_vars_in_abca4_locus(analysis_dict, out_prefix):
	combined_outfile = out_prefix + '.vars_in_abca4.xls'
	combined_avinput = out_prefix + '.vars_in_abca4.avinput'
	ann_txt = out_prefix + '.hg19_multianno.txt'
	filtered_vars = out_prefix + '.final.vars_in_abca4.xls'
	'''
	reg = 'chr1:94443600-94610800'
	for ped in analysis_dict:
		outfile = ped + '.vars_in_abca4.xls'
		with open(outfile, "w") as out_fh:
			gemini_query = subprocess.Popen([gemini, 'query', '--region', reg, '-q', "select *,(gts).(*) from variants", '--header', ped + '.gemini.db'], stdout=out_fh)
			gemini_query.wait()
	##combined the files
	with open(combined_outfile, "w") as cout_fh:
		pc = 0
		for ped in analysis_dict:
			infile = ped + '.vars_in_abca4.xls'
			pc += 1
			lc = 0
			with open(infile, "r") as in_fh:
				for line in in_fh:
					lc += 1
					if lc == 1:
						if pc == 1:
							cout_fh.write('ped' + '\t' + line)
					else:
						cout_fh.write(ped + '\t' + line)
	
	##make avinput
	with open(combined_outfile, "r") as cin_fh, open(combined_avinput, "w") as cout_fh:
		lc = 0
		for line in cin_fh:
			lc += 1
			line = line.rstrip().split(delim)
			if lc >1:
				cout_fh.write(delim.join([line[1][3:], str(int(line[2]) + 1), line[3]] + line[7:9] + [line[0], ','.join(line[152:])]) + '\n')
	##run_table_annovar(avinput, av_prefix):
	command = [table_annovar] + av_buildver + [combined_avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()
	'''
	with open(ann_txt, "r") as in_fh, open(filtered_vars, "w") as out_fh:
		for line in in_fh:
			out_fh.write(line)



def run_manta_svs(bam_files, project_name, working_dir):
	bams = []
	for bam_file in bam_files:
		in_bam = ['--bam', bam_file]
		bams.extend(in_bam)
	##call svs
	run_manta_config = subprocess.Popen([manta_config] + bams + ['--referenceFasta', fasta, '--runDir', working_dir + '/' + project_name ])
	run_manta_config.wait()
	manta_run = subprocess.Popen([working_dir + '/' + project_name + '/runWorkflow.py'])
	manta_run.wait()


def get_manta_result_around_abca4(peds, out_file, vcf_suffix):
	with open(out_file, 'w') as out_fh:
		header = ['PED','CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GTs', '\n']
		out_fh.write(delim.join(header))
		for ped in peds:
			print ped
			in_vcf = 'manta_' + ped + vcf_suffix
			with gzip.open(in_vcf, 'rb') as vcf_fh:
				for line in vcf_fh:
					line = line.split(delim)
					##abca4 gene +- 100k
					if line[0] == '1' and int(line[1]) >= 94358393 and int(line[1]) <= 94686705:
						out_fh.write(delim.join([ped] + line))








##master method
def run_analysis(working_dir, analysis_file, combined_vcf, reg_bed, bam_dict ):
	os.chdir(working_dir)
	project = analysis_file.split('.')[0]
	manta_results_file = project + '.manta.ABCA4.xls'
	##read in analysis file and make 2 dicts
	ped_file_dict, analysis_dict = make_analysis_dicts(analysis_file)
	##check...
	for p in analysis_dict:
		print(p, analysis_dict[p])
	'''
	##make ped files from first dict
	make_ped_files(ped_file_dict)
	##make vcf per ped
	# get_vcf_per_ped(analysis_dict,combined_vcf)
	##analyze data
	for ped in analysis_dict:
		standard_gemini_protocol(ped, analysis_dict[ped][1])
	get_vars_in_candidate_regions(analysis_dict, reg_bed, project)

	get_vars_in_abca4_locus(analysis_dict, project)
	get_vars_in_candidate_regions(analysis_dict, reg_bed, project)
	
	for ped in bam_dict:
		run_manta_svs(bam_dict[ped], 'manta_' + ped, working_dir)
	'''
	manta_vcf_suffix = '/results/variants/diploidSV.vcf.gz'
	get_manta_result_around_abca4(bam_dict, manta_results_file, manta_vcf_suffix)
##run methods

##run all at once
work_dir = '/home/atimms/ngs_data/genomes/cherry_genomes_1119'
# input_file = 'stargardt_genomes_test_1119.txt'
# input_file = 'stargardt_genomes_2_1119.txt'
input_file = 'stargardt_genomes_1119.txt'
uw_vcf = 'cherry_uwcmg_stgd_1.HF.final.vcf.gz'
region_bed = 'ABCA4_PROM1_PRPH2_ELOVL4_Ret_RPE_cCREs_hg19.bed'
ped_bam_dict = {'506': ['3853.328830.bam', '3854.328831.bam', '3855.328832.bam'], '529': ['3923.328834.bam', '3933.328835.bam','3905.328833.bam'],
		'568': ['3983.328837.bam', '3984.328838.bam', '3982.328836.bam'], '595': ['3150.328839.bam','9035.328840.bam', '9036.328841.bam'],
		'586': ['3939.328828.bam', '9024.328829.bam'], '82': ['3073.328842.bam'],'147': ['3212.328826.bam'], '228': ['6793.328885.bam'], '328853': ['3621.328853.bam'],
		'456': ['4429.328871.bam'], '457': ['4326.328867.bam'], '459': ['3733.328856.bam'], '461': ['3119.328844.bam'], '565': ['4565.337259.bam'], 
		'599': ['4638.328878.bam'], '636': ['4744.328880.bam'], '3118': ['3118.328843.bam'], '3133': ['3133.328845.bam'], '3197': ['3197.328847.bam'], 
		'3207': ['3207.337257.bam'], '3497': ['3497.328850.bam'], '3579': ['3579.328851.bam'], '3603': ['3603.328852.bam'], '3636': ['3636.328854.bam'], 
		'3706': ['3706.328855.bam'], '3898': ['3898.328857.bam'], '3901': ['3901.328858.bam'], '3943': ['3943.328859.bam'], '3967': ['3967.328860.bam'], 
		'3980': ['3980.328861.bam'], '3992': ['3992.328862.bam'], '4120': ['4120.328863.bam'], '4178': ['4178.328864.bam'], '4203': ['4203.328865.bam'], 
		'4280': ['4280.328866.bam'], '4367': ['4367.328868.bam'], '4410': ['4410.328869.bam'], '4422': ['4422.328870.bam'], '4487': ['4487.328872.bam'], 
		'4498': ['4498.328873.bam'], '4516': ['4516.328874.bam'], '4541': ['4541.328875.bam'], '4560': ['4560.328876.bam'], '4720': ['4720.328879.bam'], 
		'5043': ['5043.328881.bam'], '5077': ['5077.328882.bam'], '5222': ['5222.328883.bam'], '5268': ['5268.328884.bam'], '9045': ['9045.328886.bam'], 
		'9071': ['9071.328887.bam'], '9093': ['9093.328888.bam'], '9113': ['9113.328889.bam'], 'Y256': ['Y256.328890.bam']}


# bams_for_sv_analysis = ['3073.328842.bam', '3923.328834.bam']
run_analysis(work_dir, input_file, uw_vcf, region_bed, ped_bam_dict)



