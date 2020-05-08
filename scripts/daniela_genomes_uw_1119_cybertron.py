#!/usr/bin/env python
import sys
import subprocess
import os
import glob


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
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
ref_dir = '/home/atimms/ngs_data/references/hg19/'
plink = '/home/atimms/programs/plink'

##unique for this project
fasta = ref_dir + 'hs37d5.fa'
gemini_ref_dir = '/home/atimms/ngs_data/references/gemini/'
rvis_data = gemini_ref_dir + 'RVIS_ExAC_4KW.txt'
gdi_data = gemini_ref_dir + 'GDI_full_10282015.txt'
mgi_data = gemini_ref_dir + 'mgi.abnormal_outer_ear_morhology.txt'
hpo_data = gemini_ref_dir + 'genes_for_HP_0000356.csv'
sbse_rnaseq_data = gemini_ref_dir + 'esra.gene_exp.diff'
hmx1_rnaseq_data = gemini_ref_dir + 'jess.gene_exp.diff'
human_ear_rnaseq_data = gemini_ref_dir + 'human_ear.gene_exp.diff'
human_syndrome_data = gemini_ref_dir + 'microtia_human_syndromes.txt'
hoxa2_data = gemini_ref_dir + 'hoxa2.expression_data.txt'
ba_exp_data = gemini_ref_dir + 'ba_expression_papers_0816.txt'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene']
av_operation = ['-operation', 'g']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-arg', '-splicing 10 ']

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
	# gemini_load = subprocess.Popen([gemini, 'load', '--cores', '20', '-t', 'snpEff', '-p', ped_file, '--tempdir', gemini_temp_dir,  '-v', normalized_vcf + '.gz', db_name])
	gemini_load = subprocess.Popen([gemini, 'load', '--cores', '5', '-t', 'snpEff', '-p', ped_file, '--tempdir', gemini_temp_dir,  '-v', normalized_vcf + '.gz', db_name])
	gemini_load.wait()
	##remove intermediate files
	os.remove(temp_vcf)

def gemini_comp_het(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.comp_hets.xls'
	with open (outfile, 'w') as out_fh:
		##compound het: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'comp_hets', '--max-priority', '2', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter',  "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf) , '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_recessive(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.autosomal_recessive.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
# 		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', gem_db_name], stdout=out_fh)

		gemini_query.wait()
	return outfile

def gemini_denovo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

##not used atm
def gemini_dominant(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.autosomal_dominant.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_inherited(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.inher.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile


def gemini_xlinked(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.x_linked_recessive.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'x_linked_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_xlinked_de_novo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.x_linked_de_novo.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'x_linked_de_novo', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_de_novo_syn(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo_syn.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity == 'LOW' and is_coding == 1 and gerp_bp_score >= 2 AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def gemini_recessive_syn(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.auto_rec_syn.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity == 'LOW' and is_coding == 1 and gerp_bp_score >= 2 AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

def make_var_dict_from_gemini_result(gemini_file):
	var_dict = {}
	line_count = 0
	with open(gemini_file, 'r') as in_fh:
		for line in in_fh:
			line = line.split(delim)
			line_count += 1
			if line_count > 1:
				var = '_'.join(line[:3] + line[6:8])
				# print var
				var_dict[var] = 1
	return var_dict

def gemini_potential_cnv(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.potential_cnv.xls'
	all_auto_rec_file = 'all_auto_rec.temp.xls'
	correct_auto_rec_file = 'correct_auto_rec.temp.xls'
	with open (correct_auto_rec_file, 'w') as car_fh:
		gemini_query_a = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=car_fh)
		gemini_query_a.wait()
	with open (all_auto_rec_file, 'w') as aar_fh:
		gemini_query_b = subprocess.Popen([gemini, 'autosomal_recessive', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=aar_fh)
		gemini_query_b.wait()
	correct_cnv_var_dict = make_var_dict_from_gemini_result(correct_auto_rec_file)
	with open (outfile, 'w') as out_fh:
		with open(all_auto_rec_file, 'r') as all_ar_fh:
			line_count = 0
			for line in all_ar_fh:
				line = line.split(delim)
				line_count += 1
				if line_count == 1:
					out_fh.write(delim.join(line))
				else:
					var = '_'.join(line[:3] + line[6:8])
					chrom  = line[0]
					if chrom != 'chrX':
						if var not in correct_cnv_var_dict:
							out_fh.write(delim.join(line))
	##remove intermediate files
	os.remove(all_auto_rec_file)
	os.remove(correct_auto_rec_file)
	return outfile

##make dict using exac rvis file
def make_dict_from_exac_intolerence(intol_file):
	intol_dict = {}
	with open(intol_file, "U") as intol:
		line_count = 0
		for line in intol:
			line = line.strip('\n').split(delim)
			line_count += 1
			gene_list = line[0:5]
			genes = list(set(gene_list))
			# print genes, len(genes)
			rvis = line[5]
			percentile = line[6]
			# print genes, rvis, percentile
			for gene in genes:
				if gene not in intol_dict:
					intol_dict[gene] = [rvis, percentile]
					# print gene, rvis, percentile
				else:
					print 'gene %s seen multiple times'% gene
	return intol_dict

##make dict using gdi file
def make_dict_from_gdi(gdi_file):
	gdi_dict = {}
	with open(gdi_file, "U") as gdi_fh:
		for line in gdi_fh:
			line = line.strip('\n').split(delim)
			# print line
			gene = line[0]
			gdi = line[1]
			gdi_phred = line[2]
			if gene not in gdi_dict:
				gdi_dict[gene] = [gdi, gdi_phred]
				# print gene, gdi, gdi_phred
			else:
				print 'gene %s seen multiple times'% gene
	return gdi_dict

def add_rvis_and_GDI(in_file, out_file, dicts_to_add, gene_col, position_to_insert):
	with open(in_file, "U") as inf, open(out_file, "w") as outf:
		line_count, genes_not_found = 0, 0
		for line in inf:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				extra_header = ['RVIS exac', 'Percentile exac', 'GDI', 'GDI phred', 'mgi phenotype', 'hpo_syndrome', 'sbse_wt_fpkm', 'sbse_log2_fc', 'sbse_q_value', 
						'hmx1_wt_fpkm', 'hmx1_log2_fc', 'hmx1_q_value', 'ear_57d_fpkm', 'ear_59d_fpkm', 'human_syndrome', 'hoxa2_wt1', 'hoxa2_wt3', 'ba_expression_papers_0816']
				outf.write(delim.join(line[:position_to_insert] + extra_header + line[position_to_insert:] + ['\n']))
			else:
				gene = line[gene_col]
				# print gene
				##get values for all gene names and print minimum value
				extra_stuff = []
				for data_dict in dicts_to_add:
					if gene in data_dict:
						extra_stuff.extend(data_dict[gene])
					else:
						extra_stuff.extend(['na']*len(data_dict[data_dict.keys()[0]]))

				# print gene, extra_stuff
				outf.write(delim.join(line[:position_to_insert] + extra_stuff + line[position_to_insert:] + ['\n']))

##make dict using mgi file
def make_dict_from_mgi(mgi_file):
	mgi_dict = {}
	mgi_dict2 = {}
	with open(mgi_file, "U") as mgi_fh:
		for line in mgi_fh:
			line = line.strip('\n').split(delim)
			# print line
			gene = line[2]
			pheno = line[3]
			if gene not in mgi_dict:
				mgi_dict[gene] = [pheno]
				# print gene, inheritance, pheno, syndrome, loe
			else:
				mgi_dict[gene].append(pheno)
	for g in mgi_dict:
		mgi_dict2[g] = [', '.join(mgi_dict[g])]
	# for g in mgi_dict:
	# 	print g, mgi_dict[g], mgi_dict2[g]
	return mgi_dict2

def make_dict_from_hpo(hpo_file):
	hpo_dict = {}
	with open(hpo_file, "U") as hpo_fh:
		line_count = 0
		for line in hpo_fh:
			line_count += 1
			if line_count > 2:
				line = line.strip('\n').split(',',1)

				gene = line[0].strip('"').split(' ')[0]
				syndrome = line[1].strip('"')
				# print gene, syndrome
				if gene not in hpo_dict:
					hpo_dict[gene] = [syndrome]
					# print gene, gdi, gdi_phred
				else:
					print 'gene %s seen multiple times'% gene
	return hpo_dict

def make_dict_from_cufflinks_rnaseq(gene_diff_file, cols_req):
	rnaseq_dict = {}
	with open(gene_diff_file, "U") as gdf_fh:
		for line in gdf_fh:
			line = line.strip('\n').split(delim)
			# print line
			gene = line[2]
			data = []
			for col in cols_req:
				d = line[col]
				data.append(d)
			# print gene, data
			if gene != '-':
				if gene not in rnaseq_dict:
					rnaseq_dict[gene] = data
					# print gene, gdi, gdi_phred
				else:
					# print 'gene %s seen multiple times'% gene
					pass
		return rnaseq_dict

def make_dict_from_human_syndrome(human_syndrome_file):
	hum_syn_dict = {}
	with open(human_syndrome_file, "U") as hs_fh:
		line_count = 0
		for line in hs_fh:
			line_count += 1
			if line_count > 1:
				line = line.strip('\n').split(delim)
				gene = line[0]
				syndrome = line[1]
				# print gene, syndrome
				if gene not in hum_syn_dict:
					hum_syn_dict[gene] = [syndrome]
				else:
					# print 'gene %s seen multiple times'% gene
					pass
	return hum_syn_dict

def make_dict_from_hoxa2(hoxa2_file):
	hoxa2_dict = {}
	with open(hoxa2_file, "U") as hoxa2_fh:
		line_count = 0
		for line in hoxa2_fh:
			line_count += 1
			if line_count > 1:
				line = line.strip('\n').split(delim)
				gene = line[1]
				data = line[2:4]
				# print gene, data
				if gene not in hoxa2_dict:
					hoxa2_dict[gene] = data
				else:
					# print 'gene %s seen multiple times'% gene
					pass
	return hoxa2_dict

def make_dict_from_ba_expression_papers(ba_exp_file):
	ba_dict = {}
	ba_dict2 = {}
	with open(ba_exp_file, "U") as ba_fh:
		for line in ba_fh:
			line = line.strip('\n').split(delim)
			# print line
			gene = line[0]
			paper = line[1]
			if gene not in ba_dict:
				ba_dict[gene] = [paper]
				# print gene, inheritance, pheno, syndrome, loe
			else:
				ba_dict[gene].append(paper)
	for g in ba_dict:
		ba_dict2[g] = [', '.join(ba_dict[g])]
	# for g in mgi_dict:
	# 	print g, mgi_dict[g], mgi_dict2[g]
	return ba_dict2

def add_columns_to_file(in_file, out_file, dicts_to_add, gene_col, position_to_insert):
	with open(in_file, "U") as inf, open(out_file, "w") as outf:
		line_count, genes_not_found = 0, 0
		for line in inf:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				extra_header = ['RVIS exac', 'Percentile exac', 'GDI', 'GDI phred', 'mgi phenotype', 'hpo_syndrome', 'sbse_wt_fpkm', 'sbse_log2_fc', 'sbse_q_value', 
						'hmx1_wt_fpkm', 'hmx1_log2_fc', 'hmx1_q_value', 'ear_57d_fpkm', 'ear_59d_fpkm', 'human_syndrome', 'hoxa2_wt1', 'hoxa2_wt3', 'ba_expression_papers_0816']
				outf.write(delim.join(line[:position_to_insert] + extra_header + line[position_to_insert:] + ['\n']))
			else:
				gene = line[gene_col]
				# print gene
				##get values for all gene names and print minimum value
				extra_stuff = []
				for data_dict in dicts_to_add:
					if gene in data_dict:
						extra_stuff.extend(data_dict[gene])
					else:
						extra_stuff.extend(['na']*len(data_dict[data_dict.keys()[0]]))

				# print gene, extra_stuff
				outf.write(delim.join(line[:position_to_insert] + extra_stuff + line[position_to_insert:] + ['\n']))


def combine_gemini_results(ped, file_list, gene_col, position_to_insert, position_to_finish, add_annovar_pos, out_suffix):
	##file names
	outfile = ped + out_suffix
	temp_file = 'combined.temp.xls'
	temp2_file = 'combined.temp2.xls'
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
							header = delim.join(line[:position_to_insert] + line[position_to_insert:position_to_finish] + ['analysis', '\n'])
							temp_fh.write(header)
					else:
						line_out = delim.join(line[:position_to_insert] + line[position_to_insert:position_to_finish] + [analysis_type, '\n'])
						temp_fh.write(line_out)
	if total_line_count >=2:
		##add dbdb and mgi data
		##add rvis, gdi, mgi data
		rvis_data_dict = make_dict_from_exac_intolerence(rvis_data)
		gdi_data_dict = make_dict_from_gdi(gdi_data)
		mgi_data_dict = make_dict_from_mgi(mgi_data)
		hpo_data_dict = make_dict_from_hpo(hpo_data)
		sbse_rnaseq_dict = make_dict_from_cufflinks_rnaseq(sbse_rnaseq_data, [8,9,11])
		hmx1_rnaseq_dict = make_dict_from_cufflinks_rnaseq(hmx1_rnaseq_data, [8,9,11])
		human_ear_rnaseq_dict = make_dict_from_cufflinks_rnaseq(human_ear_rnaseq_data, [7,8])
		human_syndrome_dict = make_dict_from_human_syndrome(human_syndrome_data)
		hoxa2_dict = make_dict_from_hoxa2(hoxa2_data)
		ba_paper_dict = make_dict_from_ba_expression_papers(ba_exp_data)
		add_columns_to_file(temp_file, temp2_file, [rvis_data_dict, gdi_data_dict, mgi_data_dict, hpo_data_dict, sbse_rnaseq_dict, hmx1_rnaseq_dict, human_ear_rnaseq_dict, human_syndrome_dict, hoxa2_dict, ba_paper_dict], gene_col, position_to_insert)
		# add_columns_to_file(temp_file, temp2_file, dbdb_data_dict, mgi_data_dict, omim_data_dict, gene_col, position_to_insert)
		add_annovar_ts_info(ped, temp2_file, outfile, add_annovar_pos)
		##remove intermediate files
		os.remove(temp_file)
		os.remove(temp2_file)

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

def format_autosomal_dominant_analysis(ped, file_list, gene_col, position_to_insert, position_to_finish, add_annovar_pos, add_aaf_pos):
	##file names
	outfile = ped + '.auto_dominant_analysis.xls'
	temp_file = 'combined.temp.xls'
	temp2_file = 'combined.temp2.xls'
	temp3_file = 'combined.temp3.xls'
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
							header = delim.join(line[:position_to_insert] + line[position_to_insert:position_to_finish] + ['\n'])
							temp_fh.write(header)
					else:
						line_out = delim.join(line[:position_to_insert] + line[position_to_insert:position_to_finish] + ['\n'])
						temp_fh.write(line_out)

	##add rvis, gdi, mgi data
	rvis_data_dict = make_dict_from_exac_intolerence(rvis_data)
	gdi_data_dict = make_dict_from_gdi(gdi_data)
	mgi_data_dict = make_dict_from_mgi(mgi_data)
	hpo_data_dict = make_dict_from_hpo(hpo_data)
	sbse_rnaseq_dict = make_dict_from_cufflinks_rnaseq(sbse_rnaseq_data, [8,9,11])
	hmx1_rnaseq_dict = make_dict_from_cufflinks_rnaseq(hmx1_rnaseq_data, [8,9,11])
	human_ear_rnaseq_dict = make_dict_from_cufflinks_rnaseq(human_ear_rnaseq_data, [7,8])
	human_syndrome_dict = make_dict_from_human_syndrome(human_syndrome_data)
	hoxa2_dict = make_dict_from_hoxa2(hoxa2_data)
	ba_paper_dict = make_dict_from_ba_expression_papers(ba_exp_data)
	add_columns_to_file(temp_file, temp2_file, [rvis_data_dict, gdi_data_dict, mgi_data_dict, hpo_data_dict, sbse_rnaseq_dict, hmx1_rnaseq_dict, human_ear_rnaseq_dict, human_syndrome_dict, hoxa2_dict, ba_paper_dict], gene_col, position_to_insert)
	# add_columns_to_file(temp_file, temp2_file, dbdb_data_dict, mgi_data_dict, omim_data_dict, gene_col, position_to_insert)
	add_annovar_ts_info(ped, temp2_file, outfile, add_annovar_pos)
	##remove intermediate files
	os.remove(temp_file)
	os.remove(temp2_file)

##run all other methods
def standard_gemini_protocol(pedigree, ped_type):
	gemini_db = pedigree + '.gemini.db'
	ad_db = pedigree + '.gemini_ad.db'
	in_vcf = pedigree + '.vcf.gz'
	std_ped_file = pedigree + '.ped'
	ad_ped_file = pedigree + '_ad.ped'
	freq_req = 0.01
	freq_req_recessive = 0.05
	coverage_req = 10
	##load vcf file into db
	load_single_vcf_into_gemini(pedigree, gemini_db, in_vcf, std_ped_file)
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
		files_to_combine.append(gemini_de_novo_syn(pedigree, gemini_db, freq_req, coverage_req))
		files_to_combine.append(gemini_recessive_syn(pedigree, gemini_db, freq_req_recessive, coverage_req))
		files_to_combine.append(gemini_potential_cnv(pedigree, gemini_db, freq_req_recessive, coverage_req))
		auto_dom_file = gemini_dominant(pedigree, ad_db, freq_req, coverage_req)
		format_autosomal_dominant_analysis(pedigree, [auto_dom_file], 5, 13, 22, 13, 31)
	else:
		print 'ped type %s not recognized'%ped_type
	##combine results
	combine_gemini_results(pedigree, files_to_combine, 5, 13, 22, 6, '.std_analysis.xls')
	# '''



##master method
def run_analysis(working_dir, analysis_file, combined_vcf):
	os.chdir(working_dir)
	project = analysis_file.split('.')[0]
	##read in analysis file and make 2 dicts
	ped_file_dict, analysis_dict = make_analysis_dicts(analysis_file)
	##check...
	for p in analysis_dict:
		print(p, analysis_dict[p])
	# '''
	##make ped files from first dict
	make_ped_files(ped_file_dict)
	##make vcf per ped
	get_vcf_per_ped(analysis_dict,combined_vcf)
	# '''
	##analyze data
	for ped in analysis_dict:
		standard_gemini_protocol(ped, analysis_dict[ped][1])


def combine_files(peds, out_file, infile_suffix):
	with open(out_file, "w") as out_fh:
		print('making file %s by combining:'%out_file)
		file_count = 0
		for ped in peds:
			infile = ped + infile_suffix
			print(infile)
			# file_count += 1
			##add header from first file
			with open(infile, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line_count +=1
					if line_count == 1:
						file_count += 1
						if file_count == 1:
							out_fh.write(line)
					else:
						out_fh.write(line)

def combine_results_files(peds_to_analyze, work_dir):
	os.chdir(work_dir)
	combined_trio_std = 'combined.trio.std_analysis.xls'
	combined_trio_dom = 'combined.trio.autosomal_dominant.xls'
	combined_duo_single = 'combined.singleton_duo.std_analysis.xls'
	trio_peds = []
	single_duo_peds = []
	for ped in peds_to_analyze:
		if ped in duo_list or ped in single_list:
			single_duo_peds.append(ped)
		else:
			trio_peds.append(ped)
	print(trio_peds)
	print(single_duo_peds)
	combine_files(trio_peds, combined_trio_std, '.std_analysis.xls')
	combine_files(trio_peds, combined_trio_dom, '.auto_dominant_analysis.xls')
	combine_files(single_duo_peds, combined_duo_single, '.std_analysis.xls')

def filter_artifacts_from_bed(infiles, bed_file, work_dir):
	os.chdir(work_dir)
	for infile in infiles:
		outfile = infile.rsplit('.', 1)[0] + '.artifact.xls'
		bed1 = infile.rsplit('.', 1)[0] + '.temp1.bed'
		bed2 = infile.rsplit('.', 1)[0] + '.temp2.bed'
		bed3 = infile.rsplit('.', 1)[0] + '.temp3.bed'
		##write bed to query
		with open(infile, "r") as in_fh, open(bed1, "w") as b1_fh:
			lc = 0
			for line in in_fh:
				line = line.split(delim)
				lc += 1
				if lc >1:
					b1_fh.write(delim.join(line[:3]) + '\n')
		##bedtools sort then intersect
		with open(bed2, "w") as b2_fh:
			bt_sort = subprocess.Popen(['bedtools','sort', '-i', bed1], stdout=b2_fh)
			bt_sort.wait()
		with open(bed3, "w") as b3_fh:
			bt_int = subprocess.Popen(['bedtools', 'intersect', '-a', bed2, '-b', bed_file, '-sorted', '-wa'], stdout=b3_fh)
			bt_int.wait()
		vars_in_intervals = {}
		with open(bed3, "r") as b3r_fh:
			for line in b3r_fh:
				line = line.rstrip().split(delim)
				var = '_'.join(line)
				vars_in_intervals[var] = line
		##open infile and write final file
		with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
			lc = 0
			for line in in_fh:
				line = line.rstrip().split(delim)
				lc += 1
				if lc  == 1:
					out_fh.write(delim.join(line + ['artifact']) + '\n')
				else:
					var = '_'.join(line[:3])
					if var in vars_in_intervals:
						out_fh.write(delim.join(line + ['yes']) + '\n')
					else:
						out_fh.write(delim.join(line + ['no']) + '\n')




##run methods

##run all at once
work_dir = '/home/atimms/ngs_data/genomes/daniela_uw_1119'
input_file = 'daniela_genomes_1119.txt'
# input_file = 'test_dg0819.txt'
uw_vcf = 'luquetti_external_uwcmg_cfm_1.HF.final.vcf.gz'
# run_analysis(work_dir, input_file, uw_vcf)

'''
##for gemini part, takes a while, so split peds
input_file = 'daniela_genomes_1119_1.txt'
# run_analysis(work_dir, input_file, uw_vcf)
work_dir = '/home/atimms/ngs_data/genomes/daniela_uw_1119/b2'
input_file = 'daniela_genomes_1119_2.txt'
# run_analysis(work_dir, input_file, uw_vcf)
work_dir = '/home/atimms/ngs_data/genomes/daniela_uw_1119/b3'
input_file = 'daniela_genomes_1119_3.txt'
# run_analysis(work_dir, input_file, uw_vcf)
work_dir = '/home/atimms/ngs_data/genomes/daniela_uw_1119/b4'
input_file = 'daniela_genomes_1119_4.txt'
# run_analysis(work_dir, input_file, uw_vcf)
'''

##combined files, three ways all trios (std), all trios (auto dom), and duos/singelton
peds_to_analyze = ['CFM-MOS-11', 'CFM-MOS-12', 'CL1003', 'CL1016', 'CL1020', 'CL3003', 'CL3023', 'CL3033', 'CL3044', 'CL3045', 'CL3051', 'CL4001', 'CL4004', 'CL4005', 'CL4007', 'F003', 'F005', 'F007', 'F008', 'F009', 'F011', 'F015', 'F016', 'F017', 'F0200003', 'F0200013', 'F0200022', 'F020', 'F022', 'F027', 'F028', 'F0400012', 'F0400020', 'F040', 'F041', 'F042', 'F045', 'F047', 'F050', 'F093', 'F149', 'F150', 'F152', 'F154', 'F168', 'F200006', 'F300031', 'F300059', 'F300073', 'F300097', 'F300111', 'F30028', 'F400017', 'F400025', 'PM1010002', 'PM1010004', 'PM1010005', 'PM1010007', 'PM1010010', 'PM1010017', 'PM1010036', 'PM1010037', 'PM1010039', 'PM1010040', 'PM1010042', 'PM1010043', 'PM1010044', 'PM1010045', 'PM1030019', 'PM1030021', 'PM1030022', 'PM1030025', 'PM1030028', 'PM1030029', 'PM1030030', 'PM1030031', 'PM1060018', 'PM1060019', 'PM1060028', 'PM1060029', 'PM1060030', 'PM1060031', 'PM1060032', 'PM1070005', 'PM1070017', 'PM1070019', 'PM1070020', 'PM4010001', 'PM4010003', 'PM4010004', 'PM4010005']
# combine_results_files(peds_to_analyze, work_dir)
##add column if in artifact bed
combined_trio_std = 'combined.trio.std_analysis.xls'
combined_trio_dom = 'combined.trio.autosomal_dominant.xls'
combined_duo_single = 'combined.singleton_duo.std_analysis.xls'
files_to_filter = [combined_trio_std, combined_trio_dom, combined_duo_single]
# files_to_filter = [combined_trio_std]
##bed from danila, had to refmat, add chr and remove a bunch of weird lines
artifact_bed = 'luquetti_external_cfm_1_filter.sorted_merge.bed'
##filter
filter_artifacts_from_bed(files_to_filter, artifact_bed, work_dir)


