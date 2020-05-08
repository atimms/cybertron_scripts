#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

##set input variables and parameters
delim = '\t'

##programs and ref files
vt = 'vt'
# gemini = '/home/atimms/programs/gemini/data/anaconda/bin/gemini'
gemini = '/home/atimms/scripts/gemini'
snpeff_jar = '/home/atimms/programs/snpEff/snpEff.jar'
bcftools_12 = 'bcftools'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'

ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
dbdb_data = ref_dir + 'gemini/dbdb.gene_association.txt'
mgi_data = ref_dir + 'gemini/mgi.abnormal_brain.txt'
omim_data = ref_dir + 'gemini/omim_genemap2_0816.txt'
previous_exomes_data = ref_dir + 'gemini/all_exome_data_std_pipeline.0918.xls'
gnomad_exome_data = '/home/atimms/ngs_data/references/gnomad_data/gnomad.exomes.r2.1.1.avinput'
gnomad_genome_data = '/home/atimms/ngs_data/references/gnomad_data/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.avinput'
# gnomad_genome_data = '/home/atimms/ngs_data/references/gnomad_data/gnomad.genomes.r2.1.avinput'
##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_operation = ['-operation', 'g,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-arg', '-splicing 10 ,,,']


def load_single_vcf_into_gemini(ped, db_name, input_vcf, ped_file):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	temp_vcf = 'temp1.vcf'
	normalized_vcf = ped + '.int.norm.vcf'
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
	gemini_load.wait()
	##remove intermediate files
	os.remove(temp_vcf)

def cp_gemini_db(in_db, out_db):
	cp_command = subprocess.Popen(['cp', in_db, out_db])
	cp_command.wait()


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
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.auto_dom.xls'
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




##make dict using dbdb file
def make_dict_from_dbdb(dbdb_file):
	dbdb_dict = {}
	dbdb_dict2 = {}
	with open(dbdb_file, "U") as dbdb_fh:
		for line in dbdb_fh:
			line = line.strip('\n').split(delim)
			# print line
			gene = line[1]
			inheritance = line[2]
			pheno = line[3]
			syndrome = line[5]
			loe = line[7]
			if gene not in dbdb_dict:
				dbdb_dict[gene] = [[inheritance],[pheno],[syndrome],[loe]]
				# print gene, inheritance, pheno, syndrome, loe
			else:
				dbdb_dict[gene][0].append(inheritance)
				dbdb_dict[gene][1].append(pheno)
				dbdb_dict[gene][2].append(syndrome)
				dbdb_dict[gene][3].append(loe)
	for g in dbdb_dict:
		dbdb_dict2[g] = []
		for info in dbdb_dict[g]:
			info_comma = ','.join(info)
			dbdb_dict2[g].append(info_comma)
	# for g in dbdb_dict:
	# 	print g, dbdb_dict[g], dbdb_dict2[g]
	return dbdb_dict2

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

##make dict using omim file
def make_dict_from_omim(omim_file):
	omim_dict, omim_dict2 = {}, {}
	with open(omim_file, "U") as omim_fh:
		for line in omim_fh:
			if line[0] != '#':
				line = line.strip('\n').split(delim)
				# print line
				gene = line[8]
				pheno = line[12]
				if pheno == '':
					pheno = 'None'
				if gene != '':
					# print gene, pheno
					if gene not in omim_dict:
						omim_dict[gene] = [pheno]
						# print gene, pheno
					else:
						omim_dict[gene].append(pheno)
						# print 'gene %s occurs multiple times in omim file'%gene
	for g in omim_dict:
		omim_dict2[g] = [', '.join(omim_dict[g])]
	return omim_dict2

def make_dict_from_previous_exomes(exome_file):
	ex_dict, ex_dict2 = {}, {}
	with open(exome_file, "U") as ex_fh:
		line_count = 0
		for line in ex_fh:
			line = line.strip('\n').split(delim)
			line_count += 1
			if line_count >1:
				gene = line[2]
				ped = line[0]
				# print gene, ped
				if gene in ex_dict:
					ex_dict[gene].append(ped)
				else:
					ex_dict[gene] = [ped]
	for g in ex_dict:
		peds = list(set(ex_dict[g]))
		ex_dict2[g] = [', '.join(peds)]
		# print g, ex_dict[g], peds, ', '.join(peds)
	return ex_dict2

def add_columns_to_file(in_file, out_file, dbdb_dict, mgi_dict, omim_dict, exomes_data_dict, gene_column, pos_to_insert):
	genelist = []
	with open(in_file, "U") as inf, open(out_file, "w") as outf:
		line_count, genes_not_found = 0, 0
		for line in inf:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				extra_header = ['dbdb inheritance', 'dbdb phenotype', 'dbdb syndrome', 'dbdb loe', 'mgi phenotype', 'omim phenotype', 'previous exomes']
				outf.write(delim.join(line[:pos_to_insert] + extra_header + line[pos_to_insert:] + ['\n']))
			else:
				gene = line[gene_column]
				genelist.append(gene)
				##get values for all gene names and print minimum value
				extra_stuff = []
				if gene in dbdb_dict:
					extra_stuff.extend(dbdb_dict[gene])
				else:
					extra_stuff.extend(['na', 'na', 'na', 'na'])
				if gene in mgi_dict:
					extra_stuff.extend(mgi_dict[gene])
				else:
					extra_stuff.extend(['na'])
				if gene in omim_dict:
					extra_stuff.extend(omim_dict[gene])
				else:
					extra_stuff.extend(['na'])
				if gene in exomes_data_dict:
					extra_stuff.extend(exomes_data_dict[gene])
				else:
					extra_stuff.extend(['na'])
				# print gene, extra_stuff
				outf.write(delim.join(line[:pos_to_insert] + extra_stuff + line[pos_to_insert:] + ['\n']))
	return genelist

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
			extra_ann = line[10:17]
			# print chrom, pos, ts_info
			ann_dict[ch_pos] = [ts_info] + extra_ann

	with open(in_file, 'r') as in_fh, open(out_file, 'w') as out_fh:
		line_count = 0
		for line2 in in_fh:
			line2 = line2.split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join(line2[:pos_to_insert] + ['annovar_ts_info','CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG', 'cosmic90_coding', 'cosmic90_noncoding'] + line2[pos_to_insert:]))
			else:
				chrom2 = line2[0][3:]
				pos2 = line2[2]
				ch_pos2 = '_'.join([chrom2,pos2])
				if ch_pos2 in ann_dict:
					# print line2[:5],ch_pos2, ann_dict[ch_pos2]
					out_fh.write(delim.join(line2[:pos_to_insert] + ann_dict[ch_pos2] + line2[pos_to_insert:]))
				else:
					print 'not found:', line2[:5], ch_pos2
	##remove intermediate files
	# '''
	for files_to_go in [temp_vcf, av_file, multianno,region_file]:
		os.remove(files_to_go)
	# '''

def add_aaf_info(pedigree, in_file, out_file, pos_to_insert):
	##files
	normalized_vcf = pedigree + '.int.norm.vcf.gz'
	region_file = pedigree + '_temp.regions'
	temp_vcf = pedigree + 'temp0.vcf'
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
	bcftools_filter = subprocess.Popen([bcftools_12, 'view', '-R', region_file, '-o', temp_vcf, '-O', 'v', normalized_vcf])
	bcftools_filter.wait()
	##make dict with vcf data
	vcf_dict = {}
	with open(temp_vcf, 'r') as vcf_fh:
		for line in vcf_fh:
			if line[:2] != '##':
				if line[0] == '#':
					line = line.rstrip().split(delim)
					samples = line[9:]
					sample_head = [i + '_ref/alt' for i in samples]
				else:
					line = line.rstrip().split(delim)
					chrom = line[0]
					pos = line[1]
					ch_pos = '_'.join([chrom,pos])
					ref = line[3]
					alt = line[4]
					g_info = line[9:]
					aaf_info = [i.split(':')[1] for i in g_info]
					# print chrom, pos, g_info, aaf_info
					vcf_dict[ch_pos] = aaf_info

	with open(in_file, 'r') as in_fh, open(out_file, 'w') as out_fh:
		line_count = 0
		for line2 in in_fh:
			line2 = line2.split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join(line2[:pos_to_insert] + sample_head + line2[pos_to_insert:]))
			else:
				chrom2 = line2[0][3:]
				pos2 = str(int(line2[1]) + 1)
				ch_pos2 = '_'.join([chrom2,pos2])
				if ch_pos2 in vcf_dict:
					# print line2[:5],ch_pos2, vcf_dict[ch_pos2]
					out_fh.write(delim.join(line2[:pos_to_insert] + vcf_dict[ch_pos2] + line2[pos_to_insert:]))
				else:
					print 'not found:', line2[:5], ch_pos2
					nas = ['na'] * len(sample_head)
					out_fh.write(delim.join(line2[:pos_to_insert] + nas + line2[pos_to_insert:]))
	for files_to_go in [temp_vcf,region_file]:
		os.remove(files_to_go)

def make_gnomad_dicts():
	exome_dict = {}
	genome_dict = {}
	with open(gnomad_exome_data, "r") as exome_fh:
		for line in exome_fh:
			line = line.strip('\n').split(delim)
			# print line
			var = '_'.join(line[:5])
			info = line[5:]
			exome_dict[var] = info
	with open(gnomad_genome_data, "r") as genome_fh:
		for line in genome_fh:
			line = line.strip('\n').split(delim)
			# print line
			var = '_'.join(line[:5])
			info = line[5:]
			genome_dict[var] = info
	return exome_dict, genome_dict

def add_gnomad_info(in_file, out_file, pos_to_insert, exome_dict, genome_dict):
	with open(in_file, "r") as inf, open(out_file, "w") as outf:
		line_count = 0
		for iline in inf:
			iline = iline.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				extra_header = ['gnomMD_exome_AC', 'gnomMD_exome_AD', 'gnomMD_exome_AF', 'gnomMD_exome_nhom','gnomMD_exome_AC_non_neuro', 'gnomMD_exome_AD_non_neuro', 
				'gnomMD_exome_AF_non_neuro', 'gnomMD_exome_nhom_non_neuro', 'gnomMD_genome_AC', 'gnomMD_genome_AD', 'gnomMD_genome_AF', 'gnomMD_genome_nhom',
				'gnomMD_genome_AC_non_neuro', 'gnomMD_genome_AD_non_neuro', 'gnomMD_genome_AF_non_neuro', 'gnomMD_genome_nhom_non_neuro']
				outf.write(delim.join(iline[:pos_to_insert] + extra_header + iline[pos_to_insert:] + ['\n']))
			else:
				ivar = '_'.join([iline[0][3:], str(int(iline[1]) + 1)] + iline[2:5])
				extra_stuff = []
				if ivar in exome_dict:
					extra_stuff.extend(exome_dict[ivar])
				else:
					extra_stuff.extend(['na', 'na', 'na', 'na', 'na', 'na', 'na', 'na'])
				if ivar in genome_dict:
					extra_stuff.extend(genome_dict[ivar])
				else:
					extra_stuff.extend(['na', 'na', 'na', 'na', 'na', 'na', 'na', 'na'])				
				outf.write(delim.join(iline[:pos_to_insert] + extra_stuff + iline[pos_to_insert:] + ['\n']))




def combine_gemini_results(ped, file_list, gene_col, position_to_insert, position_to_finish, add_annovar_pos, add_aaf_pos,gnomad_ex_dict, gnomad_gen_dict):
	##file names
	outfile = ped + '.std_analysis.xls'
	temp_file = 'combined.temp.xls'
	temp2_file = 'combined.temp2.xls'
	temp3_file = 'combined.temp3.xls'
	temp4_file = 'combined.temp4.xls'
	print 'combining files:', file_list
	##combine all files, making sure col lengths are all the same
	with open(temp_file, "w") as temp_fh:
		file_count = 0
		for filename in file_list:
			# print filename
			analysis_type = filename.split('.')[-2]
			# print analysis_type, file_count
			##add header from first file
			with open(filename, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line_count +=1
					line = line.rstrip().split(delim)
					if line_count == 1:
						file_count += 1
						if file_count == 1:
							header = delim.join(line[:position_to_insert] + line[position_to_insert:position_to_finish] + ['analysis', '\n'])
							temp_fh.write(header)
					else:
						line_out = delim.join(line[:position_to_insert] + line[position_to_insert:position_to_finish] + [analysis_type, '\n'])
						temp_fh.write(line_out)

	##add dbdb and mgi data
	dbdb_data_dict = make_dict_from_dbdb(dbdb_data)
	mgi_data_dict = make_dict_from_mgi(mgi_data)
	omim_data_dict = make_dict_from_omim(omim_data)
	previous_exomes_data_dict = make_dict_from_previous_exomes(previous_exomes_data)
	genenames = add_columns_to_file(temp_file, temp2_file, dbdb_data_dict, mgi_data_dict, omim_data_dict, previous_exomes_data_dict, gene_col, position_to_insert)
	add_annovar_ts_info(ped, temp2_file, temp3_file, add_annovar_pos)
	add_aaf_info(ped, temp3_file, temp4_file, add_aaf_pos)
	add_gnomad_info(temp4_file, outfile, add_aaf_pos, gnomad_ex_dict, gnomad_gen_dict)
	##remove intermediate files
	os.remove(temp_file)
	os.remove(temp2_file)
	os.remove(temp3_file)
	os.remove(temp4_file)
	return genenames

def get_previous_vars(genes, outfile):
	print genes, outfile
	with open(outfile, "w") as out_fh, open(previous_exomes_data, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line = line.split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join(line))
			else:
				gene = line[2]
				if gene in genes:
					out_fh.write(delim.join(line))


def format_autosomal_dominant_analysis(ped, file_list, gene_col, position_to_insert, position_to_finish, add_annovar_pos, add_aaf_pos, gnomad_ex_dict, gnomad_gen_dict):
	##file names
	outfile = ped + '.auto_dominant_analysis.xls'
	temp_file = 'combined.temp.xls'
	temp2_file = 'combined.temp2.xls'
	temp3_file = 'combined.temp3.xls'
	temp4_file = 'combined.temp4.xls'
	print 'combining files:', file_list
	# '''
	##combine all files, making sure col lengths are all the same
	with open(temp_file, "w") as temp_fh:
		file_count = 0
		for filename in file_list:
			# print filename
			# analysis_type = filename.split('.')[-2]
			# file_count += 1
			# print analysis_type, file_count
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

	##add dbdb and mgi data
	dbdb_data_dict = make_dict_from_dbdb(dbdb_data)
	mgi_data_dict = make_dict_from_mgi(mgi_data)
	omim_data_dict = make_dict_from_omim(omim_data)
	previous_exomes_data_dict = make_dict_from_previous_exomes(previous_exomes_data)
	genenames = add_columns_to_file(temp_file, temp2_file, dbdb_data_dict, mgi_data_dict, omim_data_dict, previous_exomes_data_dict, gene_col, position_to_insert)
	add_annovar_ts_info(ped, temp2_file, temp3_file, add_annovar_pos)
	add_aaf_info(ped, temp3_file, temp4_file, add_aaf_pos)
	# '''
	add_gnomad_info(temp4_file, outfile, add_aaf_pos, gnomad_ex_dict, gnomad_gen_dict)
	##remove intermediate files
	os.remove(temp_file)
	os.remove(temp2_file)
	os.remove(temp3_file)
	os.remove(temp4_file)


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



# def gemini_comp_het(prefix, gem_db_name, test_name):
# 	outfile = prefix + '.' + test_name + '.xls'
# 	with open (outfile, 'w') as out_fh:
# 		##compound het: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
# 		gemini_query = subprocess.Popen([gemini, 'comp_hets', '--max-priority', '2', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter',  "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', gem_db_name], stdout=out_fh)
# 		gemini_query.wait()

# def gemini_recessive(prefix, gem_db_name, test_name):
# 	outfile = prefix + '.' + test_name + '.xls'
# 	with open (outfile, 'w') as out_fh:
# 		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
# 		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', gem_db_name], stdout=out_fh)
# 		gemini_query.wait()

# def gemini_dom(prefix, gem_db_name, test_name):
# 	outfile = prefix + '.' + test_name + '.xls'
# 	with open (outfile, 'w') as out_fh:
# 		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
# 		gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', gem_db_name], stdout=out_fh)
# 		gemini_query.wait()


##run all other methods
def standard_gemini_protocol(working_directory, pedigree, ped_type):
	os.chdir(working_directory)
	gemini_db = pedigree + '.gemini.db'
	ad_db = pedigree + '.gemini_ad.db'
	in_vcf = pedigree + '.intersected_vcfs/0002.vcf.gz'
	std_ped_file = pedigree + '.ped'
	ad_ped_file = pedigree + '_ad.ped'
	freq_req = 0.01
	freq_req_recessive = 0.05
	coverage_req = 10
	##load vcf file into db
	# '''
	load_single_vcf_into_gemini(pedigree, gemini_db, in_vcf, std_ped_file)
	
	##make auto dom ped and gemini db for duo and trio analysis
	if ped_type != 'singleton':
		make_ad_files(std_ped_file, ad_ped_file, gemini_db, ad_db)
	##get dicts from gnomad data
	ex_dict, genome_dict = make_gnomad_dicts()
	# '''
	##run queries and get files to combine
	files_to_combine = []
	# '''
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
		format_autosomal_dominant_analysis(pedigree, [auto_dom_file], 5, 13, 22, 13, 31, ex_dict, genome_dict)
	else:
		print 'ped type %s not recognized'%ped_type
	# '''
	##combine results
	genes_in_file = combine_gemini_results(pedigree, files_to_combine, 5, 13, 22, 13, 31, ex_dict, genome_dict)
	##get all variants seen 
	get_previous_vars(genes_in_file, pedigree + '.previous_vars_in_same_genes.xls')

	# '''


##tests
##parameters and to run
# ##add ped file which is just ped name plus .ped
# working_dir = '/data/atimms/LR12-121'
# ##xlined ped list [[probands],[parents]], if proband is female, or don't have both parents just add 'na'
# x_lists = [['LR15-121'], ['LR15-121f', 'LR15-121m']]
# standard_gemini_protocol(working_dir, 'LR15-121', x_lists)
# working_dir = '/home/atimms/ngs_data/exomes/working/ghayda_0618'
# standard_gemini_protocol(working_dir, "LR18-339", 'singleton')
# working_dir = '/home/atimms/ngs_data/exomes/working/temp'
# standard_gemini_protocol(working_dir, "LR14-053", 'trio')
















