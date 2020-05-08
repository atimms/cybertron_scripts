#!/usr/bin/python
import os
import subprocess
import glob
import shutil

##set input variables and parameters
delim = '\t'

'''
load these:
module load java/1.8.0_121 
module load biobuilds/2016.11
module load vt/0.5772-60f436c3 
'''



##programs and ref files
vt = 'vt'
gemini = '/home/atimms/scripts/gemini'
snpeff_jar = '/home/atimms/programs/snpEff/snpEff.jar'
bcftools_12 = 'bcftools'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'

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
domino_data = gemini_ref_dir + 'score_all_final_19.02.19.txt'

##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene']
av_operation = ['-operation', 'g']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-arg', '-splicing 10 ']


##methods
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
	gemini_load.wait()
	##remove intermediate files
	os.remove(temp_vcf)

def gemini_exonic_dn(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo_exonic.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL and is_coding == 1 and max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
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

def format_dn_file(ped, infile, gene_col, position_to_insert, position_to_finish, add_annovar_pos, add_aaf_pos):
	##file names
	outfile = ped + '.exonic_denovo.xls'
	temp_file = 'combined.temp.xls'
	temp2_file = 'combined.temp2.xls'
	temp3_file = 'combined.temp3.xls'
	##combine all files, making sure col lengths are all the same
	with open(temp_file, "w") as temp_fh:
		print infile
		##add header from first file
		with open(infile[0], "r") as in_fh:
			line_count = 0
			for line in in_fh:
				line_count +=1
				line = line.rstrip().split(delim)
				if line_count == 1:
					header = delim.join(line[:position_to_insert] + line[position_to_insert:position_to_finish] + ['\n'])
					temp_fh.write(header)
				else:
					line_out = delim.join(line[:position_to_insert] + line[position_to_insert:position_to_finish] + ['\n'])
					temp_fh.write(line_out)
	if line_count >=2:
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

def combine_dn_files(samples, out_file):
	with open(out_file, "w") as out_fh:
		fc = 0
		for sample in samples:
			in_file = sample + '.exonic_denovo.xls'
			if os.path.exists(in_file):
				fc += 1
				lc = 0
				with open(in_file, "r") as in_fh:
					for line in in_fh:
						lc += 1
						# line = line.rstrip().split(delim)
						if lc == 1:
							if fc == 1:
								out_fh.write(line)
						else:
							out_fh.write(line)
			else:
				print "File {} was not found".format(in_file)



##run all other methods for getting exonic de novos
def get_exonic_denovos_from_trios(working_directory, trio_pedigrees, all_outfile):
	os.chdir(working_directory)
	for pedigree in trio_pedigrees:
		gemini_db = pedigree + '.gemini.db'
		in_vcf = pedigree + '.intersected_vcfs/0002.vcf.gz'
		std_ped_file = pedigree + '.ped'
		freq_req = 0.01
		coverage_req = 10
		##load vcf file into db, and get exonic denovos
		# load_single_vcf_into_gemini(pedigree, gemini_db, in_vcf, std_ped_file)
		# exonic_dn_file = gemini_exonic_dn(pedigree, gemini_db, freq_req, coverage_req)
		# format_dn_file(pedigree, [exonic_dn_file], 5, 13, 22, 13, 31)
	combine_dn_files(trio_pedigrees, all_outfile)

def make_dict_domino(domino_file):
	dom_dict = {}
	with open(domino_file, "U") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			if line_count > 1:
				line = line.strip('\n').split(delim)

				gene = line[0]
				score = line[1]
				# print gene, syndrome
				if gene not in dom_dict:
					dom_dict[gene] = score
					# print gene, gdi, gdi_phred
				else:
					print 'gene %s seen multiple times'% gene
	return dom_dict


def filter_by_repeats_annovar_add_domino(working_directory, infile, outfile):
	os.chdir(working_directory)
	##make dict from domino file
	domino_dict = make_dict_domino(domino_data)
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)

			if line_count == 1:
				header = line[:13] + ['domnio score'] + line[13:]
				out_fh.write(delim.join(header) + '\n')
			else:
				rmsk = line[33]
				seqdups = line[34]
				annovar_ts = line[6]
				if annovar_ts == '.' or annovar_ts.split(':')[0].startswith('NM_'):
					gene = line[5]
				else:
					gene = annovar_ts.split(':')[0]
				print line[5], line[6], gene
				if gene in domino_dict:
					dom_score = domino_dict[gene]
				else:
					dom_score = 'na'
				if rmsk == 'None' and seqdups == '0' and annovar_ts != '.':
					if not annovar_ts.startswith('dist='):
						line_out = line[:13] + [dom_score] + line[13:]
						out_fh.write(delim.join(line_out) + '\n')



##run methods
work_dir = '/home/atimms/ngs_data/exomes/working/daniela_microtia_exomes_0919'

##get exonic de novos
trio_peds = ['1010002', '1010004', '1010005', '1010007', '1030001', '1010008', '1010010', '1010013', '1010019', '1010020', '1010021', '1010022', '1010023', '1010024', '1010025', '1010026', 
		'1010028', '1030002', '1030006', '1060004', '1070001', '1070002', '1070003', '1070011', '1010031', '1010032', '1030013', '1030014', '1030015', '1060020', '1070012', '1070013', '1070015', 
		'1030011', '1060023', '4010002']
# trio_peds = ['1010002'] #test
combined_exonc_dn = 'microtia_exomes.exonic_denovo.0919.xls'
##run
get_exonic_denovos_from_trios(work_dir, trio_peds, combined_exonc_dn)


##filter inherited variants from previous analysis
auto_dom_file = 'combined.auto_dom.txt'
filtered_inh_file = 'combined.auto_dom.filtered_0919.xls'
# filter_by_repeats_annovar_add_domino(work_dir, auto_dom_file, filtered_inh_file)



