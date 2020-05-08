#!/usr/bin/python
import os
import subprocess
import glob


'''

after starting interactive session run these:
module load java/1.8.0_121
module load local_python/2.7.15
source activate hg38_genomes


##installing
works on cybertron, and on ewrlnxrg26
install mini/conda2 and loading local python 2.7.14/python 2.7.15:
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
##environment to use, create and activate
conda create --name hg38_genomes
conda activate hg38_genomes
source activate hg38_genomes
##install vcfanno
conda install -c bioconda vcfanno
##install vcf2db
in programs folder...
git clone https://github.com/quinlan-lab/vcf2db
cd vcf2db
conda install -y snappy
conda install -y -c conda-forge python-snappy
conda install -y -c bioconda cyvcf2 peddy
pip install -r requirements.txt
##other tools for analysis in geneomes
conda install -c bioconda bcftools (library issue, install indepentatly)
conda install -c bioconda tabix
conda install -c bioconda vt

'''

##set input variables and parameters
delim = '\t'
##programs
snpeff_jar = '/home/atimms/programs/snpEff_201711/snpEff.jar'
vcf2db = '/home/atimms/programs/vcf2db/vcf2db.py'
bcftools = '/home/atimms/programs/bcftools-1.9/bcftools'
convert_2_annovar = '/home/atimms/programs/annovar_0618/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar_0618/table_annovar.pl'
gemini = 'gemini'
plink = '/home/atimms/programs/plink'
# fasta = '/home/atimms/ngs_data/references/hg38/hg38.fa'
fasta = '/home/atimms/ngs_data/references/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta'
##vcfanno_files
ref_dir = '/home/atimms/ngs_data/references/hg38_gemini/'
conf_file = ref_dir + 'hg38_genomes_0219.conf'
lua_file = ref_dir + 'hg38_genomes_0219.lua'
##ref data for gene based annotation
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
rho_genes = gemini_ref_dir + 'RhoGTPAse_genes.txt'

##annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,dbnsfp35a,rmsk,genomicSuperDups,exac03,gnomad_genome,gnomad_exome,kaviar_20150923,clinvar_20180603']
av_operation = ['-operation', 'g,f,r,r,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-arg', '-splicing 10 ,,,,,,,,']




##methods
def merge_vcf_files(in_vcfs, out_vcf):
	bcftools_merge = subprocess.Popen([bcftools, 'merge', '--threads', '20', '-O', 'z', '-o', out_vcf, '-f', 'PASS'] + in_vcfs)
	bcftools_merge.wait()

def get_sample_file_from_ped_file(infile, outfile):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.split(delim)
			sample = line[1]
			if line_count >1:
				out_fh.write(sample + '\n')



def vt_tidy_vcf(in_dir, in_vcfs):
	os.chdir(in_dir)
	for in_vcf in in_vcfs:
		input_vcf = in_vcf
		out_vcf = in_vcf.rsplit('.', 2)[0] + '.tidy.vcf'
		# tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', input_vcf])
		# tabix_vcf.wait()
		with open (out_vcf, 'w') as ovcf_fh:
			zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
			sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
			vt_decompose = subprocess.Popen(['vt', 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
			vt_normalize = subprocess.Popen(['vt', 'normalize', '-r', fasta, '-n', '-'], stdin=vt_decompose.stdout, stdout=ovcf_fh)
			vt_normalize.wait()
			bgzip_vcf = subprocess.Popen(['bgzip', out_vcf])
			bgzip_vcf.wait()
			tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', out_vcf + '.gz'])
			tabix_vcf.wait()

def load_vcf_into_gemini(ped, db_name, input_vcf, ped_file):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	temp_vcf = 'temp1.vcf'
	normalized_vcf = ped + '.int.norm.vcf'
	vcf_anno_vcf = ped + '.vcfanno.vcf'
	# '''
	##normalize vcf with vt
	with open (temp_vcf, 'w') as tvcf_fh:
		zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		vt_decompose = subprocess.Popen(['vt', 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
		vt_normalize = subprocess.Popen(['vt', 'normalize', '-r', fasta, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
		vt_normalize.wait()
	##annotate with snpeff and compress and index
	with open (normalized_vcf, 'w') as nvcf_fh:
		snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'GRCh38.86', '-v', '-formatEff', '-classic', temp_vcf], stdout=nvcf_fh)
		snpeff_vcf.wait()
		bgzip_vcf = subprocess.Popen(['bgzip', normalized_vcf])
		bgzip_vcf.wait()
		tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', normalized_vcf + '.gz'])
		tabix_vcf.wait()
	# '''
	
	'''
	##annotate with vcfanno --- not used
	with open(vcf_anno_vcf, 'w') as va_vcf_fh:
		run_vcfanno = subprocess.Popen(['vcfanno', '-p', '10', '-base-path', ref_dir, conf_file, normalized_vcf + '.gz'], stdout=va_vcf_fh)
		run_vcfanno.wait()
		bgzip_vcf = subprocess.Popen(['bgzip', vcf_anno_vcf])
		bgzip_vcf.wait()
		tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', vcf_anno_vcf + '.gz'])
		tabix_vcf.wait()
	'''
	##load into db
	# gemini_load = subprocess.Popen(['python', vcf2db, vcf_anno_vcf + '.gz', ped_file, db_name])
	gemini_load = subprocess.Popen(['python', vcf2db, normalized_vcf + '.gz', '-e', 'PG', ped_file, db_name])
	gemini_load.wait()

def cp_gemini_db(in_db, out_db):
	cp_command = subprocess.Popen(['cp', in_db, out_db])
	cp_command.wait()

def add_new_ped_file_to_gemini_db(pedfile, geminidb):
	gemini_amend = subprocess.Popen(['gemini', 'amend', '--sample', pedfile, geminidb])
	gemini_amend.wait()


def make_ad_files(ped, in_ped, out_ped, in_db, out_db):
	normalized_vcf = ped + '.int.norm.vcf'
	##make ad ped, so convert 1's to 0's in ped
	with open(out_ped, "w") as out_fh, open(in_ped, "r") as in_fh:
		for line in in_fh:
			line = line.strip('\n').split(delim)
			pheno = line[5]
			the_rest = line[:5]
			if pheno == '1':
				pheno = '0'
			out_fh.write(delim.join(the_rest + [pheno]) + '\n')
	##cp gemini database then add ad ped file -- doesn't work with vcf2db dbds
	# cp_gemini_db(in_db, out_db)
	# add_new_ped_file_to_gemini_db(out_ped, out_db)
	##load into db
	gemini_load = subprocess.Popen(['python', vcf2db, normalized_vcf + '.gz', out_ped, out_db])
	gemini_load.wait()

def add_annovar_info(pedigree, in_file, out_file, pos_to_insert):
	##files
	normalized_vcf = pedigree + '.int.norm.vcf.gz'
	file_prefix = in_file.rsplit('.', 1)[0]
	region_file = file_prefix + '.regions'
	temp_vcf = file_prefix + '.vcf.gz'
	av_file = file_prefix + '.avinput'
	multianno = pedigree + '.hg38_multianno.txt'
	## make_list_of_regions(var_file, region_file):
	with open(in_file, 'r') as var_fh, open(region_file, 'w') as reg_fh:
		line_count = 0
		for line in var_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count >1:
				chrom = line[0]
				# start = int(line[1]) - 2
				# end = int(line[2]) + 2
				start = line[1]
				end = line[2]
				reg_fh.write(delim.join([chrom, str(start), str(end), '\n']))
				# print chrom, start_end
	## get_var_from_regions_in_vcf(in_vcf, regions, out_vcf):
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-R', region_file, '-o', temp_vcf, '-O', 'z', normalized_vcf])
	bcftools_filter.wait()
	##convert vcf file to individual avinput file
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', temp_vcf, '-allsample', '--withfreq', '-outfile', av_file])
	con_ann.wait()
	##run_table_annovar(avinput, av_prefix):
	command = [table_annovar] + av_buildver + [av_file] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', pedigree]
	annovar = subprocess.Popen(command)
	annovar.wait()
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
			extra_info = [line[13], line[16], line[51], line[65]] + line[80:115]
			# print chrom, pos, ts_info
			ann_dict[ch_pos] = [ts_info] + extra_info

	with open(in_file, 'r') as in_fh, open(out_file, 'w') as out_fh:
		line_count = 0
		for line2 in in_fh:
			line2 = line2.split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join(line2[:pos_to_insert] + ['annovar_ts_info', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 'CADD_phred', 'GERP++_RS', 'rmsk', 
					'genomicSuperDups', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 
					'ExAC_SAS', 'gnomAD_genome_ALL', 'gnomAD_genome_AFR', 'gnomAD_genome_AMR', 'gnomAD_genome_ASJ', 
					'gnomAD_genome_EAS', 'gnomAD_genome_FIN', 'gnomAD_genome_NFE', 'gnomAD_genome_OTH', 'gnomAD_exome_ALL', 
					'gnomAD_exome_AFR', 'gnomAD_exome_AMR', 'gnomAD_exome_ASJ', 'gnomAD_exome_EAS', 'gnomAD_exome_FIN', 
					'gnomAD_exome_NFE', 'gnomAD_exome_OTH', 'gnomAD_exome_SAS', 'Kaviar_AF', 'Kaviar_AC', 'Kaviar_AN', 
					'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG'] + line2[pos_to_insert:]))
			else:
				chrom2 = line2[0]
				pos2 = line2[2]
				ch_pos2 = '_'.join([chrom2,pos2])
				if ch_pos2 in ann_dict:
					# print line2[:5],ch_pos2, ann_dict[ch_pos2]
					out_fh.write(delim.join(line2[:pos_to_insert] + ann_dict[ch_pos2] + line2[pos_to_insert:]))
				else:
					print 'not found:', line2[:5], ch_pos2

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

def make_dict_from_rho_genes(rho_file):
	rho_genes_dict = {}
	with open(rho_file, "U") as rg_fh:
		for line in rg_fh:
			gene = line.strip('\n')
			if gene not in rho_genes_dict:
				rho_genes_dict[gene] = ['RhoGTPAse gene']
			else:
				print 'gene %s seen multiple times'% gene
				# pass
	return rho_genes_dict


def add_columns_to_file(in_file, out_file, dicts_to_add, gene_col, position_to_insert):
	with open(in_file, "U") as inf, open(out_file, "w") as outf:
		line_count, genes_not_found = 0, 0
		for line in inf:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				extra_header = ['RVIS exac', 'Percentile exac', 'GDI', 'GDI phred', 'mgi phenotype', 'hpo_syndrome', 'sbse_wt_fpkm', 'sbse_log2_fc', 'sbse_q_value', 
						'hmx1_wt_fpkm', 'hmx1_log2_fc', 'hmx1_q_value', 'ear_57d_fpkm', 'ear_59d_fpkm', 'human_syndrome', 'hoxa2_wt1', 'hoxa2_wt3', 'ba_expression_papers_0816', 'RhoGTPAse_gene']
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
						extra_stuff.extend(['.']*len(data_dict[data_dict.keys()[0]]))

				# print gene, extra_stuff
				outf.write(delim.join(line[:position_to_insert] + extra_stuff + line[position_to_insert:] + ['\n']))


def add_gene_based_annotation(infile, outfile, gene_col, position_to_insert):
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
	rho_data_dict = make_dict_from_rho_genes(rho_genes)
	add_columns_to_file(infile, outfile, [rvis_data_dict, gdi_data_dict, mgi_data_dict, hpo_data_dict, sbse_rnaseq_dict, hmx1_rnaseq_dict, human_ear_rnaseq_dict, human_syndrome_dict, hoxa2_dict, ba_paper_dict, rho_data_dict], gene_col, position_to_insert)

def filter_by_popfreq(in_file, out_file, aaf_req, aaf_cols):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.split(delim)
			if line_count == 1:
				out_fh.write(delim.join(line[:74]) + '\n')
			else:
				popfreqs = line[aaf_cols[0]:aaf_cols[1]]
				# print(popfreqs)
				##convert . to 0
				new_popfreqs = [x if x != '.' else '0' for x in popfreqs]
				# print(new_popfreqs)	
				##convert string to int and get max
				new_popfreqs = [float(i) for i in new_popfreqs]
				max_popfreq = max(new_popfreqs)
				# print(new_popfreqs, max_popfreq)
				if max_popfreq <= aaf_req:
					out_fh.write(delim.join(line[:74]) + '\n')

def filter_by_popfreq_comphet(in_file, out_file, aaf_req, aaf_cols, genename_col):
	##go through vars and 
	temp_file = in_file.rsplit('.', 1)[0] + '.temp.xls'
	passed_gene_list, passed_var_list = [], []
	with open(in_file, "r") as in_fh, open(temp_file, "w") as temp_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				temp_fh.write(delim.join(line[:74]) + '\n')
			else:
				popfreqs = line[aaf_cols[0]:aaf_cols[1]]
				gene = line[genename_col]
				var = '_'.join(line[:5])
				# print(popfreqs)
				##convert . to 0
				new_popfreqs = [x if x != '.' else '0' for x in popfreqs]
				# print(new_popfreqs)	
				##convert string to int and get max
				new_popfreqs = [float(i) for i in new_popfreqs]
				max_popfreq = max(new_popfreqs)
				# print(new_popfreqs, max_popfreq)
				if max_popfreq <= aaf_req:
					if var not in passed_var_list:
						passed_var_list.append(var)
						passed_gene_list.append(gene)
						temp_fh.write(delim.join(line[:74]) + '\n')
	with open(temp_file, "r") as in_fh, open(out_file, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join(line[:74]) + '\n')
			else:
				gene = line[genename_col]
				passed_vars_in_gene = passed_gene_list.count(gene)
				if passed_vars_in_gene > 1:
					out_fh.write(delim.join(line[:74]) + '\n')


def gemini_denovo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo.xls'
	temp1_file = ped + '.de_novo_temp1.xls'
	temp2_file = ped + '.de_novo_temp2.xls'
	temp3_file = ped + '.de_novo_temp3.xls'
	with open(temp1_file, 'w') as tout_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact', '--filter', "filter IS NULL AND impact_severity != 'LOW'", '-d', str(coverage), gem_db_name], stdout=tout_fh)
		gemini_query.wait()
	add_annovar_info(ped, temp1_file, temp2_file, 6)
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [13,38]
	filter_by_popfreq(temp3_file, outfile, max_aaf, popfreq_cols)
	return outfile

def gemini_comp_het(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.comp_hets.xls'
	temp1_file = ped + '.comp_hets_temp1.xls'
	temp2_file = ped + '.comp_hets_temp2.xls'
	temp3_file = ped + '.comp_hets_temp3.xls'
	with open(temp1_file, 'w') as tout_fh:
		##compound het: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'comp_hets', '--max-priority', '2', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact', '--filter',  "filter IS NULL AND impact_severity != 'LOW' " , '-d', str(coverage), gem_db_name], stdout=tout_fh)
		gemini_query.wait()
	add_annovar_info(ped, temp1_file, temp2_file, 6)
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [13,38]
	filter_by_popfreq_comphet(temp3_file, outfile, max_aaf, popfreq_cols,5)
	return outfile

def gemini_recessive(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.autosomal_recessive.xls'
	temp1_file = ped + 'autosomal_recessive_temp1.xls'
	temp2_file = ped + '.autosomal_recessive_temp2.xls'
	temp3_file = ped + '.autosomal_recessive_temp3.xls'
	with open(temp1_file, 'w') as tout_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact', '--filter', "filter IS NULL AND impact_severity != 'LOW'", '-d', str(coverage), gem_db_name], stdout=tout_fh)
# 		gemini_query = subprocess.Popen([gemini, 'autosomal_recessive', '--columns', 'chrom,start,end,ref,alt,gene,biotype,impact,aa_change,gerp_bp_score,cadd_scaled,max_aaf_all', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= 0.01", '-d', '10', gem_db_name], stdout=out_fh)
		gemini_query.wait()
	add_annovar_info(ped, temp1_file, temp2_file, 6)
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [13,38]
	filter_by_popfreq(temp3_file, outfile, max_aaf, popfreq_cols)
	return outfile

def gemini_xlinked(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.x_linked_recessive.xls'
	temp1_file = ped + '.x_linked_recessive_temp1.xls'
	temp2_file = ped + '.x_linked_recessive_temp2.xls'
	temp3_file = ped + '.x_linked_recessive_temp3.xls'
	with open(temp1_file, 'w') as tout_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'x_linked_recessive', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact', '--filter', "filter IS NULL AND impact_severity != 'LOW' ", '-d', str(coverage), gem_db_name], stdout=tout_fh)
		gemini_query.wait()
	add_annovar_info(ped, temp1_file, temp2_file, 6)
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [13,38]
	filter_by_popfreq(temp3_file, outfile, max_aaf, popfreq_cols)
	return outfile

def gemini_xlinked_de_novo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.x_linked_de_novo.xls'
	temp1_file = ped + '.x_linked_de_novo_temp1.xls'
	temp2_file = ped + '.x_linked_de_novo_temp2.xls'
	temp3_file = ped + '.x_linked_de_novo_temp3.xls'
	with open(temp1_file, 'w') as tout_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'x_linked_de_novo', '--columns','chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact', '--filter', "filter IS NULL AND impact_severity != 'LOW' " , '-d', str(coverage), gem_db_name], stdout=tout_fh)
		gemini_query.wait()
	add_annovar_info(ped, temp1_file, temp2_file, 6)
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [13,38]
	filter_by_popfreq(temp3_file, outfile, max_aaf, popfreq_cols)
	return outfile

def gemini_dominant(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.autosomal_dominant.xls'
	temp1_file = ped + '.autosomal_dominant_temp1.xls'
	temp2_file = ped + '.autosomal_dominant_temp2.xls'
	temp3_file = ped + '.autosomal_dominant_temp3.xls'
	with open(temp1_file, 'w') as tout_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'autosomal_dominant', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact', '--filter', "filter IS NULL AND impact_severity != 'LOW' ", '-d', str(coverage), gem_db_name], stdout=tout_fh)
		gemini_query.wait()
	add_annovar_info(ped, temp1_file, temp2_file, 6)
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [13,38]
	filter_by_popfreq(temp3_file, outfile, max_aaf, popfreq_cols)
	return outfile

def combine_gemini_results(infiles, outfile):
	print 'combining files:', infiles
	##combine all files, making sure col lengths are all the same
	with open(outfile, "w") as out_fh:
		file_count, total_line_count = 0, 0
		for infile in infiles:
			analysis_type = infile.split('.')[-2]
			# file_count += 1
			print analysis_type, file_count
			##add header from first file
			with open(infile, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line_count +=1
					total_line_count += 1
					line = line.rstrip().split(delim)
					if line_count == 1:
						file_count += 1
						if file_count == 1:
							header = delim.join(line + ['analysis', '\n'])
							out_fh.write(header)
					else:
						line_out = delim.join(line + [analysis_type, '\n'])
						out_fh.write(line_out)

def change_sample_names_in_vcf(vcf_in, vcf_out, sample_file):
	temp_vcf = vcf_in.rsplit('.', 1)[0] + '.temp.vcf'
	vcf_out_uncompressed = vcf_out.rsplit('.', 1)[0]
	print(vcf_in, vcf_out, sample_file, temp_vcf, vcf_out_uncompressed)
	bcftools_rehead = subprocess.Popen([bcftools, 'reheader', '--samples', sample_file, '-o', temp_vcf + '.gz', vcf_in])
	bcftools_rehead.wait()
	bgzip_d_vcf = subprocess.Popen(['bgzip', '-d', temp_vcf + '.gz'])
	bgzip_d_vcf.wait()
	with open(vcf_out_uncompressed, "w") as out_fh, open(temp_vcf, "r") as in_fh:
		for line in in_fh:
			if line[0] == '#':
				out_fh.write(line)
			else:
				line = line.split(delim)
				info = line[7]
				if info.rsplit(';', 1)[1].startswith('ANN'):
					# print('info', info)
					info_no_ann = info.rsplit(';', 1)[0]
					line_out = line[:7] + [info_no_ann] + line[8:]
					out_fh.write(delim.join(line_out))
					# print('1', line)
					# print('2', line_out)
				else:
					print('issue with info in this line', line)

	bgzip_vcf = subprocess.Popen(['bgzip', vcf_out_uncompressed])
	bgzip_vcf.wait()

#run all other methods
def gemini_protocol(pedigrees, work_dir):
	os.chdir(work_dir)
	for pedigree in pedigrees:
		gemini_db = pedigree + '.gemini.db'
		ad_db = pedigree + '.gemini_ad.db'
		in_vcf = pedigree + '.vcf.gz'
		sample_names_file = pedigree + '.sample_names.txt'
		renamed_vcf = pedigree + '.renamed.vcf.gz'
		std_ped_file = pedigree + '.ped'
		ad_ped_file = pedigree + '_ad.ped'
		freq_req = 0.01
		freq_req_recessive = 0.05
		coverage_req = 10
		std_out_file = pedigree + '.std_analysis.xls'
		
		##change sample name
		change_sample_names_in_vcf(in_vcf, renamed_vcf, sample_names_file)
		##load vcf file into db
		load_vcf_into_gemini(pedigree, gemini_db, renamed_vcf, std_ped_file)
		##make auto dom ped and gemini db for duo and trio analysis
		# '''
		##need to remake db
		if pedigree not in single_list:
			make_ad_files(pedigree, std_ped_file, ad_ped_file, gemini_db, ad_db)
		##run queries and get files to combine
		files_to_combine = []
		if pedigree in single_list:
			print('pedigree %s is a singleton'%pedigree)
			files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
			files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
			files_to_combine.append(gemini_dominant(pedigree, gemini_db, freq_req, coverage_req))
		elif pedigree in duo_list:
			print('pedigree %s is a duo'%pedigree)
			files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
			files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
			files_to_combine.append(gemini_dominant(pedigree, ad_db, freq_req, coverage_req))
		else:
			print('pedigree %s is a trio'%pedigree)
			files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
			files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
			files_to_combine.append(gemini_denovo(pedigree, gemini_db, freq_req, coverage_req))
			files_to_combine.append(gemini_xlinked(pedigree, gemini_db, freq_req_recessive, coverage_req))
			files_to_combine.append(gemini_xlinked_de_novo(pedigree, gemini_db, freq_req, coverage_req))
			gemini_dominant(pedigree, ad_db, freq_req, coverage_req)
		##combine results
		combine_gemini_results(files_to_combine, std_out_file)
		# '''


def plink_relatadness_check(vcf, file_prefix, work_dir):
	os.chdir(work_dir)
	##correct filtering??
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(INFO/DP)>30", '-V', 'indels', '-o', 'temp_plink.vcf.gz', '-O', 'z', vcf])
	bcftools_filter.wait()
	##generate plink file from vcf
	make_plink = subprocess.Popen([plink, '--vcf', 'temp_plink.vcf.gz', '--vcf-filter', '--const-fid', '--out', 'temp.pass_q50_dp50'])
	make_plink.wait()
	##check sex -results in .sexcheck file
	plink_sex = subprocess.Popen([plink, '--bfile', 'temp.pass_q50_dp50', '--check-sex', '--out', file_prefix + '.pass_q50_dp50'])
	plink_sex.wait()
	##ibd check
	plink_ibd = subprocess.Popen([plink,  '--bfile', 'temp.pass_q50_dp50', '--genome', '--out', file_prefix + '.pass_q50_dp50'])
	plink_ibd.wait()

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





##run methods

duo_list = ['CL3003', 'CL3023', 'CL3051', 'F040', 'PM1010037']
single_list = ['CFM-MOS-11', 'F168']
# peds_to_analyze = ['CFM-MOS-12', 'PM1060018', 'CL3003', 'CL3023', 'CL3051', 'F040', 'CFM-MOS-11', 'F168', 'CL1003', 'CL1016', 
# 		'CL1020', 'CL3033', 'CL3044', 'CL3045', 'CL4001', 'CL4004', 'CL4005', 'CL4007', 'F003', 'F005', 'F007', 'F008', 'F009', 
# 		'F011', 'F015', 'F016', 'F017', 'F020', 'F0200003', 'F0200013', 'F0200022', 'F022', 'F027', 'F028', 'F0400012', 'F0400020', 
# 		'F041', 'F042', 'F045', 'F047', 'F050', 'F093', 'F149', 'F150', 'F152', 'F154', 'F200006', 'F300031', 'F300059', 'F300073', 
# 		'F300097', 'F300111', 'F30028', 'F400017', 'F400025', 'PM1010002', 'PM1010004', 'PM1010005', 'PM1010007', 'PM1010010', 'PM1010017', 
# 		'PM1010036', 'PM1010037', 'PM1010039', 'PM1010040', 'PM1010042', 'PM1010043', 'PM1010044', 'PM1010045', 'PM1030019', 'PM1030021', 
# 		'PM1030022', 'PM1030025', 'PM1030028', 'PM1030029', 'PM1030030', 'PM1030031', 'PM1060019', 'PM1060028', 'PM1060029', 'PM1060030', 
# 		'PM1060031', 'PM1060032', 'PM1070005', 'PM1070017', 'PM1070019', 'PM1070020', 'PM4010001', 'PM4010003', 'PM4010004', 'PM4010005']

peds_to_analyze = ['F003', 'F005']
# peds_to_analyze = ['F003']




##get samples from vcf, annotate, load into gemnini and then query
##need sample name in sepreate file (in same order as vcf), and ped file
working_dir = '/home/atimms/ngs_data/genomes/daniela_0519'
gemini_protocol(peds_to_analyze, working_dir)


