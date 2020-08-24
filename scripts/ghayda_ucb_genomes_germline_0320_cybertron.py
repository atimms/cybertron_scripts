#!/usr/bin/python
import os
import subprocess
import glob


'''
after starting interactive session run these:
module load java/1.8.0_121
module load local_python/2.7.14
source activate hg38_genomes


##installing again 0320
install mini/conda3 and loading local python 2.7.14:
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
##environment to use, create and activate
conda create --name hg38_genomes
conda activate hg38_genomes #or
source activate hg38_genomes
##install vcfanno #not used
conda install -c bioconda vcfanno
##install vcf2db
in programs folder...
git clone https://github.com/quinlan-lab/vcf2db
##instally dependancies, beforehand.. make file in 
##/programs/anaconda2/envs/hg38_genomes/conda-meta called 
##pinned and add python 2.7.* (keep python version 2.7)
cd vcf2db
conda install -y snappy
conda install -y -c conda-forge python-snappy
conda install -y -c bioconda cyvcf2 peddy
pip install -r requirements.txt
##issue with this so instead of pip requirements (i think due to )
conda install -c anaconda numpy
conda install -c anaconda sqlalchemy
conda install -y peddy ##conda didn't have version for py3 so had to install via github
conda install -c bioconda geneimpacts
##other tools for analysis in geneomes
conda install -c bioconda bcftools 
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

##annovar parameters --- add cosmic
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp35a,dbnsfp31a_interpro,avsnp150,gnomad211_genome,gnomad211_exome,gnomad30_genome,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.']

##methods
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
				id_name = line[6]
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
				if ped_name in analysis_dict:
					if id_name in analysis_dict[ped_name][0]:
						print('issue with sample', id_name, sample_name)
					else:
						analysis_dict[ped_name][0][id_name] = sample_name
				else:
					analysis_dict[ped_name] = [{id_name:sample_name}, anal_type]
	return ped_file_dict, analysis_dict


def make_ped_files(input_dict):
	for ped in input_dict:
		# print ped, input_dict[ped]
		outfile = ped + '.ped'
		header = ['#family_id', 'name', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
		with open(outfile, "w") as out_fh:
			out_fh.write(delim.join(header) + '\n')
			for outline in input_dict[ped]:
				out_fh.write(delim.join(outline) + '\n')

def make_sample_file(samples, outfile):
	with open(outfile, "w") as out_fh:
		for sample in samples:
			out_fh.write(sample + '\n')

def make_sample_rehead_file(sample_dict, outfile):
	with open(outfile, "w") as out_fh:
		for sample in sample_dict:
			new_id = sample_dict[sample]
			out_fh.write(sample + ' ' + new_id + '\n')

def get_vcf_per_ped(pedigree_dict, in_vcf):
	for pedigree in pedigree_dict:
		sample_list = pedigree_dict[pedigree][0].keys()
		sample_dict = pedigree_dict[pedigree][0]
		sample_file = pedigree + '.samples_temp.txt'
		rehead_sample_file = pedigree + '.samples_rehead_temp.txt'
		ped_temp_vcf = pedigree + '.temp.vcf.gz'
		ped_vcf = pedigree + '.vcf.gz'
		##make file to get samples from combined file
		make_sample_file(sample_list, sample_file)
		##get the samples we want, and remove when we don't see a call and only keep passed
		##seperate version - not used
		# bcftools_view = subprocess.Popen([bcftools, 'view', '-a', '--threads', '20', '-O', 'z', '-o', 'temp.vcf.gz', '-S', sample_file, in_vcf])
		# bcftools_view.wait()
		# bcftools_view2 = subprocess.Popen([bcftools, 'view', '-m', '2', '--threads', '20', '-O', 'z', '-o', ped_vcf, 'temp.vcf.gz'])
		# bcftools_view2.wait()
		# '''
		##combined version
		bcftools_view = subprocess.Popen(['bcftools', 'view', '-a', '--threads', '10', '-Ou', '-S', sample_file, in_vcf], stdout=subprocess.PIPE)
		bcftools_view2 = subprocess.Popen(['bcftools', 'view', '-m', '2', '-f', 'PASS', '--threads', '10', '-O', 'z', '-o', ped_temp_vcf, '-'], stdin=bcftools_view.stdout)
		bcftools_view2.wait()
		# '''
		##reheader - get new sample 
		make_sample_rehead_file(sample_dict, rehead_sample_file)
		bcftools_rehead = subprocess.Popen(['bcftools', 'reheader', '-s', rehead_sample_file, '-o', ped_vcf, ped_temp_vcf])
		bcftools_rehead.wait()


def load_vcf_into_gemini(ped, db_name, input_vcf, ped_file):
	##dcompress, change a gatk header thing, and decompse and normalize
	temp_vcf = ped + 'temp1.vcf'
	normalized_vcf = ped + '.int.norm.vcf'
	# vcf_anno_vcf = ped + '.vcfanno.vcf'
	# '''
	##normalize vcf with vt
	with open (temp_vcf, 'w') as tvcf_fh:
		zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
		sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
		vt_decompose = subprocess.Popen(['vt', 'decompose', '-s', '-'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
		vt_normalize = subprocess.Popen(['vt', 'normalize', '-r', fasta, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
		vt_normalize.wait()
	# '''
	##normalize vcf with vt etc  -with additional sed (not used)
	# with open (temp_vcf, 'w') as tvcf_fh:
	# 	zless_vcf = subprocess.Popen(['zless', input_vcf], stdout=subprocess.PIPE)
	# 	##recomendation from gemini
	# 	sed_vcf = subprocess.Popen(['sed', 's/ID=AD,Number=./ID=AD,Number=R/'], stdin=zless_vcf.stdout, stdout=subprocess.PIPE)
	# 	##new sed sed -e 's/Number=A/Number=1/g' from https://github.com/bcbio/bcbio-nextgen/issues/1775
	# 	sed_vcf2 = subprocess.Popen(['sed', '-e', 's/Number=A/Number=1/g'], stdin=sed_vcf.stdout, stdout=subprocess.PIPE)
	# 	vt_decompose = subprocess.Popen(['vt', 'decompose', '-s', '-'], stdin=sed_vcf2.stdout, stdout=subprocess.PIPE)
	# 	vt_normalize = subprocess.Popen(['vt', 'normalize', '-r', fasta, '-'], stdin=vt_decompose.stdout, stdout=tvcf_fh)
	# 	vt_normalize.wait()
	# '''
	##annotate with snpeff and compress and index
	with open (normalized_vcf, 'w') as nvcf_fh:
		snpeff_vcf = subprocess.Popen(['java', '-Xmx20G', '-jar', snpeff_jar, 'GRCh38.86', '-v', '-formatEff', '-classic', temp_vcf], stdout=nvcf_fh)
		snpeff_vcf.wait()
		bgzip_vcf = subprocess.Popen(['bgzip', normalized_vcf])
		bgzip_vcf.wait()
		tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', normalized_vcf + '.gz'])
		tabix_vcf.wait()
	# '''
	# ##annotate with vcfanno --- not used
	# with open(vcf_anno_vcf, 'w') as va_vcf_fh:
	# 	run_vcfanno = subprocess.Popen(['vcfanno', '-p', '10', '-base-path', ref_dir, conf_file, normalized_vcf + '.gz'], stdout=va_vcf_fh)
	# 	run_vcfanno.wait()
	# 	bgzip_vcf = subprocess.Popen(['bgzip', vcf_anno_vcf])
	# 	bgzip_vcf.wait()
	# 	tabix_vcf = subprocess.Popen(['tabix', '-p', 'vcf', vcf_anno_vcf + '.gz'])
	# 	tabix_vcf.wait()
	##load into db
	# gemini_load = subprocess.Popen(['python', vcf2db, vcf_anno_vcf + '.gz', ped_file, db_name])
	gemini_load = subprocess.Popen(['python', vcf2db, normalized_vcf + '.gz', ped_file, db_name])
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
			# extra_info = [line[13], line[16], line[51], line[65]] + line[80:115]
			extra_info = line[10:12] + [line[15], line[18], line[53], line[67]] + line[83:138]
			# print chrom, pos, ts_info
			ann_dict[ch_pos] = [ts_info] + extra_info

	with open(in_file, 'r') as in_fh, open(out_file, 'w') as out_fh:
		line_count = 0
		for line2 in in_fh:
			line2 = line2.split(delim)
			line_count += 1
			if line_count == 1:
				out_fh.write(delim.join(line2[:pos_to_insert] + ['annovar_ts_info', 'rmsk', 'genomicSuperDups', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 
					'CADD_phred', 'GERP++_RS', 'avsnp150', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 
					'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 
					'non_cancer_AF_popmax', 'controls_AF_popmax', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 
					'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 
					'controls_AF_popmax', 'AF', 'AF_raw', 'AF_male', 'AF_female', 'AF_afr', 'AF_ami', 'AF_amr', 'AF_asj', 'AF_eas', 'AF_fin', 'AF_nfe', 
					'AF_oth', 'AF_sas', 'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG', 'cosmic90_coding', 'cosmic90_noncoding'] + line2[pos_to_insert:]))
			else:
				chrom2 = line2[0]
				pos2 = line2[2]
				ch_pos2 = '_'.join([chrom2,pos2])
				if ch_pos2 in ann_dict:
					# print line2[:5],ch_pos2, ann_dict[ch_pos2]
					out_fh.write(delim.join(line2[:pos_to_insert] + ann_dict[ch_pos2] + line2[pos_to_insert:]))
				else:
					print 'not found:', line2[:5], ch_pos2

def filter_by_popfreq(in_file, out_file, aaf_req, aaf_cols):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.split(delim)
			if line_count == 1:
				out_fh.write(delim.join(line[:74]) + '\n')
			else:
				popfreqs = [ line[i] for i in aaf_cols]
				# popfreqs = line[aaf_cols[0]:aaf_cols[1]]
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

def gemini_denovo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo.xls'
	temp1_file = ped + '.de_novo_temp1.xls'
	temp2_file = ped + '.de_novo_temp2.xls'
	temp3_file = ped + '.de_novo_temp3.xls'
	with open(temp1_file, 'w') as tout_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact', '--filter', "filter IS NULL AND impact_severity != 'LOW'", '-d', str(coverage), gem_db_name], stdout=tout_fh)
		gemini_query.wait()
	add_annovar_info(ped, temp1_file, temp2_file, 10)
	'''
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [18,35,52]
	filter_by_popfreq(temp3_file, outfile, max_aaf, popfreq_cols)
	'''
	##without gene based
	popfreq_cols = [18,35,52]
	filter_by_popfreq(temp2_file, outfile, max_aaf, popfreq_cols)
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
	add_annovar_info(ped, temp1_file, temp2_file, 10)
	'''
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [18,35,52]
	filter_by_popfreq_comphet(temp3_file, outfile, max_aaf, popfreq_cols,5)
	'''
	##without gene based
	popfreq_cols = [18,35,52]
	filter_by_popfreq_comphet(temp2_file, outfile, max_aaf, popfreq_cols,5)
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
	add_annovar_info(ped, temp1_file, temp2_file, 10)
	'''
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [13,38]
	filter_by_popfreq(temp3_file, outfile, max_aaf, popfreq_cols)
	'''
	##without gene based
	popfreq_cols = [18,35,52]
	filter_by_popfreq(temp2_file, outfile, max_aaf, popfreq_cols)
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
	add_annovar_info(ped, temp1_file, temp2_file, 10)
	'''
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [13,38]
	filter_by_popfreq(temp3_file, outfile, max_aaf, popfreq_cols)
	'''
	##without gene based
	popfreq_cols = [18,35,52]
	filter_by_popfreq(temp2_file, outfile, max_aaf, popfreq_cols)
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
	add_annovar_info(ped, temp1_file, temp2_file, 10)
	'''
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [13,38]
	filter_by_popfreq(temp3_file, outfile, max_aaf, popfreq_cols)
	'''
	##without gene based
	popfreq_cols = [18,35,52]
	filter_by_popfreq(temp2_file, outfile, max_aaf, popfreq_cols)
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
	add_annovar_info(ped, temp1_file, temp2_file, 10)
	'''
	##add gene based, include gene column and col to insert data
	add_gene_based_annotation(temp2_file, temp3_file, 5, 50)
	popfreq_cols = [13,38]
	filter_by_popfreq(temp3_file, outfile, max_aaf, popfreq_cols)
	'''
	##without gene based
	popfreq_cols = [18,35,52]
	filter_by_popfreq(temp2_file, outfile, max_aaf, popfreq_cols)
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


#run all gemini methods
def gemini_protocol(pedigree_dict):
	for pedigree in pedigree_dict:
		gemini_db = pedigree + '.gemini.db'
		ad_db = pedigree + '.gemini_ad.db'
		in_vcf = pedigree + '.vcf.gz'
		std_ped_file = pedigree + '.ped'
		ad_ped_file = pedigree + '_ad.ped'
		freq_req = 0.01
		freq_req_recessive = 0.05
		coverage_req = 10
		std_out_file = pedigree + '.std_analysis.xls'
		ped_type = pedigree_dict[pedigree][1]
		
		##load vcf file into db
		'''
		load_vcf_into_gemini(pedigree, gemini_db, in_vcf, std_ped_file)
		'''
		##make auto dom ped and gemini db for duo and trio analysis
		make_ad_files(pedigree, std_ped_file, ad_ped_file, gemini_db, ad_db)
		##run queries and get files to combine
		files_to_combine = []
		if ped_type == 'singleton':
			print('pedigree %s is a singleton'%pedigree)
			# files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
			# files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
			# files_to_combine.append(gemini_dominant(pedigree, gemini_db, freq_req, coverage_req))
		elif ped_type == 'duo':
			print('pedigree %s is a duo'%pedigree)
			# files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
			# files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
			# files_to_combine.append(gemini_dominant(pedigree, ad_db, freq_req, coverage_req))
		elif ped_type == 'trio' or ped_type == 'quad':
			print('pedigree %s is a trio or quad'%pedigree)
			# files_to_combine.append(gemini_comp_het(pedigree, gemini_db, freq_req_recessive, coverage_req))
			# files_to_combine.append(gemini_recessive(pedigree, gemini_db, freq_req_recessive, coverage_req))
			# files_to_combine.append(gemini_denovo(pedigree, gemini_db, freq_req, coverage_req))
			# files_to_combine.append(gemini_xlinked(pedigree, gemini_db, freq_req_recessive, coverage_req))
			# files_to_combine.append(gemini_xlinked_de_novo(pedigree, gemini_db, freq_req, coverage_req))
			##not used right now
			gemini_dominant(pedigree, ad_db, freq_req, coverage_req)
		else:
			print('ped type %s not recognized'%ped_type)
		##combine results
		# combine_gemini_results(files_to_combine, std_out_file)
		


##master method
def run_analysis(working_dir, analysis_file, combined_vcf):
	os.chdir(working_dir)
	project = analysis_file.split('.')[0]
	##read in analysis file and make 2 dicts
	ped_file_dict, analysis_dict = make_analysis_dicts(analysis_file)
	##check...
	for p in analysis_dict:
		print(p, analysis_dict[p])
	##make ped files from first dict
	make_ped_files(ped_file_dict)
	##split vcf and change ids
	# get_vcf_per_ped(analysis_dict,combined_vcf)
	##gemini_protocol
	gemini_protocol(analysis_dict)

##run methods
work_dir = '/home/atimms/ngs_data/genomes/ghayda_ucb_0320'
# work_dir = '/home/atimms/ngs_data/genomes/test'
broad_vcf = 'GhaydaMirzaa_SeattleChildrens_WGS_Joint_VCF_all_samples.vcf.gz'
# sample_file = 'ghayda_ucb_0320_1.txt'
sample_file = 'ghayda_ucb_0320_all.txt'
run_analysis(work_dir, sample_file, broad_vcf)












