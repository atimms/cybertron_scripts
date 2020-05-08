#!/usr/bin/python
import os
import subprocess
import glob

##modules 
'''
module load biobuilds/2017.11
'''


##set input variables and parameters
delim = '\t'


##programs
convert_2_annovar = '/home/atimms/programs/annovar_0618/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar_0618/table_annovar.pl'
# bcftools = '/home/atimms/programs/bcftools-1.9/bcftools'

##files
fasta = '/home/atimms/ngs_data/references/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta'
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

##get family ids
proband_list = ['CFM-MOS-11-01', 'CFM-MOS-12-01', 'CL100301', 'CL101601', 'CL102001', 'CL300301', 'CL302301', 'CL303301', 'CL304401', 'CL304501', 
		'CL305101', 'CL400101', 'CL400401', 'CL400501', 'CL400701', 'F03-00006', 'F03-00011', 'F03-00015', 'F03-00018', 'F03-00020', 'F03-00025', 
		'F03-00034', 'F03-00037', 'F03-00040', 'F03-00045', 'F02-00003', 'F02-00013', 'F02-00022', 'F03-00051', 'F03-00062', 'F03-00065', 'F3-0021-01', 
		'F04-00012', 'F04-00020', 'F03-00091', 'F03-00094', 'F03-00129', 'F03-00106', 'F03-00114', 'PE2D1012', 'F3-0002-01', 'F3-0003-01', 'F3-0006-01', 
		'F3-0010-01', 'F3-0017-01', 'F200006', 'F300031', 'F300059', 'F300073', 'F300097', 'F300111', 'PE2D1004', 'F400017', 'F400025', 'PM101000201', 
		'PM101000401', 'PM101000501', 'PM101000701', 'PM101001001', 'PM101001701', 'PM101003601', 'PM101003701', 'PM101003901', 'PM101004001', 
		'PM101004201', 'PM101004301', 'PM101004401', 'PM101004501', 'PM103001901', 'PM103002101', 'PM103002202', 'PM103002501', 'PM103002801', 'PM103002901', 
		'PM103003001', 'PM103003101', 'PM106001801', 'PM106001809', 'PM106001901', 'PM106002801', 'PM106002901', 'PM106002907', 'PM106003001', 'PM106003101', 
		'PM106003201', 'PM107000501', 'PM107001701', 'PM107001901', 'PM107002001', 'PM401000101', 'PM401000301', 'PM401000401', 'PM401000501']
ped_id_list = ['CFM-MOS-11', 'CFM-MOS-12', 'CL1003', 'CL1016', 'CL1020', 'CL3003', 'CL3023', 'CL3033', 'CL3044', 'CL3045', 'CL3051', 'CL4001', 'CL4004', 
		'CL4005', 'CL4007', 'F003', 'F005', 'F007', 'F008', 'F009', 'F011', 'F015', 'F016', 'F017', 'F020', 'F0200003', 'F0200013', 'F0200022', 'F022', 
		'F027', 'F028', 'F040', 'F0400012', 'F0400020', 'F041', 'F042', 'F045', 'F047', 'F050', 'F093', 'F149', 'F150', 'F152', 'F154', 'F168', 'F200006', 
		'F300031', 'F300059', 'F300073', 'F300097', 'F300111', 'F30028', 'F400017', 'F400025', 'PM1010002', 'PM1010004', 'PM1010005', 'PM1010007', 'PM1010010', 
		'PM1010017', 'PM1010036', 'PM1010037', 'PM1010039', 'PM1010040', 'PM1010042', 'PM1010043', 'PM1010044', 'PM1010045', 'PM1030019', 'PM1030021', 
		'PM1030022', 'PM1030025', 'PM1030028', 'PM1030029', 'PM1030030', 'PM1030031', 'PM1060018', 'PM1060018', 'PM1060019', 'PM1060028', 'PM1060029', 
		'PM1060029', 'PM1060030', 'PM1060031', 'PM1060032', 'PM1070005', 'PM1070017', 'PM1070019', 'PM1070020', 'PM4010001', 'PM4010003', 'PM4010004', 'PM4010005']
pro_ped_dict = dict(zip(proband_list, ped_id_list))


##annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,dbnsfp35a,rmsk,genomicSuperDups,exac03,gnomad_genome,gnomad_exome,kaviar_20150923,clinvar_20180603']
av_operation = ['-operation', 'g,f,r,r,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput', '-arg', '-splicing 10 ,,,,,,,,']




##methods
def get_vcf_per_ped(sample_file, in_vcf, ped_vcf):
	temp1_vcf = ped_vcf.split('.')[0] + '.temp1.vcf.gz'
	temp2_vcf = ped_vcf.split('.')[0] + '.temp2.vcf.gz'
	# temp3_vcf = ped_vcf.split('.')[0] + '.temp3.vcf.gz'
	##get the samples we want, and remove when we don't see a call
	bcftools_view = subprocess.Popen(['bcftools', 'view', '-a', '--threads', '10', '-Ou', '-S', sample_file, in_vcf], stdout=subprocess.PIPE)
	bcftools_view2 = subprocess.Popen(['bcftools', 'view', '-m', '2', '--threads', '10', '-f', 'PASS', '-O', 'z', '-o', temp1_vcf, '-'], stdin=bcftools_view.stdout)
	bcftools_view2.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', temp2_vcf, temp1_vcf])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fasta, '-O', 'z', '-o', ped_vcf, temp2_vcf])
	bcf_norm2.wait()
	# bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	# bcf_index.wait()

def make_sample_file(outfile, sample_names):
	with open(outfile, "w") as out_fh:
		for sample_name in sample_names:
			out_fh.write(sample_name + '\n')

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


def filter_by_popfreq_genotype_reformat(in_file, all_outfile, dn_outfile, aaf_req, aaf_cols, sample_names):
	with open(in_file, "r") as in_fh, open(all_outfile, "w") as aout_fh, open(dn_outfile, "w") as dout_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				extra_header = ['rs_id', 'vcf_format', 'proband_info', 'father_info', 'mother_info']
				aout_fh.write(delim.join(line[:133] + extra_header ) + '\n')
				dout_fh.write(delim.join(line[:133] + extra_header ) + '\n')
			else:
				##get metrics
				popfreqs = line[aaf_cols[0]:aaf_cols[1]]
				# print(popfreqs)
				##convert . to 0
				new_popfreqs = [x if x != '.' else '0' for x in popfreqs]
				# print(new_popfreqs)	
				##convert string to int and get max
				new_popfreqs = [float(i) for i in new_popfreqs]
				max_popfreq = max(new_popfreqs)
				# print(new_popfreqs, max_popfreq)
				rmsk = line[80]
				seqdups = line[81]
				genotypes = [i.split(':')[0] for i in line[146:]]
				# print(genotypes)
				func = line[5]
				ex_func = line[8]
				line_out = line[:133] + [line[139]] + line[145:] + [', '.join(sample_names)]
				if max_popfreq <= aaf_req and genotypes == ['0/1', '0/0', '0/0']:
					if rmsk == '.' and seqdups == '.':
						aout_fh.write(delim.join(line_out) + '\n')
					if (func == 'exonic' or func == 'splicing') and ex_func != 'synonymous SNV':
						dout_fh.write(delim.join(line_out) + '\n')

def annotate_vcf(in_vcf, out_prefix, sample_ids):
	multi_anno = out_prefix + '.hg38_multianno.txt'
	temp_ann = out_prefix + '.anno_temp.xls'
	all_denovo_file = out_prefix + '.filtered.all.xls'
	exonic_denovo_file = out_prefix + '.filtered.exonic.xls'
	# '''
	##run_table_annovar
	command = [table_annovar] + av_buildver + [in_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	run_annovar = subprocess.Popen(command)
	run_annovar.wait()
	##add other data i.e. based on genes
	add_gene_based_annotation(multi_anno, temp_ann, 6, 115)
	##filter and reformat ann files and get exonic files
	max_aaf = 0.01
	popfreq_cols = [90,106] ##all column between these numbers
	filter_by_popfreq_genotype_reformat(temp_ann, all_denovo_file, exonic_denovo_file, max_aaf, popfreq_cols, sample_ids)
	# '''
	return [exonic_denovo_file, all_denovo_file]


def combine_ann_files(infiles, out_file):
	with open(out_file, "w") as out_fh:
		print('making file %s by combining:'%out_file)
		file_count = 0
		for infile in infiles:
			print(infile)
			ped = infile.split('.')[0]
			# file_count += 1
			##add header from first file
			with open(infile, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line = line.split(delim)
					line_count +=1
					if line_count == 1:
						file_count += 1
						if file_count == 1:
							header = ['ped'] + line
							out_fh.write(delim.join(header))
					else:
						line_out = [ped] + line
						out_fh.write(delim.join(line_out))


#run all other methods
def master_results_files_from_vcf_fams(ped_files, work_dir, no_pcr_vcf, pcr_vcf, out_prefix):
	os.chdir(work_dir)
	exonic_dn_files = []
	all_dn_files = []
	for ped_file in ped_files:
		with open(ped_file, "r") as in_fh:
			for line in in_fh:
				line = line.split(delim)
				ped_id = pro_ped_dict[line[1]]
				samples = line[1:4]
				sample_file = ped_id + '.samples_temp.txt'
				pedigree_vcf = ped_id + '.vcf.gz'
				if ped_file.startswith('pcr_free'):
					combined_vcf = no_pcr_vcf
				else:
					combined_vcf = pcr_vcf
				print(ped_id, samples, combined_vcf)
				# '''
				##make file with list of samples
				make_sample_file(sample_file, samples)
				##get vcf with just a single ped
				get_vcf_per_ped(sample_file, combined_vcf, pedigree_vcf)
				# '''
				##annoate with annovar
				ann_files = annotate_vcf(pedigree_vcf, ped_id, samples)
				##get list of all files
				exonic_dn_files.append(ann_files[0])
				all_dn_files.append(ann_files[1])
	print(exonic_dn_files, all_dn_files)
	combine_ann_files(exonic_dn_files, out_prefix + '.filtered.exonic.xls')
	combine_ann_files(all_dn_files, out_prefix + '.filtered.all.xls')





##run methods
working_dir = '/home/atimms/ngs_data/genomes/daniela_broad_dn_vcfs_0819'
project_name = 'broad_dn_analysis_0919'
# fam_files = ['pcr_free.test.fam']
fam_files = ['pcr_free.pedigree.fam', 'pcr_plus.pedigree.fam']
pcr_free_vcf = 'GMKF_Gabriel_Luquetti_Craniofacial_WGS_v2.HIGH_confidence_denovos.vcf.gz'
pcr_plus_vcf = 'GMKF_Gabriel_Luquetti_Craniofacial_WGS_PCR_Plus_v2.HIGH_confidence_denovos.vcf.gz'

master_results_files_from_vcf_fams(fam_files, working_dir, pcr_free_vcf, pcr_plus_vcf, project_name)


