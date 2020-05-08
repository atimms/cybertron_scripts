#!/usr/bin/env python
import sys
import subprocess
import os
import glob



'''
##load modules required for analysis
module load java/1.8.0_121 
module load biobuilds/2016.11
module load vt/0.5772-60f436c3 
'''

##parameters
delim = '\t'
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
##programs
snpeff_jar = '/home/atimms/programs/snpEff/snpEff.jar'
bcftools = '/home/atimms/programs/bcftools-1.9/bcftools'
vt = 'vt'
gemini = '/home/atimms/scripts/gemini'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'
table_annovar = '/home/atimms/programs/annovar/table_annovar.pl'
plink = '/home/atimms/programs/plink'
##annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,clinvar_20170905,cosmic83_coding,cosmic83_noncoding']
av_operation = ['-operation', 'g,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-arg', '-splicing 10 ,,,']

##methods
def make_ped_files(input_file):
	ped_file_dict ={}
	with open(input_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				ped_name = line[0]
				print ped_name, line
				##add data to ped file dict
				ped_file_info = line
				# print ped_file_info
				if ped_name in ped_file_dict:
					ped_file_dict[ped_name].append(ped_file_info)
				else:
					ped_file_dict[ped_name] = [ped_file_info]
	for ped in ped_file_dict:
		# print ped, ped_file_dict[ped]
		outfile = ped + '.ped'
		header = ['#family_id', 'name', 'paternal_id', 'maternal_id', 'sex', 'phenotype', '\n']
		with open(outfile, "w") as out_fh:
			out_fh.write(delim.join(header))
			for outline in ped_file_dict[ped]:
				out_fh.write(delim.join(outline) + '\n')
	return ped_file_dict.keys()

def get_sample_file_from_ped_file(infile, outfile):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.split(delim)
			sample = line[1]
			if line_count >1:
				out_fh.write(sample + '\n')

def get_vcf_per_ped(pedigrees, in_vcf):
	for pedigree in pedigrees:
		ped_file = pedigree + '.ped'
		sample_file = pedigree + '.samples_temp.txt'
		ped_vcf = pedigree + '.vcf.gz'
		get_sample_file_from_ped_file(ped_file, sample_file)
		##get the samples we want, and remove when we don't see a call
		'''
		bcftools_view = subprocess.Popen([bcftools, 'view', '-a', '--threads', '20', '-O', 'z', '-o', 'temp.vcf.gz', '-S', sample_file, in_vcf])
		bcftools_view.wait()
		bcftools_view2 = subprocess.Popen([bcftools, 'view', '-m', '2', '--threads', '20', '-O', 'z', '-o', ped_vcf, 'temp.vcf.gz'])
		bcftools_view2.wait()
		'''
		bcftools_view = subprocess.Popen([bcftools, 'view', '-a', '--threads', '5', '-Ou', '-S', sample_file, in_vcf], stdout=subprocess.PIPE)
		bcftools_view2 = subprocess.Popen([bcftools, 'view', '-m', '2', '--threads', '5', '-O', 'z', '-o', ped_vcf, '-'], stdin=bcftools_view.stdout)
		bcftools_view2.wait()
		# '''

def load_single_vcf_into_gemini(ped, db_name, input_vcf, ped_file):
	##dcompress, change a gatk header thing, and decompse and normalize (should be redundant)
	temp_vcf = ped + '.temp_vt.vcf'
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

def gemini_denovo(ped, gem_db_name, max_aaf, coverage):
	outfile = ped + '.cov' + str(coverage) + '.maf'+ str(max_aaf) + '.de_novo.xls'
	with open (outfile, 'w') as out_fh:
		##auto recessive: passed by gatk, impact med or high, aaf <=1% and coverage >=10 in all family members
		gemini_query = subprocess.Popen([gemini, 'de_novo', '--lenient', '--columns', 'chrom,start,end,ref,alt,gene,transcript,biotype,aa_change,impact,gerp_bp_score,cadd_scaled,max_aaf_all,clinvar_gene_phenotype,rmsk,in_segdup', '--filter', "filter IS NULL AND impact_severity != 'LOW' AND max_aaf_all <= " + str(max_aaf), '-d', str(coverage), gem_db_name], stdout=out_fh)
		gemini_query.wait()
	return outfile

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
	print(line_count)
	##if have any variants
	if line_count >= 1:
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
					out_fh.write(delim.join(line2[:pos_to_insert] + ['annovar_ts_info','CLINSIG', 'CLNDBN', 'CLNACC', 'CLNDSDB', 'CLNDSDBID', 'cosmic83_coding', 'cosmic83_noncoding'] + line2[pos_to_insert:]))
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

def plink_relatadness_check(vcf, file_prefix):
	##correct filtering??
	# bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(INFO/DP)>30", '-V', 'indels', '-o', 'temp_plink.vcf.gz', '-O', 'z', vcf])
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-i', "QUAL>50", '-V', 'indels', '-o', 'temp_plink.vcf.gz', '-O', 'z', vcf])

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

##master method
def run_exome_analysis(working_dir, analysis_file, all_sample_vcf):
	os.chdir(working_dir)
	# project = analysis_file.split('.')[0]
	freq_req = 0.01
	coverage_req = 10
	add_annovar_pos = 6
	##read in analysis file and make 2 pedigree files and get name of all pedigrees
	peds = make_ped_files(analysis_file)
	print(peds)
	##split combined vcf
	get_vcf_per_ped(peds, all_sample_vcf)
	##load into gemini and 
	for pedigree in peds:
		gemini_db = pedigree + '.gemini.db'
		in_vcf = pedigree + '.vcf.gz'
		std_ped_file = pedigree + '.ped'
		gem_denovo = pedigree + '.cov' + str(coverage_req) + '.maf'+ str(freq_req) + '.de_novo.xls'
		format_denovo =  pedigree + '.cov' + str(coverage_req) + '.maf'+ str(freq_req) + '.de_novo.formatted.xls'
		# '''
		load_single_vcf_into_gemini(pedigree, gemini_db, in_vcf, std_ped_file)
		gemini_denovo(pedigree, gemini_db, freq_req, coverage_req)
		add_annovar_ts_info(pedigree, gem_denovo, format_denovo, add_annovar_pos)
		# '''
		plink_relatadness_check(in_vcf, pedigree)


def combine_results(infiles, outfile):
	print 'combining files:', infiles
	##combine all files, making sure col lengths are all the same
	with open(outfile, "w") as out_fh:
		file_count, total_line_count = 0, 0
		for infile in infiles:
			# file_count += 1
			##add header from first file
			with open(infile, "r") as in_fh:
				line_count = 0
				for line in in_fh:
					line_count +=1
					total_line_count += 1
					if line_count == 1:
						file_count += 1
						if file_count == 1:
							out_fh.write(line)
					else:
						out_fh.write(line)


##run methods
work_dir = '/home/atimms/ngs_data/exomes/working/dave_cili_exomes_0319'
# anal_file = 'cili_exomes_0319.txt'
# anal_file = 'cili_exomes_0319_test.txt'
##issue with NDD-238, so repeat
anal_file = 'cili_exomes_0319_rpt.txt'
##had to decompress, bgzip and tabix
combined_vcf = 'ciliopathies_exomes_2569.vcf.gz'
##run analysi
run_exome_analysis(work_dir, anal_file, combined_vcf)
##combine results file
results_suffix = 'de_novo.formatted.xls'
files_to_combine = glob.glob('*' + results_suffix)
combined_results = 'combined.' + results_suffix
combine_results(files_to_combine, combined_results)
