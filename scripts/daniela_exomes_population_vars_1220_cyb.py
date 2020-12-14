#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil


'''
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 9aa67182-09d6-405c-b7f7-93b0320828c1
module load biobuilds

'''


##parameters
delim = '\t'
star_threads = '10'

##programs
table_annovar = '/home/atimms/programs/annovar_0618/table_annovar.pl'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'


#annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,gnomad211_exome,generic', '-genericdbfile', 'microtia_exomes.non_sa_samples_1220.freq.avinput']
av_operation = ['-operation', 'g,r,r,f,f']
# av_options = ['-otherinfo', '-remove', '-nastring', '.','-vcfinput']
av_options = ['-otherinfo', '-remove', '-nastring', '.']


##methods
def filter_vcf_pass_samples(in_vcf, sa_sample_file, not_sa_sample_file, sa_vcf_file, not_sa_vcf_file):
	##get the samples we want, and remove when we don't see a call
	bcftools_view = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-a', '--threads', '2', '-Ou', '-S', sa_sample_file, in_vcf], stdout=subprocess.PIPE)
	bcftools_view2 = subprocess.Popen(['bcftools', 'view', '-m', '2', '--threads', '2', '-O', 'z', '-o', sa_vcf_file, '-'], stdin=bcftools_view.stdout)
	bcftools_view2.wait()
	##get the samples we want, and remove when we don't see a call
	bcftools_view = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-a', '--threads', '2', '-Ou', '-S', not_sa_sample_file, in_vcf], stdout=subprocess.PIPE)
	bcftools_view2 = subprocess.Popen(['bcftools', 'view', '-m', '2', '--threads', '2', '-O', 'z', '-o', not_sa_vcf_file, '-'], stdin=bcftools_view.stdout)
	bcftools_view2.wait()

def annotate_vcf(in_vcf, ref_vcf):
	'''
	##make avinput files withfreq command and then copy reference aninout to ref dir
	ref_avinput = ref_vcf.rsplit('.', 3)[0] + '.freq.avinput'
	# con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', ref_vcf, '-includeinfo', '-withfreq', '-allsample', '-outfile', ref_avinput])
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', ref_vcf, '-withfreq', '-allsample', '-outfile', ref_avinput])
	con_ann.wait()
	shutil.copy(ref_avinput, str(av_ref_dir[0]))
	'''
	in_avinput = in_vcf.rsplit('.', 3)[0] + '.freq.avinput'
	# con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', ref_vcf, '-includeinfo', '-withfreq', '-allsample', '-outfile', ref_avinput])
	# con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', in_vcf, '-includeinfo', '-withfreq', '-allsample', '-outfile', in_avinput])
	# con_ann.wait()
	out_prefix = in_vcf.rsplit('.', 3)[0]
	##annotate vcf file
	command = [table_annovar] + av_buildver + [in_avinput] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()	

def counts_and_make_ann_txt(infile, exonic_outfile, final_outfile, header, afs_req):
	with open(infile, "r") as in_fh, open(exonic_outfile, "w") as out_fh, open(final_outfile, "w") as fout_fh:
		lc, ec, sac, gnc, nsac = 0, 0,0,0,0
		out_fh.write(delim.join(header) + '\n')
		fout_fh.write(delim.join(header) + '\n')
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.strip('\n').split(delim)
				func = line[5]
				gnomad_AF = line[12]
				non_SA_AF = line[29]
				SA_AF = float(line[30])
				rmsk = line[10]
				segdup = line[11]
				if func == 'exonic' and rmsk == '.' and segdup == '.':
					ec += 1
					line_out = line[:31] + line[40:]
					out_fh.write(delim.join(line_out) + '\n')
					if gnomad_AF == '.':
						gnomad_AF = 0
					else:
						gnomad_AF = float(gnomad_AF)
					if non_SA_AF == '.':
						non_SA_AF = 0
					else:
						non_SA_AF = float(non_SA_AF)
					if SA_AF >= afs_req[0]:
						sac += 1
						if gnomad_AF <=  afs_req[1]:
							gnc += 1
							if non_SA_AF <=  afs_req[2]:
								nsac += 1
								fout_fh.write(delim.join(line_out) + '\n')
	print(lc, ec, sac, gnc, nsac)









##run methods
working_dir = '/home/atimms/ngs_data/exomes/working/daniela_exomes_recall_0720'
os.chdir(working_dir)

##files
recall_vcf = 'microtia_exomes_recall_0720.gatkHC.vcf.gz'
sa_samples = 'sa_samples.txt'
non_sa_samples = 'non_sa_samples.txt'
sa_vcf = 'microtia_exomes.sa_samples_1220.gatkHC.vcf.gz'
not_sa_vcf = 'microtia_exomes.non_sa_samples_1220.gatkHC.vcf.gz'
multi_anno = 'microtia_exomes.sa_samples_1220.hg19_multianno.txt'
exonic_ann_txt = 'microtia_exomes.exonic_vars_1220.annotated.xls'
filtered_ann_txt = 'microtia_exomes.filtered_vars_1220.annotated.xls'
aaf_values = [0.1, 0.01, 0.05] #SA, gnomad and non SA AFs 
aaf_values = [0.05, 0.01, 0.01] #SA, gnomad and non SA AFs 


##filter vcf for passed vars and south american samples
# filter_vcf_pass_samples(recall_vcf, sa_samples, non_sa_samples, sa_vcf, not_sa_vcf)

##annotate
# annotate_vcf(sa_vcf, not_sa_vcf)

##get var counts and filter vars
final_header = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 
		'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 
		'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 
		'controls_AF_popmax', 'non_SA_AF', 'SA_AF', 'info', 'format', '101000201', '101000202', '101000203', '101000301', '101000302', 
		'101000401', '101000402', '101000403', '101000501', '101000502', '101000503', '101000601', '101000602', '101000701', '101000702', 
		'101000703', '101000801', '101000802', '101000803', '101001001', '101001002', '101001003', '101001301', '101001302', '101001303', 
		'101001901', '101001902', '101001903', '101002001', '101002002', '101002003', '101002101', '101002102', '101002103', '101002201', 
		'101002202', '101002203', '101002301', '101002302', '101002303', '101002401', '101002402', '101002403', '101002501', '101002502', 
		'101002503', '101002601', '101002602', '101002603', '101002801', '101002802', '101002803', '101003101', '101003102', '101003103', 
		'101003201', '101003202', '101003203', '103000101', '103000102', '103000103', '103000201', '103000202', '103000203', '103000601', 
		'103000602', '103000603', '103000701', '103000702', '103001101', '103001102', '103001103', '103001301', '103001302', '103001303', 
		'103001401', '103001402', '103001403', '103001501', '103001502', '103001503', '106000401', '106000402', '106000403', '106001201', 
		'106001401', '106001501', '106001502', '106001503', '106001801', '106001805', '106001806a', '106001806b', '106001809', '106002001', 
		'106002002', '106002003', '106002301', '106002302', '106002303', '106002901', '106002902', '106002903', '106002907', '107000101', 
		'107000102', '107000103', '107000201', '107000202', '107000203', '107000301', '107000302', '107000303', '107001101', '107001102', 
		'107001103', '107001201', '107001202', '107001203', '107001301', '107001302', '107001303', '107001501', '107001502', '107001503']


counts_and_make_ann_txt(multi_anno, exonic_ann_txt, filtered_ann_txt, final_header, aaf_values)









