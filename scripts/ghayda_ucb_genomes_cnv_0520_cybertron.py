#!/usr/bin/env python
import os
import subprocess
import glob
import shutil


'''
load these for canvas:
module load dotnet/2.1.805
for formatting manta:
module load biobuilds
'''

##set input variables and parameters
delim = '\t'

##programs
canvas_dll = '/home/atimms/programs/Canvas-1.40.0.1613+master_x64/Canvas.dll'
manta_config = '/home/atimms/programs/manta-1.6.0.centos6_x86_64/bin/configManta.py'
manta_denovo = '/home/atimms/programs/manta-1.6.0.centos6_x86_64/libexec/denovo_scoring.py'
vcf_bed = '/home/atimms/programs/svtools/vcfToBedpe'
##files etc
fasta = '/home/atimms/ngs_data/references/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta'
hg38_chr_bed = '/home/atimms/ngs_data/references/hg38/hg38.main_chr.bed.gz'
template_ploidy_vcf = '/home/atimms/ngs_data/references/canvas/GRCh38/template.ploidy.vcf'
hg38_refgene_exons = '/home/atimms/ngs_data/references/hg38/hg38_refFlat_coding_exons.dup_removed.bed'

##sample dict
sample_dict = {'UC1402089': 'LR13-282_UC1402089', 'UC1309030': 'LR13-282_UC1309030', 'UC1403115': 'LR13-282f_UC1403115', 
		'UC1309021': 'LR13-282m_UC1309021', 'UC1405044': 'LR14-155_UC1405044', 'UC1406081': 'LR14-155_UC1406081', 
		'UC1406080': 'LR14-155f_UC1406080', 'UC1406082': 'LR14-155m_UC1406082', 'UC1711041': 'LR17-408_UC1711041', 
		'UC1711046': 'LR17-408_UC1711046', 'UC1711032': 'LR17-408f_UC1711032', 'UC1711033': 'LR17-408m_UC1711033', 
		'UC1210009': 'LR12-123_UC1210009', 'UC1401004': 'LR12-123_UC1401004', 'UC1410166': 'LR12-123f_UC1410166', 
		'UC1205034': 'LR12-123m_UC1205034', 'UC1401017': 'LR12-269_UC1401017', 'UC1109008': 'LR12-269_UC1109008', 
		'UC1305096': 'LR12-269f_UC1305096', 'UC1305097': 'LR12-269m_UC1305097', 'UC1401005': 'LR12-257_UC1401005', 
		'UC1906025': 'LR12-257_UC1906025', 'UC1906053': 'LR12-257f_UC1906053', 'UC1906054': 'LR12-257m_UC1906054', 
		'UC1812014': 'LR16-432a1_UC1812014', 'UC1611017': 'LR16-432a1_UC1611017', 'UC1908020': 'LR16-432f_UC1908020', 
		'UC1908021': 'LR16-432m_UC1908021', 'UC1401007': 'LR12-259_UC1401007', 'UC1305007': 'LR12-259_UC1305007', 'UC1305008': 'LR12-259f_UC1305008', 'UC1305009': 'LR12-259m_UC1305009', 'UC1401013_2': 'LR12-265_UC1401013', 'UC1311084': 'LR12-265_UC1311084', 'UC1302036': 'LR12-265f_UC1302036', 'UC1302037': 'LR12-265m_UC1302037', 'UC1701038-2': 'LR12-101_UC1701038', 'UC1206006': 'LR12-101_UC1206006', 'UC1206140': 'LR12-101f_UC1206140', 'UC1206141': 'LR12-101m_UC1206141', 'UC1401039': 'LR12-249_UC1401039', 'UC1908029': 'LR12-249_UC1908029', 'UC1908030': 'LR12-249f_UC1908030', 'UC1908031': 'LR12-249m_UC1908031', 'UC1708005': 'LR17-337_UC1708005', 'UC1906037': 'LR17-337_UC1906037', 'UC1906038': 'LR17-337f_UC1906038', 'UC1906039': 'LR17-337m_UC1906039', 'UC1308074': 'LR10-246_UC1308074', 'UC1504075': 'LR10-246f_UC1504075', 'UC1308075': 'LR10-246m_UC1308075', 'UC1604151': 'LR16-214_UC1604151', 'UC1605024': 'LR16-214f_UC1605024', 'UC1604152': 'LR16-214m_UC1604152'}


##methods
def make_analysis_dict(input_file):
	analysis_dict = {}
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
				affection = line[5]
				sex = line[4]
			##add data to analysis_dict
				if ped_name in analysis_dict:
					if id_name in analysis_dict[ped_name][0]:
						print('issue with sample', id_name, sample_name)
					else:
						analysis_dict[ped_name][0][id_name] = sample_name
						if affection == '2':
							analysis_dict[ped_name][1].append(id_name)
							analysis_dict[ped_name][4].append(sex)
						elif affection == '1':
							analysis_dict[ped_name][2].append(id_name)
				else:
					##if affected sample
					if affection == '2':
						analysis_dict[ped_name] = [{id_name:sample_name}, [id_name], [], anal_type, [sex]]
					##or parent
					elif affection == '1':
						analysis_dict[ped_name] = [{id_name:sample_name}, [], [id_name], anal_type]
	return analysis_dict

def make_plody_vcf(out_vcf, proband, parents, proband_sex):
	with open(template_ploidy_vcf, "r") as in_fh, open(out_vcf, "w") as out_fh:
		for line in in_fh:
			if line.startswith('##'):
				out_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				print(line)
				if line[0] == '#CHROM':
					out_fh.write(delim.join(line + [proband] + parents) + '\n')

				elif line[0] == 'chrX':
					if proband_sex == '1':
						out_fh.write(delim.join(line + ['1', '1', '2']) + '\n')
					if proband_sex == '2':
						out_fh.write(delim.join(line + ['2', '1', '2']) + '\n')						
				elif line[0] == 'chrY':
					if proband_sex == '1':
						out_fh.write(delim.join(line + ['1', '1', '0']) + '\n')
					if proband_sex == '2':
						out_fh.write(delim.join(line + ['0', '1', '0']) + '\n')		





def run_canvas_sp_wgs(ped_dict, c_ref_dir):
	##reference files
	kmer_fa = c_ref_dir + 'kmer.fa'
	dbsnp_vcf = c_ref_dir + 'dbsnp.vcf'
	filter13_bed = c_ref_dir + 'filter13.bed'

	##go through sample by sample
	for ped in ped_dict:
		##call variants individually
		affected_samples = ped_dict[ped][1]
		parent_ids = ped_dict[ped][2]
		new_id_dict = ped_dict[ped][0]
		pro_sex = ped_dict[ped][4][0]
		##go sample by sample
		# '''
		for affected_sample in affected_samples:
			outdir = affected_sample + '.canvas_sp_wgs'
			ploidy_vcf = c_ref_dir + affected_sample + '.ploidy.vcf'
			# print(affected_sample, parent_ids, pro_sex, ped)
			make_plody_vcf(ploidy_vcf, affected_sample, parent_ids, pro_sex)
			##run canvas SP-WGS
			#dotnet Canvas.dll SmallPedigree-WGS --bam=father.bam --bam=mother.bam --bam=child1.bam --mother=mother --father=father --proband=child1 -r kmer.fa -g WholeGenomeFasta --sample-b-allele-vcf Pedigree.vcf.gz -f filter13.bed -o /tmp/gHapMixDemo --ploidy-vcf=MultiSamplePloidy.vcf
			canvas_run = subprocess.Popen(['dotnet', canvas_dll, 'SmallPedigree-WGS', '-b', parent_ids[0] + '.bam', 'father', parent_ids[0], 
				'-b', parent_ids[1] + '.bam', 'mother', parent_ids[1], '-b', affected_sample + '.bam', 'proband', affected_sample, 
				'-r', kmer_fa, '-g', c_ref_dir, '--population-b-allele-vcf=' + dbsnp_vcf,  '-f', filter13_bed, '-o', outdir, 
				'--ploidy-vcf=' + ploidy_vcf])
			canvas_run.wait()
		# '''

def filter_denovo_manta_vcf(in_vcf, all_file, exonic_file, no_of_as):
	temp_vcf = 'temp.vcf'
	temp_bed = 'temp.bed'
	temp2_bed = 'temp2.bed'
	##filter for de novo using dq ==60
	with open(temp_vcf, "w") as out_fh:
		bcftools_filter = subprocess.Popen(['bcftools', 'filter', "-e", "FORMAT/DQ<50", in_vcf], stdout=out_fh)
		bcftools_filter.wait()

	##make vcf into bed
	manta_config = subprocess.Popen([vcf_bed, '-i', temp_vcf, '-o', temp_bed])
	manta_config.wait()
	##get header from bed file
	with open(temp_bed, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc == 1:
				header = line.rstrip().split(delim) + ['chr_exon', 'exon_start', 'exon_end', 'exon_name']

	##use bedtools to add genename
	with open(temp2_bed, "w") as out_fh:
		manta_config = subprocess.Popen(['bedtools', 'intersect', '-a', temp_bed, '-b', hg38_refgene_exons, '-wao'], stdout=out_fh)
		manta_config.wait()		

	with open(temp2_bed, "r") as in_fh, open(all_file, "w") as all_fh, open(exonic_file, "w") as ex_fh:
		all_fh.write(delim.join(header) + '\n')
		ex_fh.write(delim.join(header) + '\n')
		for line in in_fh:
			line = line.split(delim)
			pass_filter = line[11]
			if pass_filter == 'PASS':
				if no_of_as == 1:
					exon_name = line[20]
					all_fh.write(delim.join(line[:21]) + '\n')
					if exon_name != '.':
						ex_fh.write(delim.join(line[:21]) + '\n')
				elif no_of_as == 2:
					exon_name = line[21]
					all_fh.write(delim.join(line[:22]) + '\n')
					if exon_name != '.':
						ex_fh.write(delim.join(line[:22]) + '\n')

def run_manta_looking_for_denovo_svs(ped_dict):
	##go through ped by ped
	for ped in ped_dict:
		##call variants individually
		affected_samples = ped_dict[ped][1]
		parent_ids = ped_dict[ped][2]
		##get all bams
		bams = []
		for sample in affected_samples + parent_ids:
			bam = ['--bam', sample + '.bam']
			bams.extend(bam)
		outdir = ped + '.manta_joint'
		##call svs
		'''
		make_manta_config = subprocess.Popen([manta_config] + bams + ['--referenceFasta', fasta, '--runDir', outdir, '--callRegions', hg38_chr_bed])
		make_manta_config.wait()
		manta_run = subprocess.Popen([outdir + '/runWorkflow.py', '-j', '20'])
		manta_run.wait()
		'''
		##get de novos
		original_vcf_gz = outdir + '/results/variants/diploidSV.vcf.gz'
		original_vcf = outdir + '/results/variants/diploidSV.vcf'
		##gzip so can be read by denovo command
		'''
		gzip_d_vcf = subprocess.Popen(['gzip', '-d', original_vcf_gz])
		gzip_d_vcf.wait()
		'''
		##run de novo scripts and change to specific file name
		de_novo_vcfs = []	
		for affected_sample in affected_samples:
			# print(affected_sample, parent_ids)
			'''
			#denovo_scoring.py <vcf file> <proband sample ID> <father sample ID> <mother sample ID>
			manta_dn_run = subprocess.Popen([manta_denovo, original_vcf, affected_sample] + parent_ids)
			manta_dn_run.wait()
			'''
			##rename/mv results files
			in_vcf = outdir + '/results/variants/diploidSV.de_novo.vcf'
			in_stats = outdir + '/results/variants/diploidSV.de_novo.stats.txt'
			out_vcf = outdir + '/' + affected_sample + '.de_novo.vcf'
			out_stats = outdir + '/' + affected_sample + '.de_novo.stats.txt'
			all_denovos = sample_dict[affected_sample] + '.manta_sv.de_novo.xls'
			exonic_denovos = sample_dict[affected_sample] + '.manta_sv.de_novo_exonic.xls'

			'''
			shutil.move(in_vcf, out_vcf)
			shutil.move(in_stats, out_stats)
			'''
			##filter and annotate denovo vcf
			filter_denovo_manta_vcf(out_vcf, all_denovos,exonic_denovos, len(affected_samples))

		##gzip original vcf
		'''
		gzip_vcf = subprocess.Popen(['gzip', original_vcf])
		gzip_vcf.wait()
		'''


def filter_somatic_manta_vcf(sample_name, mom_vcf, dad_vcf, all_file, exonic_file, no_of_as):
	temp_vcf = 'temp.vcf'
	mom_temp_bed = sample_name + 'mom_temp.bed'
	dad_temp_bed = sample_name + 'dad_temp.bed'
	mom_temp_bed2 = sample_name + 'mom_temp2.bed'
	dad_temp_bed2 = sample_name + 'dad_temp2.bed'
	intersect_temp_bed = sample_name + 'intersect_temp.bed'
	# all_temp = sample_name + 'temp_all.txt'
	# ex_temp = sample_name + 'temp_ex.txt'
	##make vcfs into bed
	mom_uncomp = subprocess.Popen(['gunzip', '-d', mom_vcf])
	mom_uncomp.wait()
	dad_uncomp = subprocess.Popen(['gunzip', '-d', dad_vcf])
	dad_uncomp.wait()
	mom_vcf_bed = subprocess.Popen([vcf_bed, '-i', mom_vcf.rsplit('.',1)[0], '-o', mom_temp_bed])
	mom_vcf_bed.wait()
	dad_vcf_bed = subprocess.Popen([vcf_bed, '-i', dad_vcf.rsplit('.',1)[0], '-o', dad_temp_bed])
	dad_vcf_bed.wait()
	##use bedtools to add genename
	with open(mom_temp_bed2, "w") as out_fh:
		manta_config = subprocess.Popen(['bedtools', 'intersect', '-a', mom_temp_bed, '-b', hg38_refgene_exons, '-wao', '-header'], stdout=out_fh)
		manta_config.wait()		
	with open(dad_temp_bed2, "w") as out_fh:
		manta_config = subprocess.Popen(['bedtools', 'intersect', '-a', dad_temp_bed, '-b', hg38_refgene_exons, '-wao', '-header'], stdout=out_fh)
		manta_config.wait()	

	##use bedtools to SVs in mom and dad
	with open(intersect_temp_bed, "w") as out_fh:
		manta_config = subprocess.Popen(['bedtools', 'intersect', '-a', mom_temp_bed2, '-b', dad_temp_bed2, '-wo'], stdout=out_fh)
		manta_config.wait()		


	##get header from bed files
	header = []
	with open(mom_temp_bed2, "r") as m_fh, open(dad_temp_bed2, "r") as d_fh:
		mlc, dlc = 0, 0
		for line in m_fh:
			mlc += 1
			if mlc == 1:
				h1 = line.rstrip().split(delim) + ['chr_exon', 'exon_start', 'exon_end', 'exon_name']
				header.extend(h1)
		for line in d_fh:
			dlc += 1
			if dlc == 1:
				h2 = line.rstrip().split(delim) + ['chr_exon', 'exon_start', 'exon_end', 'exon_name']
				header.extend(h2)
	header = header + ['bp_intesect']


	##filter etc
	with open(intersect_temp_bed, "r") as in_fh, open(all_file, "w") as all_fh, open(exonic_file, "w") as ex_fh:
		all_fh.write(delim.join(header) + '\n')
		ex_fh.write(delim.join(header) + '\n')
		for line in in_fh:
			line = line.rstrip().split(delim)
			pass_filter = [line[11], line[34]]
			if pass_filter[0] == 'PASS' and pass_filter[1] == 'PASS':
				exon_name = [line[19], line[42]]
				line_out = line[:20] + line[23:43] + [line[46]]
				# print(intersect_temp_bed, exon_name)
				all_fh.write(delim.join(line) + '\n')
				if exon_name[0] != '.' or exon_name[1] != '.':
					ex_fh.write(delim.join(line) + '\n')
	##remove duplicate lines
	# with open(all_file, "w") as all_fh:
	# 	#sort garbage.txt | uniq -u
	# 	rm_dup1 = subprocess.Popen(['sort', all_temp], stdout=subprocess.PIPE)
	# 	rm_dup2 = subprocess.Popen(['uniq', '-u'], stdin=rm_dup1.stdout, stdout=all_fh)
	# 	rm_dup2.wait()
	# with open(exonic_file, "w") as ex_fh:
	# 	#sort garbage.txt | uniq -u
	# 	rm_dup1 = subprocess.Popen(['sort', ex_temp], stdout=subprocess.PIPE)
	# 	rm_dup2 = subprocess.Popen(['uniq', '-u'], stdin=rm_dup1.stdout, stdout=ex_fh)
	# 	rm_dup2.wait()		

def run_manta_looking_for_somatic_svs(ped_dict):
	##go through ped by ped
	for ped in ped_dict:
		##call variants individually
		affected_samples = ped_dict[ped][1]
		parent_ids = ped_dict[ped][2]
		##get all proband bams
		for sample in affected_samples:
			tbam = ['--tumorBam', sample + '.bam']
			sample_out_dir = sample_dict[sample] + '.manta_somatic/'
			##call SVs in dad
			outdir = sample_out_dir + 'dad_analysis'
			'''
			make_manta_config = subprocess.Popen([manta_config, '--normalBam', parent_ids[0] + '.bam'] + tbam + ['--referenceFasta', fasta, '--runDir', outdir, '--callRegions', hg38_chr_bed])
			make_manta_config.wait()
			manta_run = subprocess.Popen([outdir + '/runWorkflow.py', '-j', '20'])
			manta_run.wait()
			##call SVs in mom
			outdir = sample_out_dir + 'mom_analysis'
			make_manta_config = subprocess.Popen([manta_config, '--normalBam', parent_ids[1] + '.bam'] + tbam + ['--referenceFasta', fasta, '--runDir', outdir, '--callRegions', hg38_chr_bed])
			make_manta_config.wait()
			manta_run = subprocess.Popen([outdir + '/runWorkflow.py', '-j', '20'])
			manta_run.wait()
			'''
			##filter and add genes/exons to somatic vcf
			mom_manta_vcf = sample_out_dir + 'mom_analysis/results/variants/somaticSV.vcf.gz'
			dad_manta_vcf = sample_out_dir + 'dad_analysis/results/variants/somaticSV.vcf.gz'
			all_somatic = sample_dict[sample] + '.manta_sv.somatic.xls'
			exonic_somatic = sample_dict[sample] + '.manta_sv.somatic_exonic.xls'
			filter_somatic_manta_vcf(sample_dict[sample], mom_manta_vcf, dad_manta_vcf, all_somatic, exonic_somatic, len(affected_samples))


##master method
def run_analysis(working_dir, analysis_file, canvas_ref_dir):
	os.chdir(working_dir)
	project = analysis_file.split('.')[0]
	##read in analysis file and make 2 dicts
	analysis_dict = make_analysis_dict(analysis_file)
	##check...
	for p in analysis_dict:
		print(p, analysis_dict[p])
	##run canvas small pedigree (de novo cnvs)
	# run_canvas_sp_wgs(analysis_dict, canvas_ref_dir)
	##run manta for de novo analysis
	# run_manta_looking_for_denovo_svs(analysis_dict)
	##run manta for somatic analysis
	run_manta_looking_for_somatic_svs(analysis_dict)









##run methods
work_dir = '/home/atimms/ngs_data/genomes/ghayda_ucb_0320'
ref_dir = '/home/atimms/ngs_data/references/canvas/GRCh38/'
##test on two peds
# sample_file = 'ghayda_ucb_0320_1.txt'
##the rest of the peds
# sample_file = 'ghayda_ucb_0320_2.txt'
##all samples
sample_file = 'ghayda_ucb_0320_all.txt'
##run master method
run_analysis(work_dir, sample_file, ref_dir)


