#!/usr/bin/python
import os
import subprocess
import glob
import shutil

##set input variables and parameters
delim = '\t'

'''

need 20 threads for bwa...
module load java/1.8.0_121 ##for gatk and picard


qsub -Iq cdbrmq -l mem=50gb,ncpus=2 -P 0f241fcd-0dd3-491c-acac-15c6cff1c880
##for MF -- set up env
conda activate MF

'''

##programs
bwa = '/home/atimms/programs/bwa-0.7.17/bwa'
samtools = '/home/atimms/programs/samtools-1.11/bin/samtools'
picard = '/home/atimms/programs/picard_2.19/picard.jar'
gatk = '/home/atimms/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
gatk4 = '/home/atimms/programs/gatk-4.1.3.0/gatk'
##added path to bedtools
mt2_pon_filter = '/home/atimms/programs/MosaicForecast/MuTect2-PoN_filter_AT.py'
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
ann_var = '/home/atimms/programs/annovar_1019/annotate_variation.pl'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
readlevel_fe = '/home/atimms/programs/MosaicForecast/ReadLevel_Features_extraction.py'
prediction_R = '/home/atimms/programs/MosaicForecast/Prediction.R'


##files
ref_dir = '/home/atimms/ngs_data/references/MosaicForecast_hg19/'
# fasta = ref_dir + 'hs37d5.fa'
hg19_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = hg19_dir + 'human_g1k_v37.fasta'
indels_mills = hg19_dir + 'Mills_and_1000G_gold_standard.indels.b37.vcf'
indels_1000g = hg19_dir + '1000G_phase1.indels.b37.vcf'
dbsnp = hg19_dir + 'dbsnp_138.b37.vcf'
mutect_gnomad_vcf = hg19_dir + 'af-only-gnomad.raw.sites.b37.vcf.gz'
rtc_intervals = hg19_dir + 'Mills_1000G.intervals'
final_bam_suffix = '.bwa_gatk.bam'
bamlist = 'bams.list'
mt2_pon_vcf = 'mt2_pon.vcf.gz'
mf_segdup_bed = '/home/atimms/programs/MosaicForecast/resources/SegDup_and_clustered.GRCh37.bed'
mf_indel_bed = '/home/atimms/programs/MosaicForecast/resources/allrepeats_forindel.GRCh37.bed'
k24_bw = '/home/atimms/ngs_data/references/MosaicForecast_hg19/hg19/k24.umap.wg.bw'
model_snv_50x = '/home/atimms/programs/MosaicForecast/models_trained/50xRFmodel_addRMSK_Refine.rds'
model_snv_250x = '/home/atimms/programs/MosaicForecast/models_trained/250xRFmodel_addRMSK_Refine.rds'
model_ins = '/home/atimms/programs/MosaicForecast/models_trained/insertions_250x.RF.rds'
model_del = '/home/atimms/programs/MosaicForecast/models_trained/deletions_250x.RF.rds'

##methods

def make_mutect2_pons(work_dir, samples):
	os.chdir(work_dir)
	pon_vcfs = []
	for sample in samples:
		bam = sample + final_bam_suffix
		raw_vcf = sample + '.mt2_pon.vcf.gz'
		'''
		if os.path.isfile(bam + '.bai'):
			print 'bam %s alreaded indexed'%bam
		else:
			print 'indexing bam file:', bam
			st_index = subprocess.Popen([samtools, 'index', bam])
			st_index.wait()
		##run each sample
		#gatk Mutect2 -R ref_fasta.fa -I normal1.bam -tumor normal1_sample_name --germline-resource af-only-gnomad.vcf.gz -L intervals.list -O normal1_for_pon.vcf.gz
		mt2_cmd = [gatk4, 'Mutect2', '-R', fasta, '-I', bam, '--max-mnp-distance', '0', '--tmp-dir', work_dir + '/tmp', '-O', raw_vcf]
		run_mt2 = subprocess.Popen(mt2_cmd)
		run_mt2.wait()
		'''
		pon_vcfs.append(raw_vcf)
	# '''
	##gatk GenomicsDBImport -R reference.fasta -L intervals.interval_list --genomicsdb-workspace-path pon_db -V normal1.vcf.gz \
	add_vcfs = []
	for pon_vcf in pon_vcfs:
		pvcf = ['-V', pon_vcf]
		add_vcfs.extend(pvcf)
	##intervals.list is just a file with all the chromosomes listed
	gdb_cmd = [gatk4, 'GenomicsDBImport', '-R', fasta, '-L', 'intervals.list', '--genomicsdb-workspace-path', 'pon_db', '--tmp-dir', work_dir + '/tmp'] + add_vcfs
	run_gdb = subprocess.Popen(gdb_cmd)
	run_gdb.wait()

	#gatk CreateSomaticPanelOfNormals -R reference.fasta --germline-resource af-only-gnomad.vcf.gz -V gendb://pon_db -O pon.vcf.gz
	pon_cmd = [gatk4, 'CreateSomaticPanelOfNormals', '-R', fasta, '--germline-resource', mutect_gnomad_vcf, '-V', 'gendb://pon_db', '-O', mt2_pon_vcf, '--tmp-dir', work_dir + '/tmp']
	run_pon = subprocess.Popen(pon_cmd)
	run_pon.wait()
	# '''

def align_with_bwa_one_at_time(prefix, samples):
	sample_names = [prefix + '.' + i for i in samples]
	for sample in sample_names:
		r1_fq = sample + '.R1.fastq.gz'
		r2_fq = sample + '.R2.fastq.gz'
		rg = '@RG\\tID:' + sample + '_RG\\tSM:' + sample + "\\tPL:ILLUMINA"
		post_bwa_bam = sample + '.bwa.bam'
		sort_bam = sample + '.bwa_sort.bam'
		mkdup_bam = sample + '.bwa_mkdup.bam'
		realigned_bam = sample + '.bwa_religned.bam'
		gatk_bam = sample + final_bam_suffix
		# mkdup_bam = sample + '.bwa_mkdup.bam'
		##bwa and convert to bam
		'''
		bwa_pe = subprocess.Popen([bwa, 'mem', '-M', '-t', '15', '-R', rg, fasta, r1_fq, r2_fq], stdout=subprocess.PIPE)
		st_sam_bam_pe = subprocess.Popen([samtools, 'view', '-b', '-@', '5', '-o', post_bwa_bam, '-'], stdin=bwa_pe.stdout)
		st_sam_bam_pe.wait()
		##sort bam
		st_sort_pe = subprocess.Popen([samtools, 'sort', '-O', 'bam', '-T', sample, '-o', sort_bam, '-@', '10', '-m', '10G', post_bwa_bam])
		st_sort_pe.wait()
		##mark duplicates
		picard_md = subprocess.Popen(['java', '-Xmx100g', '-jar', picard, 'MarkDuplicates', 'VALIDATION_STRINGENCY=LENIENT', 'CREATE_INDEX=true', 'INPUT=' + sort_bam, 'OUTPUT=' + mkdup_bam, 'METRICS_FILE=' + sample + '.metrics'])
		picard_md.wait()
		
		##realign around indels
		gatk_ir = subprocess.Popen(['java', '-Xmx80g', '-jar', gatk, '-T', 'IndelRealigner', '-R', fasta, '-I', mkdup_bam, '-targetIntervals', rtc_intervals, '--consensusDeterminationModel', 'KNOWNS_ONLY', '-LOD', '0.4', '-o', realigned_bam, '-known', indels_mills, '-known', indels_1000g])
		gatk_ir.wait()
		'''
		##bqsr
		gatk_br = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'BaseRecalibrator', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-knownSites', indels_mills, '-knownSites', dbsnp, '-knownSites', indels_1000g, '-o', sample + '.recal_data.table'])
		gatk_br.wait()
		gatk_pr = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'PrintReads', '-nct', '8', '-R', fasta, '-I', realigned_bam, '-BQSR', sample + '.recal_data.table', '-o', gatk_bam])
		gatk_pr.wait()


def run_mutect2_tumor_samples(prefix, samples, work_dir):
	sample_names = [prefix + '.' + i for i in samples]
	for sample in sample_names:
		bam = sample + final_bam_suffix
		raw_vcf = sample + '.mt2_raw.vcf.gz'
		temp_vcf = sample + '.mt2_temp.vcf.gz'
		filtered_vcf = sample + '.mt2_filtered.vcf.gz'
		if os.path.isfile(bam + '.bai'):
			print 'bam %s alreaded indexed'%bam
		else:
			print 'indexing bam file:', bam
			st_index = subprocess.Popen([samtools, 'index', bam])
			st_index.wait()
		##run each sample
		#gatk Mutect2 -R ref_fasta.fa  -I tumor.bam -tumor tumor_sample_name --germline-resource af-only-gnomad.vcf.gz --pon pon.vcf.gz -L intervals.list --interval-padding 100 -O tumor_unmatched_m2_snvs_indels.vcf.gz
		# mt2_cmd = [gatk4, 'Mutect2', '-R', fasta, '-I', bam, '-tumor', sample , '--germline-resource', mutect_gnomad_vcf, '--pon', mt2_pon_vcf, '--tmp-dir', work_dir + '/tmp', '-O', raw_vcf]
		# run_mt2 = subprocess.Popen(mt2_cmd)
		# run_mt2.wait()
		##gatk FilterMutectCalls -R ref_fasta.fa -v tumor_unmatched_m2_snvs_indels.vcf.gz -O tumor_unmatched_m2_snvs_indels_filtered.vcf.gz 
		mt2_filter = [gatk4, 'FilterMutectCalls', '-R', fasta, '-V', raw_vcf, '-O', temp_vcf]
		run_mt2_filter = subprocess.Popen(mt2_filter)
		run_mt2_filter.wait()
		##using single or 10 threads with no sed command as ad good in vcf
		bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '-w', '10000', '-f', fasta, '-o', filtered_vcf, '-O', 'z', temp_vcf])
		# bcftools_norm = subprocess.Popen([bcftools, 'norm', '-m', '-', '--threads', '10', '-w', '10000', '-f', fasta, '-o', temp_vcf, '-O', 'z', input_vcf])
		bcftools_norm.wait()

def filter_annovar_data(ann_filtered, ann_dropped, out_maf0, out_maf1):
	temp_snv_bed = ann_filtered.rsplit('.',2)[0] + '.snv.bed'
	temp_indel_bed1 = ann_filtered.rsplit('.',2)[0] + '.indel1.bed'
	temp_indel_bed2 = ann_filtered.rsplit('.',2)[0] + '.indel2.bed'
	with open(temp_snv_bed, "w") as tsb_fh, open(temp_indel_bed1, "w") as tib1_fh:
		##get the filtered results i.e. not in gnomad
		with open(ann_filtered, "r") as af_fh:
			for line in af_fh:
				##remove unwanted chromosome
				if line[:2] != 'GL':
					line = line.rstrip().split(delim)
					bed_line = [line[0], str(int(line[1]) -1)] + line[2:] + ['.']
					##get SNVs
					if len(line[3]) == 1 and len(line[4]) == 1:
						tsb_fh.write(delim.join(bed_line) + '\n')
					else:
						tib1_fh.write(delim.join(bed_line) + '\n')
		##get the dropped results i.e. filter now for af
		with open(ann_dropped, "r") as ad_fh:
			lc = 0
			for line in ad_fh:
				line = line.rstrip().split(delim)
				##remove unwanted chromosome
				if not line[2].startswith('GL'):
					maf = line[1]
					if maf == '.':
						maf = '0'
					if float(maf) <= 0.001:
						if maf == '0':
							bed_line = [line[2], str(int(line[3]) - 1)] + line[4:] + ['0']
						else:
							bed_line = [line[2], str(int(line[3]) - 1)] + line[4:] + ['1']
						##get SNVs
						if len(line[5]) == 1 and len(line[6]) == 1:
							tsb_fh.write(delim.join(bed_line) + '\n')
						else:
							tib1_fh.write(delim.join(bed_line) + '\n')
	##remove indels in known regions
	with open(temp_indel_bed2, "w") as out_fh:
		manta_config = subprocess.Popen([bedtools, 'subtract', '-a', temp_indel_bed1, '-b', mf_indel_bed], stdout=out_fh)
		manta_config.wait()
	##combine indels and snvs into final files
	with open(out_maf0, 'w') as m0_fh, open(out_maf1, 'w') as m1_fh, open(temp_snv_bed, "r") as tsb_fh, open(temp_indel_bed2, "r") as tib2_fh:
		for line in tsb_fh:
			m1_fh.write(line)
			if line.rstrip().endswith('0') or line.rstrip().endswith('.'):
				m0_fh.write(line)
		for line in tib2_fh:
			m1_fh.write(line)
			if line.rstrip().endswith('0') or line.rstrip().endswith('.'):
				m0_fh.write(line)

def filter_mt2_calls(prefix, samples):
	sample_names = [prefix + '.' + i for i in samples]
	for sample in sample_names:
		filtered_vcf = sample + '.mt2_filtered.vcf.gz'
		mt2_filtered_bed = sample + '.bed'
		annovar_file = sample + '.temp.avinput'
		annovar_filtered = sample + '.temp.avinput.hg19_gnomad_genome_filtered'
		annovar_dropped = sample + '.temp.avinput.hg19_gnomad_genome_dropped'
		annovar_maf0 = sample + ".MAF0.bed"
		annovar_maf1 = sample + ".MAF1.bed"
		# '''
		##filter mt2 calls
		#python MuTect2-PoN_filter.py test demo/test.Mutect2.vcf resources/SegDup_and_clustered.GRCh37.bed
		run_mt2 = subprocess.Popen(['python', mt2_pon_filter, sample, filtered_vcf, mf_segdup_bed])
		run_mt2.wait()
		##annovar i.e. make file and use annotate variants
		make_annovar="cat "+mt2_filtered_bed+"|awk '\''{{OFS=\"\\t\";len=length($4)-length($5);if(len<=0){{print $1,$3,$3,$4,$5,$6}}if(len>0){{print $1,$3,$3+len,$4,$5,$6}}}}'\''>"+annovar_file
		subprocess.Popen(make_annovar, shell=True, stdout=subprocess.PIPE).communicate()
		##"annotate_variation.pl -filter -buildver hg19 -dbtype gnomad_genome {input} /n/data1/hms/dbmi/park/yanmei/tools/annovar/humandb/ --outfile MF/{wildcards.sample}.mt2pon.AF0.02.ANNOVAR.list"
		annovar_vars = subprocess.Popen([ann_var, '-filter', '-buildver', 'hg19', '-dbtype', 'gnomad_genome', annovar_file, '/home/atimms/ngs_data/references/annovar/hg19'])
		annovar_vars.wait()
		# '''
		##filter annovar calls
		filter_annovar_data(annovar_filtered, annovar_dropped, annovar_maf0, annovar_maf1)

def extract_read_level_features(prefix, samples):
	sample_names = [prefix + '.' + i for i in samples]
	for sample in sample_names:
		annovar_maf0 = sample + ".MAF0.bed"
		annovar_maf1 = sample + ".MAF1.bed"
		features_maf0 = sample + ".MAF0.features"
		features_maf1 = sample + ".MAF1.features"
		bam_dir = sample + '_bamdir'
		# '''
		##cp bam files
		make_bamdir = subprocess.Popen(['mkdir', bam_dir])
		make_bamdir.wait()
		cp_bam = subprocess.Popen(['cp', sample + '.bwa_gatk.bam', bam_dir + '/' + sample + '.bam'])
		cp_bam.wait()
		cp_bai = subprocess.Popen(['cp', sample + '.bwa_gatk.bam.bai', bam_dir + '/' + sample + '.bam.bai'])
		cp_bai.wait()
		# '''
		##using 2 threads
		#python ReadLevel_Features_extraction.py demo/test.input demo/test.features demo ${ref.fa} ${k24.umap.wg.bw} 2 bam  
		run_mt2 = subprocess.Popen(['python', readlevel_fe, annovar_maf0, features_maf0, bam_dir, fasta, k24_bw, '2', 'bam'])
		run_mt2.wait()
		run_mt2_2 = subprocess.Popen(['python', readlevel_fe, annovar_maf1, features_maf1, bam_dir, fasta, k24_bw, '2', 'bam'])
		run_mt2_2.wait()

def split_feature_file(infile):
	snv_outfile = infile + '.snv.temp'
	ins_outfile = infile + '.ins.temp'
	del_outfile = infile + '.del.temp'
	with open(infile, "r") as in_fh, open(snv_outfile, "w") as so_fh, open(del_outfile, "w") as do_fh, open(ins_outfile, "w") as io_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			line = line.split(delim)
			if lc == 1:
				so_fh.write(delim.join(line))
				io_fh.write(delim.join(line))
				do_fh.write(delim.join(line))
			else:
				mut_type = line[3]
				if mut_type == 'SNP':
					so_fh.write(delim.join(line))
				elif mut_type == 'DEL':
					do_fh.write(delim.join(line))
				elif mut_type == 'INS':
					io_fh.write(delim.join(line))
				else:
					print(mut_type + ' not recognized')

def combine_pred_files(file_prefix):
	outfile = file_prefix + '.all.xls'
	with open(outfile, "w") as out_fh:
		fc = 0
		for infile in [file_prefix + '.snv.temp', file_prefix + '.ins.temp', file_prefix + '.del.temp']:
			with open(infile, "r") as in_fh:
				fc += 1
				lc = 0
				for line in in_fh:
					lc += 1
					print(fc, lc)
					if lc == 1:
						if fc == 1:
							out_fh.write(line)
					else:
						out_fh.write(line)


def genotype_prediction(prefix, samples):
	sample_names = [prefix + '.' + i for i in samples]
	for sample in sample_names:
		features_maf0 = sample + ".MAF0.features"
		features_maf1 = sample + ".MAF1.features"
		predictions_maf0_prefix = sample + ".MAF0.predictions"
		predictions_maf1_prefix = sample + ".MAF1.predictions"
		# '''
		##split features by type i.e snv, del and ins
		split_feature_file(features_maf0)
		split_feature_file(features_maf1)
		##predict SNVs
		##Rscript Prediction.R demo/test.SNP.features models_trained/250xRFmodel_addRMSK_Refine.rds Refine demo/test.SNP.predictions   
		if 'comb' in sample:
			snv_model = model_snv_250x
		else:
			snv_model = model_snv_50x
		run_gp_snv0 = subprocess.Popen(['Rscript', prediction_R, features_maf0 + '.snv.temp', snv_model, 'Refine', predictions_maf0_prefix + '.snv.temp'])
		run_gp_snv0.wait()
		run_gp_snv1 = subprocess.Popen(['Rscript', prediction_R, features_maf1 + '.snv.temp', snv_model, 'Refine', predictions_maf1_prefix + '.snv.temp'])
		run_gp_snv1.wait()
		##predict ins/deletions
		##Rscript Prediction.R demo/test.DEL.features models_trained/deletions_250x.RF.rds Phase demo/test.DEL.predictions
		run_gp_ins0 = subprocess.Popen(['Rscript', prediction_R, features_maf0 + '.ins.temp', model_ins, 'Phase', predictions_maf0_prefix + '.ins.temp'])
		run_gp_ins0.wait()
		run_gp_ins1 = subprocess.Popen(['Rscript', prediction_R, features_maf1 + '.ins.temp', model_ins, 'Phase', predictions_maf1_prefix + '.ins.temp'])
		run_gp_ins1.wait()
		run_gp_del0 = subprocess.Popen(['Rscript', prediction_R, features_maf0 + '.del.temp', model_ins, 'Phase', predictions_maf0_prefix + '.del.temp'])
		run_gp_del0.wait()
		run_gp_del1 = subprocess.Popen(['Rscript', prediction_R, features_maf1 + '.del.temp', model_ins, 'Phase', predictions_maf1_prefix + '.del.temp'])
		run_gp_del1.wait()
		# '''
		##combine prediction files
		combine_pred_files(predictions_maf0_prefix)
		combine_pred_files(predictions_maf1_prefix)



def run_all_methods(work_dir, samples, project_prefix):
	os.chdir(work_dir)
	##map fastqs -- not used, just use original bams
	# align_with_bwa_one_at_time(project_prefix, samples)
	##run mutect2
	# run_mutect2_tumor_samples(project_prefix, samples, work_dir)
	##filter mutect2 calls
	# filter_mt2_calls(project_prefix, samples)
	##extract read-level features
	extract_read_level_features(project_prefix, samples)
	##Genotype Prediction
	genotype_prediction(project_prefix, samples)
	##phase -- used to generate model, start by using their models
	# phase_calls()



##run methods
working_dir = '/home/atimms/ngs_data/genomes/jimmy_fetal_genomes_mf_0321'


##MF details
#https://github.com/parklab/MosaicForecast
#https://github.com/parklab/MosaicForecast/blob/master/FAQ.md

##make panel of normals from mom and dad data
##mutect details: https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
normal_samples = ['genome_0819.Mom', 'genome_0819.Dad', 'LR16-302', 'LR16-302f', 'LR16-302m']
# normal_samples = ['genome_0819.Mom']
# normal_samples = ['genome_0819.Dad']
# normal_samples = ['LR16-302']
# normal_samples = ['LR16-302f']
# normal_samples = ['LR16-302m']
# make_mutect2_pons(working_dir, normal_samples)


##for genomes 0819
project_name_0819 = 'genome_0819'
# samples_0819 = ['Brain', 'Cerebellum', 'Heart', 'Kidney', 'Lung', 'Mom', 'Dad', 'Child_combined']
samples_0819 = ['Brain', 'Cerebellum', 'Heart', 'Kidney', 'Lung', 'Child_combined']
# samples_0819 = ['Brain', 'Cerebellum', 'Heart']
# samples_0819 = ['Kidney', 'Lung']
# samples_0819 = ['Child_combined']
# samples_0819 = ['Brain']
# run_all_methods(working_dir, samples_0819, project_name_0819)


##for genomes 0816
project_name_0816 = 'genome_0816'
samples_0816 = ['Brain', 'Liver', 'Heart', 'Lung', 'Muscle', 'Placenta', 'fg_combined']
# samples_0816 = ['Brain', 'Liver', 'Heart']
# samples_0816 = ['Muscle', 'Placenta', 'Lung']
# samples_0816 = ['fg_combined']
# run_all_methods(working_dir, samples_0816, project_name_0816)


def combine_vars(infiles, outfile):
	# print(infiles)
	var_dict = {}
	samples = []
	for infile in infiles:
		sample = infile.split('.')[1]
		samples.append(sample)
		with open(infile, "r") as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				if lc >1:
					line = line.rstrip().split(delim)
					var = '_'.join(line[0].split('~')[1:])
					var_type = line[34]
					# print(var, var_type)
					if var_type == 'mosaic':
						if var in var_dict:
							var_dict[var].append(sample)
						else:
							var_dict[var] = [sample]
	print(samples)
	with open(outfile, "w") as out_fh:
		header = ['var', 'sample_count', 'samples']
		out_fh.write(delim.join(header) + '\n')
		for v in var_dict:
			print(v, len(var_dict[v]))
			line_out = [v, str(len(var_dict[v])), ', '.join(var_dict[v])]
			out_fh.write(delim.join(line_out) + '\n')

##compare variants in all samples
os.chdir(working_dir)
prediction_files = glob.glob('genome_0816*MAF0.predictions.all.xls')
combined_var_file = 'genome_0816_MAF0_variants.xls'
combine_vars(prediction_files, combined_var_file)

prediction_files = glob.glob('genome_0819*MAF0.predictions.all.xls')
combined_var_file = 'genome_0819_MAF0_variants.xls'
combine_vars(prediction_files, combined_var_file)



