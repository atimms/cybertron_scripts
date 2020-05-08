#!/usr/bin/env python
import subprocess
import os


##parameters
delim = '\t'
thread_number = '20'

##setup 
'''
so need to install miso:
module load local_python/2.7.14
python -m pip install misopy ##dependancies already satisfied except bedtools
conda install -c bioconda bedtools

some files missing for pip version, so used git to get repo, and transferred from /home/atimms/programs/MISO/misopy
cp -r settings /home/atimms/programs/anaconda2/lib/python2.7/site-packages/misopy

##also download gff file and move to cluster:
curl -o ensGene.gff3 http://hollywood.mit.edu/burgelab/miso/annotations/ucsc_tables/mm10/ensGene.gff3
or: http://hollywood.mit.edu/burgelab/miso/annotations/gene-models/Mus_musculus.NCBIM37.65.gff.zip


'''

##working dir
# working_dir = '/data/atimms/timon_0317'
working_dir = '/home/atimms/ngs_data/rnaseq/sophie_carm1_1119'
os.chdir(working_dir)

##files
# mm10_gff = '/home/atimms/ngs_data/references/mm10/Mus_musculus.NCBIM37.65.gff'
# mm10_gff_index = '/home/atimms/ngs_data/references/mm10/miso_index_Mus_musculus_NCBIM37_65'
mm10_gff = '/home/atimms/ngs_data/references/mm10/ensGene.gff3'
mm10_gff_index = '/home/atimms/ngs_data/references/mm10/miso_index_mm10_ensGene_gff3'
const_exons = working_dir + '/const_exons'
const_exon_mm10 = working_dir + '/const_exons/ensGene.min_1000.const_exons.gff'
insert_sizes = working_dir + '/insert-dist'
fq_read_length = '150'


##methods
def get_const_exons(gff_file, out_dir):
	#exon_utils --get-const-exons Mus_musculus.NCBIM37.65.gff --min-exon-size 1000 --output-dir exons/
	miso_ce = subprocess.Popen(['exon_utils', '--get-const-exons', gff_file, '--min-exon-size', '1000', '--output-dir', out_dir])
	miso_ce.wait()

def get_insert_stats(bam_files, out_dir, ce_gff):
	#pe_utils --compute-insert-len sample.bam Mus_musculus.NCBIM37.65.min_1000.const_exons.gff --output-dir insert-dist/
	bam_file_comma = ','.join(bam_files)
	# miso_is = subprocess.Popen(['pe_utils', '--compute-insert-len', bam_file_comma, ce_gff, '--output-dir', out_dir, '--no-bam-filter'])
	miso_is = subprocess.Popen(['pe_utils', '--compute-insert-len', bam_file_comma, ce_gff, '--output-dir', out_dir])
	miso_is.wait()


def index_gff_file(gff, gff_index_dir):
	#index_gff --index SE.gff3 indexed_SE_events/
	miso_index = subprocess.Popen(['index_gff', '--index', gff, gff_index_dir])
	miso_index.wait()


def run_miso(sample_dict, read_length, gff_index):
	#miso --run mm9/pickled/SE data/control.bam --output-dir SE/control/ --read-len 35 --paired-end 250 15 --use-cluster
	#summarize_miso --summarize-samples SE/control/ SE/control/
	for sample in sample_dict:
		bam = sample + '.Aligned.sortedByCoord.out.bam'
		out_dir = working_dir + '/' + sample + '_miso_output'
		insert_size = sample_dict[sample][0]
		is_sd = sample_dict[sample][1]
		print(sample, insert_size, is_sd)
		# miso_run = subprocess.Popen(['miso', '--run', gff_index, bam, '--output-dir', out_dir, '--read-len', read_length, '--paired-end', insert_size, is_sd, '--use-cluster'])
		miso_run = subprocess.Popen(['miso', '--run', gff_index, bam, '--output-dir', out_dir, '--read-len', read_length, '--paired-end', insert_size, is_sd])
		miso_run.wait()
		miso_sum = subprocess.Popen(['summarize_miso', '--summarize-samples', out_dir, out_dir])
		miso_sum.wait()

def compare_miso_runs(wt_samples, carm1_samples):
	#compare_miso --compare-samples SE/control/ SE/knockdown/ SE/comparisons/
	for wt_sample in wt_samples:
		carm1_sample = carm1_samples[wt_samples.index(wt_sample)]
		out_dir = wt_sample + '_' + carm1_sample + '_miso'
		miso_compare = subprocess.Popen(['compare_miso', '--compare-samples', wt_sample + '_miso_output', carm1_sample + '_miso_output', out_dir])
		miso_compare.wait()




##run methods
##all bams
# bams = ['cko13581.Aligned.sortedByCoord.out.bam', 'cko13584.Aligned.sortedByCoord.out.bam', 'cko1374.Aligned.sortedByCoord.out.bam', 'cko1378.Aligned.sortedByCoord.out.bam', 
# 		'cko14102.Aligned.sortedByCoord.out.bam', 'cko14117.Aligned.sortedByCoord.out.bam', 'cko14271.Aligned.sortedByCoord.out.bam', 'cko14272.Aligned.sortedByCoord.out.bam', 
# 		'cko14273.Aligned.sortedByCoord.out.bam', 'cko14275.Aligned.sortedByCoord.out.bam', 'cko14277.Aligned.sortedByCoord.out.bam', 'cko14278.Aligned.sortedByCoord.out.bam', 
# 		'cko15511.Aligned.sortedByCoord.out.bam', 'cko15512.Aligned.sortedByCoord.out.bam', 'cko15514.Aligned.sortedByCoord.out.bam', 'cko15517.Aligned.sortedByCoord.out.bam', 
# 		'cko15518.Aligned.sortedByCoord.out.bam', 'cko15561.Aligned.sortedByCoord.out.bam']
##just e12
bams = ['cko14271.Aligned.sortedByCoord.out.bam', 'cko14272.Aligned.sortedByCoord.out.bam', 'cko14273.Aligned.sortedByCoord.out.bam', 'cko14275.Aligned.sortedByCoord.out.bam', 
	'cko14277.Aligned.sortedByCoord.out.bam', 'cko14278.Aligned.sortedByCoord.out.bam']
# bams = ['cko15561.Aligned.sortedByCoord.out.bam']


e12_sample_dict = {'cko14271': ['289', '58'], 'cko14272': ['291', '59'], 'cko14273': ['291', '59'], 
		'cko14275': ['287', '59'], 'cko14277': ['293', '60'], 'cko14278': ['290', '60']}
e12_carm1_samples = ['cko14271', 'cko14273','cko14275']
e12_wt_samples = ['cko14272', 'cko14277','cko14278']
##get all constitutive exons from a gene models GF
# get_const_exons(mm10_gff, const_exons)

##get stats for all samples
# get_insert_stats(bams, insert_sizes, const_exon_mm10)

##index gff file
# index_gff_file(mm10_gff, mm10_gff_index)

##run miso on all e12 samples
run_miso(e12_sample_dict, fq_read_length, mm10_gff_index)

##compare samples
compare_miso_runs(e12_wt_samples, e12_carm1_samples)



