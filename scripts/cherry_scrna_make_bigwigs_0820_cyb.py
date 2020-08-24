#!/usr/bin/env python
import subprocess
import os
import glob

'''
info....

##for sinto part i.e. getting bam files 
##go to worker node and load biobuilds and python3 and open env 
qsub -Iq cdbrmq -l mem=200gb,ncpus=20,walltime=10:00:00 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load biobuilds
module load local_python/3.7.6
source activate sinto
##program is now at:
/home/atimms/programs/sinto/scripts/sinto



##for getting stranded use bamCoverage from deepTools:
##to setup:
module load local_python/3.7.6
conda create --name deeptools
source activate deeptools
conda install -c bioconda deeptools

##from now on:
qsub -Iq cdbrmq -l mem=200gb,ncpus=20,walltime=10:00:00 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load local_python/3.7.6
source activate deeptools
module load biobuilds
'''

##parameters
delim = '\t'
thread_number = '20'
##programs
sinto = '/home/atimms/programs/sinto/scripts/sinto'
bedtools = '/home/atimms/programs/bedtools2.28/bin/bedtools'
bg_to_bw = '/home/atimms/programs/bedGraphToBigWig'

##methods
def run_sinto_method(sample_name, bam_file, cells_file, work_dir):
	directory = './' + sample_name
	if not os.path.exists(sample_name):
		os.makedirs(sample_name)
	os.chdir(directory)
	run_sinto = subprocess.Popen([sinto, 'filterbarcodes', '-b', bam_file,  '-c', cells_file, '-p', thread_number])
	run_sinto.wait()
	os.chdir(work_dir)

def make_bigwig_files(samples, tissues):
	for tissue in tissues:
		for sample in samples:
			bam = sample + '/' + tissue + '.bam'
			bedgraph1 = sample + '_' + tissue + '.temp1.bedgraph'
			bedgraph2 = sample + '_' + tissue + '.temp2.bedgraph'
			out_for_bw = sample + '_' + tissue  + '.fwd.bigwig'
			out_rev_bw = sample + '_' + tissue  + '.rev.bigwig'
			if os.path.exists(bam):
				##index file
				samtools_index = subprocess.Popen(['samtools', 'index', bam])
				samtools_index.wait()
				##get bigwigs
				bc_f = subprocess.Popen(['bamCoverage', '-b', bam, '-o', out_for_bw, '--filterRNAstrand', 'forward', '-p', '18', '--normalizeUsing', 'RPKM'])
				bc_f.wait()
				bc_r = subprocess.Popen(['bamCoverage', '-b', bam, '-o', out_rev_bw, '--filterRNAstrand', 'reverse', '-p', '18', '--normalizeUsing', 'RPKM'])
				bc_r.wait()
			else:
				print("bam file doesn't exist", bam)



##run methods

##human adult
working_dir = '/home/atimms/ngs_data/misc/cherry_scRNA_make_bigwigs_0820/human_adult_0820'
os.chdir(working_dir)
##info
##bam and barcode file
hu37_cell_info = working_dir + '/' + 'human_adult_scrnaseq.hu37.cell_types.txt'
hu5_cell_info = working_dir + '/' + 'human_adult_scrnaseq.hu5.cell_types.txt'
hu7_cell_info = working_dir + '/' + 'human_adult_scrnaseq.hu7.cell_types.txt'
hu37_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/H37/outs/possorted_genome_bam.bam'
hu5_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/Hu5/outs/possorted_genome_bam.bam'
hu7_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/Hu7/outs/possorted_genome_bam.bam'

##samples and tissue
sample_dict = {'hu37':[hu37_bam, hu37_cell_info], 'hu5':[hu5_bam, hu5_cell_info], 'hu7':[hu7_bam, hu7_cell_info]}
sample_names = ['hu37', 'hu5', 'hu7']
cell_type_names = ['Amacrines', 'Astrocytes', 'Bipolars', 'Cones', 'Ganglions', 'Horizontals', 'Microglia', 'Mulllers', 'Rods', 'Rods_with_MT_reads']

##make bams per tissue type
##need python3 and sinto conda env
# for sample in sample_dict:
# 	run_sinto_method(sample, sample_dict[sample][0], sample_dict[sample][1], working_dir)
##make bigwigs
# make_bigwig_files(sample_names, cell_type_names)


##human adult
working_dir = '/home/atimms/ngs_data/misc/cherry_scRNA_make_bigwigs_0820/human_all_0820'
os.chdir(working_dir)
##info
##bam and barcode file
d113_cell_info = working_dir + '/' + 'human_all_scrnaseq.d113.cell_types.txt'
d132_cell_info = working_dir + '/' + 'human_all_scrnaseq.d132.cell_types.txt'
d53_cell_info = working_dir + '/' + 'human_all_scrnaseq.d53.cell_types.txt'
d59_cell_info = working_dir + '/' + 'human_all_scrnaseq.d59.cell_types.txt'
d74_cell_info = working_dir + '/' + 'human_all_scrnaseq.d74.cell_types.txt'
d78_cell_info = working_dir + '/' + 'human_all_scrnaseq.d78.cell_types.txt'
hu37_cell_info = working_dir + '/' + 'human_all_scrnaseq.hu37.cell_types.txt'
hu5_cell_info = working_dir + '/' + 'human_all_scrnaseq.hu5.cell_types.txt'
hu7_cell_info = working_dir + '/' + 'human_all_scrnaseq.hu7.cell_types.txt'
d113_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/d113/outs/possorted_genome_bam.bam'
d132_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/d132/outs/possorted_genome_bam.bam'
d53_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/d53/outs/possorted_genome_bam.bam'
d59_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/d59/outs/possorted_genome_bam.bam'
d74_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/d74/outs/possorted_genome_bam.bam'
d78_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/d78/outs/possorted_genome_bam.bam'
hu37_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/H37/outs/possorted_genome_bam.bam'
hu5_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/Hu5/outs/possorted_genome_bam.bam'
hu7_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scRNA/Hu7/outs/possorted_genome_bam.bam'

##samples and tissue
sample_dict = {'d113':[d113_bam, d113_cell_info], 'd132':[d132_bam, d132_cell_info], 'd53':[d53_bam, d53_cell_info],
		'd59':[d59_bam, d59_cell_info], 'd74':[d74_bam, d74_cell_info], 'd78':[d78_bam, d78_cell_info],
		'hu37':[hu37_bam, hu37_cell_info], 'hu5':[hu5_bam, hu5_cell_info], 'hu7':[hu7_bam, hu7_cell_info]}
sample_dict = {'hu5':[hu5_bam, hu5_cell_info], 'hu7':[hu7_bam, hu7_cell_info]}
sample_names = sample_dict.keys()
cell_type_names = ['Amacrines', 'Astrocytes', 'bipolar_photoreceptor_cell_precursors', 'Bipolars', 'Cones', 'Cycling_progenitors', 
		'Early_Primary_Progenitors', 'Ganglions', 'high_mito_rods', 'Horizontals', 'Microglia', 'Mullers', 'Rods', 'Transitional_1']


##make bams per tissue type
##need python3 and sinto conda env
# for sample in sample_dict:
# 	run_sinto_method(sample, sample_dict[sample][0], sample_dict[sample][1], working_dir)
##make bigwigs
make_bigwig_files(sample_names, cell_type_names)





