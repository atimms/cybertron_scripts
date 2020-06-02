#!/usr/bin/env python
import subprocess
import os

##parameters
delim = '\t'
cellranger = '/home/atimms/programs/cellranger-3.1.0/cellranger'
cellranger_atac = '/home/atimms/programs/cellranger-atac-1.2.0/cellranger-atac'
grc38_prerna_ref = '/home/atimms/ngs_data/references/10x/GRCh38-3.0.0_premrna_new'
mm10_prerna_ref = '/gpfs/home/atimms/ngs_data/references/10x/mm10-3.0.0_premrna'
grc38_atac_ref = '/home/atimms/ngs_data/references/10x/refdata-cellranger-atac-GRCh38-1.2.0'


##methods

def make_dict_from_info_file(in_file):
	idict = {}
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				line = line.rstrip().split(delim)
				sample = line[0]
				fqs = line[1]
				expected_cells = line[2]
				group = line[3]
				idict[sample] = [fqs, expected_cells, group]
	return(idict)

def run_cellranger_count(infodict, refdir):
	for sample in infodict:
		#cellranger count --id=Mut2-6 --transcriptome=/gpfs/home/atimms/ngs_data/references/10x/refdata-cellranger-GRCh38-3.0.0 --fastqs=/home/atimms/ngs_data/misc/cherry_10x_organoid_0120/MacTel_organoid_wk5_rerun_done/outs/fastq_path/HHT75BGXC/Mut2-6,/home/atimms/ngs_data/misc/cherry_10x_organoid_0120/MacTel_Organoid_5wk_done/outs/fastq_path/HFLF2BGXC/Mut2-6 --sample=Mut2-6 --expect-cells=1000
		fqs = infodict[sample][0]
		ec = infodict[sample][1]
		print(cellranger, 'count', '--id=', sample, '--transcriptome=', refdir, '--fastqs=', fqs, '--sample=', sample, '--expect-cells=', ec)
		cr_count = subprocess.Popen([cellranger, 'count', '--id=' + sample, '--transcriptome=' + refdir, '--fastqs=' + fqs, '--sample=' + sample, '--expect-cells=' + ec])
		cr_count.wait()

def make_aggr_csv(infodict, out_suffix, w_dir):
	outfile = out_suffix + '.csv'
	with open(outfile, "w") as out_fh:
		out_fh.write(','.join(['library_id','molecule_h5','group']) + '\n')
		for sample in infodict:
			mol_info = w_dir + '/' + sample + '/outs/molecule_info.h5'
			group = infodict[sample][2]
			out_fh.write(','.join([sample,mol_info,group]) + '\n')

def run_cellranger_aggr(comb_name):
	#/home/atimms/programs/cellranger-3.1.0/cellranger aggr --id=organoid_5wk_combined_0120 --csv=organoid_5wk_combined_0120.csv --normalize=mapped
	cr_aggr = subprocess.Popen([cellranger, 'aggr', '--id=' + comb_name, '--csv=' + comb_name + '.csv', '--normalize=mapped'])
	cr_aggr.wait()

def cellranger_scrnaseq_master(work_dir, infile, suffix_for_aggr, ref_dir):
	os.chdir(work_dir)
	info_dict = make_dict_from_info_file(infile)
	for s in info_dict:
		print(s, info_dict[s])
	run_cellranger_count(info_dict, ref_dir)
	make_aggr_csv(info_dict, suffix_for_aggr, work_dir)
	run_cellranger_aggr(suffix_for_aggr)

def cellranger_scrnaseq_master_no_aggr(work_dir, infile, ref_dir):
	os.chdir(work_dir)
	info_dict = make_dict_from_info_file(infile)
	for s in info_dict:
		print(s, info_dict[s])
	run_cellranger_count(info_dict, ref_dir)


def run_cellranger_atac_count(infodict, refdir):
	for sample in infodict:
		fqs = infodict[sample][0]
		##/home/atimms/programs/cellranger-atac-1.2.0/cellranger-atac count --id=Hu8 --reference=/home/atimms/ngs_data/references/10x/refdata-cellranger-atac-GRCh38-1.2.0 --fastqs=/home/atimms/ngs_data/misc/cherry_10x_atac_1119/organoid_timecourse_scATAC_rd2_cherry_done/outs/fastq_path/HHCLGDRXX/Hu8,/home/atimms/ngs_data/misc/cherry_10x_atac_1119/organoid_timecourse_scATAC_rd2_run2_done/outs/fastq_path/HHF2MDRXX/Hu8 --sample=Hu8
		# print(cellranger_atac, 'count', '--id=', sample, '--transcriptome=', refdir, '--fastqs=', fqs, '--sample=', sample)
		cr_count = subprocess.Popen([cellranger_atac, 'count', '--id=' + sample, '--reference=' + refdir, '--fastqs=' + fqs, '--sample=' + sample])
		cr_count.wait()

def make_aggr_atac_csv(infodict, out_suffix, w_dir):
	outfile = out_suffix + '.csv'
	with open(outfile, "w") as out_fh:
		# out_fh.write(','.join(['library_id','molecule_h5','group']) + '\n')
		out_fh.write(','.join(['library_id','fragments','cells','group']) + '\n')
		for sample in infodict:
			fragments = w_dir + '/' + sample + '/outs/fragments.tsv.gz'
			cells = w_dir + '/' + sample + '/outs/singlecell.csv'
			group = infodict[sample][2]
			out_fh.write(','.join([sample,fragments,cells,group]) + '\n')

def run_cellranger_atac_aggr(comb_name, refdir):
	#/home/atimms/programs/cellranger-atac-1.2.0/cellranger-atac aggr --id=organoid_atac_combined_cr12_1119 --csv=organoid_atac_combined_1119.csv --reference=/home/atimms/ngs_data/references/10x/refdata-cellranger-atac-GRCh38-1.2.0
	cr_aggr = subprocess.Popen([cellranger_atac, 'aggr', '--id=' + comb_name, '--csv=' + comb_name + '.csv', '--reference=' + refdir])
	cr_aggr.wait()

def cellranger_scatac_master(work_dir, infile, suffix_for_aggr, ref_dir):
	os.chdir(work_dir)
	info_dict = make_dict_from_info_file(infile)
	for s in info_dict:
		print(s, info_dict[s])
	run_cellranger_atac_count(info_dict, ref_dir)
	make_aggr_atac_csv(info_dict, suffix_for_aggr, work_dir)
	run_cellranger_atac_aggr(suffix_for_aggr, ref_dir)

def cellranger_scatac_just_aggr(work_dir, infile, suffix_for_aggr, ref_dir):
	os.chdir(work_dir)
	info_dict = make_dict_from_info_file(infile)
	for s in info_dict:
		print(s, info_dict[s])
	make_aggr_atac_csv(info_dict, suffix_for_aggr, work_dir)
	run_cellranger_atac_aggr(suffix_for_aggr, ref_dir)

##run methods

##cherry organoid analysis 0320
working_dir = '/home/atimms/ngs_data/cellranger/cherry_10x_organoid_0120'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
info_file = 'cherry_10x_organoid_5wk_0320.txt'
combined_suffix = 'organoid_5wk_combined_0320'
transciptome_ref = grc38_prerna_ref
# cellranger_scrnaseq_master(working_dir, info_file, combined_suffix, transciptome_ref)

##cherry organoid analysis 0320 with 1 extra sample - new samples bad, so don't use
working_dir = '/home/atimms/ngs_data/cellranger/cherry_10x_organoid_0320'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
##for the cell ranger just the new samples
# info_file = 'cherry_10x_organoid_5wk_extras_0320.txt'
##for aggr so combining 6 samples
info_file = 'cherry_10x_organoid_5wk_6samples_0320.txt'
combined_suffix = 'organoid_5wk_6samples_0320'
transciptome_ref = grc38_prerna_ref
# cellranger_scrnaseq_master(working_dir, info_file, combined_suffix, transciptome_ref)


##make premrna ref for mm10 - https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#header
##acomy scRNASeq
working_dir = '/home/atimms/ngs_data/cellranger/acomy_10x_0320'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
info_file = 'acomy_10x_0320.txt'
combined_suffix = 'acomy_10x_aggr_0320'
transciptome_ref = mm10_prerna_ref
# cellranger_scrnaseq_master(working_dir, info_file, combined_suffix, transciptome_ref)

##mouse kidney scRNASeq for Dave
working_dir = '/home/atimms/ngs_data/cellranger/beier_kidney_10x_0320'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
info_file = 'beier_kidney_10x_0320.txt'
combined_suffix = 'kidney_10x_aggr_0320'
transciptome_ref = mm10_prerna_ref
# cellranger_scrnaseq_master(working_dir, info_file, combined_suffix, transciptome_ref)

##cherry 3 sets of data 0320 - 2x rnaseq 1x scATAC *******
working_dir = '/home/atimms/ngs_data/cellranger/cherry_human_retina_scRNAseq_0320'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
info_file = 'human_retina_scRNAseq_0320.txt'
combined_suffix = 'human_retina_scRNAseq_combined_0320'
transciptome_ref = grc38_prerna_ref
# cellranger_scrnaseq_master(working_dir, info_file, combined_suffix, transciptome_ref)

working_dir = '/home/atimms/ngs_data/cellranger/cherry_MacTel_20wkorg_scRNAseq_0320'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
info_file = '20wkorg_scRNAseq_0320.txt'
combined_suffix = '20wkorg_scRNAseq_combined_0320'
transciptome_ref = grc38_prerna_ref
# cellranger_scrnaseq_master(working_dir, info_file, combined_suffix, transciptome_ref)
##using the new data and also re analyzing the adult human data
working_dir = '/home/atimms/ngs_data/cellranger/cherry_embryonic_retina_scATACseq_0320'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group (expect cells not actually used)
info_file = 'human_retina_scATACseq_0320.txt'
combined_suffix = 'human_retina_scATACseq_combined_0320'
transciptome_ref = grc38_atac_ref
# cellranger_scatac_master(working_dir, info_file, combined_suffix, transciptome_ref)
##issue with mwthods just run aggr
# cellranger_scatac_just_aggr(working_dir, info_file, combined_suffix, transciptome_ref)
##redo an analysis looking at all the organoid scATAC
working_dir = '/home/atimms/ngs_data/cellranger/cherry_org_scATAC_reanalysis_0320'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group (expect cells not actually used)
info_file = 'organoid_scATACseq_0320.txt'
combined_suffix = 'organoid_scATACseq_combined_0320'
transciptome_ref = grc38_atac_ref
# cellranger_scatac_master(working_dir, info_file, combined_suffix, transciptome_ref)
##issue with mwthods just run aggr
# cellranger_scatac_just_aggr(working_dir, info_file, combined_suffix, transciptome_ref)

##new data 0620
working_dir = '/home/atimms/ngs_data/cellranger/cherry_scRNASeq_0620'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
info_file = '28wkorg_scRNAseq_0620.txt'
combined_suffix = 'organoid_28wk_scATACseq_0620'
transciptome_ref = grc38_prerna_ref
# cellranger_scrnaseq_master(working_dir, info_file, combined_suffix, transciptome_ref)

info_file = 'extra_human_scRNAseq_0620.txt'
transciptome_ref = grc38_prerna_ref
cellranger_scrnaseq_master_no_aggr(working_dir, info_file, transciptome_ref)

