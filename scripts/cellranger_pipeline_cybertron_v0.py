#!/usr/bin/env python
import subprocess
import os

##parameters
delim = '\t'
# cellranger = '/home/atimms/programs/cellranger-3.1.0/cellranger'
# cellranger = '/home/atimms/programs/cellranger-4.0.0/cellranger'
# cellranger = '/home/atimms/programs/cellranger-5.0.0/cellranger' ##has include-introns 
# cellranger = '/home/atimms/programs/cellranger-5.0.1/cellranger' ##has include-introns 
cellranger = '/home/atimms/programs/cellranger-6.1.1/cellranger' 
cellranger_atac = '/home/atimms/programs/cellranger-atac-1.2.0/cellranger-atac'
cellranger_atac2 = '/home/atimms/programs/cellranger-atac-2.0.0/cellranger-atac'
# cellranger_arc = '/home/atimms/programs/cellranger-arc-1.0.1/cellranger-arc'
cellranger_arc = '/home/atimms/programs/cellranger-arc-2.0.0/cellranger-arc'
grc38_prerna_ref = '/home/atimms/ngs_data/references/10x/GRCh38-3.0.0_premrna_new'
grc38_ref = '/home/atimms/ngs_data/references/10x/refdata-gex-GRCh38-2020-A'
mm10_prerna_ref = '/gpfs/home/atimms/ngs_data/references/10x/mm10-3.0.0_premrna'
grc38_atac_ref = '/home/atimms/ngs_data/references/10x/refdata-cellranger-atac-GRCh38-1.2.0'
mm10_ref = '/home/atimms/ngs_data/references/10x/refdata-gex-mm10-2020-A'
# grc38_arc_ref = '/home/atimms/ngs_data/references/10x/refdata-cellranger-arc-GRCh38-2020-A'
grc38_arc_ref = '/home/atimms/ngs_data/references/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0'
mm10_arc_ref = '/home/atimms/ngs_data/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0'
GRCz11_arc_ref = '/home/atimms/ngs_data/references/10x/GRCz11'

##methods
def make_dict_from_info_file(in_file):
	# print(in_file)
	idict = {}
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				# print(line)
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

def run_cellranger_count_5prime(infodict, refdir):
	for sample in infodict:
		#cellranger count --id=Mut2-6 --transcriptome=/gpfs/home/atimms/ngs_data/references/10x/refdata-cellranger-GRCh38-3.0.0 --fastqs=/home/atimms/ngs_data/misc/cherry_10x_organoid_0120/MacTel_organoid_wk5_rerun_done/outs/fastq_path/HHT75BGXC/Mut2-6,/home/atimms/ngs_data/misc/cherry_10x_organoid_0120/MacTel_Organoid_5wk_done/outs/fastq_path/HFLF2BGXC/Mut2-6 --sample=Mut2-6 --expect-cells=1000
		fqs = infodict[sample][0]
		ec = infodict[sample][1]
		print(cellranger, 'count', '--id=', sample, '--transcriptome=', refdir, '--fastqs=', fqs, '--sample=', sample, '--expect-cells=', ec)
		cr_count = subprocess.Popen([cellranger, 'count', '--id=' + sample, '--transcriptome=' + refdir, '--fastqs=' + fqs, '--sample=' + sample, '--expect-cells=' + ec, '--chemistry=fiveprime' ])
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

def cellranger_scrnaseq_5prime_master(work_dir, infile, suffix_for_aggr, ref_dir):
	os.chdir(work_dir)
	info_dict = make_dict_from_info_file(infile)
	for s in info_dict:
		print(s, info_dict[s])
	run_cellranger_count_5prime(info_dict, ref_dir)
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

def cellranger_scatac_no_aggr_master(work_dir, infile, ref_dir):
	os.chdir(work_dir)
	info_dict = make_dict_from_info_file(infile)
	for s in info_dict:
		print(s, info_dict[s])
	run_cellranger_atac_count(info_dict, ref_dir)


def run_cellranger5_count(infodict, refdir, use_intronic_reads):
	for sample in infodict:
		#cellranger count --id=Mut2-6 --transcriptome=/gpfs/home/atimms/ngs_data/references/10x/refdata-cellranger-GRCh38-3.0.0 --fastqs=/home/atimms/ngs_data/misc/cherry_10x_organoid_0120/MacTel_organoid_wk5_rerun_done/outs/fastq_path/HHT75BGXC/Mut2-6,/home/atimms/ngs_data/misc/cherry_10x_organoid_0120/MacTel_Organoid_5wk_done/outs/fastq_path/HFLF2BGXC/Mut2-6 --sample=Mut2-6 --expect-cells=1000
		fqs = infodict[sample][0]
		print(cellranger, 'count', '--id=', sample, '--transcriptome=', refdir, '--fastqs=', fqs)
		##run using just exonic and both exonic and intronic reads
		if use_intronic_reads == 'no':
			cr_count1 = subprocess.Popen([cellranger, 'count', '--id=' + sample, '--transcriptome=' + refdir, '--fastqs=' + fqs, '--sample=' + sample])
			cr_count1.wait()
		elif use_intronic_reads == 'yes':
			cr_count2 = subprocess.Popen([cellranger, 'count', '--id=' + sample, '--transcriptome=' + refdir, '--fastqs=' + fqs, '--sample=' + sample, '--include-introns'])
			cr_count2.wait()
		else:
			print('issues....')

def run_cellranger5_count_force_cells(infodict, refdir, use_intronic_reads, cell_number_wanted):
	for sample in infodict:
		#cellranger count --id=Mut2-6 --transcriptome=/gpfs/home/atimms/ngs_data/references/10x/refdata-cellranger-GRCh38-3.0.0 --fastqs=/home/atimms/ngs_data/misc/cherry_10x_organoid_0120/MacTel_organoid_wk5_rerun_done/outs/fastq_path/HHT75BGXC/Mut2-6,/home/atimms/ngs_data/misc/cherry_10x_organoid_0120/MacTel_Organoid_5wk_done/outs/fastq_path/HFLF2BGXC/Mut2-6 --sample=Mut2-6 --expect-cells=1000
		fqs = infodict[sample][0]
		print(cellranger, 'count', '--id=', sample, '--transcriptome=', refdir, '--fastqs=', fqs)
		##run using just exonic and both exonic and intronic reads
		if use_intronic_reads == 'no':
			cr_count1 = subprocess.Popen([cellranger, 'count', '--id=' + sample, '--transcriptome=' + refdir, '--fastqs=' + fqs, '--sample=' + sample, '--force-cells', cell_number_wanted])
			cr_count1.wait()
		elif use_intronic_reads == 'yes':
			cr_count2 = subprocess.Popen([cellranger, 'count', '--id=' + sample, '--transcriptome=' + refdir, '--fastqs=' + fqs, '--sample=' + sample, '--include-introns', '--force-cells', cell_number_wanted])
			cr_count2.wait()
		else:
			print('issues....')

def cellranger5_scrnaseq_master_no_aggr(work_dir, infile, ref_dir, use_intronic_reads):
	os.chdir(work_dir)
	info_dict = make_dict_from_info_file(infile)
	for s in info_dict:
		print(s, info_dict[s])
	run_cellranger5_count(info_dict, ref_dir, use_intronic_reads)

def cellranger5_scrnaseq_master_force_cells(work_dir, infile, ref_dir, use_intronic_reads, cells_wanted):
	os.chdir(work_dir)
	info_dict = make_dict_from_info_file(infile)
	for s in info_dict:
		print(s, info_dict[s])
	run_cellranger5_count_force_cells(info_dict, ref_dir, use_intronic_reads, cells_wanted)

def make_arc_library_file_from_dict(sample_name, fq_dicts, out_file):
	# print(in_file)
	with open(out_file, "w") as out_fh:
		out_fh.write(','.join(['fastqs','sample','library_type']) + '\n')
		# print(line)
		ge_fqs = fq_dicts[0]
		ca_fqs = fq_dicts[1]
		out_fh.write(','.join([ge_fqs,sample_name,'Gene Expression']) + '\n')
		out_fh.write(','.join([ca_fqs,sample_name,'Chromatin Accessibility']) + '\n')

def make_arc_dict_from_info_file(in_file):
	# print(in_file)
	i_dict = {}
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 1:
				# print(line)
				line = line.rstrip().split(delim)
				sample = line[0]
				ge_fqs = line[1]
				ca_fqs = line[2]
				i_dict[sample] = [ge_fqs, ca_fqs]
	return(i_dict)

def run_cellranger_arc_count(sample_dict, refdir, lib_file_suffix):
	for sample in sample_dict:
		fq_list = sample_dict[sample]
		lib_file = sample + lib_file_suffix
		make_arc_library_file_from_dict(sample, fq_list, lib_file)
		##cellranger-arc count --id=sample345 --reference=/opt/refdata-cellranger-arc-GRCh38-2020-A --libraries=/home/jdoe/runs/libraries.csv
		cr_count = subprocess.Popen([cellranger_arc, 'count', '--id=' + sample, '--reference=' + refdir, '--libraries=' + lib_file])
		cr_count.wait()

def cellranger_arc_master(work_dir, infile, ref_dir):
	os.chdir(work_dir)
	library_file_suffix = '.library.csv'
	sample_dict = make_arc_dict_from_info_file(infile)
	print(sample_dict)
	run_cellranger_arc_count(sample_dict, ref_dir, library_file_suffix)

def run_cellranger_arc_count_modified(sample_dict, refdir, lib_file_suffix, gex_umi, atac_umi):
	for sample in sample_dict:
		fq_list = sample_dict[sample]
		lib_file = sample + lib_file_suffix
		make_arc_library_file_from_dict(sample, fq_list, lib_file)
		##cellranger-arc count --id=sample345 --reference=/opt/refdata-cellranger-arc-GRCh38-2020-A --libraries=/home/jdoe/runs/libraries.csv
		cr_count = subprocess.Popen([cellranger_arc, 'count', '--id=' + sample, '--reference=' + refdir, '--libraries=' + lib_file, '--min-atac-count=' + atac_umi, '--min-gex-count=' + gex_umi])
		cr_count.wait()

def cellranger_arc_master_modified(work_dir, infile, ref_dir, gex_umi, atac_umi):
	os.chdir(work_dir)
	library_file_suffix = '.library.csv'
	sample_dict = make_arc_dict_from_info_file(infile)
	print(sample_dict)
	run_cellranger_arc_count_modified(sample_dict, ref_dir, library_file_suffix, gex_umi, atac_umi)

def run_cellranger_atac2_count(sample, fq_dir, refdir):
	cr_count = subprocess.Popen([cellranger_atac2, 'count', '--id=' + sample, '--fastqs=' + fq_dir, '--reference=' + refdir])
	cr_count.wait()


def cellranger_atac2_master(work_dir, in_file):
	os.chdir(work_dir)
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			sample = line[0]
			fqs = line[1]
			ref = line[2]
			run_cellranger_atac2_count(sample,fqs,ref)

def run_cellranger_count_v1(sample, fqs, refdir, use_intronic_reads):
	if use_intronic_reads == 'no':
		cr_count1 = subprocess.Popen([cellranger, 'count', '--id=' + sample, '--transcriptome=' + refdir, '--fastqs=' + fqs, '--sample=' + sample])
		cr_count1.wait()
	elif use_intronic_reads == 'yes':
		cr_count2 = subprocess.Popen([cellranger, 'count', '--id=' + sample, '--transcriptome=' + refdir, '--fastqs=' + fqs, '--sample=' + sample, '--include-introns'])
		cr_count2.wait()
	else:
		print('issues....')

def cellranger_scrnaseq_master_no_aggr_v1(work_dir, in_file):
	os.chdir(work_dir)
	with open(in_file, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			sample = line[0]
			fqs = line[1]
			ref = line[2]
			intronic_reads = line[3]
			run_cellranger_count_v1(sample,fqs,ref,intronic_reads)


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
# cellranger_scrnaseq_master_no_aggr(working_dir, info_file, transciptome_ref)
##even more data 12 wk orgaoid rna and 3 embryonic atac
working_dir = '/home/atimms/ngs_data/cellranger/cherry_extra_atac_rna_0620/scRNA_timcherry_done'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
info_file = '12wkorg_scRNAseq_0620.txt'
combined_suffix = '12wkorg_scRNAseq_0620'
transciptome_ref = grc38_prerna_ref
# cellranger_scrnaseq_master(working_dir, info_file, combined_suffix, transciptome_ref)

working_dir = '/home/atimms/ngs_data/cellranger/cherry_extra_atac_rna_0620/scATAC_timcherry_done'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group (expect cells not actually used)
info_file = 'extra_human_scATACseq_0620.txt'
transciptome_ref = grc38_atac_ref
# cellranger_scatac_no_aggr_master(working_dir, info_file, transciptome_ref)


##new data 0720 - Jenn Chao Tims collabarator
working_dir = '/home/atimms/ngs_data/cellranger/jchao_rna_0720'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
info_file = 'SFD_CRISPR_SFD_snRNAseq_0720.txt'
combined_suffix = 'SFD_CRISPR_SFD_snRNAseq_0720'
transciptome_ref = grc38_prerna_ref
# cellranger_scrnaseq_master(working_dir, info_file, combined_suffix, transciptome_ref)


##kim data 0920
working_dir = '/archive/millen_k/kims_data/kim_10Xv3_RNA_DSbrain'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group
info_file = 'kim_10Xv3_RNA_DSbrain_0920.txt'
combined_suffix = 'kim_10Xv3_RNA_DSbrain_0920'
transciptome_ref = grc38_ref
# cellranger_scrnaseq_master(working_dir, info_file, combined_suffix, transciptome_ref)


##dana/jimmy data 1220
##do 2x one with std hg38 and one with premrna
##issue with count, so had to specify chemistry type
working_dir = '/home/atimms/ngs_data/cellranger/dana_scrna_1220/hg38'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group (expect cells not actually used)
info_file = '201112_Bennett_scRNAseq.txt'
combined_suffix = '201112_Bennett_scRNAseq_hg38'
transciptome_ref = grc38_ref
# cellranger_scrnaseq_5prime_master(working_dir, info_file, combined_suffix, transciptome_ref)
working_dir = '/home/atimms/ngs_data/cellranger/dana_scrna_1220/hg38_premrna'
##info on the analysis, 4x columns with a header: sample fqs expected_cells group (expect cells not actually used)
info_file = '201112_Bennett_scRNAseq.txt'
combined_suffix = '201112_Bennett_scRNAseq_hg38_premrna'
transciptome_ref = grc38_prerna_ref
# cellranger_scrnaseq_5prime_master(working_dir, info_file, combined_suffix, transciptome_ref)

##data for michael/gus
working_dir = '/home/atimms/ngs_data/cellranger/cunn_mouse_scrna_0121/std_ref'
##info on the analysis, 4x columns with a header: sample fqs na na
info_file = 'cunn_scRNA_0121.txt'
transciptome_ref = mm10_ref
intronic_ref = 'no'
# cellranger5_scrnaseq_master_no_aggr(working_dir, info_file, transciptome_ref, intronic_ref)
working_dir = '/home/atimms/ngs_data/cellranger/cunn_mouse_scrna_0121/intronic_ref'
intronic_ref = 'yes'
# cellranger5_scrnaseq_master_no_aggr(working_dir, info_file, transciptome_ref, intronic_ref)

##data for tim/eric --- run atac and rnaseq both ways
working_dir = '/home/atimms/ngs_data/cellranger/cherry_10x_0121/'
##info on the analysis, 4x columns with a header: sample fq_ge fq_atac
info_file = '12wk_snporg_0121.txt'
transciptome_ref = grc38_arc_ref
# cellranger_arc_master(working_dir, info_file, transciptome_ref)

##data for tim/eric, rnaseq on mouse and human (organoid) data
working_dir = '/home/atimms/ngs_data/cellranger/cherry_scrnaseq_0221/mactel_mouse_scRNAseq'
info_file = 'mactel_mouse_scRNAseq_0221.txt'
transciptome_ref = mm10_ref
intronic_ref = 'yes'
# cellranger5_scrnaseq_master_no_aggr(working_dir, info_file, transciptome_ref, intronic_ref)
working_dir = '/home/atimms/ngs_data/cellranger/cherry_scrnaseq_0221/organoid_microglia_scRNAseq'
info_file = 'organoid_microglia_scRNAseq_0221.txt'
transciptome_ref = grc38_ref
intronic_ref = 'yes'
# cellranger5_scrnaseq_master_no_aggr(working_dir, info_file, transciptome_ref, intronic_ref)
##try using force cells to get extra data
working_dir = '/home/atimms/ngs_data/cellranger/cherry_scrnaseq_0221/mactel_mouse_scRNAseq_force_cells'
info_file = 'mactel_mouse_scRNAseq_0221.txt'
transciptome_ref = mm10_ref
intronic_ref = 'yes'
# cellranger5_scrnaseq_master_force_cells(working_dir, info_file, transciptome_ref, intronic_ref, '8000')


working_dir = '/home/atimms/ngs_data/cellranger/cherry_human_scATAC_periph_vs_macula_0521'
##info on the analysis, 3x columns with a header: sample/id fqs ref_file (atac now uses ARC ref)
info_file = 'peripheral_vs_macula_huret.txt'
# cellranger_atac2_master(working_dir, info_file)

##data for tim/eric --- run atac and rnaseq both ways
working_dir = '/home/atimms/ngs_data/cellranger/cherry_10x_multiome_0621/'
##info on the analysis, 4x columns with a header: sample fq_ge fq_atac
info_file = '20wk_snporg_0621.txt'
transciptome_ref = grc38_arc_ref
# cellranger_arc_master(working_dir, info_file, transciptome_ref)

##mouse multiome data for kim 
working_dir = '/home/atimms/ngs_data/cellranger/kim_mouse_multiome_0721/'
##info on the analysis, 4x columns with a header: sample fq_ge fq_atac
info_file = 'kim_mouse_multiome_0721.txt'
transciptome_ref = mm10_arc_ref
# cellranger_arc_master(working_dir, info_file, transciptome_ref)

##kim atac data 0821
working_dir = '/archive/millen_k/kims_data/kim_cbl_10X_atac_0821'
##info on the analysis, 3x columns with no header: sample/id fqs ref_file (atac now uses ARC ref)
info_file = 'kim_cbl_10X_atac.txt'
# cellranger_atac2_master(working_dir, info_file)

##dana atac data test 0821
working_dir = '/home/atimms/ngs_data/cellranger/dana_atac_10x_0821'
##info on the analysis, 3x columns with a header: sample/id fqs ref_file (atac now uses ARC ref)
info_file = 'dana_atac_10x_0821.txt'
# cellranger_atac2_master(working_dir, info_file)


##tim/eric scrna
working_dir = '/active/cherry_t/OrgManuscript_SingleCell_Data/org_snp_multiome/28wk/'
##4x columns: sample, fqs, transciptome_ref, and if to include intronic reads (no header)
info_file = '28wk_snporg_scRNAseq.txt'
transciptome_ref = grc38_ref
intronic_ref = 'yes'
# cellranger_scrnaseq_master_no_aggr_v1(working_dir, info_file)

##tim/eric atac data 0921
working_dir = '/active/cherry_t/OrgManuscript_SingleCell_Data/org_scATAC2/wk12_wt_vs_ko'
##info on the analysis, 3x columns with no header: sample/id fqs ref_file (atac now uses ARC ref)
info_file = 'wk12_wt_vs_ko_atac_0921.txt'
# cellranger_atac2_master(working_dir, info_file)

##dana atac data test 0821
working_dir = '/home/atimms/ngs_data/cellranger/dana_atac_10x_0921'
##info on the analysis, 3x columns with a header: sample/id fqs ref_file (atac now uses ARC ref)
info_file = 'dana_atac_10x_0921.txt'
# cellranger_atac2_master(working_dir, info_file)

##zebrafish multiome data for lisa/hank
working_dir = '/home/atimms/ngs_data/cellranger/lisa_zf_multiome_1021/'
##info on the analysis, 4x columns with a header: sample fq_ge fq_atac
info_file = 'lisa_zf_multiome_1021.txt'
transciptome_ref = GRCz11_arc_ref
min_gex_umi = '500'
min_atac_umi = '1000'
# cellranger_arc_master(working_dir, info_file, transciptome_ref)
# cellranger_arc_master_modified(working_dir, info_file, transciptome_ref, min_gex_umi, min_atac_umi)

##dana atac data test 0821
working_dir = '/archive/bennett_j/210930_Bennett_ATAC/'
##info on the analysis, 3x columns without a header: sample/id fqs ref_file (atac now uses ARC ref)
info_file = 'dana_atac_10x_1021.txt'
# cellranger_atac2_master(working_dir, info_file)

##4 human multiome samples for kim
working_dir = '/archive/millen_k/kims_data/kim_10x_multi_0921/'
##info on the analysis, 4x columns with a header: sample fq_ge fq_atac
info_file = 'kim_10x_multiome_0921.txt'
transciptome_ref = grc38_arc_ref
# cellranger_arc_master(working_dir, info_file, transciptome_ref)

##4 human multiome samples for kim
working_dir = '/archive/millen_k/kims_data/kim_10x_multi_1021/'
##info on the analysis, 4x columns with a header: sample fq_ge fq_atac
info_file = 'kim_10x_multiome_1021.txt'
transciptome_ref = grc38_arc_ref
# cellranger_arc_master(working_dir, info_file, transciptome_ref)

##dana atac data from miseq 1221
working_dir = '/home/atimms/ngs_data/cellranger/dana_atac_10x_1221'
##info on the analysis, 3x columns without a header: sample/id fqs ref_file (atac now uses ARC ref)
info_file = 'dana_mouse_kidney_atac_10x_1221.txt'
# cellranger_atac2_master(working_dir, info_file)
info_file = 'dana_avm_10x_atac_1221.txt'
# cellranger_atac2_master(working_dir, info_file)

##extra sequenceing for previous human multiome samples for kim
working_dir = '/archive/millen_k/kims_data/kim_10x_multi_1021_embredo'
##info on the analysis, 4x columns with a header: sample fq_ge fq_atac
info_file = 'kim_10x_multi_1021_embredo.txt'
transciptome_ref = grc38_arc_ref
cellranger_arc_master(working_dir, info_file, transciptome_ref)

