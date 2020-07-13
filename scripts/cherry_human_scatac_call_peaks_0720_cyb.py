#!/usr/bin/env python
import subprocess
import os

'''
info....

##for sinto part i.e. getting bam files 
##go to worker node and load biobuilds and python3 and open env 
qsub -Iq longq -l mem=200gb,ncpus=20
module load biobuilds
module load local_python/3.7.6
source activate sinto
##program is now at:
/home/atimms/programs/sinto/scripts/sinto

'''

##parameters
delim = '\t'
thread_number = '20'

sinto = '/home/atimms/programs/sinto/scripts/sinto'
genrich = '/home/atimms/programs/Genrich/Genrich'



##methods
def run_sinto_method(sample_name, bam_file, cells_file):
	directory = './' + sample_name
	if not os.path.exists(sample_name):
		os.makedirs(sample_name)
	os.chdir(directory)
	run_sinto = subprocess.Popen([sinto, 'filterbarcodes', '-b', bam_file,  '-c', cells_file, '-p', thread_number])
	run_sinto.wait()

def call_peaks_on_bams(samples, tissues):
	chrs_to_exclude = 'chrY,chrM,chr1_KI270706v1_random,chr1_KI270707v1_random,chr1_KI270708v1_random,chr1_KI270709v1_random,chr1_KI270710v1_random,chr1_KI270711v1_random,chr1_KI270712v1_random,chr1_KI270713v1_random,chr1_KI270714v1_random,chr2_KI270715v1_random,chr2_KI270716v1_random,chr3_GL000221v1_random,chr4_GL000008v2_random,chr5_GL000208v1_random,chr9_KI270717v1_random,chr9_KI270718v1_random,chr9_KI270719v1_random,chr9_KI270720v1_random,chr11_KI270721v1_random,chr14_GL000009v2_random,chr14_GL000225v1_random,chr14_KI270722v1_random,chr14_GL000194v1_random,chr14_KI270723v1_random,chr14_KI270724v1_random,chr14_KI270725v1_random,chr14_KI270726v1_random,chr15_KI270727v1_random,chr16_KI270728v1_random,chr17_GL000205v2_random,chr17_KI270729v1_random,chr17_KI270730v1_random,chr22_KI270731v1_random,chr22_KI270732v1_random,chr22_KI270733v1_random,chr22_KI270734v1_random,chr22_KI270735v1_random,chr22_KI270736v1_random,chr22_KI270737v1_random,chr22_KI270738v1_random,chr22_KI270739v1_random,chrY_KI270740v1_random,chrUn_KI270302v1,chrUn_KI270304v1,chrUn_KI270303v1,chrUn_KI270305v1,chrUn_KI270322v1,chrUn_KI270320v1,chrUn_KI270310v1,chrUn_KI270316v1,chrUn_KI270315v1,chrUn_KI270312v1,chrUn_KI270311v1,chrUn_KI270317v1,chrUn_KI270412v1,chrUn_KI270411v1,chrUn_KI270414v1,chrUn_KI270419v1,chrUn_KI270418v1,chrUn_KI270420v1,chrUn_KI270424v1,chrUn_KI270417v1,chrUn_KI270422v1,chrUn_KI270423v1,chrUn_KI270425v1,chrUn_KI270429v1,chrUn_KI270442v1,chrUn_KI270466v1,chrUn_KI270465v1,chrUn_KI270467v1,chrUn_KI270435v1,chrUn_KI270438v1,chrUn_KI270468v1,chrUn_KI270510v1,chrUn_KI270509v1,chrUn_KI270518v1,chrUn_KI270508v1,chrUn_KI270516v1,chrUn_KI270512v1,chrUn_KI270519v1,chrUn_KI270522v1,chrUn_KI270511v1,chrUn_KI270515v1,chrUn_KI270507v1,chrUn_KI270517v1,chrUn_KI270529v1,chrUn_KI270528v1,chrUn_KI270530v1,chrUn_KI270539v1,chrUn_KI270538v1,chrUn_KI270544v1,chrUn_KI270548v1,chrUn_KI270583v1,chrUn_KI270587v1,chrUn_KI270580v1,chrUn_KI270581v1,chrUn_KI270579v1,chrUn_KI270589v1,chrUn_KI270590v1,chrUn_KI270584v1,chrUn_KI270582v1,chrUn_KI270588v1,chrUn_KI270593v1,chrUn_KI270591v1,chrUn_KI270330v1,chrUn_KI270329v1,chrUn_KI270334v1,chrUn_KI270333v1,chrUn_KI270335v1,chrUn_KI270338v1,chrUn_KI270340v1,chrUn_KI270336v1,chrUn_KI270337v1,chrUn_KI270363v1,chrUn_KI270364v1,chrUn_KI270362v1,chrUn_KI270366v1,chrUn_KI270378v1,chrUn_KI270379v1,chrUn_KI270389v1,chrUn_KI270390v1,chrUn_KI270387v1,chrUn_KI270395v1,chrUn_KI270396v1,chrUn_KI270388v1,chrUn_KI270394v1,chrUn_KI270386v1,chrUn_KI270391v1,chrUn_KI270383v1,chrUn_KI270393v1,chrUn_KI270384v1,chrUn_KI270392v1,chrUn_KI270381v1,chrUn_KI270385v1,chrUn_KI270382v1,chrUn_KI270376v1,chrUn_KI270374v1,chrUn_KI270372v1,chrUn_KI270373v1,chrUn_KI270375v1,chrUn_KI270371v1,chrUn_KI270448v1,chrUn_KI270521v1,chrUn_GL000195v1,chrUn_GL000219v1,chrUn_GL000220v1,chrUn_GL000224v1,chrUn_KI270741v1,chrUn_GL000226v1,chrUn_GL000213v1,chrUn_KI270743v1,chrUn_KI270744v1,chrUn_KI270745v1,chrUn_KI270746v1,chrUn_KI270747v1,chrUn_KI270748v1,chrUn_KI270749v1,chrUn_KI270750v1,chrUn_KI270751v1,chrUn_KI270752v1,chrUn_KI270753v1,chrUn_KI270754v1,chrUn_KI270755v1,chrUn_KI270756v1,chrUn_KI270757v1,chrUn_GL000214v1,chrUn_KI270742v1,chrUn_GL000216v2,chrUn_GL000218v1,chrEBV'
	for tissue in tissues:
		bams = []
		np_default_file = tissue + '.default.narrowPeak'
		np_pcrdup_file = tissue + '.pcr_dup.narrowPeak'
		np_d200_file = tissue + '.d200.narrowPeak'
		np_p1_default_file = tissue + '.default.p0.001.narrowPeak'
		np_p1_pcrdup_file = tissue + '.pcr_dup.p0.001.narrowPeak'
		np_p1_d200_file = tissue + '.d200.p0.001.narrowPeak'
		np_p2_default_file = tissue + '.default.p0.0001.narrowPeak'
		np_p2_pcrdup_file = tissue + '.pcr_dup.p0.0001.narrowPeak'
		np_p2_d200_file = tissue + '.d200.p0.0001.narrowPeak'
		default_log_file = tissue + '.default.log'
		pcrdup_log_file = tissue + '.pcrdup.log'
		d200_log_file = tissue + '.d200.log'

		##get sorted bams
		for sample in samples:
			bam = sample + '/' + tissue + '.bam'
			nsort_bam = tissue + '.' + sample + '.nsorted_temp.bam'
			name = tissue + '_' + sample
			# st_sort_pe = subprocess.Popen(['samtools', 'sort', '-nO', 'bam', '-T', name, '-o', nsort_bam, '-@', '10', '-m', '10G', bam])
			# st_sort_pe.wait()
			bams.append(nsort_bam)
		print(bams)
		##run genrich
		## -j = atac mode, ? use -d cut sites (default =100), or -y so use unpaired reads, -r duplicate reads ?
		##this produces a log file i.e. -f so can run with different parameters using the log, so do 3x and then change p values
		run_genrich = subprocess.Popen([genrich, '-t', ','.join(bams), '-j', '-e', chrs_to_exclude, '-f', default_log_file, '-X'])
		run_genrich.wait()
		run_genrich = subprocess.Popen([genrich, '-t', ','.join(bams), '-j', '-e', chrs_to_exclude, '-f', pcrdup_log_file, '-X', '-r'])
		run_genrich.wait()
		run_genrich = subprocess.Popen([genrich, '-t', ','.join(bams), '-j', '-e', chrs_to_exclude, '-f', d200_log_file, '-X', '-d', '200'])
		run_genrich.wait()
		##then use log to generate results
		run_genrich2 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_default_file, '-e', chrs_to_exclude])
		run_genrich2.wait()
		run_genrich3 = subprocess.Popen([genrich, '-P', '-f', pcrdup_log_file, '-j', '-o', np_pcrdup_file, '-e', chrs_to_exclude])
		run_genrich3.wait()
		run_genrich4 = subprocess.Popen([genrich, '-P', '-f', d200_log_file, '-j', '-o', np_d200_file, '-e', chrs_to_exclude])
		run_genrich4.wait()
		run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p1_default_file, '-e', chrs_to_exclude, '-p', '0.001'])
		run_genrich5.wait()
		run_genrich6 = subprocess.Popen([genrich, '-P', '-f', pcrdup_log_file, '-j', '-o', np_p1_pcrdup_file, '-e', chrs_to_exclude, '-p', '0.001'])
		run_genrich6.wait()
		run_genrich7 = subprocess.Popen([genrich, '-P', '-f', d200_log_file, '-j', '-o', np_p1_d200_file, '-e', chrs_to_exclude, '-p', '0.001'])
		run_genrich7.wait()
		run_genrich8 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p2_default_file, '-e', chrs_to_exclude, '-p', '0.0001'])
		run_genrich8.wait()
		run_genrich9 = subprocess.Popen([genrich, '-P', '-f', pcrdup_log_file, '-j', '-o', np_p2_pcrdup_file, '-e', chrs_to_exclude, '-p', '0.0001'])
		run_genrich9.wait()
		run_genrich10 = subprocess.Popen([genrich, '-P', '-f', d200_log_file, '-j', '-o', np_p2_d200_file, '-e', chrs_to_exclude, '-p', '0.0001'])
		run_genrich10.wait()
		##egs without using log
		# run_genrich2 = subprocess.Popen([genrich, '-t', ','.join(bams), '-j', '-o', ea, '-e', chrs_to_exclude])
		# run_genrich2.wait()
		# run_genrich3 = subprocess.Popen([genrich, '-t', ','.join(bams), '-j', '-r', '-o', np_pcrdup_file, '-e', chrs_to_exclude])
		# run_genrich3.wait()
		# run_genrich4 = subprocess.Popen([genrich, '-t', ','.join(bams), '-j', '-d', '200', '-o', np_d200_file, '-e', chrs_to_exclude])
		# run_genrich4.wait()

##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_human_scatac_call_peaks_0720'
os.chdir(working_dir)
##info
##bam and barcode file
hu5_cell_info = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/call_peaks_macs2/human_adult.hu5.cell_types.txt'
hu7_cell_info = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/call_peaks_macs2/human_adult.hu7.cell_types.txt'
hu8_cell_info = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/call_peaks_macs2/human_adult.hu8.cell_types.txt'
hu5_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu5/outs/possorted_bam.bam'
hu7_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu7/outs/possorted_bam.bam'
hu8_bam = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/Hu8/outs/possorted_bam.bam'
##samples and tissue
sample_names = ['hu5', 'hu7', 'hu8']
tissue_names = ['Amacrines', 'Bipolars', 'Cones', 'Ganglions', 'Horizontals', 'Mullers', 'Rods']

##make bams per tissue type
##need python3 and sinto conda env
# run_sinto_method('hu5', hu5_bam, hu5_cell_info)
# run_sinto_method('hu7', hu7_bam, hu7_cell_info)
# run_sinto_method('hu8', hu8_bam, hu8_cell_info)

##run genrich and macs2 to call peaks on bams
call_peaks_on_bams(sample_names, tissue_names)




