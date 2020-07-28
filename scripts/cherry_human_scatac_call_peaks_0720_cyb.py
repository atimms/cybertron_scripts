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

##for macs2
qsub -Iq cdbrmq -l mem=100gb,ncpus=5,walltime=10:00:00 -P 19833a08-f6fb-4bea-8526-8a79069da878
module load local_python/3.7.6 
source activate macs2

##idr 
can just load module i.e.
module load idr/2.0.4
https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html

##homer
load modules:
qsub -Iq cdbrmq -l mem=100gb,ncpus=5,walltime=10:00:00 -P 19833a08-f6fb-4bea-8526-8a79069da87
module load biobuilds
module load homer/4.9.1

'''

##parameters
delim = '\t'
thread_number = '20'
##programs
sinto = '/home/atimms/programs/sinto/scripts/sinto'
genrich = '/home/atimms/programs/Genrich/Genrich'
bedtools = '/home/atimms/programs/bedtools2.28/bin/bedtools'
bg_to_bw = '/home/atimms/programs/bedGraphToBigWig'

##files
##copied from /active/cherry_t/NGS_Data/ChIP_input_for_bkgrnd/Hu1-ret-input-peakcalling-R1_trimmed-fixchr.bed.gz
control_bed = 'Hu1-ret-input-peakcalling-R1_trimmed-fixchr.bed'
chrom_sizes = '/home/atimms/ngs_data/references/10x/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.chrom_sizes'

##methods
def run_sinto_method(sample_name, bam_file, cells_file):
	directory = './' + sample_name
	if not os.path.exists(sample_name):
		os.makedirs(sample_name)
	os.chdir(directory)
	run_sinto = subprocess.Popen([sinto, 'filterbarcodes', '-b', bam_file,  '-c', cells_file, '-p', thread_number])
	run_sinto.wait()

def call_genrich_on_bams_first(samples, tissues):
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


def call_genrich_on_bams_second(samples, tissues):
	chrs_to_exclude = 'chrY,chrM,chr1_KI270706v1_random,chr1_KI270707v1_random,chr1_KI270708v1_random,chr1_KI270709v1_random,chr1_KI270710v1_random,chr1_KI270711v1_random,chr1_KI270712v1_random,chr1_KI270713v1_random,chr1_KI270714v1_random,chr2_KI270715v1_random,chr2_KI270716v1_random,chr3_GL000221v1_random,chr4_GL000008v2_random,chr5_GL000208v1_random,chr9_KI270717v1_random,chr9_KI270718v1_random,chr9_KI270719v1_random,chr9_KI270720v1_random,chr11_KI270721v1_random,chr14_GL000009v2_random,chr14_GL000225v1_random,chr14_KI270722v1_random,chr14_GL000194v1_random,chr14_KI270723v1_random,chr14_KI270724v1_random,chr14_KI270725v1_random,chr14_KI270726v1_random,chr15_KI270727v1_random,chr16_KI270728v1_random,chr17_GL000205v2_random,chr17_KI270729v1_random,chr17_KI270730v1_random,chr22_KI270731v1_random,chr22_KI270732v1_random,chr22_KI270733v1_random,chr22_KI270734v1_random,chr22_KI270735v1_random,chr22_KI270736v1_random,chr22_KI270737v1_random,chr22_KI270738v1_random,chr22_KI270739v1_random,chrY_KI270740v1_random,chrUn_KI270302v1,chrUn_KI270304v1,chrUn_KI270303v1,chrUn_KI270305v1,chrUn_KI270322v1,chrUn_KI270320v1,chrUn_KI270310v1,chrUn_KI270316v1,chrUn_KI270315v1,chrUn_KI270312v1,chrUn_KI270311v1,chrUn_KI270317v1,chrUn_KI270412v1,chrUn_KI270411v1,chrUn_KI270414v1,chrUn_KI270419v1,chrUn_KI270418v1,chrUn_KI270420v1,chrUn_KI270424v1,chrUn_KI270417v1,chrUn_KI270422v1,chrUn_KI270423v1,chrUn_KI270425v1,chrUn_KI270429v1,chrUn_KI270442v1,chrUn_KI270466v1,chrUn_KI270465v1,chrUn_KI270467v1,chrUn_KI270435v1,chrUn_KI270438v1,chrUn_KI270468v1,chrUn_KI270510v1,chrUn_KI270509v1,chrUn_KI270518v1,chrUn_KI270508v1,chrUn_KI270516v1,chrUn_KI270512v1,chrUn_KI270519v1,chrUn_KI270522v1,chrUn_KI270511v1,chrUn_KI270515v1,chrUn_KI270507v1,chrUn_KI270517v1,chrUn_KI270529v1,chrUn_KI270528v1,chrUn_KI270530v1,chrUn_KI270539v1,chrUn_KI270538v1,chrUn_KI270544v1,chrUn_KI270548v1,chrUn_KI270583v1,chrUn_KI270587v1,chrUn_KI270580v1,chrUn_KI270581v1,chrUn_KI270579v1,chrUn_KI270589v1,chrUn_KI270590v1,chrUn_KI270584v1,chrUn_KI270582v1,chrUn_KI270588v1,chrUn_KI270593v1,chrUn_KI270591v1,chrUn_KI270330v1,chrUn_KI270329v1,chrUn_KI270334v1,chrUn_KI270333v1,chrUn_KI270335v1,chrUn_KI270338v1,chrUn_KI270340v1,chrUn_KI270336v1,chrUn_KI270337v1,chrUn_KI270363v1,chrUn_KI270364v1,chrUn_KI270362v1,chrUn_KI270366v1,chrUn_KI270378v1,chrUn_KI270379v1,chrUn_KI270389v1,chrUn_KI270390v1,chrUn_KI270387v1,chrUn_KI270395v1,chrUn_KI270396v1,chrUn_KI270388v1,chrUn_KI270394v1,chrUn_KI270386v1,chrUn_KI270391v1,chrUn_KI270383v1,chrUn_KI270393v1,chrUn_KI270384v1,chrUn_KI270392v1,chrUn_KI270381v1,chrUn_KI270385v1,chrUn_KI270382v1,chrUn_KI270376v1,chrUn_KI270374v1,chrUn_KI270372v1,chrUn_KI270373v1,chrUn_KI270375v1,chrUn_KI270371v1,chrUn_KI270448v1,chrUn_KI270521v1,chrUn_GL000195v1,chrUn_GL000219v1,chrUn_GL000220v1,chrUn_GL000224v1,chrUn_KI270741v1,chrUn_GL000226v1,chrUn_GL000213v1,chrUn_KI270743v1,chrUn_KI270744v1,chrUn_KI270745v1,chrUn_KI270746v1,chrUn_KI270747v1,chrUn_KI270748v1,chrUn_KI270749v1,chrUn_KI270750v1,chrUn_KI270751v1,chrUn_KI270752v1,chrUn_KI270753v1,chrUn_KI270754v1,chrUn_KI270755v1,chrUn_KI270756v1,chrUn_KI270757v1,chrUn_GL000214v1,chrUn_KI270742v1,chrUn_GL000216v2,chrUn_GL000218v1,chrEBV'
	for tissue in tissues:
		bams = []
		default_log_file = tissue + '.default.log'
		np_p1_default_file = tissue + '.gr.default.p0.0001.narrowPeak'
		np_p2_default_file = tissue + '.gr.default.p0.00001.narrowPeak'
		np_p3_default_file = tissue + '.gr.default.p0.000001.narrowPeak'
		np_p4_default_file = tissue + '.gr.default.p0.0000001.narrowPeak'
		np_p5_default_file = tissue + '.gr.default.p0.00000001.narrowPeak'
		np_p6_default_file = tissue + '.gr.default.p0.000000001.narrowPeak'
		np_p7_default_file = tissue + '.gr.default.p0.0000000001.narrowPeak'
		np_p8_default_file = tissue + '.gr.default.p0.00000000001.narrowPeak'
		np_p9_default_file = tissue + '.gr.default.p0.000000000001.narrowPeak'
		np_p10_default_file = tissue + '.gr.default.p0.0000000000001.narrowPeak'
		np_p11_default_file = tissue + '.gr.default.p0.00000000000001.narrowPeak'
		np_p12_default_file = tissue + '.gr.default.p0.000000000000001.narrowPeak'
		np_p13_default_file = tissue + '.gr.default.q0.05.narrowPeak'
		np_p14_default_file = tissue + '.gr.default.q0.01.narrowPeak'
		# np_p15_default_file = tissue + '.gr.default.q0.05_a20.narrowPeak'
		# np_p16_default_file = tissue + '.gr.default.q0.05_a100.narrowPeak'		
		np_pileup = tissue + '.gr.default.pileup.bed'
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
		# run_genrich = subprocess.Popen([genrich, '-t', ','.join(bams), '-j', '-e', chrs_to_exclude, '-f', default_log_file, '-X', '-k', np_pileup])
		# run_genrich.wait()

		##then use log to generate results
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p1_default_file, '-e', chrs_to_exclude, '-p', '0.00001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p2_default_file, '-e', chrs_to_exclude, '-p', '0.000001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p3_default_file, '-e', chrs_to_exclude, '-p', '0.0000001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p4_default_file, '-e', chrs_to_exclude, '-p', '0.00000001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p5_default_file, '-e', chrs_to_exclude, '-p', '0.000000001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p6_default_file, '-e', chrs_to_exclude, '-p', '0.0000000001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p7_default_file, '-e', chrs_to_exclude, '-p', '0.00000000001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p8_default_file, '-e', chrs_to_exclude, '-p', '0.000000000001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p9_default_file, '-e', chrs_to_exclude, '-p', '0.0000000000001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p10_default_file, '-e', chrs_to_exclude, '-p', '0.00000000000001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p11_default_file, '-e', chrs_to_exclude, '-p', '0.000000000000001'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p12_default_file, '-e', chrs_to_exclude, '-p', '0.0000000000000001'])
		# run_genrich5.wait()
		##didn't work
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p13_default_file, '-e', chrs_to_exclude, '-q', '0.05'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p14_default_file, '-e', chrs_to_exclude, '-q', '0.01'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p15_default_file, '-e', chrs_to_exclude, '-q', '0.05', '-a', '20.0'])
		# run_genrich5.wait()
		# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p16_default_file, '-e', chrs_to_exclude, '-q', '0.05', '-a', '100.0'])
		# run_genrich5.wait()
		##egs without using log
		run_genrich2 = subprocess.Popen([genrich, '-t', ','.join(bams), '-j', '-o', np_p13_default_file, '-e', chrs_to_exclude, '-q', '0.05'])
		run_genrich2.wait()
		run_genrich2 = subprocess.Popen([genrich, '-t', ','.join(bams), '-j', '-o', np_p14_default_file, '-e', chrs_to_exclude, '-q', '0.01'])
		run_genrich2.wait()

def call_genrich_on_individual_samples(samples, tissues):
	chrs_to_exclude = 'chrY,chrM,chr1_KI270706v1_random,chr1_KI270707v1_random,chr1_KI270708v1_random,chr1_KI270709v1_random,chr1_KI270710v1_random,chr1_KI270711v1_random,chr1_KI270712v1_random,chr1_KI270713v1_random,chr1_KI270714v1_random,chr2_KI270715v1_random,chr2_KI270716v1_random,chr3_GL000221v1_random,chr4_GL000008v2_random,chr5_GL000208v1_random,chr9_KI270717v1_random,chr9_KI270718v1_random,chr9_KI270719v1_random,chr9_KI270720v1_random,chr11_KI270721v1_random,chr14_GL000009v2_random,chr14_GL000225v1_random,chr14_KI270722v1_random,chr14_GL000194v1_random,chr14_KI270723v1_random,chr14_KI270724v1_random,chr14_KI270725v1_random,chr14_KI270726v1_random,chr15_KI270727v1_random,chr16_KI270728v1_random,chr17_GL000205v2_random,chr17_KI270729v1_random,chr17_KI270730v1_random,chr22_KI270731v1_random,chr22_KI270732v1_random,chr22_KI270733v1_random,chr22_KI270734v1_random,chr22_KI270735v1_random,chr22_KI270736v1_random,chr22_KI270737v1_random,chr22_KI270738v1_random,chr22_KI270739v1_random,chrY_KI270740v1_random,chrUn_KI270302v1,chrUn_KI270304v1,chrUn_KI270303v1,chrUn_KI270305v1,chrUn_KI270322v1,chrUn_KI270320v1,chrUn_KI270310v1,chrUn_KI270316v1,chrUn_KI270315v1,chrUn_KI270312v1,chrUn_KI270311v1,chrUn_KI270317v1,chrUn_KI270412v1,chrUn_KI270411v1,chrUn_KI270414v1,chrUn_KI270419v1,chrUn_KI270418v1,chrUn_KI270420v1,chrUn_KI270424v1,chrUn_KI270417v1,chrUn_KI270422v1,chrUn_KI270423v1,chrUn_KI270425v1,chrUn_KI270429v1,chrUn_KI270442v1,chrUn_KI270466v1,chrUn_KI270465v1,chrUn_KI270467v1,chrUn_KI270435v1,chrUn_KI270438v1,chrUn_KI270468v1,chrUn_KI270510v1,chrUn_KI270509v1,chrUn_KI270518v1,chrUn_KI270508v1,chrUn_KI270516v1,chrUn_KI270512v1,chrUn_KI270519v1,chrUn_KI270522v1,chrUn_KI270511v1,chrUn_KI270515v1,chrUn_KI270507v1,chrUn_KI270517v1,chrUn_KI270529v1,chrUn_KI270528v1,chrUn_KI270530v1,chrUn_KI270539v1,chrUn_KI270538v1,chrUn_KI270544v1,chrUn_KI270548v1,chrUn_KI270583v1,chrUn_KI270587v1,chrUn_KI270580v1,chrUn_KI270581v1,chrUn_KI270579v1,chrUn_KI270589v1,chrUn_KI270590v1,chrUn_KI270584v1,chrUn_KI270582v1,chrUn_KI270588v1,chrUn_KI270593v1,chrUn_KI270591v1,chrUn_KI270330v1,chrUn_KI270329v1,chrUn_KI270334v1,chrUn_KI270333v1,chrUn_KI270335v1,chrUn_KI270338v1,chrUn_KI270340v1,chrUn_KI270336v1,chrUn_KI270337v1,chrUn_KI270363v1,chrUn_KI270364v1,chrUn_KI270362v1,chrUn_KI270366v1,chrUn_KI270378v1,chrUn_KI270379v1,chrUn_KI270389v1,chrUn_KI270390v1,chrUn_KI270387v1,chrUn_KI270395v1,chrUn_KI270396v1,chrUn_KI270388v1,chrUn_KI270394v1,chrUn_KI270386v1,chrUn_KI270391v1,chrUn_KI270383v1,chrUn_KI270393v1,chrUn_KI270384v1,chrUn_KI270392v1,chrUn_KI270381v1,chrUn_KI270385v1,chrUn_KI270382v1,chrUn_KI270376v1,chrUn_KI270374v1,chrUn_KI270372v1,chrUn_KI270373v1,chrUn_KI270375v1,chrUn_KI270371v1,chrUn_KI270448v1,chrUn_KI270521v1,chrUn_GL000195v1,chrUn_GL000219v1,chrUn_GL000220v1,chrUn_GL000224v1,chrUn_KI270741v1,chrUn_GL000226v1,chrUn_GL000213v1,chrUn_KI270743v1,chrUn_KI270744v1,chrUn_KI270745v1,chrUn_KI270746v1,chrUn_KI270747v1,chrUn_KI270748v1,chrUn_KI270749v1,chrUn_KI270750v1,chrUn_KI270751v1,chrUn_KI270752v1,chrUn_KI270753v1,chrUn_KI270754v1,chrUn_KI270755v1,chrUn_KI270756v1,chrUn_KI270757v1,chrUn_GL000214v1,chrUn_KI270742v1,chrUn_GL000216v2,chrUn_GL000218v1,chrEBV'
	for tissue in tissues:

		##get sorted bams
		for sample in samples:
			bam = sample + '/' + tissue + '.bam'
			nsort_bam = tissue + '.' + sample + '.nsorted_temp.bam'
			name = tissue + '_' + sample
			np_p8_default_file = name + '.gr.default.p1e11.narrowPeak'
			np_p9_default_file = name + '.gr.default.p1e12.narrowPeak'
			np_1_default_file = name + '.gr.default.p1e5.narrowPeak'
			np_2_default_file = name + '.gr.default.p1e6.narrowPeak'
			np_3_default_file = name + '.gr.default.p1e4.narrowPeak'
			np_4_default_file = name + '.gr.default.p1e3.narrowPeak'
			default_log_file = name + '.default.log'
			##sort bams
			# st_sort_pe = subprocess.Popen(['samtools', 'sort', '-nO', 'bam', '-T', name, '-o', nsort_bam, '-@', '10', '-m', '10G', bam])
			# st_sort_pe.wait()
			##run genrich
			## -j = atac mode, ? use -d cut sites (default =100), or -y so use unpaired reads, -r duplicate reads ?
			##this produces a log file i.e. -f so can run with different parameters using the log, so do 3x and then change p values
			'''
			run_genrich = subprocess.Popen([genrich, '-t', nsort_bam, '-j', '-e', chrs_to_exclude, '-f', default_log_file, '-X'])
			run_genrich.wait()
			'''
			##then use log to generate results
			# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p8_default_file, '-e', chrs_to_exclude, '-p', '0.000000000001'])
			# run_genrich5.wait()
			# run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p9_default_file, '-e', chrs_to_exclude, '-p', '0.0000000000001'])
			# run_genrich5.wait()
			run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_1_default_file, '-e', chrs_to_exclude, '-p', '0.00001'])
			run_genrich5.wait()
			run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_2_default_file, '-e', chrs_to_exclude, '-p', '0.000001'])
			run_genrich5.wait()
			run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_3_default_file, '-e', chrs_to_exclude, '-p', '0.0001'])
			run_genrich5.wait()
			run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_4_default_file, '-e', chrs_to_exclude, '-p', '0.001'])
			run_genrich5.wait()

def call_genrich_on_bams_third(samples, tissues):
	chrs_to_exclude = 'chrY,chrM,chr1_KI270706v1_random,chr1_KI270707v1_random,chr1_KI270708v1_random,chr1_KI270709v1_random,chr1_KI270710v1_random,chr1_KI270711v1_random,chr1_KI270712v1_random,chr1_KI270713v1_random,chr1_KI270714v1_random,chr2_KI270715v1_random,chr2_KI270716v1_random,chr3_GL000221v1_random,chr4_GL000008v2_random,chr5_GL000208v1_random,chr9_KI270717v1_random,chr9_KI270718v1_random,chr9_KI270719v1_random,chr9_KI270720v1_random,chr11_KI270721v1_random,chr14_GL000009v2_random,chr14_GL000225v1_random,chr14_KI270722v1_random,chr14_GL000194v1_random,chr14_KI270723v1_random,chr14_KI270724v1_random,chr14_KI270725v1_random,chr14_KI270726v1_random,chr15_KI270727v1_random,chr16_KI270728v1_random,chr17_GL000205v2_random,chr17_KI270729v1_random,chr17_KI270730v1_random,chr22_KI270731v1_random,chr22_KI270732v1_random,chr22_KI270733v1_random,chr22_KI270734v1_random,chr22_KI270735v1_random,chr22_KI270736v1_random,chr22_KI270737v1_random,chr22_KI270738v1_random,chr22_KI270739v1_random,chrY_KI270740v1_random,chrUn_KI270302v1,chrUn_KI270304v1,chrUn_KI270303v1,chrUn_KI270305v1,chrUn_KI270322v1,chrUn_KI270320v1,chrUn_KI270310v1,chrUn_KI270316v1,chrUn_KI270315v1,chrUn_KI270312v1,chrUn_KI270311v1,chrUn_KI270317v1,chrUn_KI270412v1,chrUn_KI270411v1,chrUn_KI270414v1,chrUn_KI270419v1,chrUn_KI270418v1,chrUn_KI270420v1,chrUn_KI270424v1,chrUn_KI270417v1,chrUn_KI270422v1,chrUn_KI270423v1,chrUn_KI270425v1,chrUn_KI270429v1,chrUn_KI270442v1,chrUn_KI270466v1,chrUn_KI270465v1,chrUn_KI270467v1,chrUn_KI270435v1,chrUn_KI270438v1,chrUn_KI270468v1,chrUn_KI270510v1,chrUn_KI270509v1,chrUn_KI270518v1,chrUn_KI270508v1,chrUn_KI270516v1,chrUn_KI270512v1,chrUn_KI270519v1,chrUn_KI270522v1,chrUn_KI270511v1,chrUn_KI270515v1,chrUn_KI270507v1,chrUn_KI270517v1,chrUn_KI270529v1,chrUn_KI270528v1,chrUn_KI270530v1,chrUn_KI270539v1,chrUn_KI270538v1,chrUn_KI270544v1,chrUn_KI270548v1,chrUn_KI270583v1,chrUn_KI270587v1,chrUn_KI270580v1,chrUn_KI270581v1,chrUn_KI270579v1,chrUn_KI270589v1,chrUn_KI270590v1,chrUn_KI270584v1,chrUn_KI270582v1,chrUn_KI270588v1,chrUn_KI270593v1,chrUn_KI270591v1,chrUn_KI270330v1,chrUn_KI270329v1,chrUn_KI270334v1,chrUn_KI270333v1,chrUn_KI270335v1,chrUn_KI270338v1,chrUn_KI270340v1,chrUn_KI270336v1,chrUn_KI270337v1,chrUn_KI270363v1,chrUn_KI270364v1,chrUn_KI270362v1,chrUn_KI270366v1,chrUn_KI270378v1,chrUn_KI270379v1,chrUn_KI270389v1,chrUn_KI270390v1,chrUn_KI270387v1,chrUn_KI270395v1,chrUn_KI270396v1,chrUn_KI270388v1,chrUn_KI270394v1,chrUn_KI270386v1,chrUn_KI270391v1,chrUn_KI270383v1,chrUn_KI270393v1,chrUn_KI270384v1,chrUn_KI270392v1,chrUn_KI270381v1,chrUn_KI270385v1,chrUn_KI270382v1,chrUn_KI270376v1,chrUn_KI270374v1,chrUn_KI270372v1,chrUn_KI270373v1,chrUn_KI270375v1,chrUn_KI270371v1,chrUn_KI270448v1,chrUn_KI270521v1,chrUn_GL000195v1,chrUn_GL000219v1,chrUn_GL000220v1,chrUn_GL000224v1,chrUn_KI270741v1,chrUn_GL000226v1,chrUn_GL000213v1,chrUn_KI270743v1,chrUn_KI270744v1,chrUn_KI270745v1,chrUn_KI270746v1,chrUn_KI270747v1,chrUn_KI270748v1,chrUn_KI270749v1,chrUn_KI270750v1,chrUn_KI270751v1,chrUn_KI270752v1,chrUn_KI270753v1,chrUn_KI270754v1,chrUn_KI270755v1,chrUn_KI270756v1,chrUn_KI270757v1,chrUn_GL000214v1,chrUn_KI270742v1,chrUn_GL000216v2,chrUn_GL000218v1,chrEBV'
	for tissue in tissues:
		bams = []
		default_log_file = tissue + '.default.log'
		np_p8_default_file = tissue + '.gr.default.p1e11.narrowPeak'
		np_p9_default_file = tissue + '.gr.default.p1e12.narrowPeak'		
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
		##then use log to generate results
		run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p8_default_file, '-e', chrs_to_exclude, '-p', '0.000000000001'])
		run_genrich5.wait()
		run_genrich5 = subprocess.Popen([genrich, '-P', '-f', default_log_file, '-j', '-o', np_p9_default_file, '-e', chrs_to_exclude, '-p', '0.0000000000001'])
		run_genrich5.wait()



def call_macs2_on_bams(samples, tissues):
	'''
	##run macs2
	for tissue in tissues:
		for sample in samples:
			bam = sample + '/' + tissue + '.bam'
			bed = sample + '_' + tissue + '.temp.bed'
			name = tissue + '_' + sample
			##outfile names
			out_bampe_def = name + '.macs2.bampe_p0.01'
			out_bampe_keepdups =  name + '.macs2.bampe_p0.01_keepdups'
			out_bed_def = name + '.macs2.bed_p0.01'
			out_bed_keepdups =  name + '.macs2.bed_p0.01_keepdups'	
			out_bed_control =  name + '.macs2.bed_p0.01_control'	
			##run macs2
			##? bam vs bed, paired vs not, dups vs not, p value cut off?, control
			##from bam files
			run_macs2_1 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_def, '-g', 'hs', '-p', '1e-2'])
			run_macs2_1.wait()
			run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups, '-g', 'hs', '-p', '1e-2', '--keep-dup', 'all'])
			run_macs2_2.wait()
			##make bed and then run macs2
			with open(bed, "w") as bedfh:
				bedtools_run = subprocess.Popen([bedtools, 'bamtobed', '-i', bam], stdout=bedfh)
				bedtools_run.wait()
			run_macs2_3 = subprocess.Popen(['macs2', 'callpeak', '-t', bed, '-f', 'BED', '-n', out_bed_def, '-g', 'hs', '-p', '1e-2', '--nomodel', '--extsize', '200'])
			run_macs2_3.wait()
			run_macs2_4 = subprocess.Popen(['macs2', 'callpeak', '-t', bed, '-f', 'BED', '-n', out_bed_keepdups, '-g', 'hs', '-p', '1e-2', '--nomodel', '--extsize', '200', '--keep-dup', 'all'])
			run_macs2_4.wait()
			run_macs2_5 = subprocess.Popen(['macs2', 'callpeak', '-t', bed, '-c', control_bed, '-f', 'BED', '-n', out_bed_control, '-g', 'hs', '-p', '1e-2', '--nomodel', '--extsize', '200'])
			run_macs2_5.wait()
	'''
	##run idr on the 3 replicates
	for tissue in tissues:
		# analysis_prefix_type = ['.macs2.bampe_p0.01', '.macs2.bampe_p0.01_keepdups', '.macs2.bed_p0.01', '.macs2.bed_p0.01_keepdups', '.macs2.bed_p0.01_control']
		analysis_prefix_type = ['.macs2.bampe_p0.01_keepdups']

		for file_prefix in analysis_prefix_type:
			sorted_np_files = []
			for sample in samples:
				##sort each narrowpeak file by the -log10(p-value)
				#sort -k8,8nr NAME_OF_INPUT_peaks.narrowPeak > macs/NAME_FOR_OUPUT_peaks.narrowPeak
				np_file = tissue + '_' + sample + file_prefix + '_peaks.narrowPeak'
				sorted_np_file = tissue + '_' + sample + file_prefix + '_peaks_sorted.narrowPeak'
				sorted_np_files.append(sorted_np_file)
				# with open(sorted_np_file, "w") as snp_fh:
				# 	sort_np = subprocess.Popen(['sort', '-k8,8nr', np_file], stdout=snp_fh)
				# 	sort_np.wait()
			print(sorted_np_files)
			##run idr on 3 replicates in pairwise 
			idr_files = []
			for pair in ['hu5_hu7', 'hu5_hu8', 'hu7_hu8']:
				idr_out_file =  tissue + '.' + pair + file_prefix + '.idr_all.narrowPeak'
				idr_files.append(idr_out_file)
				#idr -s ${outPrefix}_combined_pseudorep1.regionPeak ${outPrefix}_combined_pseudorep2.regionPeak -o ${outPrefix}_combined_rep1_vs_rep2.txt 
				if pair == 'hu5_hu7':
					sort_np = subprocess.Popen(['idr', '-s', sorted_np_files[0], sorted_np_files[1], '--input-file-type', 'narrowPeak', '-o', idr_out_file, '--output-file-type', 'narrowPeak'])
					sort_np.wait()			
				elif pair == 'hu5_hu8':
					sort_np = subprocess.Popen(['idr', '-s', sorted_np_files[0], sorted_np_files[2], '--input-file-type', 'narrowPeak', '-o', idr_out_file, '--output-file-type', 'narrowPeak'])
					sort_np.wait()	
				elif pair == 'hu7_hu8':
					sort_np = subprocess.Popen(['idr', '-s', sorted_np_files[1], sorted_np_files[2], '--input-file-type', 'narrowPeak', '-o', idr_out_file, '--output-file-type', 'narrowPeak'])
					sort_np.wait()	
			##get peaks that pass all 3 pairwise test
			idr_filter_file =  tissue + file_prefix + '.idr_p0.005.narrowPeak'
			print(idr_files)

def make_norm_bedgraph(infile, outfile, mapped_reads):
	with open(infile, "r") as in_fh, open(outfile, "w") as out_fh:
		line_count = 0
		for line in in_fh:
			line = line.strip().split(delim)
			norm_count = int(line[3]) / (10000000.0/mapped_reads)
			chrom = line[0]
			##remove weird chrs
			if '_' not in chrom:
				out_fh.write(delim.join(line[:3] + [str(norm_count)]) + '\n')

def make_bigwig_files(samples, tissues):
	for tissue in tissues:
		for sample in samples:
			bam = sample + '/' + tissue + '.bam'
			bedgraph1 = sample + '_' + tissue + '.temp1.bedgraph'
			bedgraph2 = sample + '_' + tissue + '.temp2.bedgraph'
			out_bw = sample + '_' + tissue + '.bigwig'
			##make bedgraph from bam
			# '''
			with open(bedgraph1, 'w') as out_fh:
				bedtools_gc = subprocess.Popen(['bedtools', 'genomecov', '-g', chrom_sizes,  '-bg', '-ibam', bam], stdout=subprocess.PIPE)
				bedtools_sort = subprocess.Popen(['bedtools', 'sort', '-i', 'stdin'], stdin=bedtools_gc.stdout, stdout=out_fh)
				bedtools_sort.wait()
			# '''
			##get mapped and total read calculation
			# total_reads = subprocess.check_output(['samtools', 'view', '-c', bam])
			# total_reads.wait()
			# total_reads = total_reads.rstrip()
			# '''
			mapped_reads = subprocess.check_output(['samtools', 'view', '-c', '-F', '4', bam])
			mapped_reads = int(mapped_reads.rstrip())
			print(bam, mapped_reads)
			
			##normalize read count
			# mapped_reads = 25000000
			make_norm_bedgraph(bedgraph1, bedgraph2, mapped_reads)
			# '''
			##make bigwig
			bgbw = subprocess.Popen([bg_to_bw, bedgraph2, chrom_sizes, out_bw])
			bgbw.wait()			
			

def call_macs2_on_bams_second(samples, tissues):
	'''
	##run macs2
	for tissue in tissues:
		for sample in samples:
			bam = sample + '/' + tissue + '.bam'
			bed = sample + '_' + tissue + '.temp.bed'
			name = tissue +'_' + sample
			##outfile names
			out_bampe_keepdups =  name + '.macs2.bampe_p0.01_keepdups'
			##run macs2
			##from bam files
			run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups, '-g', 'hs', '-p', '1e-2', '--keep-dup', 'all'])
			run_macs2_2.wait()

	'''
	##run idr on the 3 replicates
	for tissue in tissues:
		analysis_prefix_type = ['.macs2.bampe_p0.01_keepdups']
		for file_prefix in analysis_prefix_type:
			sorted_np_files = []
			for sample in samples:
				##sort each narrowpeak file by the -log10(p-value)
				#sort -k8,8nr NAME_OF_INPUT_peaks.narrowPeak > macs/NAME_FOR_OUPUT_peaks.narrowPeak
				np_file = tissue + '_' + sample + file_prefix + '_peaks.narrowPeak'
				sorted_np_file = tissue + '_' + sample + file_prefix + '_peaks_sorted.narrowPeak'
				sorted_np_files.append(sorted_np_file)
				'''
				with open(sorted_np_file, "w") as snp_fh:
					sort_np = subprocess.Popen(['sort', '-k8,8nr', np_file], stdout=snp_fh)
					sort_np.wait()
				'''
			# print(sorted_np_files)
			##run idr on 3 replicates in pairwise 
			idr_files = []
			for pair in ['hu5_hu7', 'hu5_hu8', 'hu7_hu8']:
				idr_out_file =  tissue + '.' + pair + file_prefix + '.idr_all.narrowPeak'
				idr_files.append(idr_out_file)
				'''
				#idr -s ${outPrefix}_combined_pseudorep1.regionPeak ${outPrefix}_combined_pseudorep2.regionPeak -o ${outPrefix}_combined_rep1_vs_rep2.txt 
				if pair == 'hu5_hu7':
					sort_np = subprocess.Popen(['idr', '-s', sorted_np_files[0], sorted_np_files[1], '--input-file-type', 'narrowPeak', '-o', idr_out_file, '--output-file-type', 'narrowPeak'])
					sort_np.wait()			
				elif pair == 'hu5_hu8':
					sort_np = subprocess.Popen(['idr', '-s', sorted_np_files[0], sorted_np_files[2], '--input-file-type', 'narrowPeak', '-o', idr_out_file, '--output-file-type', 'narrowPeak'])
					sort_np.wait()	
				elif pair == 'hu7_hu8':
					sort_np = subprocess.Popen(['idr', '-s', sorted_np_files[1], sorted_np_files[2], '--input-file-type', 'narrowPeak', '-o', idr_out_file, '--output-file-type', 'narrowPeak'])
					sort_np.wait()
				'''

			##get peaks that pass all 3 pairwise test
			print(idr_files)
			temp_bt_file = tissue + file_prefix + '.idr_temp.txt'
			combined_idr_bed = tissue + file_prefix + '.idr_combined.bed'
			##sort idr files
			sorted_beds = []
			for idr_file in idr_files:
				out_bed = idr_file.rsplit('.',1)[0] + 'sorted_temp.bed'
				sorted_beds.append(out_bed)
				'''
				with open(out_bed, 'w') as out_fh:
					bedtools_sort = subprocess.Popen(['bedtools', 'sort', '-i', idr_file], stdout=out_fh)
					bedtools_sort.wait()
				'''
			##run bt mutiint on all three file
			'''
			with open(temp_bt_file, 'w') as out_fh:
				bt_mi = subprocess.Popen(['bedtools', 'multiinter', '-header', '-i'] + sorted_beds, stdout=out_fh)
				bt_mi.wait()
			'''
			with open(combined_idr_bed, 'w') as out_fh, open(temp_bt_file, 'r') as in_fh:
				lc = 0
				for line in in_fh:
					line = line.split(delim)
					lc += 1
					if lc > 1:
						##if in all 3 files
						if line[3] == '3':
							out_fh.write(delim.join(line[:3]) + '\n')

def call_macs2_on_bams_third(samples, tissues):
	# '''
	##run macs2
	for tissue in tissues:
		for sample in samples:
			bam = sample + '/' + tissue + '.bam'
			bed = sample + '_' + tissue + '.temp.bed'
			name = tissue +'_' + sample
			##outfile names
			out_bampe_keepdups_1 =  name + '.macs2.bampe_p1e-2_keepdups'
			out_bampe_keepdups_2 =  name + '.macs2.bampe_p1e-5_keepdups'
			out_bampe_keepdups_3 =  name + '.macs2.bampe_p1e-10_keepdups'
			out_bampe_keepdups_4 =  name + '.macs2.bampe_q0.01_keepdups'
			##run macs2
			##from bam files
			run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_1, '-g', 'hs', '-p', '1e-2', '--keep-dup', 'all'])
			run_macs2_2.wait()
			run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_2, '-g', 'hs', '-p', '1e-5', '--keep-dup', 'all'])
			run_macs2_2.wait()			
			run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_3, '-g', 'hs', '-p', '1e-10', '--keep-dup', 'all'])
			run_macs2_2.wait()
			run_macs2_2 = subprocess.Popen(['macs2', 'callpeak', '-t', bam, '-f', 'BAMPE', '-n', out_bampe_keepdups_4, '-g', 'hs', '-q', '0.01', '--keep-dup', 'all'])
			run_macs2_2.wait()
	# '''

def make_tag_dirs(name, bam):
	outdir = name + '.tag_dir'
	mk_tag_dir = subprocess.Popen(['makeTagDirectory', outdir, bam])
	mk_tag_dir.wait()

def merge_peaks(out_prefix, d_value, peak_files, out_suffix):
	out_file = out_prefix + d_value + 'd' + out_suffix
	print(len(peak_files), out_file)
	with open(out_file, 'w') as out_fh:
		# mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value, '-matrix', out_prefix + d_value] + peak_files, stdout=out_fh)
		mk_tag_dir = subprocess.Popen(['mergePeaks', '-d', d_value] + peak_files, stdout=out_fh)
		mk_tag_dir.wait()

def make_bed_from_narrowpeak(samples, tissues, in_suffices):
	for tissue in tissues:
		for sample in samples:
			for in_suffix in in_suffices:
				np_file =  tissue + '_' + sample + in_suffix
				out_bed = np_file.rsplit('.', 1)[0] + '.bed'
				with open(out_bed, 'w') as out_fh, open(np_file, 'r') as in_fh:
					for line in in_fh:
						line = line.split(delim)
						out_fh.write(delim.join(line[:6]) + '\n')

def annotate_peaks(samples, tissues, peak_prefix, peak_suffix, d_value, tag_suffix, out_prefix, sizes_req):
	peak_file = peak_prefix + d_value + 'd' + peak_suffix
	tag_dirs = []
	for tissue in tissues:
		for sample in samples:
			tag_dir = sample + '_' + tissue + tag_suffix
			tag_dirs.append(tag_dir)
	out_file_no_size = out_prefix + d_value + 'd' + peak_suffix
	print(len(tag_dirs), peak_file, out_file_no_size)
	with open(out_file_no_size, 'w') as out_ns_fh:
		mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-d'] + tag_dirs, stdout=out_ns_fh)
		mk_tag_dir.wait()
	for size_req in sizes_req:
		out_file_using_size = out_prefix + d_value + 'd.' + size_req + 'size' + peak_suffix
		with open(out_file_using_size, 'w') as out_fh:
			#annotatePeaks.pl pu1peaks.txt mm8 -size 400 -d Macrophage-PU.1/ Bcell-PU.1/ > output.txt
			##alternatives
			# mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg19', '-size', size_req], stdout=out_fh)
			# mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg19', '-noann', '-size', size_req, '-norm', '10000000', '-d'] + tag_dirs, stdout=out_fh)
			mk_tag_dir = subprocess.Popen(['annotatePeaks.pl', peak_file, 'hg38', '-size', size_req, '-d'] + tag_dirs, stdout=out_fh)
			mk_tag_dir.wait()

def compute_correlation(annotate_prefix, annotate_suffix, d_value, sizes_req, out_prefix):
	##get list of infile
	infiles = []
	no_size_file = annotate_prefix + d_value + 'd' + annotate_suffix
	infiles.append(no_size_file)
	for size_req in sizes_req:
		out_file_using_size = annotate_prefix + d_value + 'd.' + size_req + 'size' + annotate_suffix
		infiles.append(out_file_using_size)
	# print(len(infiles), infiles)
	##format file by file
	for infile in infiles:
		out_file = out_prefix  + infile.split('.', 1)[1]
		print(infile, out_file)
		with open(out_file, 'w') as out_fh, open(infile, 'r') as in_fh:
			lc = 0
			for line in in_fh:
				lc += 1
				line = line.strip('\n').split(delim)
				if lc == 1:
					sample_info = line[19:]
					sample_info = [s.split('.')[0] for s in sample_info]
					out_fh.write(delim.join(['peak_id'] + sample_info + ['\n']))
					print(sample_info)
				else:
					line_out = [line[0]] + line[19:]
					out_fh.write(delim.join(line_out + ['\n']))



def get_correlation_info_homer(samples, tissues, d_values, size_values, bed_suffixes):
	##make tag dirs
	for tissue in tissues:
		for sample in samples:
			name = sample + '_' + tissue
			bam = sample + '/' + tissue + '.bam'
			##make tag dirs so can analyze
			# make_tag_dirs(name, bam)
	##for different parameters merge peaks and then annoate
	tag_dir_suffix = '.tag_dir'
	merge_prefix = 'merged.'
	annotate_prefix = 'annotated.'
	correlation_prefix = 'correlation.'
	for d in d_values:
		for bed_suffix in bed_suffixes:
			##merge peak files for annotating the peaks
			peak_beds = sorted(glob.glob('*' + bed_suffix))
			merge_file_suffix = bed_suffix.rsplit('.', 1)[0] + '.txt'
			# merge_peaks(merge_prefix, d, peak_beds, merge_file_suffix)
			##annotate those peaks, and then make a correlation file
			annotate_peaks(samples, tissues, merge_prefix, merge_file_suffix, d, tag_dir_suffix, annotate_prefix, size_values)
			##make file to use in r for heatmaps etc
			compute_correlation(annotate_prefix, merge_file_suffix, d, size_values, correlation_prefix)

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
##homer values
d_values = ['100', '200', '500']
size_values = ['500', '2000']
ind_gr_narrowpeak_suffices = ['.gr.default.p1e3.narrowPeak', '.gr.default.p1e4.narrowPeak',
		'.gr.default.p1e5.narrowPeak', '.gr.default.p1e6.narrowPeak']
peak_bed_suffices = ['.macs2.bampe_p1e-10_keepdups_summits.bed', '.macs2.bampe_p1e-2_keepdups_summits.bed', 
	'.macs2.bampe_p1e-5_keepdups_summits.bed', '.macs2.bampe_q0.01_keepdups_summits.bed',
	'.gr.default.p1e3.bed', '.gr.default.p1e4.bed', '.gr.default.p1e5.bed', '.gr.default.p1e6.bed']


##temp
# sample_names = ['hu5']
# tissue_names = ['Amacrines']

##make bams per tissue type
##need python3 and sinto conda env
# run_sinto_method('hu5', hu5_bam, hu5_cell_info)
# run_sinto_method('hu7', hu7_bam, hu7_cell_info)
# run_sinto_method('hu8', hu8_bam, hu8_cell_info)

##make bigwig files from each bam
# make_bigwig_files(sample_names, tissue_names)

##run tests on genrich and macs2 to call peaks on bams
# call_genrich_on_bams_first(sample_names, tissue_names)
# call_macs2_on_bams(sample_names, tissue_names)
# call_genrich_on_bams_second(sample_names, tissue_names)

##decision make use p-value cut off for genrich
# call_genrich_on_bams_third(sample_names, tissue_names)
# call_genrich_on_individual_samples(sample_names, tissue_names)
# make_bed_from_narrowpeak(sample_names, tissue_names, ind_gr_narrowpeak_suffices)
##use bampe keepdups from macs2
# call_macs2_on_bams_second(sample_names, tissue_names)

##use different parameters for macs2 for entry to homer
# call_macs2_on_bams_third(sample_names, tissue_names)

##use homer to look for correlation
get_correlation_info_homer(sample_names, tissue_names, d_values, size_values, peak_bed_suffices)

