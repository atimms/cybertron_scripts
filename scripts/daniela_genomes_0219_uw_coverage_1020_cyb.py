#!/usr/bin/python
import sys
import subprocess
import os
import glob
# import dobyns_gemini_pipeline_cybertron_v8

##parameters
delim = '\t'


##file names etc
ref_dir = '/home/atimms/ngs_data/references/hg19/'
fasta = ref_dir + 'human_g1k_v37.fasta'
fasta_fai = ref_dir + 'human_g1k_v37.genome.txt'
exome_bed_for_coverage = ref_dir + 'hg19_refGene_coding_exons.nochr_sorted_merged.bed'



##get exome coverage for bam files
def calculate_exome_coverage(samples, bam_suffix, exome_bed, output_file):
	##parameters and header
	coverage_required = [1,5,10,20,50,100,200]
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
	print output_file
	with open(output_file, "w") as outfh:
		outfh.write(header)
		##go through bam files and get coverage histogram
		for sample in samples:
			bam = sample + bam_suffix
			print 'calculating coverge for bam file', bam
			##get temp coverage files
			coverage_histogram = sample + '.hist.temp'
			hist_fh = open(coverage_histogram, 'w')
			# bedtools_cov = subprocess.Popen([bedtools, 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen(['bedtools', 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted'], stdout=subprocess.PIPE)
			# bedtools_cov = subprocess.Popen(['bedtools', 'coverage', '-a', exome_bed, '-b', bam, '-hist', '-sorted', '-g', fasta_fai], stdout=subprocess.PIPE)

			grep_cmd = subprocess.Popen(['grep', '^all'], stdin=bedtools_cov.stdout, stdout=hist_fh)
			grep_cmd.wait()
			hist_fh.close()
			##calculate numbers from histogram
			with open(coverage_histogram, "r") as cov_temp:
				seq_total = 0
				for line in cov_temp:
					line = line.strip('\n').split(delim)
					target_size = float(line[3])
					if line[1] != '0':
						seq_total += int(line[1]) *  int(line[2])
			cov_stats = []
			for cr in coverage_required:
				coverage_bp = 0
				with open(coverage_histogram, "r") as cov_temp:
					for line in cov_temp:
						line = line.strip('\n').split(delim)
						target_size = float(line[3])
						if int(line[1]) >= cr:
							coverage_bp += int(line[2])
					cov_stats.extend([str(coverage_bp), str((coverage_bp / target_size) * 100)])
			line_out = delim.join([bam, str(target_size), str(seq_total), str(seq_total / target_size)] + cov_stats + ['\n'])
			outfh.write(line_out)

def calculate_genome_coverage(samples, bam_suffix, output_file):
	##parameters and header
	coverage_required = [1,5,10,20,50,100,200]
	cov_head = []
	for cr in coverage_required:
		cr1 = 'total bp >=' + str(cr)
		cr2 = 'percentage covered >=' + str(cr)
		cov_head.extend([cr1, cr2])
	header = delim.join(['bam file', 'target size', 'total bp', 'coverage'] + cov_head + ['\n'])
	##open results file and write header
	print output_file
	with open(output_file, "w") as outfh:
		outfh.write(header)
		##go through bam files and get coverage histogram
		for sample in samples:
			bam = sample + bam_suffix
			print 'calculating coverge for bam file', bam
			##get temp coverage files
			coverage_histogram = sample + '.hist.temp'
			'''
			hist_fh = open(coverage_histogram, 'w')
			##just the file
			# bedtools_cov = subprocess.Popen([bedtools, 'genomecov', '-ibam', bam, '-g', genome_file], stdout=hist_fh)
			# bedtools_cov.wait()
			bedtools_cov = subprocess.Popen(['bedtools', 'genomecov', '-ibam', bam, '-g', fasta_fai], stdout=subprocess.PIPE)
			grep_cmd = subprocess.Popen(['grep', '^genome'], stdin=bedtools_cov.stdout, stdout=hist_fh)
			grep_cmd.wait()
			hist_fh.close()
			'''
			##calculate numbers from histogram
			with open(coverage_histogram, "r") as cov_temp:
				seq_total = 0
				for line in cov_temp:
					line = line.strip('\n').split(delim)
					target_size = float(line[3])
					if line[1] != '0':
						seq_total += int(line[1]) *  int(line[2])
			cov_stats = []
			for cr in coverage_required:
				coverage_bp = 0
				with open(coverage_histogram, "r") as cov_temp:
					for line in cov_temp:
						line = line.strip('\n').split(delim)
						target_size = float(line[3])
						if int(line[1]) >= cr:
							coverage_bp += int(line[2])
					cov_stats.extend([str(coverage_bp), str((coverage_bp / target_size) * 100)])
			line_out = delim.join([bam, str(target_size), str(seq_total), str(seq_total / target_size)] + cov_stats + ['\n'])
			outfh.write(line_out)


##
sample_names = ['CFM-MOS-11-01.338865', 'CFM-MOS-12-01.338868', 'CFM-MOS-12-11.338869', 'CFM-MOS-12-21.338870', 'CL100301.338615', 
		'CL100302.338616', 'CL100303.338617', 'CL101302.366650', 'CL101303.366651', 'CL101601.366652', 'CL101602.366653', 'CL101603.366654', 
		'CL101902.366656', 'CL102001.338618', 'CL102002.338619', 'CL102003.338620', 'CL300301.338621', 'CL300302.338622', 'CL302301.338624', 
		'CL302303.338626', 'CL303301.338627', 'CL303302.338628', 'CL303303.338629', 'CL304401.338630', 'CL304402.338631', 'CL304403.338632', 
		'CL304501.338633', 'CL304502.338634', 'CL304503.338635', 'CL305101.366658', 'CL305102.366659', 'CL308102.366662', 'CL400101.338636', 
		'CL400102.338637', 'CL400103.338638', 'CL400401.338639', 'CL400402.338640', 'CL400403.338641', 'CL400501.338642', 'CL400502.338643', 
		'CL400503.338644', 'CL400701.338645', 'CL400702.338646', 'CL400703.338647', 'F02-00003.338648', 'F02-00004.338649', 'F02-00005.338650', 
		'F02-00013.338651', 'F02-00014.338652', 'F02-00015.338653', 'F02-00020.338654', 'F02-00021.338655', 'F02-00022.338656', 'F03-00006.338657', 
		'F03-00007.338658', 'F03-00008.338659', 'F03-00011.338660', 'F03-00012.338661', 'F03-00015.338662', 'F03-00016.338663', 'F03-00017.338664', 
		'F03-00018.338665', 'F03-00019.338666', 'F03-00020.338667', 'F03-00021.338674', 'F03-00022.338668', 'F03-00024.338669', 'F03-00025.338670', 
		'F03-00026.338671', 'F03-00027.338672', 'F03-00034.338673', 'F03-00035.338675', 'F03-00036.338676', 'F03-00037.338677', 'F03-00038.338678', 
		'F03-00039.338679', 'F03-00040.338680', 'F03-00041.338681', 'F03-00045.338682', 'F03-00046.338683', 'F03-00047.338684', 'F03-00050.338685', 
		'F03-00051.338686', 'F03-00052.338687', 'F03-00053.338688', 'F03-00062.338689', 'F03-00063.338690', 'F03-00064.338691', 'F03-00065.338692', 
		'F03-00066.338693', 'F03-00067.338694', 'F03-00091.338695', 'F03-00092.338696', 'F03-00093.338697', 'F03-00094.338698', 'F03-00095.338699', 
		'F03-00096.338700', 'F03-00106.338701', 'F03-00107.338702', 'F03-00108.338703', 'F03-00114.338704', 'F03-00115.338705', 'F03-00116.338747', 
		'F03-00129.338706', 'F03-00130.338707', 'F03-00131.338708', 'F03-00133.338710', 'F03-00134.338711', 'F04-00012.338712', 'F04-00013.338713', 
		'F04-00014.338714', 'F04-00020.338715', 'F04-00021.338716', 'F04-00022.338717', 'F200006.366625', 'F200007.366626', 'F200008.366627', 
		'F3-0002-01.338718', 'F3-0002-11.338719', 'F3-0002-21.338720', 'F3-0003-01.338721', 'F3-0003-11.338722', 'F300031.366628', 'F3-0003-21.338723', 
		'F300032.366629', 'F300033.366630', 'F300059.366631', 'F3-0006-01.338724', 'F300060.366632', 'F3-0006-11.338725', 'F300061.366633', 
		'F3-0006-21.338726', 'F3-0007-11.338727', 'F300073.366634', 'F300074.366635', 'F300075.366636', 'F3-0009-11.338730', 'F3-0009-21.338731', 
		'F300097.366637', 'F300098.366638', 'F300099.366639', 'F3-0010-01.338732', 'F3-0010-11.338733', 'F3-0010-21.338734', 'F300111.366640', 
		'F300112.366641', 'F300113.366642', 'F3-0013-11.338736', 'F3-0017-01.338738', 'F3-0021-01.338744', 'F3-0021-21.338746', 'F3002811.366669', 
		'F3002821.366668', 'F400017.366643', 'F400018.366644', 'F400019.366645', 'F400025.366646', 'F400026.366647', 'F400027.366648', 'PE2D1004.366667', 
		'PE2D1012.338729', 'PM101000201.366673', 'PM101000202.366674', 'PM101000203.366675', 'PM101000401.366676', 'PM101000402.366677', 
		'PM101000403.366678', 'PM101000501.366679', 'PM101000502.366680', 'PM101000503.366681', 'PM101000701.366682', 'PM101000702.366683', 
		'PM101000703.366684', 'PM101001001.366685', 'PM101001002.366686', 'PM101001003.366687', 'PM101001701.338748', 'PM101001702.338749', 
		'PM101001703.338750', 'PM101003601.338751', 'PM101003602.338752', 'PM101003603.338753', 'PM101003701.338754', 'PM101003702.338755', 
		'PM101003703.338756', 'PM101003901.338760', 'PM101003902.338761', 'PM101003903.338762', 'PM101004001.338763', 'PM101004002.338764', 
		'PM101004003.338765', 'PM101004201.338766', 'PM101004202.338767', 'PM101004203.338768', 'PM101004301.338769', 'PM101004302.338770', 
		'PM101004303.338771', 'PM101004401.338772', 'PM101004402.338773', 'PM101004403.338774', 'PM101004501.338775', 'PM101004502.338776', 
		'PM101004503.338777', 'PM103001901.338781', 'PM103001902.338782', 'PM103001903.338783', 'PM103002101.338787', 'PM103002102.338788', 
		'PM103002103.338789', 'PM103002201.338790', 'PM103002202.338791', 'PM103002203.338792', 'PM103002501.338795', 'PM103002502.338796', 
		'PM103002503.338797', 'PM103002801.338798', 'PM103002802.338799', 'PM103002803.338800', 'PM103002901.338801', 'PM103002902.338802', 
		'PM103002903.338803', 'PM103003001.338804', 'PM103003002.338805', 'PM103003003.338806', 'PM103003101.338807', 'PM103003102.338808', 
		'PM103003103.338809', 'PM106001801.338810', 'PM106001802.338811', 'PM106001803.338812', 'PM106001806a.338813', 'PM106001806b.338814', 
		'PM106001807a.338815', 'PM106001807b.338816', 'PM106001809.338817', 'PM106001901.338818', 'PM106001902.338819', 'PM106001903.338820', 
		'PM106002801.338829', 'PM106002802.338830', 'PM106002803.338831', 'PM106002901.338832', 'PM106002902.338833', 'PM106002903.338834', 
		'PM106002907.338835', 'PM106003001.338836', 'PM106003002.338873', 'PM106003003.338837', 'PM106003101.366670', 'PM106003102.366671', 
		'PM106003103.366672', 'PM106003201.338838', 'PM106003202.338839', 'PM106003203.338840', 'PM107000501.338841', 'PM107000502.338842', 
		'PM107000503.338843', 'PM107001701.338844', 'PM107001702.338845', 'PM107001703.338846', 'PM107001901.338847', 'PM107001902.338848', 
		'PM107001903.338849', 'PM107002001.338850', 'PM107002002.338851', 'PM107002003.338852', 'PM401000101.338853', 'PM401000102.338854', 
		'PM401000103.338855', 'PM401000301.338856', 'PM401000302.338857', 'PM401000303.338858', 'PM401000401.338859', 'PM401000402.338860', 
		'PM401000403.338861', 'PM401000501.338862', 'PM401000502.338863', 'PM401000503.338864']

# sample_names = ['CFM-MOS-11-01.338865']

##params etc
working_dir = '/home/atimms/ngs_data/genomes/daniela_genomes_0219/daniela_svs_1020'
os.chdir(working_dir)

# calculate_exome_coverage(sample_names, '.bam', exome_bed_for_coverage, 'genomes_0219.refGene_coding_exon_coverage.xls')
# sample_names = ['CFM-MOS-11-01.338865']
calculate_genome_coverage(sample_names, '.bam', 'genomes_0219.genome_coverage.xls')



