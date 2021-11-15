#!/usr/bin/env python
import subprocess
import os

##paramters
delim = '\t'


def reformat_gt_file(in_file):
	baf_out_file = in_file.rsplit('_', 1)[0] + '.baf.txt'
	cnv_out_file = in_file.rsplit('_', 1)[0] + '.cnv.txt'
	logr_out_file = in_file.rsplit('_', 1)[0] + '.lrr.txt'
	print(in_file, baf_out_file, cnv_out_file, logr_out_file)
	sample_list = []
	gt_dict = {}

	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			lc += 1
			if lc > 10:
				line = line.rstrip().split(delim)
				snp = line[0]
				snp_info = line[2:6]
				cnv_value = line[7]
				cnv_conf = line[8]
				baf = line[9]
				lrr = line[10]
				sample = line[1]
				if sample not in sample_list:
					sample_list.append(sample)
				if snp in gt_dict:
					gt_dict[snp][1].append(cnv_value)
					gt_dict[snp][2].append(cnv_conf)
					gt_dict[snp][3].append(baf)
					gt_dict[snp][4].append(lrr)
				else:
					gt_dict[snp] = [snp_info, [cnv_value], [cnv_conf], [baf], [lrr]]
	##check have all samples for all snps
	for gt in gt_dict:
	# 	print(gt, gt_dict[gt])
		if len(gt_dict[gt][1]) != len(sample_list):
			print('bugger.......', len(gt_dict[gt][1]), len(sample_list))
	with open(baf_out_file, "w") as bout_fh, open(cnv_out_file, "w") as cout_fh, open(logr_out_file, "w") as lout_fh:
		header = ['SNP Name', 'Chr', 'Position', 'Allele1 - Forward', 'Allele2 - Forward']
		baf_head = ['BAF_' + s for s in sample_list]
		cnvv_head = ['CNV_value_' + s for s in sample_list]
		cnvc_head = ['CNV_conf_' + s for s in sample_list]
		lrr_head = ['LRR_' + s for s in sample_list]
		bout_fh.write(delim.join(header + baf_head) + '\n')
		lout_fh.write(delim.join(header + lrr_head) + '\n')
		cout_fh.write(delim.join(header + cnvv_head + cnvc_head) + '\n')
		for g in gt_dict:
			bout_fh.write(delim.join([g] + gt_dict[g][0] + gt_dict[g][3]) + '\n')
			lout_fh.write(delim.join([g] + gt_dict[g][0] + gt_dict[g][4]) + '\n')
			cout_fh.write(delim.join([g] + gt_dict[g][0] + gt_dict[g][1] + gt_dict[g][2]) + '\n')


def get_cnv_region_by_sample(infile, sum_file):
	print(infile, sum_file)
	infile_sorted = infile.rsplit('.', 1)[0] + '.sorted.xls'
	infile_filtered = infile.rsplit('.', 1)[0] + '.filtered.xls'
	infile_filtered_sorted = infile.rsplit('.', 1)[0] + '.filtered_sorted.xls'
	'''
	with open(infile_sorted, "w") as ifs_fh:
		#sort -k2,2n -k3,3n test.cnv.txt > test.cnv_sorted.txt 
		if_sort = subprocess.Popen(['sort','-k2,2n', '-k3,3n', infile], stdout=ifs_fh)
		if_sort.wait()
	with open(infile_sorted, "r") as ins_fh, open(infile_filtered, "w") as inf_fh:
		lc = 0
		for line in ins_fh:
			line = line.rstrip("\n").split(delim)
			# line = line.split(delim)
			lc += 1
			if lc == 1:
				header = ['sample_id'] + line[1:5] + ['CNV_value', 'CNV_conf', 'var_number']
				inf_fh.write(delim.join(header) + '\n')
				value_conf_header = line[5:]
				samples = []
				for info in value_conf_header:
					sample_name = info.split('_')[2]
					if sample_name not in samples:
						samples.append(sample_name)

			# print(samples)
			else:
				info = line[1:5]
				value_conf = line[5:]
				# print(value_conf, len(value_conf))
				for sample in samples:
					v_index = samples.index(sample)
					c_index = samples.index(sample) + len(samples)
					value = value_conf[v_index]
					conf = value_conf[c_index]
					# print(sample, v_index, c_index, value, conf, lc)
					if conf != '':
						if float(conf) >= 35:
							line_out = [sample] + info + [value, conf, str(lc)]
							inf_fh.write(delim.join(line_out) + '\n')
	'''
	# '''
	##the sort the filtered file
	with open(infile_filtered_sorted, "w") as iffs_fh:
		#sort -k1,1n -k8,8n test.cnv.filtered.xls
		if_sort = subprocess.Popen(['sort','-k1,1n', '-k8,8n', infile_filtered], stdout=iffs_fh)
		if_sort.wait()
	# '''
	with open(infile_filtered_sorted, "r") as infs_fh, open(sum_file, "w") as out_fh:
		lc = 0
		var_numbers, chroms, positions, vals, confs, samples = [], [], [], [], [], []
		for line in infs_fh:
			line = line.rstrip("\n").split(delim)
			lc += 1
			if lc == 1:
				header = ['Sample ID', 'Chr', 'Start', 'Stop', 'Size', 'Probes >=3', 'CNV_value', 'CNV_conf >35']
				out_fh.write(delim.join(header) + '\n')
			else:
				sample = line[0]
				chrom = line[1]
				pos = line[2]
				val = line[5]
				con = line[6]
				var_number = line[7]
				# print(var_number, var_numbers, len(var_numbers))
				if len(var_numbers) == 0:
					var_numbers.append(var_number)
					chroms.append(chrom)
					positions.append(pos)
					vals.append(val)
					confs.append(con)
					samples.append(sample)
				else:
					if int(var_number) == int(var_numbers[-1]) + 1:
						var_numbers.append(var_number)
						chroms.append(chrom)
						positions.append(pos)
						vals.append(val)
						confs.append(con)
						samples.append(sample)
					else:
						if len(var_numbers) >= 3:
							size = str(int(positions[-1]) - int(positions[0]))
							# print(var_numbers)
							print(samples[0], chroms[0], positions[0], positions[-1], size, len(var_numbers), vals[0], confs[0])
							line_out = [samples[0], chroms[0], positions[0], positions[-1], str(size), str(len(var_numbers)), vals[0], confs[0]]
							out_fh.write(delim.join(line_out) + '\n')
						var_numbers, chroms, positions, vals, confs, samples = [var_number], [chrom], [pos], [val], [con], [sample]
		##print last line
		size = str(int(positions[-1]) - int(positions[0]))
		print(samples[0], chroms[0], positions[0], positions[-1], size, len(var_numbers), vals[0], confs[0])
		line_out = [samples[0], chroms[0], positions[0], positions[-1], str(size), str(len(var_numbers)), vals[0], confs[0], samples[0]]
		out_fh.write(delim.join(line_out) + '\n')





##run methods
# gt_files = ['glass_grc_genotyping_1_181029_FinalReport.txt', 'glass_grc_genotyping_2_190730_FinalReport.txt', 
# 		'glass_grc_genotyping_3_190813_FinalReport.txt', 'glass_grc_genotyping_4_200116_FinalReport.txt']
# gt_files = ['glass_grc_genotyping_1_181029_FinalReport.txt']
# gt_files = ['glass_grc_genotyping_5_200219_FinalReport.txt']
gt_files = ['glass_grc_genotyping_1_181029_FinalReport.txt', 'glass_grc_genotyping_2_190730_FinalReport.txt', 
		'glass_grc_genotyping_3_190813_FinalReport.txt', 'glass_grc_genotyping_4_200116_FinalReport.txt',
		'glass_grc_genotyping_5_200219_FinalReport.txt']

# gt_files = ['test_FinalReport.txt']

##working dir
# working_dir = '/home/atimms/ngs_data/misc/kim_genotyping_0320'
# os.chdir(working_dir)

##make baf and cnv files
# for gt_file in gt_files:
# 	reformat_gt_file(gt_file)

##filter the cnv file and make table of when >=3 probes with conf >=35
# for gt_file in gt_files:
# 	cnv_file = gt_file.rsplit('_', 1)[0] + '.cnv.txt'
# 	cnv_summary_file = gt_file.rsplit('_', 1)[0] + '.cnv_summary.xls'
# 	get_cnv_region_by_sample(cnv_file, cnv_summary_file)

##redo 1021 with previous 5 and 1 new, also extract lrr info
working_dir = '/archive/millen_k/kims_data/kim_bdrlgenotypes'
os.chdir(working_dir)
gt_files = ['glass_grc_genotyping_1_181029_FinalReport.txt', 'glass_grc_genotyping_2_190730_FinalReport.txt', 
		'glass_grc_genotyping_3_190813_FinalReport.txt', 'glass_grc_genotyping_4_200116_FinalReport.txt',
		'glass_grc_genotyping_5_200219_FinalReport.txt', 'glass_grc_genotyping_6_211010_FinalReport.txt']
##make baf and cnv files
for gt_file in gt_files:
	reformat_gt_file(gt_file)




