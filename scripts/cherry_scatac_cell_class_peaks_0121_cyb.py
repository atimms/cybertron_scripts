#!/usr/bin/env python
import subprocess
import os
import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# import glob



'''
info:
##need to use correct env
qsub -Iq cdbrmq -l mem=100gb,ncpus=5 -P 19833a08-f6fb-4bea-8526-8a79069da878
conda activate pd_np_plt_etc
'''
##parameters
delim = '\t'


##methods
def filter_annotated_for_human_data(infiles):
	for infile in infiles:
		outfile = infile.rsplit('.', 3)[0] + '.human.txt'
		with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
			line_count = 0
			for line in in_fh:
				line_count += 1
				line = line.rstrip().split(delim)
				if line_count == 1:
					cell_class_names = line[19:]
					cell_class_names = [i.split('.')[0] + '.' + i.split('.')[1] for i in cell_class_names]
					human_cc_names = [i.split('.')[1] for i in cell_class_names if i.startswith('human')]
					# print(cell_class_names, human_cc_names)
					out_fh.write(delim.join(line[:19] + human_cc_names) + '\n')
				else:
					cell_classes = line[19:]
					human_indices = [i for i, x in enumerate(cell_class_names) if x.startswith('human')]
					human_values_wanted = [cell_classes[i] for i in human_indices]
					out_fh.write(delim.join(line[:19] + human_values_wanted) + '\n')

def peaks_for_cell_class(infiles, cc_dict, nums_req, out_prefix):
	for infile in infiles:
		outfile = out_prefix + 'fc' + str(nums_req[0]) + '_gt' + str(nums_req[1]) + '.' + infile.split('.')[2] + '.' + infile.split('.')[3] + '.txt'
		outbed = out_prefix + 'fc' + str(nums_req[0]) + '_gt' + str(nums_req[1]) + '.' + infile.split('.')[2] + '.' + infile.split('.')[3] + '.bed'
		with open(outfile, "w") as out_fh, open(infile, "r") as in_fh, open(outbed, "w") as outb_fh:
			line_count = 0
			for line in in_fh:
				line_count += 1
				line = line.rstrip().split(delim)
				if line_count == 1:
					cell_classes = line[19:]
					header = line + ['cell class name', 'min tag count', 'fold change (min cell class / max others)']
					out_fh.write(delim.join(header) + '\n')
				else:
					peak = line[0]
					counts = [float(i) for i in line[19:]]
					##for each different cell class combination
					for cc in cc_dict:
						counts_for_cc = []
						counts_not_for_cc = []
						for cell_class in cell_classes:
							cell_class_count = counts[cell_classes.index(cell_class)]
							if cell_class in cc_dict[cc]:
								counts_for_cc.append(cell_class_count)
							else:
								counts_not_for_cc.append(cell_class_count)
						# print(peak, cc, cell_classes, counts)
						# print(counts_for_cc)
						# print(counts_not_for_cc)
						# print(type(counts_not_for_cc[0]))
						# print(cc, len(counts_for_cc), len(counts_not_for_cc), len(cell_classes))
						##filter for specific peaks
						min_cc = min(counts_for_cc)
						##if cell class peak all over certain threshold
						if min_cc >= nums_req[1]:
							max_ncc = max(counts_not_for_cc)
							cc_fc = min_cc/max_ncc
							if cc_fc >= nums_req[0]:
								line_out = line + [cc, str(min_cc), str(cc_fc)]
								line_out_bed = line[1:4] + [cc, str(min_cc), str(cc_fc)]
								out_fh.write(delim.join(line_out) + '\n')
								outb_fh.write(delim.join(line_out_bed) + '\n')



##run methods
working_dir = '/home/atimms/ngs_data/misc/cherry_sc_org_project_20/cell_class_specific_peaks_0121'
os.chdir(working_dir)

##params
annotated_cc_peaks = ['annotated.collapsed_cc_only.100d.2000size.macs2_q0.01_summits.txt', 
	'annotated.collapsed_cc_only.100d.500size.macs2_q0.01_summits.txt', 
	'annotated.collapsed_cc_only.500d.2000size.macs2_q0.01_summits.txt', 
	'annotated.collapsed_cc_only.500d.500size.macs2_q0.01_summits.txt']
human_cc_peaks = ['annotated.collapsed_cc_only.100d.2000size.human.txt', 'annotated.collapsed_cc_only.100d.500size.human.txt', 
		'annotated.collapsed_cc_only.500d.2000size.human.txt', 'annotated.collapsed_cc_only.500d.500size.human.txt']
# annotated_cc_peaks = ['annotated.collapsed_cc_only.100d.2000size.macs2_q0.01_summits.temp.txt']

cell_class_dict = {'rods': ['Mature_Rods', 'Developing_Rods'], 'cones': ['Mature_Cones', 'Developing_Cones'],
		'horizontals': ['Mature_Horizontals', 'Developing_Horizontals'], 'bipolar': ['Mature_Bipolars', 'Developing_Bipolars'],
		'mullers': ['Mature_Mullers'], 'amacrines': ['Mature_Amacrines', 'Developing_Amacrines'], 'ganglions': ['Developing_Ganglions']}

human_cc_peak_prefix = 'human.cell_class_peaks.'
##thresholds... fold change required, and must be above value
threshold_values = [[10, 1], [10, 10], [5, 1], [5, 10]]
# thresholds = [10, 1]

##filter annotated files for just human data
# filter_annotated_for_human_data(annotated_cc_peaks)

##get peaks for each cell class
# human_cc_peaks = ['annotated.collapsed_cc_only.100d.2000size.human.txt']
for thresholds in threshold_values:
	peaks_for_cell_class(human_cc_peaks, cell_class_dict, thresholds, human_cc_peak_prefix)






