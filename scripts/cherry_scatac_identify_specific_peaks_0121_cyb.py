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
def transpose_matrix(infile, outfile):
	data = pd.read_csv(infile, index_col=0)
	ts_data = data.transpose()
	ts_data.to_csv(outfile, sep='\t')


def peaks_for_cell_class(infile, cc_dict, nums_req, out_prefix):
	outfile = out_prefix + 'gt' + str(nums_req[0]) + '_lt' + str(nums_req[1]) + '.bed'
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		line_count = 0
		for line in in_fh:
			line_count += 1
			line = line.rstrip().split(delim)
			if line_count == 1:
				cell_classes = line[1:]
			else:
				peak = line[0]
				counts = [float(i) for i in line[1:]]
				##for each different cell class combination
				for cc in cc_dict:
					counts_for_cc = []
					counts_not_for_cc = []
					for cell_class in cell_classes:
						cell_class_count = counts[cell_classes.index(cell_class)]
						if cell_class != 'Mature Ganglions':
							if cell_class in cc_dict[cc]:
								counts_for_cc.append(cell_class_count)
							else:
								counts_not_for_cc.append(cell_class_count)
					# print(peak, cc, cell_classes, counts)
					# print(counts_for_cc)
					# print(counts_not_for_cc)
					if min(counts_for_cc) >= nums_req[0] and max(counts_not_for_cc) <= nums_req[1]:
						peak_list = peak.replace(':', '-').split('-')
						out_fh.write(delim.join(peak_list + [cc]) + '\n')



				# 
				# 	cc_req_counts = []
				# 	cc_not_req_counts = []
				# 	for ind_cc in cc_dict[cc]:
				# 		ind_cc_count = counts[cell_classes.index(ind_cc)]
				# 		cc_req_counts.append(ind_cc_count)
				# 		counts.pop(cell_classes.index(ind_cc))



##run methods
working_dir = '/active/cherry_t/OrgManuscript_SingleCell_Data/human_scATAC/archr_analysis/all_third'
os.chdir(working_dir)

##params
heatmap_matrix = 'heatmap_cell_class_marker_peaks.matrix.csv'
heatmap_matrix_ts = 'heatmap_cell_class_marker_peaks.matrix.transposed.txt'
cell_class_dict = {'rods': ['Mature Rods', 'Developing Rods'], 'cone': ['Mature Cones', 'Developing Cones'],
		'horizontals': ['Mature Horizontals', 'Developing Horizontals'], 'bipolar': ['Mature Bipolars', 'Developing Bipolars'],
		'muller': ['Mature Mullers'], 'amacrine': ['Mature Amacrines', 'Developing Amacrines'], 'ganglion': ['Developing Ganglions']}
thresholds = [1, 0.5]
cc_specific_peak_prefix = 'cell_class_peaks.'



##transpose the matrix file so we can investigate
# transpose_matrix(heatmap_matrix, heatmap_matrix_ts)

# heatmap_matrix_ts = 'temp.txt'
##get peaks for each cell class
peaks_for_cell_class(heatmap_matrix_ts, cell_class_dict, thresholds, cc_specific_peak_prefix)






