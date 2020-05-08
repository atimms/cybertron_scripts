#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import shutil

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/exomes/working/kim_0220'
os.chdir(working_dir)

##methods
def combine_case_std_anal_add_counts(cases, controls, infile_suffix, out_file):
	case_dict, control_dict = {}, {}
	for control in controls:
		ctl_infile = control + infile_suffix
		with open(ctl_infile, 'r') as ctin_fh:
			lc1 = 0
			for line in ctin_fh:
				lc1 += 1
				if lc1 > 1:
					line = line.rstrip().split(delim)
					var = '_'.join(line[:5])
					if var in control_dict:
						if control not in control_dict[var]:
							control_dict[var].append(control)
					else:
						control_dict[var] = [control]
	with open(out_file, 'w') as out_fh:
		for case in cases:
			case_infile = case + infile_suffix
			fc = 0
			with open(case_infile, 'r') as cain_fh:
				fc += 1
				lc2 = 0
				for line in cain_fh:
					lc2 += 1
					line = line.rstrip().split(delim)
					if lc2 == 1:
						if fc == 1:
							header = ['sample', 'ctl_count', 'ctl_ids'] + line
							out_fh.write(delim.join(header) + '\n')
					else:
						var2 = '_'.join(line[:5])
						if var2 in control_dict:
							ctl_count = len(control_dict[var2])
							ctl_ids = ','.join(control_dict[var2])
						else:
							ctl_count = 0
							ctl_ids = ''
						line_out = [case, str(ctl_count), ctl_ids] + line
						out_fh.write(delim.join(line_out) + '\n')



##run methods
autism_samples = ['B03_12', 'B05_13', 'B40_17', 'N04_15', 'NBB_13646', 'NBB_5294', 'NBB_5419', 'NBB_5531', 
		'NBB_5565', 'NBB_5841', 'NBB_5864', 'NBB_5978', 'NBB_6033']
control_samples = ['B28_17', 'B36_17', 'M23_17', 'NBB_1347', 'NBB_1714', 'NBB_4599', 'NBB_4787', 'NBB_4916', 
		'NBB_5242', 'NBB_5334', 'NBB_5936', 'U30_17', 'U38_18']
anal_file_suffix = '.std_analysis.xls'
combined_file = 'combined_cases' + anal_file_suffix



combine_case_std_anal_add_counts(autism_samples, control_samples, anal_file_suffix, combined_file)