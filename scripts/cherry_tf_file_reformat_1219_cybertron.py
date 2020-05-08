#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import itertools
import csv
import pandas as pd

'''
module load python/3.6.6
'''

##parameters
delim = ','
working_dir = '/home/atimms/ngs_data/misc/cherry_tf_file_reformat_1219'
os.chdir(working_dir)


infiles = ['180513_TF_ATAC_boxplot.csv', '180513_TF_H3K27ac_boxplot.csv', '180517_TF_RNA_boxplot.csv']
# print(infiles)

# infiles = ['test.csv']
def transpose_file(infile, outfile):
	with open(infile, 'r') as file, open(outfile, 'w') as final:
		#make list for each line
		a = [x.split(',') for x in file]
		#print a
		#turn into array
		b = np.array(a)
		#print b
		#transpose array
		c = b.T
		#convert array back to list of list
		d = c.tolist()
		#print d
		#convert each list back to 
		for e in d:
			line = delim.join(e) + '\n'
			final.write(line)

for infile in infiles:
	print(infile)
	temp_file = infile.split('.')[0] + '_temp.csv'
	outfile = infile.split('.')[0] + '_reformatted.csv'
	tf_dict = {}
	lc = 0
	with open(infile, "U") as in_fh:
		for line in in_fh:
			# print(line)
			lc += 1
			if lc >1:
				line = line.rstrip().split(delim)
				tf = line[1]
				number = line[0]
				if tf in tf_dict:
					tf_dict[tf].append(number)
				else:
					tf_dict[tf] = [number]
	##add dummy values to dict keys so all max length
	lens = []
	for t in tf_dict:
		print(t, len(tf_dict[t]))
		lens.append(len(tf_dict[t]))
	max_lens = max(lens)
	for t in tf_dict:
		if not len(tf_dict[t]) == max_lens:
			tf_dict[t].extend(['']*(max_lens-len(tf_dict[t])))
	# print(tf_dict.keys())
	# csv_columns = tf_dict.keys()
	# with open(outfile, 'w') as csvfile:
	# 	writer = csv.writer(csvfile)
	# 	writer.writerow(csv_columns)
	# 	# print(list(zip(*tf_dict.values())))
	# 	print(list(itertools.zip_longest(tf_dict.values())))
	# 	# writer.writerows(zip(*tf_dict.values()))
	# 	writer.writerows(itertools.zip_longest(tf_dict.values(), fillvalue=''))
	# with open(temp_file, 'w') as temp_fh:
	# 	for t in tf_dict:
	# 		temp_fh.write(','.join([t] + tf_dict[t]) + '\n')
	# transpose_file(temp_file, outfile)

	tf_df = pd.DataFrame(tf_dict)
	# print(tf_df)
	tf_df.to_csv(outfile, index=False)