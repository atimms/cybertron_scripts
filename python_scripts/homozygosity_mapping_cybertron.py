#!/usr/bin/env python
import shutil
# import pybedtools as pbt
import os 
import subprocess
import operator

##programs
bedtools = '/home/atimms/programs/bedtools2-master/bin/bedtools'
# bedtools = 'bedtools'
#parameters
delim = '\t'
acceptable_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
				'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
				'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 
				'chr23', 'chr24', 'chr25', 'chrx', 'chry', 'chrX', 'chrY', 
				'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
				'13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', 
				'24', '25', 'x', 'y', 'X', 'Y'] 


##edit bedfile to remove non standard chromosomes
def remove_nonStd_chr(infile, outfile):
	fh = open(infile, 'r')
	outfh = open(outfile, 'w')
	for line in fh:
		line = line.split(delim)
		if line[0] in acceptable_chr:
			outfh.write(delim.join(line))
	fh.close()
	outfh.close()


##divide genome into windows
def make_windows(working_dir, genome_fai, window_size, step_size):
	os.chdir(working_dir)
	genome = genome_fai.split('/')[-1].split('.')[0]
	genome_window_bed = genome + '_' + str(window_size / 1000) + 'kb_' + str(step_size / 1000) + 'kb' + '.bed'
	print 'bed file %s created from file: %s with window size %r and step size %r' % (genome_window_bed, genome_fai, window_size, step_size)
	print ''
	##run bedtools makewindows without pybedtools installed
	command = bedtools + ' makewindows -g ' + genome_fai + ' -w ' + str(window_size) + ' -s ' + str(step_size) + ' > ' + 'bed.temp'
	subprocess.call(command, shell=True)
	##filters chromosomes if required
	remove_nonStd_chr('bed.temp', genome_window_bed)
	##if you don't require any filtering
	#shutil.copy('bed.temp', genome_window_bed)
	return genome_window_bed

##make bed file from ann.txt contains: chr,start,end,zygosity,novel_allele_freq
def make_bed_from_ann(working_dir, var_caller, filename, zygosity_col, last_col):
	os.chdir(working_dir)
	fh = open(filename, 'r')
	outfile = filename.split(".", 1)[0] + '.bed'
	outfh = open(outfile, 'w')
	line_count = 0
	kept_count = 0
	for line in fh:
		line_count += 1
		if line_count > 1:
			line = line.split(delim)
			info = line[last_col - 1].split(':')
			#print info
			if var_caller == 'gatk' or var_caller == 'gatk_hc':
				if len(info[1].split(',')) != 2: #gives snps with mutiple alleles naf of 0
					non_ref_allele, total_alleles = 0,0
				elif info[2] == '.':
					non_ref_allele, total_alleles = 0,0
				else:
					non_ref_allele = float(info[1].split(',')[1])
					total_alleles = int(info[2])


			elif var_caller == 'samtools':
				if len(info[1].split(',')) != 3: #gives snps with mutiple alleles naf of 0
					non_ref_allele, total_alleles = 0,0
				else:
					non_ref_allele = float(info[3])
					total_alleles = int(info[2])

					
			elif var_caller == 'freebayes':
				if ',' in info[4]: #gives snps with mutiple alleles naf of 0
					non_ref_allele, total_alleles = 0,0
				else:
					non_ref_allele = float(info[4])
					total_alleles = int(info[1])

			else:
				print "doofus: %s isn't a variant caller" % var_caller
			if total_alleles == 0:
				naf = 0
			else:
				naf = non_ref_allele / total_alleles
			
			#print var_caller, info, non_ref_allele, total_alleles, naf
			
			if naf > 0.1:            #prevents very low frequency SNPs and those we've filtered
				zygosity = line [zygosity_col -1]
				line_out = [line [0], str(int(line [1]) - 1), line [2], zygosity, str(naf), '\n']
				outfh.write(delim.join(line_out))
				kept_count += 1
				#print 'kept'
			#else:
				#print info, naf
	print 'bed file %s created from %s' % (outfile, filename)
	print '%r variants checked and %r kept' % (line_count -1, kept_count)
	print ''
	fh.close()
	outfh.close()

##take bedfile with zygosity in 4th col and write a het and hom temp file
def split_het_hom(working_dir, filename):
	fh = open(filename, 'r')
	het_fh = open('het.temp', 'w')
	hom_fh = open('hom.temp', 'w')
	line_count = 0
	het_count = 0
	hom_count = 0
	for line in fh:
		line_count += 1 
		line = line.split(delim)
		if line[3] == 'het':
			het_fh.write(delim.join(line))
			het_count += 1
		elif line[3] == 'hom':
			hom_fh.write(delim.join(line))
			hom_count += 1
		else:
			print 'issue with snp on line', line_count
	fh.close()
	het_fh.close()
	hom_fh.close()
	print 'for file %s we have %r snps, of which %r are het and %r are hom' % (filename, line_count, het_count, hom_count)
	print ''

##combines multiple beds
##'chr:start-end' must match and be sorted in same way
##will only add last col, so need to adjust if need more
def combine_multiple_bed(bed_list, outfile):
	results = []
	no_of_bed = 0
	##loop through bedfiles and add to lists
	for bed in bed_list:
		fh = open(bed, 'r')
		no_of_bed += 1
		line_number = -1
		for line in fh:
			line_number += 1
			line = line.strip('\n').split(delim)
			##if first file add all columns
			if no_of_bed == 1:
				results.append(line)
			##if not first file check 'chr:start-end' and add last column
			if no_of_bed >1:
				if results[line_number][:3] == line[:3]:
					results[line_number].append(line[-1])
		fh.close()
	##write to outfile
	outfh = open(outfile, 'w')
	for i in range(len(results)):
		outfh.write(delim.join(results[i]) + '\n')
	outfh.close()
	print "combined the %r bedfiles %s" % (no_of_bed, ','.join(bed_list))
	print "checked %r and printed %r lines" % (line_number + 1 , len(results))

##make het and hom count bedgraphs for specified window size
def het_and_hom_bed(working_dir, genome_and_window, filename):
	os.chdir(working_dir)
	print 'making het and hom bedfiles for:', filename
	split_het_hom(working_dir, filename)
	sample = filename.split(".", 1)[0]
	window_bed = genome_and_window + '.bed'
	hom_count_bed = sample + '_' + genome_and_window + '_hom_count.bedgraph'
	het_count_bed = sample + '_' + genome_and_window + '_het_count.bedgraph'
	##pbt version
	# a = pbt.BedTool(genome_and_window + '.bed')
	# b = pbt.BedTool('hom.temp')
	# c = pbt.BedTool('het.temp')
	# a_with_b = a.intersect(b, c=True).moveto(hom_count_bed)
	# a_with_c = a.intersect(c, c=True).moveto(het_count_bed)
	##bedtools version
	#bedtools intersect -a A.bed -b B.bed -c
	with open(hom_count_bed, "w") as hom_fh: 
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', window_bed, '-b', 'hom.temp', '-c'], stdout=hom_fh)
		hom_bt_intersect.wait()
	with open(het_count_bed, "w") as het_fh: 
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', window_bed, '-b', 'het.temp', '-c'], stdout=het_fh)
		hom_bt_intersect.wait()
	return hom_count_bed, het_count_bed

##count hom and het, and hom percentage in specified window
def count_and_percentage(working_dir, genome_and_window, filename):
	os.chdir(working_dir)
	print 'calculating hom count, het count and hom percentage for:', filename
	hom_count_bed, het_count_bed = het_and_hom_bed(working_dir, genome_and_window, filename)
	combine_multiple_bed([hom_count_bed, het_count_bed], 'hom_het.temp')
	fh = open('hom_het.temp', 'r')
	sample = filename.split(".", 1)[0]
	outfile = sample + '_' + genome_and_window + '_hom_percentage.bedgraph'
	outfh = open(outfile, 'w')
	for line in fh:
		line = line.strip('\n').split(delim)
		if line[3] == '0':
			hom_percentage = 0
		else:
			hom_percentage = float(line[3]) / (float(line[3]) + float(line[4]))
		outfh.write(delim.join(line[:3]) + delim + str(hom_percentage) + '\n')
	fh.close()
	outfh.close()
	print 'generated files:'
	print 'hom_count:', hom_count_bed
	print 'het_count:', het_count_bed
	print 'hom_percentage:', outfile


##novel allele frequency within windows
def naf_in_window(working_dir, genome_and_window, bed_file):
	os.chdir(working_dir)
	print 'calculating novel allele freq within windows for file:', bed_file
	sample = bed_file.split(".", 1)[0]
	##bedtools interect window bed and sample bed 
	# a = pbt.BedTool(genome_and_window + '.bed')
	# b = pbt.BedTool(bed_file)
	# a_with_b = a.intersect(b, wa=True, wb=True).moveto(sample + '_naf.temp')
	##bedtools version
	window_bed = genome_and_window + '.bed'
	#bedtools intersect -a A.bed -b B.bed -c
	with open(sample + '_naf.temp', "w") as naf_fh: 
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', window_bed, '-b', bed_file, '-wa', '-wb'], stdout=naf_fh)
		hom_bt_intersect.wait()
	##make dictionary: key = chr:start-end and value a list of naf values
	naf_dict = {}
	fh = open(sample + '_naf.temp', 'r')
	for line in fh:
		line = line.strip('\n').split(delim)
		chr_start_end = delim.join(line[:3])
		naf_value = float(line[-1])
		if naf_dict.has_key(chr_start_end):
			naf_dict[chr_start_end].append(naf_value)
		else:
			naf_dict[chr_start_end] = []
			naf_dict[chr_start_end].append(naf_value)
	fh.close()
	
	##make list of lists containing chr,start,end,average_naf, snp#
	list_of_lists = []
	for cse in naf_dict.keys():
		split_cse = cse.split(delim)
		average_naf = [sum(naf_dict[cse])/float(len(naf_dict[cse]))]
		total_snps = [len(naf_dict[cse])]
		outlist = split_cse + average_naf + total_snps
		outlist[1] = int(outlist[1])
		outlist[2] = int(outlist[2])
		list_of_lists.append(outlist)
	list_of_lists.sort(key = operator.itemgetter(0, 1))
	
	##print bed file with average_naf
	outfh = open(sample + '_' + genome_and_window + '_naf.bedgraph', 'w')
	for line in list_of_lists:
		for i in range(len(line)):
			line[i] = str(line[i])
		outfh.write(delim.join(line) + '\n')
	outfh.close()

##total snp number within windows (for testing)
def total_snp_in_window(working_dir, genome_and_window, bed_file):
	os.chdir(working_dir)
	print 'caluclating total snp number within windows for file:', bed_file
	sample = bed_file.split(".", 1)[0]
	##bedtools interect window bed and sample bed 
	# a = pbt.BedTool(genome_and_window + '.bed')
	# b = pbt.BedTool(bed_file)
	# a_with_b = a.intersect(b, wa=True, wb=True).moveto(sample + '_naf.temp')
	##bedtools version
	window_bed = genome_and_window + '.bed'
	#bedtools intersect -a A.bed -b B.bed -c
	with open(sample + '_naf.temp', "w") as naf_fh: 
		hom_bt_intersect = subprocess.Popen([bedtools, 'intersect', '-a', window_bed, '-b', bed_file, '-wa', '-wb'], stdout=naf_fh)
		hom_bt_intersect.wait()

	##make dictionary: key = chr:start-end and value a list of naf values
	naf_dict = {}
	fh = open(sample + '_naf.temp', 'r')
	for line in fh:
		line = line.strip('\n').split(delim)
		chr_start_end = delim.join(line[:3])
		naf_value = float(line[-1])
		if naf_dict.has_key(chr_start_end):
			naf_dict[chr_start_end].append(naf_value)
		else:
			naf_dict[chr_start_end] = []
			naf_dict[chr_start_end].append(naf_value)
	fh.close()
	
	##make list of lists containing chr,start,end,total_snp, snp#
	list_of_lists = []
	for cse in naf_dict.keys():
		split_cse = cse.split(delim)
		total_snps = [len(naf_dict[cse])]
		outlist = split_cse + total_snps
		outlist[1] = int(outlist[1])
		outlist[2] = int(outlist[2])
		list_of_lists.append(outlist)
	list_of_lists.sort(key = operator.itemgetter(0, 1))
	
	##print bed file with total snp number
	outfh = open(sample + '_' + genome_and_window + '_snp_number.bedgraph', 'w')
	for line in list_of_lists:
		for i in range(len(line)):
			line[i] = str(line[i])
		outfh.write(delim.join(line) + '\n')
	outfh.close()


##method to combine hom mapping data i.e. bedgraphs
##suply working dir, file prefix name, and genome and window information
def combine_bedgraphs_for_r(working_directory, prefix, genome_window):
	chr_wanted = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
		'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 
		'chr21', 'chr22', 'chr23', 'chr24', 'chr25','1', '2', '3', '4', '5', '6', '7', '8', '9', 
		'10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25'] 
	snp_count = prefix + '_' + genome_window + '_snp_number.bedgraph'
	het_count = prefix + '_' + genome_window + '_het_count.bedgraph'
	hom_count = prefix + '_' + genome_window + '_hom_count.bedgraph'
	hom_perc = prefix + '_' + genome_window + '_hom_percentage.bedgraph'
	naf_ph = prefix + '_' + genome_window + '_naf.bedgraph'
	combined_r_file = prefix + '_' + genome_window + '_combined_hom_mapping.txt'
	with open(combined_r_file, "w") as r_file:
		r_file.write(delim.join(['chr', 'start', 'end', 'value', 'chromosome', 'analysis', '\n']))
		with open(hom_count, "r") as homc:
			for line in homc:
				line = line.strip('\n').split(delim)
				line = [line[0].replace('chr','')] + line[1:] #removes chr from start of line
				##remove chromosomes we're not interested in
				if line[0] in chr_wanted:
					start = int(line[1])
					end = int(line[2])
					midpoint = start + ((end - start) / 2)
					r_file.write(delim.join(line[:4] + [str(midpoint), 'hom count', '\n']))
		with open(hom_perc, "r") as homp:
			for line in homp:
				line = line.strip('\n').split(delim)
				line = [line[0].replace('chr','')] + line[1:] #removes chr from start of line
				if line[0] in chr_wanted:
					start = int(line[1])
					end = int(line[2])
					midpoint = start + ((end - start) / 2)
					r_file.write(delim.join(line[:4] + [str(midpoint), 'hom percentage', '\n']))
		with open(naf_ph, "r") as naf:
			for line in naf:
				line = line.strip('\n').split(delim)
				line = [line[0].replace('chr','')] + line[1:] #removes chr from start of line
				if line[0] in chr_wanted:
					start = int(line[1])
					end = int(line[2])
					midpoint = start + ((end - start) / 2)
					r_file.write(delim.join(line[:4] + [str(midpoint), 'average naf', '\n']))