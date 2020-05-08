#!/usr/bin/env python

''' Methods for filtering annotated.txt files:'''

import csv
import operator
import os
import re

delim = '\t'

##check if 's' is a number and return boolean
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False
	
##dictionary used to convert operator to .operator command
def operator_dict(op):
	ops = {"==": operator.eq, "!=": operator.ne,">=": operator.ge,"<=": operator.le, ">": operator.gt,"<": operator.lt,"in": operator.contains}
	return ops[op]


def filter(working_dir, test_type, input_file, output_file, col_number, op, comparison):
	"""
	method for filtering files by operator and a value, need to provide:
	working directory,
	test type i,e, 'and' or 'or',
	input filename,
	output filename.
	parameters for tests i.e;
	column number
	operator (must be in dict)
	comparison (what to compare value in col number with)
	"""
	os.chdir(working_dir)
	print "filtering variants from the file '%s' in an '%s' fashion" % (input_file,test_type)
	##variables
	row_number, variants_kept = 0,0                             #keep track of number of tests and variants kept
	delim = '\t'
	##open files for reading and writing
	infile = open(input_file, 'rb')
	outfile = open(output_file, 'wb')
	writer = csv.writer(outfile, delimiter=delim)
	##check row by row
	for row in csv.reader(infile, delimiter=delim):
		row_number += 1
		if row_number == 1:
			writer.writerow(row)                                  #writes header to out file
			header = row                                          #writes header to header
		else:
			tests_correct = 0                                     #keep track of number of correct tests and reset
			for i in range(len(comparison)):
				op_func = operator_dict(op[i])                    #queries operator dictionary so can perform comparisons
				col = col_number[i] -1                            #adjust column number so starts at 1 not 0

				if row_number == 2:                               #print info on test completed
					print "test %r: '%s %s' in column '%s'" % (i, op[i], comparison[i], header[col])

				if is_number(comparison[i]):                      #if comparison is a number
					# print row[col]
					if row[col] == 'NA' or row[col] == '.' or row[col] == '':       #if using na as 'empty' if can't compare number to text so say test is true
						tests_correct += 1
					elif row[col] == '.':                         # in gatk sometime coverage = 0 but printed as .
						pass
					elif op_func(float(row[col]), comparison[i]): #proper comparison
						tests_correct += 1
				else:
					if op_func(row[col], comparison[i]):          #if comparison is not a number
						tests_correct += 1
			##looks if a 'and' a 'or' or a 'once' test type, and writes into new file
			if test_type == "or":                                 #if any test correct per row
				if tests_correct > 0:
					writer.writerow(row)
					variants_kept += 1
			elif test_type == "and":                              #if all tests correct per row
				if tests_correct == len(comparison):
					writer.writerow(row)
					variants_kept += 1
			elif test_type == "once":                             #if one test correct per row
				if tests_correct == 1:
					writer.writerow(row)
					variants_kept += 1
			else:
				print "incorrect test type!!"
	##print overall numbers
	print "%s variants kept out of %s checked" % (variants_kept, row_number - 1)
	print "filtered variants in file:", output_file 
	print "" #empty line
	infile.close()
	outfile.close()





def multiple_variants_in_gene(working_dir, input_file, output_file, genes_col):
	"""
	method for looking for compound hets i.e. multiple variants within a gene:
	-superseded by gene_in_multiple_files
	need to provide:
	working directory,
	single input file, either annotated.txt or derivative of.
	name of output file
	column gene name is located (column count start at 1 not 0)
	"""
	##change to working dir
	os.chdir(working_dir)
	
	##set up variables
	row_number, exonic_var = 0,0                                 #variables to store totals
	gene_check,multiple_var = [],[]                              #set up lists, to query and results
	var_gene_col = genes_col - 1                                 #col with gene name 
	var_location_col = genes_col - 2                             #col with gene location i.e. exonic
	delim = '\t'
	
	##open files for reading and writing
	infile = open(input_file, 'rb')
	outfile = open(output_file, 'wb')
	writer = csv.writer(outfile, delimiter=delim)
	
	##make list of all gene name, query and if another occurs put in multiple_var list
	print "checking file '%s' for genes with more than one variant" % input_file
	print "and writing to file '%s' " % output_file
	
	##check row by row
	for row in csv.reader(infile, delimiter=delim):
		#if 'exonic' in row[var_location_col] or 'splicing' in row[var_location_col]: #if exonic with this gene db
		if row[var_location_col] == 'exonic' or row[var_location_col] == 'splicing': #if exonic with this gene db
			genes = re.sub(r'\(.+?\)', '', row[var_gene_col])    #get field with gene names and remove () and anything within
			gene_list = genes.split(',')                         #splits by , and puts into list
			for gene in gene_list:
				if gene in gene_check:                           #if gene name seen before
					multiple_var.append(gene)                    #put in results list
				else:
					gene_check.append(gene)                      #if not put in query list
	infile.close()
	
	##reopen infile and check row by row
	infile = open(input_file, 'rb')
	for row in csv.reader(infile, delimiter=delim):
		row_number += 1
		if row_number == 1:
			writer.writerow(row)                                   #writes header to out file
		else:
			#if 'exonic' in row[var_location_col] or 'splicing' in row[var_location_col]: #if exonic with this gene db
			if row[var_location_col] == 'exonic' or row[var_location_col] == 'splicing': #if exonic with this gene db
				exonic_var += 1                                    #keeps count of variants checked
				genes = re.sub(r'\(.+?\)', '', row[var_gene_col])  #get field with gene names and remove () and anything within
				gene_list = genes.split(',')                       #splits by , and puts into list
				for gene in gene_list:
					if gene in multiple_var:                       #if gene name in results list
						writer.writerow(row)                       #write to file
						break
						
	print "rows checked:", row_number
	print "exonic variants checked:", exonic_var
	print "genes with multiple variants:", len(set(multiple_var))  #number of unique gene names in results list
	print ""
	infile.close()
	outfile.close()





def var_in_multiple_files(working_dir, input_filelist, output_file, sample_pos):
	'''
	method looking if the same variant is in multiple files
	need to provide:
	working directory,
	list of input files, either annotated.txt or derivative of.
	name of output file
	location of sample name in output file when split with a '.'
	'''
	##change to working dir
	os.chdir(working_dir)
	
	##set up variables
	to_check = {}                       #set up lists, to query and results
	results = []
	delim = '\t'
	row_number,file_number = 0,0                          #variables to store totals
	
	##open file for writing
	outfile = open(output_file, 'wb')
	writer = csv.writer(outfile, delimiter=delim)
	
	##go through file(s) if var not seen before but if query list, 
	##if seen multiple times put in results
	for file in input_filelist:
		row_number = 0
		infile = open(file, 'rb')                         #open file
		for row in csv.reader(infile, delimiter=delim):
			row_number += 1
			if row_number >1:                             #if not header
				var = ','.join(row[0:5]) #varaint details
				proband = file.split('.')[sample_pos]
				#print var, proband
				if var in to_check and proband not in to_check[var]:                       #if variant seen before
					results.append(var)                   #put in results list
				else:
					to_check[var] = [proband]
		print "for file %s, %s variants checked" % (file, row_number -1)
		infile.close()
	print results
	##go through file(s) if var in results list to in outfile
	for file in input_filelist:
		infile = open(file, 'rb')
		row_number = 0
		file_number += 1
		for row in csv.reader(infile, delimiter=delim):
			row_number += 1
			if row_number == 1 and file_number == 1:                    #if header and 1st file
				writer.writerow(['Sample'] + row)                       #write to file (add 'sample' to start of list)
			else:
				var = ','.join(row[0:5])
				#print var
				if var in results:                                      #if variant in results list
					row = [file.split(".")[sample_pos]] + row                 #change row to include 1st part of file name i.e. before 1st .
					#row = [file.split(".")[1]] + row                 #change row to include 1st part of file name i.e. before 1st .
					writer.writerow(row)                                #write to file 
		infile.close()
	outfile.close()


def gene_in_multiple_files(working_dir, input_filelist, output_file, sample_pos, genes_col):
	"""
	method looking for the same gene in multiple files
	will work for individual files i.e. looking for compound hets
	- set up to look for exonic/splicing variants
	need to provide:
	working directory,
	list of input files, either annotated.txt or derivative of.
	name of output file
	location of sample name in output file when split with a '.'
	column gene name is located (column count start at 1 not 0)
	"""
	##change to working dir
	os.chdir(working_dir)
	
	##set up variables
	row_number,file_number = 0,0                                 #variables to store totals
	gene_check,multiple_var = [],[]                              #set up lists, to query and results
	var_gene_col = genes_col - 1                                 #col with gene name 
	var_location_col = genes_col - 2                             #col with gene location i.e. exonic
	delim = '\t'
	
	##open files for writing
	outfile = open(output_file, 'wb')
	writer = csv.writer(outfile, delimiter=delim)
	
	##make list of all gene name, query and if another occurs put in multiple_var list
	print "checking files '%r' for genes with more than one variant" % input_filelist
	print "and writing to file '%s' " % output_file
	
	##check row by row and update lists
	for file in input_filelist:
		infile = open(file, 'rb')
		for row in csv.reader(infile, delimiter=delim):
			row_number += 1
			if 'exonic' in row[var_location_col] or 'splicing' in row[var_location_col]: #if exonic with this gene db
				genes = re.sub(r'\(.+?\)', '', row[var_gene_col])    #get field with gene names and remove () and anything within
				gene_list = genes.split(',')                         #splits by , and puts into list
				print gene_list
				for gene in gene_list:
					if gene in gene_check:                           #if gene name seen before
						multiple_var.append(gene)                    #put in results list
					else:
						gene_check.append(gene)                      #if not put in query list
		infile.close()
		print "for file %s, %s variants checked" % (file, row_number -1)

	##reopen infile and check row by row if gene in multiple list
	for file in input_filelist:
		infile = open(file, 'rb')
		row_number = 0
		file_number += 1
		for row in csv.reader(infile, delimiter=delim):
			row_number += 1
			if row_number == 1 and file_number == 1:                   #if header and 1st file
				writer.writerow(['Sample'] + row)                      #writes header to out file
			else:
				if 'exonic' in row[var_location_col] or 'splicing' in row[var_location_col]: #if exonic with this gene db
					genes = re.sub(r'\(.+?\)', '', row[var_gene_col])  #get field with gene names and remove () and anything within
					gene_list = genes.split(',')                       #splits by , and puts into list
					for gene in gene_list:
						if gene in multiple_var:                       #if gene name in results list
							row = [file.split(".")[sample_pos]] + row        #adds sample name to first column
							writer.writerow(row)                       #write to file
							break
		infile.close()
	print "genes with multiple variants:", len(set(multiple_var))  #number of unique gene names in results list
	print ""
	outfile.close()
	
def compound_het_from_both_parents(working_dir, input_file, output_file, genes_col, parents_col):
	"""
	method for looking for compound hets and makes sure a het is coming from both parents:
	need to provide:
	working directory,
	single input file, either annotated.txt or derivative of.
	name of output file
	column gene name is located (column count start at 1 not 0)
	column number of both parents (column count start at 1 not 0)
	caveats:
	var must not be het in both parents
	"""
	##change to working dir
	os.chdir(working_dir)
	
	##set up variables
	row_number, exonic_var = 0,0                                 #variables to store totals
	gene_dict = {}                              #set up gene dict
	passed_genes = []
	var_gene_col = genes_col - 1                                 #col with gene name 
	var_location_col = genes_col - 2                             #col with gene location i.e. exonic
	dad_col = parents_col[0] - 1
	mum_col = parents_col[1] - 1
	delim = '\t'
	
	##open files for reading and writing
	infile = open(input_file, 'rb')
	outfile = open(output_file, 'wb')
	writer = csv.writer(outfile, delimiter=delim)
	
	##make list of all gene name, query and if another occurs put in multiple_var list
	print "checking file '%s' for genes with more than one variant, and getting variants from both parents" % input_file
	print "and writing to file '%s' " % output_file
	
	##make dictionary with each gene and where the variant comes from
	for row in csv.reader(infile, delimiter=delim):
		#if 'exonic' in row[var_location_col] or 'splicing' in row[var_location_col]: #if exonic with this gene db
		if row[var_location_col] == 'exonic' or row[var_location_col] == 'splicing': #if exonic with this gene db
			genes = re.sub(r'\(.+?\)', '', row[var_gene_col])    #get field with gene names and remove () and anything within
			gene_list = genes.split(',')                         #splits by , and puts into list
			for gene in gene_list:
				if gene in gene_dict:                           #if gene name seen before
					if row[dad_col] == 'het':
						gene_dict[gene].append('dad')                    #put in results list
					elif row[mum_col] == 'het':
						gene_dict[gene].append('mum')                    #put in results list
					else:
						gene_dict[gene].append('none')                    #put in results list
				else:
					if row[dad_col] == 'het':
						gene_dict[gene] = ['dad']                      #if not put in query list
					elif row[mum_col] == 'het':
						gene_dict[gene] = ['mum']                      #if not put in query list
					else:
						gene_dict[gene] = ['none']                      #if not put in query list
	infile.close()
	##make list of gene dictionary that have 2 or more variants and from mum and dad
	for gene in gene_dict:
		if len(gene_dict[gene]) > 1 and 'dad' in gene_dict[gene] and 'mum' in gene_dict[gene]:
			passed_genes.append(gene)

	##reopen infile and check row by row
	infile = open(input_file, 'rb')
	for row in csv.reader(infile, delimiter=delim):
		row_number += 1
		if row_number == 1:
			writer.writerow(row)                                   #writes header to out file
		else:
			#if 'exonic' in row[var_location_col] or 'splicing' in row[var_location_col]: #if exonic with this gene db
			if row[var_location_col] == 'exonic' or row[var_location_col] == 'splicing': #if exonic with this gene db
				exonic_var += 1                                    #keeps count of variants checked
				genes = re.sub(r'\(.+?\)', '', row[var_gene_col])  #get field with gene names and remove () and anything within
				gene_list = genes.split(',')                       #splits by , and puts into list
				for gene in gene_list:
					if gene in passed_genes:                       #if gene name in results list
						writer.writerow(row)                       #write to file
						break
						
	print "rows checked:", row_number
	print "exonic variants checked:", exonic_var
	print "genes with multiple variants:", len(passed_genes)  #number of unique gene names in results list
	print ""
	infile.close()
	outfile.close()

def gene_pairs_in_multiple_samples(working_dir, sample_list, input_prefix, input_suffix, output_file, genes_col, samples_required):
	"""
	method looking for if a number of individuals share paires of genes
	- set up to look for exonic/splicing variants
	need to provide:
	working directory,
	list of input files, either annotated.txt or derivative of.
	name of output file
	column gene name is located (column count start at 1 not 0)
	number of samples seen to have passed
	"""
	##change to working dir
	os.chdir(working_dir)
	
	##set up variables
	gene_sample_var_dict = {}                                    #set up gene dict
	var_gene_col = genes_col - 1                                 #col with gene name 
	var_location_col = genes_col - 2                             #col with gene location i.e. exonic
	delim = '\t'
	comparisons_made = 0
	passed_test = 0
	passed_list = []

	##fill in dictionary for all genes with samples and with variants (no duplicates)
	for sample in sample_list:
		infile = open(input_prefix + sample + input_suffix, 'rb')
		row_count = 0
		for row in csv.reader(infile, delimiter=delim):
			row_count += 1
			if 'exonic' in row[var_location_col] or 'splicing' in row[var_location_col]: #if exonic with this gene db
				#genes = re.sub(r'\(.+?\)', '', row[var_gene_col])    #get field with gene names and remove () and anything within
				#genes = genes.split(',')                         #splits by , and puts into list
				genes = [row[var_gene_col]]      #don't seperate gene names, so don't get duplicates (cheating)
				for gene in genes:
					if gene in gene_sample_var_dict:
						if sample in gene_sample_var_dict[gene]:
							gene_sample_var_dict[gene][sample].append(row) #if already have gene and sample within gene add new var to var list with double dict
						else:
							gene_sample_var_dict[gene][sample] = [row] #if entry there but not in the same sample, make entry in dict with sample and variant info as first entry in list
							
					else:
						gene_sample_var_dict[gene] = {sample:[row]}  #if no entry for gene make dict with sample name and variant info as first entry in list
		infile.close()

	##get header for outfile
	header_file = input_prefix + sample_list[0] + input_suffix
	with open(header_file, "r") as headfh, open(output_file, "w") as outfh:
		line_count = 0
		for line in headfh:
			line_count += 1
			if line_count == 1:
				outfh.write(delim.join(['Gene Pair','Number of Samples', 'Samples', 'Sample']) + line)
	
	##compare gene in list against each other
	with open(output_file, "a") as outfh:
		for gene in gene_sample_var_dict:
			for other_gene in gene_sample_var_dict:
				samples_shared = 0
				samples_shared_list = []
				samples_shared_vars = []
				if gene != other_gene:
					comparisons_made += 1
					for sample in gene_sample_var_dict[gene]:
						if sample in gene_sample_var_dict[other_gene]:
							samples_shared += 1
							samples_shared_list.append(sample)

				##if passed test i.e. enough individuals had sample in those gene
				if samples_shared >= samples_required:
					passed_test += 1
					#print gene, other_gene
					#print samples_shared, samples_shared_list
					gene_pair = gene + '_' + other_gene
					joined_ss = '_'.join(samples_shared_list)
					for sample in samples_shared_list:
						for var in gene_sample_var_dict[gene][sample]:
							outfh.write(delim.join([gene_pair] + [str(samples_shared)] + [joined_ss] + [sample] + var + ['\n']))
							
						for var in gene_sample_var_dict[other_gene][sample]:
							outfh.write(delim.join([gene_pair] + [str(samples_shared)] + [joined_ss] + [sample] + var + ['\n']))
	print "%s comparsions made and %s gene pairs were in %i or more samples" %  (comparisons_made, passed_test, samples_required)


def variants_pairs_in_multiple_samples(working_dir, sample_list, input_prefix, input_suffix, output_file, samples_required):
	"""
	method looking for if a number of individuals share paires of genes
	- set up to look for exonic/splicing variants
	need to provide:
	working directory,
	list of input files, either annotated.txt or derivative of.
	name of output file
	column gene name is located (column count start at 1 not 0)
	number of samples seen to have passed
	"""
	##change to working dir
	os.chdir(working_dir)
	
	##set up variables
	var_sample_dict = {}                                    #set up gene dict
	delim = '\t'
	comparisons_made = 0
	passed_test = 0
	passed_list = []

	##fill in dictionary for all genes with samples and with variants (no duplicates)
	for sample in sample_list:
		infile = open(input_prefix + sample + input_suffix, 'rb')
		row_count = 0
		
		for row in csv.reader(infile, delimiter=delim):
			row_count += 1
			if row_count >1:
				variant = ','.join(row[0:5] + [row[-1]])
				if variant in var_sample_dict:
					var_sample_dict[variant][sample] = row #if entry there but not in the same sample, make entry in dict with sample and variant info as first entry in list
						
				else:
					var_sample_dict[variant] = {sample:row}  #if no entry for gene make dict with sample name and variant info as first entry in list
		infile.close()

	##get header for outfile
	header_file = input_prefix + sample_list[0] + input_suffix
	with open(header_file, "r") as headfh, open(output_file, "w") as outfh:
		line_count = 0
		for line in headfh:
			line_count += 1
			if line_count == 1:
				outfh.write(delim.join(['Variant Pair','Number of Samples', 'Samples', 'Sample'] + line.split(delim)))
	
	##compare variants in list against each other
	with open(output_file, "a") as outfh:
		for variant in var_sample_dict:
			for other_variant in var_sample_dict:
				samples_shared = 0
				samples_shared_list = []
				samples_shared_vars = []
				
				if variant != other_variant:
					comparisons_made += 1
					for sample in var_sample_dict[variant]:
						if sample in var_sample_dict[other_variant]:
							samples_shared += 1
							samples_shared_list.append(sample)
				not_in_ctl = 0   #for checking if variant in ctls


				##if passed test i.e. enough individuals had sample in those variant
				if samples_shared >= samples_required:
					
					#check in variant not in controls
					if variant.rsplit(',',1)[1] == '.':
						not_in_ctl +=1
					if other_variant.rsplit(',',1)[1] == '.':
						not_in_ctl +=1
					if not_in_ctl >= 1:
						passed_test += 1
						#print variant, other_variant
						#print samples_shared, samples_shared_list
						variant_pair = variant + '_' + other_variant
						joined_ss = '_'.join(samples_shared_list)
						#print samples_shared_list
						for sample in samples_shared_list:
							#print delim.join([variant_pair] + [str(samples_shared)] + [joined_ss] + [sample] + var_sample_dict[variant][sample])
							outfh.write(delim.join([variant_pair] + [str(samples_shared)] + [joined_ss] + [sample] + var_sample_dict[variant][sample] + ['\n']))
							outfh.write(delim.join([variant_pair] + [str(samples_shared)] + [joined_ss] + [sample] + var_sample_dict[other_variant][sample] + ['\n']))
	print "%s comparsions made and %s variant pairs were in %i or more samples" %  (comparisons_made, passed_test, samples_required)


def convert_list_to_dict_with_zero(gene_list):
	##convert list to dict to track var_count
	gene_dict = {}
	for gene in gene_list:
		gene_dict[gene] = 0
	return gene_dict

def filter_var_by_genename(genelist, input_file,output_file, gene_col):
	##convert list to dict to track var_count
	genedict = convert_list_to_dict_with_zero(genelist)
# 	print genedict
# 	print len(genedict), len(genelist)
	with open(output_file, "w") as outfh, open(input_file, "U") as infh:
		line_count = 0
		for line in infh:
			line_count += 1
			if line_count == 1:
				outfh.write(line)
			else:
				line = line.strip('\n').strip('\r').rstrip().split(delim)
				genes = re.sub(r'\(.+?\)', '', line[gene_col])
				genes = genes.split(',')
				for gene in genes:
					if gene in genedict:
						outfh.write(delim.join(line) + '\n')
						genedict[gene] += 1
						break
		#print genedict
		print 'of the variants checked %i are in the genes specified'% (sum(genedict.values()))
		print '\n'
# 		for g in genedict:
# 			print 'in file %s we have %s variants in gene %s'%(output_file, genedict[g], g)

##method to combine files
def combine_ann_txt(sample_list, prefix, suffix, outfile):
	with open(outfile, "w") as final_file:
		sample_count = 0
		for sample in sample_list:
			file = prefix + sample + suffix
			sample_count += 1
			with open(file, "r") as open_file:
				line_count = 0
				if sample_count == 1:
					for line in open_file:
						line_count += 1
						if line_count == 1:
							final_file.write('Proband' + delim + line)
						else:
							final_file.write(sample + delim + line)
				else:
					for line in open_file:
						line_count += 1
						if line_count > 1:
							final_file.write(sample + delim + line)