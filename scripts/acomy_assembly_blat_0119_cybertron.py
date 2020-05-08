#!/usr/bin/env python
import subprocess
import re
import os

##note
'''
module load local_python/3.6.4
module load biobuilds/2017.11
'''
##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/misc/acomy_assembly_0119/seq_query'
os.chdir(working_dir)
blat = 'blat'
acomy_fa = '/home/atimms/ngs_data/misc/acomy_assembly_0119/ref_files/2018_10_15_gap_filled_Acomys_v1_soft_masked.fa'
mm10_fa = '/home/atimms/ngs_data/references/mm10/mm10.fa'
def convert_dna_to_exons(fa_file, results_fa):
	seq_string = ''
	with open (fa_file, 'r') as fasta:
		for line in fasta:
			if line[0] != '>':
				line = line.strip('\n')
				seq_string = seq_string + line
	#print seq_string
	#replace lowercase letters with '#'
	p = re.compile('[a-z]+')
	new_seq = p.sub('#', seq_string)
	##add '#' to start if not there and remove 
	if new_seq[0] != '#':
		new_seq = '#' + new_seq
	if new_seq[-1] == '#':
		new_seq = new_seq[:-1]
	##make exons
	with open (results_fa, 'w') as final:
		exon_count = 0
		for l in new_seq:
			if l == '#':
				exon_count += 1
				if exon_count == 1:
					final.write('>exon' + str(exon_count) + '\n')
				else:
					final.write('\n' + '>exon' + str(exon_count) + '\n')
			else:
				final.write(l)		

##run blat
def run_blat(genome, query, type):
	outfile = query.split('.')[0] + '.pslx'
	outfile_x = query.split('.')[0] + '_x.pslx'
	if type == 'rna':
		con_ann = subprocess.Popen([blat, genome, query, '-q=rna', '-t=dna', '-out=pslx', outfile])
		con_ann = subprocess.Popen([blat, genome, query, '-q=rnax', '-t=dnax', '-out=pslx', outfile_x])
		con_ann.wait()
	elif type == 'dna':
		con_ann = subprocess.Popen([blat, genome, query, '-q=dna', '-t=dna', '-out=pslx', outfile])
		con_ann = subprocess.Popen([blat, genome, query, '-q=dnax', '-t=dnax', '-out=pslx', outfile_x])
		con_ann.wait()
	
def convert_bed_to_exon_seq(in_bed, out_fa):
	bt_gf = subprocess.Popen(['bedtools', 'getfasta', '-fi', mm10_fa, '-bed', in_bed, '-fo', out_fa])
	bt_gf.wait()

def get_exon_seq_from_ucsc_dna(gene):
	dna_seq_file = gene + '.gene_region.fa'
	exon_seq_file = gene + '.exon_seq.fa'
	convert_dna_to_exons(dna_seq_file, exon_seq_file)
	run_blat(acomy_fa, exon_seq_file, 'dna')

def get_exon_seq_from_coding_bed(gene, bed_file):
	exon_seq_file = gene + '.exon_seq.fa'
	convert_bed_to_exon_seq(bed_file, exon_seq_file)
	run_blat(acomy_fa, exon_seq_file, 'dna')

def run_blat_against_fa(query_fa, compare_type):
	run_blat(acomy_fa, query_fa, compare_type)

##extract exons (uppercase) from gene region file, and then blat the exons
# convert_dna_to_exons(gene_region_fa, gene_exons_fa)
# run_blat(acomy_fa, gene_exons_fa, 'dna')


##run methods
'''
##get sequence from exons
##get genomic sequence from ucsc dna tools (exons as caps)
gene = 'Cdh11'
get_exon_seq_from_ucsc_dna(gene)
'''

##get sequence from bed file
##do to ucsc, table browser, get refgene bed, coding exons only
'''
gene = 'Tlr4'
bed = 'mm10_tlr4.bed'
get_exon_seq_from_coding_bed(gene, bed)
'''
##get sequence from exons
##get genomic sequence from ucsc dna tools (exons as caps)
# gene = 'rosa26'
# get_exon_seq_from_ucsc_dna(gene)


##get sequence from bed file
##go to ucsc, table browser, get refgene bed, coding exons only
'''
gene = 'Faah'
bed = 'mm10_Faah.bed'
get_exon_seq_from_coding_bed(gene, bed)
'''


##get sequence from fa file
'''
fa_file = 'hamster_pro_2900_4020.fa'
run_blat_against_fa(fa_file, 'rna')
'''

##get sequence from bed file
##do to ucsc, table browser, get refgene bed, coding exons only
'''
gene = 'Cttn'
bed = 'mm10_Cttn.bed'
get_exon_seq_from_coding_bed(gene, bed)
'''
##get sequence from fa file
# '''
gene = 'mouse_insulin'
get_exon_seq_from_ucsc_dna(gene)
# '''


