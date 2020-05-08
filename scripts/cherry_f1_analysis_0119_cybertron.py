#!/usr/bin/env python
import sys
import subprocess
import os
import glob


##note
'''
analysis of f1

load modules:
module load biobuilds/2017.11
module load homer/4.9.1
module load local_python/3.6.4

for modmap need old python i.e.
module load local_python/2.7.14


'''

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/misc/cherry_f1_hybrid/spreb_atac_0119'
os.chdir(working_dir)
##programs
vcf2bed = '/home/atimms/programs/bedops_linux_x86_64-v2.4.35/vcf2bed'


def homer_find_motifs(motif_file, out_pre_suff, fa):
	#scanMotifGenomeWide.pl pu1.motif mm9 -bed > pu1.sites.mm9.bed
	out_bed = out_pre_suff[0] + fa.split('.')[0] + out_pre_suff[1]
	with open(out_bed, "w") as out_fh: 
		hom_bt_intersect = subprocess.Popen(['scanMotifGenomeWide.pl', motif_file, fa, '-bed'], stdout=out_fh)
		hom_bt_intersect.wait()


def get_passed_snps(in_vcf, out_vcf):
	##correct filtering??
	bcftools_filter = subprocess.Popen(['bcftools', 'view', '-f', 'PASS', '-o', out_vcf, '-O', 'z', in_vcf])
	bcftools_filter.wait()



def convert_vcf_to_bed(in_vcf_gz, bed):
	## vcf to bed  "/opt/bedops-2.4.14/bin/vcf2bed < ${describer}.vcf > ${describer}_pre.bed"
	##vcf to bed
	# in_vcf = in_vcf_gz.rsplit('.', 1)[0]
	# print in_vcf
	# print([vcf2bed, '<', in_vcf])
	with open(bed, 'w') as bed_fh:
		cat_vcf = subprocess.Popen(['zcat', in_vcf_gz], stdout=subprocess.PIPE)
		vcf2bed_bo = subprocess.Popen([vcf2bed, '-'], stdin=cat_vcf.stdout, stdout=bed_fh)
		vcf2bed_bo.wait()

def split_and_format_motif_bed(in_bed, out_file_suffix):
	motif_dict = {}
	with open(in_bed, "r") as in_fh:
		for line in in_fh:
			line = line.rstrip().split(delim)
			tf = line[3].split("(")[0]
			info = line[:3] + [line[5], line[4]]
			if tf in motif_dict:
				motif_dict[tf].append(info)
			else:
				motif_dict[tf] = [info]
	for trans in motif_dict:
		print(trans)
		outfile = trans + out_file_suffix
		with open(outfile, "w") as out_fh:
			# out_fh.write(delim.join(['Chr', 'Start', 'Stop', 'Strand', 'PwmScore']) + '\n')
			for inf in motif_dict[trans]:
				out_fh.write(delim.join(inf) + '\n')
		# run_gzip = subprocess.Popen(['gzip', outfile])
		# run_gzip.wait()

def get_motifs_in_peaks_bed(in_bed, out_bed, peaks_bed):
	with open(out_bed, "w") as out_fh: 
		hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', in_bed, '-b', peaks_bed, '-wa'], stdout=out_fh)
		hom_bt_intersect.wait()

def modmap_coords_to_mm10_spret(infile, outfile, changed_file, mod_file):
	#modmap -b -d '\t' SPRET_EiJ_mgp_v5_snps_indels_pass_copy9.mod MA0477.1_consensus_iupac_SPRET.bed MA0477.1_consensus_iupac_SPRET_mod_start.bed 1,2
	temp_inbed = infile.split('.')[0] + '.temp_in.bed'
	temp_start_bed = infile.split('.')[0] + '.temp_start.bed'
	temp_end_bed = infile.split('.')[0] + '.temp_end.bed'
	# '''
	##modify infile
	with open(temp_inbed, "w") as out_fh, open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.split(delim)
			line_out = [line[0].split(' ')[0]] + line[1:]
			out_fh.write(delim.join(line_out))
	##get start and end coordinates 
	run_modmap_start = subprocess.Popen(['modmap', '-b', '-d', "\t", mod_file, temp_inbed, temp_start_bed, '1,2'])
	run_modmap_start.wait()
	run_modmap_start = subprocess.Popen(['modmap', '-b', '-d', "\t", mod_file, temp_inbed, temp_end_bed, '1,3'])
	run_modmap_start.wait()
	# '''
	##combine start and end
	bed_dict = {}
	with open(temp_start_bed, "r") as start_fh:
		lc = 0
		for line in start_fh:
			line = line.split(delim)
			# start_coord = str(abs(int(line[1])))
			# line_out = [line[0]] + [start_coord] + line[2:]
			lc += 1
			bed_dict[lc] = line
	with open(temp_end_bed, "r") as end_fh:
		lc = 0
		for line in end_fh:
			line = line.split(delim)
			end_coord = line[2]
			lc += 1
			bed_dict[lc][2] = end_coord
	with open(outfile, "w") as out_fh, open(changed_file, "w") as ch_fh:
		changed_count, total_count = 0,0
		for l in bed_dict:
			line_out = bed_dict[l]
			start = line_out[1]
			end = line_out[2]
			total_count += 1
			# print(line_out, start, end)
			if start.startswith('-') or end.startswith('-'):
				changed_count += 1
				start_coord = abs(int(start))
				end_coord = abs(int(end))
				m_size = end_coord - start_coord
				if m_size > 0:
					new_lo = [line_out[0]] + [str(start_coord), str(end_coord)] + line_out[3:]
					ch_fh.write(delim.join(new_lo))
				else:
					new_lo = [line_out[0]] + [str(start_coord -1), str(end_coord)] + line_out[3:]
					ch_fh.write(delim.join(new_lo))	
				print(line_out, new_lo)
			else:
				out_fh.write(delim.join(line_out))
	print(changed_count, total_count)

def modmap_coords_to_mm10(infile, outfile, mod_file):
	#modmap -b -d '\t' SPRET_EiJ_mgp_v5_snps_indels_pass_copy9.mod MA0477.1_consensus_iupac_SPRET.bed MA0477.1_consensus_iupac_SPRET_mod_start.bed 1,2
	temp_inbed = infile.split('.')[0] + '.temp_in.bed'
	temp_start_bed = infile.split('.')[0] + '.temp_start.bed'
	temp_end_bed = infile.split('.')[0] + '.temp_end.bed'
	# '''
	##modify infile
	with open(temp_inbed, "w") as out_fh, open(infile, "r") as in_fh:
		for line in in_fh:
			line = line.split(delim)
			line_out = [line[0].split(' ')[0]] + line[1:]
			out_fh.write(delim.join(line_out))
	##get start and end coordinates 
	run_modmap_start = subprocess.Popen(['modmap', '-b', '-d', "\t", mod_file, temp_inbed, temp_start_bed, '1,2'])
	run_modmap_start.wait()
	run_modmap_start = subprocess.Popen(['modmap', '-b', '-d', "\t", mod_file, temp_inbed, temp_end_bed, '1,3'])
	run_modmap_start.wait()
	# '''
	##combine start and end
	bed_dict = {}
	with open(temp_start_bed, "r") as start_fh:
		lc = 0
		for line in start_fh:
			line = line.split(delim)
			lc += 1
			bed_dict[lc] = line
	with open(temp_end_bed, "r") as end_fh:
		lc = 0
		for line in end_fh:
			line = line.split(delim)
			end_coord = line[2]
			lc += 1
			bed_dict[lc][2] = end_coord
	with open(outfile, "w") as out_fh:
		for l in bed_dict:
			out_fh.write(delim.join(bed_dict[l]))

def add_chr_to_name_bed(infile, outfile):
	with open(outfile, "w") as out_fh, open(infile, "r") as in_fh:
		for line in in_fh:
			out_fh.write('chr' + line)

def intersect_peaks_files_with_snps(in_bed, out_bed, snp_bed):
	with open(out_bed, "w") as out_fh: 
		hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', in_bed, '-b', snp_bed, '-c'], stdout=out_fh)
		hom_bt_intersect.wait()

def combine_motifs(tf_names, spret_in_suffix, bl6_in_suffix, outfile_suffix):
	for tf in tf_names:
		print(tf)
		spret_bed = tf + spret_in_suffix
		bl6_bed = tf + bl6_in_suffix
		comb_bed = tf + outfile_suffix
		infiles = [spret_bed, bl6_bed]
		with open(comb_bed, "w") as out_fh:
			# sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n'] + infiles, stdout=out_fh)
			# sort_file.wait()

			sort_file = subprocess.Popen(['sort', '-k1,1', '-k2,2n'] + infiles, stdout=subprocess.PIPE)
			bt_merge = subprocess.Popen(['bedtools', 'merge', '-i', '-', '-c', '4,5,4,5', '-o', 'distinct,distinct,count_distinct,count_distinct'], stdin=sort_file.stdout, stdout=out_fh)
			bt_merge.wait()

def remove_duplicates_from_peaks_files(infiles):
	for infile in infiles:
		outfile = infile.split('.')[0] + '.dup_rm.bed'
		lines_seen = set()
		with open(outfile, "w") as out_fh, open(infile, "r") as infh:
			for line in infh:
				if line not in lines_seen:
					lines_seen.add(line)
					line = line.rstrip().split(delim)
					out_fh.write(delim.join(line[:4]) + '\n')

def intersect_peaks_files_with_motifs_get_counts(peak_beds, tf_names, motif_suffixes, changed_motif_suffix, snp_bed, outfile_suffix):
	##step one - get motifs with and without snps
	bed_files = []
	for tf in tf_names:
		for motif_suffix in motif_suffixes:
			print(tf, motif_suffix)
			motif_bed = tf + motif_suffix
			motif_with_snp = motif_bed.rsplit('.',1)[0] + '.snp.bed'
			motif_no_snp = motif_bed.rsplit('.',1)[0] + '.no_snp.bed'
			bed_files.extend([motif_with_snp, motif_no_snp])
			##get motifs with snps
			'''
			with open(motif_with_snp, "w") as out_fh:
				bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', motif_bed, '-b', snp_bed, '-u'], stdout=out_fh)
				bt_intersect.wait()
			with open(motif_no_snp, "w") as out_fh:
				bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-a', motif_bed, '-b', snp_bed, '-v'], stdout=out_fh)
				bt_intersect.wait()
			'''
	##add changed bed files to comparison
	changed_files = [tfac + changed_motif_suffix for tfac in tf_names]
	bed_files.extend(changed_files)
	# print(bed_files)
	##compare all peak files with the motif files
	for peak_bed in peak_beds:
		count_dict = { i : [] for i in bed_files}
		print(count_dict)
		##interect against all bed files
		peak_counts = peak_bed.rsplit('.',1)[0] + '.counts.bed'
		sum_file = peak_bed.rsplit('.',1)[0] + outfile_suffix
		# with open(peak_counts, "w") as out_fh:
		# 	hom_bt_intersect = subprocess.Popen(['bedtools', 'intersect', '-wao', '-a', peak_bed, '-b'] + bed_files + ['-filenames'], stdout=out_fh)
		# 	hom_bt_intersect.wait()
		with open(peak_counts, "r") as in_fh:
			for line in in_fh:
				line = line.split(delim)
				bed = line[4]
				peak_id = line[3]
				##if peak has a motif
				if bed != '.':
					##add to dict if not already seen
					if peak_id not in count_dict[bed]:
						count_dict[bed].append(peak_id)
		print(count_dict)
		with open(sum_file, "w") as out_fh:
			for b in count_dict:
				out_fh.write(delim.join([b.rsplit('.',1)[0], str(len(count_dict[b]))]) + '\n')





##run methods

##parameters

##working
original_peak_beds = ['f1_b6_enriched_diffpeaks.bed', 'f1_similarpeaks.bed', 'f1_spret_enriched_diffpeaks.bed', 
		'intersect_p0_f1_b6_diffpeaks.bed', 'intersect_p0_f1_similarpeaks.bed', 'intersect_p0_f1_spret_diffpeaks.bed', 
		'parent_b6_enriched_diffpeaks.bed', 'parent_similarpeaks.bed', 'parent_spret_enriched_diffpeaks.bed']
homer_db_motifs = 'homer_motifs_1018.txt'
homer_motifs_bed_pre_suff = ['homer_motifs_', '_0119.bed']
tfs = ['BORIS', 'CRE', 'CRX', 'MafA', 'Mef2d', 'Otx2', 'RORgt']
##genomes - dl from http://www.csbio.unc.edu/CCstatus/index.py?run=Pseudo
bl6_fa = 'C57BL6J_b38.fa'
spret_fa = 'SPRETEiJ_b38_f.fa'
bl6_mod = 'C57BL6J_b38.mod'
spret_mod = 'SPRETEiJ_b38_f.mod'
spret_snp_vcf = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz'
spret_passed_snp_vcf = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.passed.vcf.gz'
spret_motifs = 'homer_motifs_SPRETEiJ_b38_f_0119.bed'
b6_motifs = 'homer_motifs_C57BL6J_b38_0119.bed'
mm10_motifs = 'homer_motifs_mm10_0119.bed'
spret_motifs_mm10 = 'homer_motifs_SPRETEiJ_mm10_0119.bed'
b6_motifs_mm10 = 'homer_motifs_C57BL6J_mm10_0119.bed'
combined_motifs_mm10 = 'homer_motifs_combined_mm10_0119.bed'
spret_motifs_changed_mm10 = 'homer_motifs_SPRETEiJ_mm10_changed_0119.bed'
spret_passed_snp_bed = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.passed.bed'

##files to compare
adjusted_peak_beds = ['f1_b6_enriched_diffpeaks.dup_rm.bed', 'f1_similarpeaks.dup_rm.bed', 'f1_spret_enriched_diffpeaks.dup_rm.bed', 
		'intersect_p0_f1_b6_diffpeaks.dup_rm.bed', 'intersect_p0_f1_similarpeaks.dup_rm.bed', 'intersect_p0_f1_spret_diffpeaks.dup_rm.bed', 
		'parent_b6_enriched_diffpeaks.dup_rm.bed', 'parent_similarpeaks.dup_rm.bed', 'parent_spret_enriched_diffpeaks.dup_rm.bed']
# peak_beds = ['test.bed', 'f1_b6_enriched_diffpeaks.bed']
spret_passed_snp_bed_adj = 'SPRET_EiJ.mgp.v5.snps.dbSNP142.passed.added_chr.bed'
spret_motif_bed_suffix = '.homer_SPRETEiJ_mm10_0119.bed'
b6_motif_bed_suffix = '.homer_C57BL6J_mm10_0119.bed'
combined_motif_bed_suffix = '.homer_combined_mm10_0119.bed'
spret_changed_motif_bed_suffix = '.homer_SPRETEiJ_mm10_changed_0119.bed'
motif_file_suffixes_to_get_snps = [spret_motif_bed_suffix, b6_motif_bed_suffix, combined_motif_bed_suffix]
peaks_motif_counts_suffix = '.motif_counts_0119.xls'

###input files
'''
need 3 things:
differential and shared peaks - have from tim, generated from previous work
snps - have vcf from sanger, need to keep passed and make into bed
motifs - get peaks from homer (have from atac analysis 1018) then map to bl6 and speb genome
'''

##motifs
##get motifs for pseudogenome - motif from homer db, fa's from mod file guys
# homer_find_motifs(homer_db_motifs, homer_motifs_bed_pre_suff, bl6_fa)
# homer_find_motifs(homer_db_motifs, homer_motifs_bed_pre_suff, spret_fa)
##not used
# homer_find_motifs(homer_db_motifs, homer_motifs_bed_pre_suff, 'mm10')
##convert motifs to mm10 coordinates
# modmap_coords_to_mm10_spret(spret_motifs, spret_motifs_mm10, spret_motifs_changed_mm10,spret_mod)
# modmap_coords_to_mm10(b6_motifs, b6_motifs_mm10,bl6_mod)
##split into individual motif files
# split_and_format_motif_bed(spret_motifs_mm10, spret_motif_bed_suffix)
# split_and_format_motif_bed(b6_motifs_mm10, b6_motif_bed_suffix)
# split_and_format_motif_bed(spret_motifs_changed_mm10, spret_changed_motif_bed_suffix)
##combine motif files
# combine_motifs(tfs,spret_motif_bed_suffix, b6_motif_bed_suffix, combined_motif_bed_suffix)

##snps
##clean snps and make bed
# get_passed_snps(spret_snp_vcf, spret_passed_snp_vcf)
# convert_vcf_to_bed(spret_passed_snp_vcf, spret_passed_snp_bed)
# add_chr_to_name_bed(spret_passed_snp_bed, spret_passed_snp_bed_adj)

##peaks
##remove duplicate lines from peak files and format i.e. make all 4 columns wide
# remove_duplicates_from_peaks_files(original_peak_beds)

##intersect snps with motifs/peaks
##experiment 1. snps in peaks
# for pbed in adjusted_peak_beds:
# 	count_bed = pbed.split('.')[0] + '.snp_count.bed'
# 	intersect_peaks_files_with_snps(pbed, count_bed, spret_passed_snp_bed_adj)

##experiment 2. intersect peaks with motifs with/without a snp and the changed motifs
intersect_peaks_files_with_motifs_get_counts(adjusted_peak_beds, tfs, motif_file_suffixes_to_get_snps, spret_changed_motif_bed_suffix, spret_passed_snp_bed_adj, peaks_motif_counts_suffix)







