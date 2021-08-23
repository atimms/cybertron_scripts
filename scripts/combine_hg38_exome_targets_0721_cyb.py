#!/usr/bin/env python
import subprocess
import os

##parameters
delim = '\t'

bigbedtobed = '/home/atimms/programs/bigBedToBed'
bedtools = '/home/atimms/programs/bedtools-2.29.2/bin/bedtools'
bt_genome = '/home/atimms/ngs_data/references/10x/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.chrom_sizes'

def convert_bigbed_bed(big_beds, bed_suffix):
	for big_bed in big_beds:
		temp_bed = big_bed.rsplit('.', 1)[0] + '_temp.bed'
		temp2_bed = big_bed.rsplit('.', 1)[0] + '_temp2.bed'
		final_bed = big_bed.rsplit('.', 1)[0] + '_padded.bed'
		run_bigbedtobed = subprocess.Popen([bigbedtobed, big_bed, temp_bed])
		run_bigbedtobed.wait()
		with open(temp2_bed, "w") as t2_fh:  
			bedtools_merge = subprocess.Popen([bedtools, 'merge', '-i', temp_bed], stdout=t2_fh)
			bedtools_merge.wait()
		with open(final_bed, "w") as out_fh:  
			bt_slop = subprocess.Popen([bedtools, 'slop', '-i', temp2_bed, '-g', bt_genome, '-b', '50'], stdout=out_fh)
			bt_slop.wait()

def cat_sort_merge(big_beds, bed_suffix, comb_bed):
	beds = [b.rsplit('.', 1)[0] + '_padded.bed' for b in big_beds]
	comined_bed_t1 = comb_bed.rsplit('.', 1)[0] + '_temp1.bed'
	comined_bed_t2 = comb_bed.rsplit('.', 1)[0] + '_temp2.bed'
	with open(comined_bed_t1, 'w') as t1_fh:
		cat_files = subprocess.Popen(['cat'] + beds, stdout=t1_fh)
		cat_files.wait()
	with open(comined_bed_t2, 'w') as t2_fh:
		bedtools_sort = subprocess.Popen([bedtools, 'sort', '-i', comined_bed_t1], stdout=t2_fh)
		bedtools_sort.wait()
	with open(comb_bed, 'w') as comb_fh:
		bedtools_merge = subprocess.Popen([bedtools, 'merge', '-i', comined_bed_t2], stdout=comb_fh)
		bedtools_merge.wait()

##run methods
working_dir = '/home/atimms/ngs_data/references/exome_beds_hg38/'
os.chdir(working_dir)

big_bed_files = ['KAPA_HyperExome_hg38_capture_targets.bb', 'KAPA_HyperExome_hg38_primary_targets.bb', 'S04380110_Covered.bb', 
		'S04380110_Regions.bb', 'S04380219_Covered.bb', 'S04380219_Regions.bb', 'S07084713_Covered.bb', 'S07084713_Regions.bb', 
		'S07604514_Covered.bb', 'S07604514_Regions.bb', 'S07604624_Covered.bb', 'S07604624_Regions.bb', 'S07604715_Covered.bb', 
		'S07604715_Regions.bb', 'S30409818_Covered.bb', 'S30409818_Regions.bb', 'S31285117_Covered.bb', 'S31285117_Regions.bb', 
		'SeqCap_EZ_MedExome_hg38_capture_targets.bb', 'SeqCap_EZ_MedExome_hg38_empirical_targets.bb', 
		'SeqCap_EZ_MedExomePlusMito_hg38_capture_targets.bb', 'SeqCap_EZ_MedExomePlusMito_hg38_empirical_targets.bb', 
		'Twist_ComprehensiveExome_targets_hg38.bb', 'Twist_Exome_RefSeq_targets_hg38.bb', 'Twist_Exome_Target_hg38.bb', 
		'xgen-exome-research-panel-v2-probes-hg38.bb', 'xgen-exome-research-panel-v2-targets-hg38.bb']
bed_file_suffix = '_padded.bed'
combined_bed = 'hg38_targets_combined_padded_0721.bed'

convert_bigbed_bed(big_bed_files, bed_file_suffix)

cat_sort_merge(big_bed_files, bed_file_suffix, combined_bed)


