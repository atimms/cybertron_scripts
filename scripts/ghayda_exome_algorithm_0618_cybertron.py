#!/usr/bin/env python
import subprocess
import os

'''
setting up environment
1. interactive qsub: qsub -Iq longq -l mem=50gb,ncpus=2
2. load anaconda python: module load local_python/3.6.4 
3a. make environment, only needs to be done once: conda create --name ghayda_exome_algorithm_0618
3b. load environment: source activate ghayda_exome_algorithm_0618

'''
run_mh_single1 = subprocess.Popen(['which', 'python'])
run_mh_single1.wait()

##parameters
delim = '\t'
working_dir = '/home/atimms/ngs_data/misc/ghayda_exome_algorithm_0618'
os.chdir(working_dir)


##methods
def check_for_vcf(peds, vcf_suffix):
	for ped in peds:
		vcf = ped + vcf_suffix
		if os.path.exists(vcf):
			print(ped, 'yes')
		else:
			print(ped, 'no')


def get_hits_and_non_hits(all_vars, vars_of_interest, out_prefix):
	##make ped and var dict for vars we want
	ped_dict = {}
	lc = 0
	with open(vars_of_interest, "r") as vi_fh:
		for line in vi_fh:
			lc += 1
			if lc > 1:
				line = line.split(delim)
				ped_id = line[1]
				gene = line[7]
				analysis = line[8]
				if ped_id in ped_dict:
					ped_dict[ped_id][0].append(gene)
					ped_dict[ped_id][3].append(analysis)
				else:
					ped_dict[ped_id] = [[gene],0,0, [analysis]]
	##check ped dict
	# for i in ped_dict:
	# 	print(i, ped_dict[i][0])
	##open all vars files and get hits and no hits
	hits_file = out_prefix + '.hits.xls'
	nothits_file = out_prefix + '.notahit.xls'
	lc = 0
	with open(all_vars, "r") as av_fh, open(hits_file, "w") as hit_fh, open(nothits_file, "w") as nah_fh:
		for line in av_fh:
			lc += 1
			if lc == 1:
				##write headers to out files
				hit_fh.write(line)
				nah_fh.write(line)
			else:
				line = line.rstrip().split(delim)
				ped = line[20]
				gene = line[5]
				anal = line[24]
				# print(ped, gene, anal)
				if ped in ped_dict:
					# print(ped, gene, anal)
					# print(anal, ped_dict[ped][3])
					if gene in ped_dict[ped][0] and anal in ped_dict[ped][3]: 
					# if gene in ped_dict[ped][0]: 
						# print(type(anal), type(ped_dict[ped][3][0]))
						hit_fh.write(delim.join(line + ['\n']))
						ped_dict[ped][1] += 1
					else:
						nah_fh.write(delim.join(line+ ['\n']))
						ped_dict[ped][2] += 1
	##check ped dict
	for i in ped_dict:
		print(i, ped_dict[i][1], ped_dict[i][2], ped_dict[i][3])

###check we have all pedigrees


























##run methods

##parameters
# pedigrees = ['LR01-314', 'LR03-214', 'LR05-162', 'LR06-105', 'LR07-200', 'LR12-112', 'LR05-054', 'LR07-054', 'LR03-340', 'LR04-199', 
# 		'LR12-068', 'LR03-304', 'LR04-222', 'LR14-075', 'LR01-173', 'LR09-416', 'LR12-018', 'LP98-078', 'LR00-016', 'LR01-306', 'LR01-242', 
# 		'LR12-049', 'LR12-475', 'LP97-141', 'LR01-279', 'LR03-039', 'LR07-037', 'LR04-315', 'LR12-068', 'LR02-064', 'LR07-041', 'LR03-260', 
# 		'LR11-424', 'LR11-331', 'LR17-332', 'LR06-300', 'LR06-029a1', 'LR11-144', 'LR06-099', 'LR03-077', 'LP88-024', 'LR09-203', 'LR11-330', 
# 		'LR12-032', 'LR11-033', 'LR09-149', 'LR04-208', 'LR12-101', 'LR10-102', 'LR14-226', 'LR04-376', 'LR13-037', 'LR13-175', 'LR11-042', 
# 		'LR06-397', 'LR13-003', 'LR09-248a2', 'LR07-150', 'LR11-367', 'LR03-056a1', 'LR08-023', 'LR05-160', 'LR02-264a1', 'LR09-280', 
# 		'LR09-227', 'LR15-011', 'LR09-006', 'LR10-193', 'LR12-225', 'LR06-231', 'LR14-268', 'LP97-105', 'LR07-227', 'LR13-013', 'LR08-390', 
# 		'LR10-046', 'LR11-001', 'LR13-013', 'LR01-265a1', 'LR12-053', 'LR04-017', 'LP98-095', 'LR09-266', 'LR03-332', 'LR13-270', 'LR11-054', 
# 		'LR11-241', 'LR00-144', 'LR10-199', 'LR02-046', 'LP97-114', 'LR13-278', 'LR05-398', 'LR13-103', 'LR05-007', 'LR03-130', 'LR12-264', 
# 		'LR12-463', 'LR13-300', 'LR06-085', 'LR07-016', 'LR08-323', 'LP89-036', 'LR00-220', 'LR11-328', 'LP92-083', 'LR04-233', 'LR06-339', 
# 		'LR11-233', 'LR10-230', 'LR03-206', 'LR08-002', 'LR06-207', 'LR11-189', 'LR12-436a2', 'LP99-206', 'LR01-338', 'LR15-155', 'LR13-160', 
# 		'LR12-434', 'LR02-027', 'LR12-379a1', 'LR12-323', 'LR08-257', 'LR09-141', 'LR12-156', 'LR04-289', 'LR13-150', 'LR13-061', 'LR03-170', 
# 		'LR15-236', 'LR13-200', 'LR12-304', 'LR14-048', 'LR12-339', 'LR15-061', 'LR11-124', 'LR11-124', 'LR09-023', 'LR11-003', 'LR12-492', 
# 		'LR05-120', 'LR13-068', 'LR12-343', 'LR07-031', 'LR10-246', 'LR04-101', 'LR14-166', 'LR14-255']

pedigrees = ['LP88-024', 'LP89-036', 'LP97-114', 'LP98-095', 'LP99-206', 'LR00-144', 'LR00-220', 'LR01-265a1', 'LR01-338', 'LR02-027', 'LR02-046', 'LR03-056a1', 'LR03-077', 'LR03-130', 'LR03-170', 'LR03-206', 'LR03-332', 'LR04-017', 'LR04-101', 'LR04-233', 'LR04-289', 'LR04-376', 'LR05-007', 'LR05-120', 'LR05-160', 'LR05-398', 'LR06-029a1', 'LR06-085', 'LR06-099', 'LR06-207', 'LR06-231', 'LR06-300', 'LR06-397', 'LR07-150', 'LR07-227', 'LR08-002', 'LR08-023', 'LR08-257', 'LR08-323', 'LR08-390', 'LR09-023', 'LR09-141', 'LR09-149', 'LR09-203', 'LR09-227', 'LR09-248a2', 'LR09-266', 'LR09-280', 'LR10-046', 'LR10-102', 'LR10-199', 'LR10-230', 'LR11-001', 'LR11-003', 'LR11-033', 'LR11-042', 'LR11-054', 'LR11-124', 'LR11-144', 'LR11-189', 'LR11-241', 'LR11-328', 'LR11-330', 'LR11-331', 'LR11-367', 'LR12-032', 'LR12-049_2', 'LR12-053', 'LR12-068', 'LR12-101', 'LR12-112a1', 'LR12-156', 'LR12-264', 'LR12-323', 'LR12-339', 'LR12-343', 'LR12-379a1', 'LR12-434', 'LR12-436a2', 'LR12-463', 'LR13-003', 'LR13-013', 'LR13-037', 'LR13-068', 'LR13-103', 'LR13-200', 'LR13-278', 'LR14-048', 'LR14-268', 'LR17-332', 'LR15-004', 'LR12-354', 'LR10-064', '1090', '1533', 'LR12-225', 'LR15-011', 'LR13-279', 'LR12-230', 'LR13-270', 'LR09-266', 'LR16-034', 'LR13-300', 'LR13-199', 'LR06-339', 'LR11-233', 'LR15-081', 'LR15-155', 'LR10-084', 'LR13-160', 'LR13-321', 'LR16-170', 'LR15-097', 'LR17-038', 'LR16-079', '1270', 'LR05-069', 'LR13-150', 'LR15-003', 'LR13-061', 'LR15-236', 'LR12-304', 'LR15-061', 'LR15-338', '1588', 'LR12-492', 'LR15-272', 'LR10-246']

##parameters
project = 'exome_algorithm_0618'
all_var_file = 'all_exome_data_std_pipeline.1217.xls'
peds_of_intrest_var_file = project + '.data_std_pipeline.1217.xls'
vcf_suffix = '.intersected_vcfs/0002.vcf.gz'
vars_to_analyze = 'exome_algorithm_vars_to_analyze_0318.txt'



##check for data
##check for vcf file and then combine and annotate
check_for_vcf(pedigrees, vcf_suffix)

##seperate into hits and non hits
# get_hits_and_non_hits(all_var_file, vars_to_analyze, project)













