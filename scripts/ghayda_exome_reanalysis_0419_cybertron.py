#!/usr/bin/env python
import dobyns_gemini_pipeline_cybertron_v11
import subprocess
import glob
import os
'''
##modules 
module load java/1.8.0_121 
module load biobuilds/2016.11
module load vt/0.5772-60f436c3 
module load mono/5.10.1.47
module load Pisces/5.1.6.54
'''

##programs
bcftools = '/home/atimms/programs/bcftools-1.9/bcftools'
plink = '/home/atimms/programs/plink'


##run all other methods
def run_standard_gemini_protocol(working_directory, pedigree_dict):
	for pedigree in pedigree_dict:
		in_vcf = pedigree + '.intersected_vcfs/0002.vcf.gz'
		ped_type = pedigree_dict[pedigree]
		dobyns_gemini_pipeline_cybertron_v11.standard_gemini_protocol(working_directory, pedigree, ped_type)
		plink_relatadness_check(in_vcf, pedigree)
		delete_unwated_files(['*temp*'])


def plink_relatadness_check(vcf, file_prefix):
	##correct filtering??
	# bcftools_filter = subprocess.Popen([bcftools, 'view', '-f', 'PASS', '-i', "QUAL>50 & MIN(INFO/DP)>30", '-V', 'indels', '-o', 'temp_plink.vcf.gz', '-O', 'z', vcf])
	bcftools_filter = subprocess.Popen([bcftools, 'view', '-i', "QUAL>50", '-V', 'indels', '-o', 'temp_plink.vcf.gz', '-O', 'z', vcf])
	bcftools_filter.wait()
	##generate plink file from vcf
	make_plink = subprocess.Popen([plink, '--vcf', 'temp_plink.vcf.gz', '--vcf-filter', '--const-fid', '--out', 'temp.pass_q50_dp50'])
	make_plink.wait()
	##check sex -results in .sexcheck file
	plink_sex = subprocess.Popen([plink, '--bfile', 'temp.pass_q50_dp50', '--check-sex', '--out', file_prefix + '.pass_q50_dp50'])
	plink_sex.wait()
	##ibd check
	plink_ibd = subprocess.Popen([plink,  '--bfile', 'temp.pass_q50_dp50', '--genome', '--out', file_prefix + '.pass_q50_dp50'])
	plink_ibd.wait()

def delete_unwated_files(file_extensions):
	files_to_delete = []
	for ext in file_extensions:
		files = glob.glob(ext)
		files_to_delete.extend(files)
	for f in files_to_delete:
		os.remove(f)

##tests
##parameters and to run
##test dir
'''
# work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/1002506'
# ped_dict = {'LR04-101': 'trio', 'LR13-378': 'trio', 'LR14-268': 'trio'}
# run_standard_gemini_protocol(work_dir, ped_dict)
'''
##first batch
'''
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/acc_fp1'
ped_dict = {'LR07-200':'trio', 'LR11-330':'trio', 'LR11-331':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/broad'
ped_dict = {'LP97-114':'trio', 'LP97-141':'singleton', 'LP98-078':'singleton', 'LR01-099':'trio', 'LR01-173':'singleton', 'LR01-314':'singleton', 
		'LR02-264':'trio', 'LR04-315':'singleton', 'LR05-054':'singleton', 'LR05-160':'trio', 'LR07-016':'trio', 'LR07-227':'trio', 
		'LR09-141':'trio', 'LR10-024':'trio', 'LR10-102':'trio', 'LR10-227':'trio', 'LR10-260':'singleton', 'LR10-270':'duo', 'LR11-006':'trio', 
		'LR11-019':'duo', 'LR11-033':'trio', 'LR11-105':'trio', 'LR11-112':'trio', 'LR11-117':'trio', 'LR11-124_2':'trio', 'LR11-132':'singleton', 'LR11-144':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/dwm_rpt'
ped_dict = {'LR04-185':'duo', 'LR04-414':'duo', 'LR06-105':'duo', 'LR14-098':'duo'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/dwm_b1'
ped_dict = {'LR03-305':'trio', 'LR04-106':'trio', 'LR04-186':'trio', 'LR04-276':'trio', 'LR05-203_2':'trio', 'LR05-354':'trio', 'LR07-079':'trio', 'LR08-323':'trio', 
		'LR08-390':'trio', 'LR08-396':'trio', 'LR10-230':'trio', 'LR11-169':'trio', 'LR12-032':'trio', 'LR12-115':'trio', 'LR12-434':'trio', 'LR12-443':'trio', 
		'LR12-455':'trio', 'LR12-463':'trio', 'LR12-464':'trio', 'LR12-473':'trio', 'LR13-002':'trio', 'LR13-003':'trio', 'LR13-037':'trio', 'LR13-085':'trio', 
		'LR13-153':'trio', 'LR13-199':'trio', 'LR13-200':'trio', 'LR13-279':'trio', 'LR13-315':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/dwm_b2'
ped_dict = {'LR02-148':'trio', 'LR03-055':'trio', 'LR03-206':'trio', 'LR03-274':'trio', 'LR03-332':'trio', 'LR04-020':'trio', 'LR04-371':'trio', 'LR06-085':'trio', 
		'LR09-023_2':'trio', 'LR09-227':'trio', 'LR09-280':'trio', 'LR11-241':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/dwm_b3'
ped_dict = {'LR03-077':'trio', 'LR03-169':'trio', 'LR03-170':'trio', 'LR03-223':'trio', 'LR03-278':'trio', 'LR03-298':'trio', 'LR04-017':'trio', 'LR04-084':'trio', 
		'LR04-208':'trio', 'LR04-233':'trio', 'LR05-118':'trio', 'LR14-075':'singleton'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/dwm_b4'
ped_dict = {'LR00-225':'singleton', 'LR02-027':'trio', 'LR03-384':'duo', 'LR05-396':'trio', 'LR06-087':'trio', 'LR07-054':'singleton', 'LR12-475':'singleton'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/dwm_b5'
ped_dict = {'DWM10':'singleton', 'DWM13':'singleton', 'DWM3':'singleton', 'LR01-079':'trio', 'LR01-242':'singleton', 'LR01-306':'singleton', 'LR02-263':'duo', 
		'LR03-039':'singleton', 'LR03-304':'singleton', 'LR03-340':'singleton', 'LR04-022':'duo', 'LR04-222':'singleton', 'LR04-341':'singleton', 'LR05-162':'duo', 
		'LR05-265':'singleton', 'LR05-398':'trio', 'LR08-056':'singleton', 'LR12-313':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/e3'
ped_dict = {'LP98-041_2':'trio', 'LR10-036':'trio', 'LR10-199':'trio', 'LR11-001':'trio', 'LR11-328':'trio', 'LR12-112':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/e4'
ped_dict = {'LR03-037':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
'''
##second batch
'''
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/e5'
ped_dict = {'LR09-129':'trio', 'LR12-206':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/e6'
ped_dict = {'LR06-005':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/e7'
ped_dict = {'LR01-190_2':'trio', 'LR04-350':'trio', 'LR06-157':'trio', 'LR08-337':'trio', 'LR09-203':'trio', 'LR10-222':'trio', 'LR11-016':'trio', 
		'LR11-152':'trio', 'LR11-219':'trio', 'LR12-097':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/e8'
ped_dict = {'LR03-120':'trio', 'LR03-419':'trio', 'LR04-376':'trio', 'LR05-007_2':'trio', 'LR05-120':'trio', 'LR06-278':'trio', 'LR06-430':'trio', 
		'LR07-096':'trio', 'LR07-150':'trio', 'LR08-002':'trio', 'LR08-360':'trio', 'LR09-407':'trio', 'LR10-046':'trio', 'LR10-276':'trio', 
		'LR11-061':'trio', 'LR11-343':'trio', 'LR12-049':'trio', 'LR12-068':'trio', 'LR12-156':'trio', 'LR12-323':'trio', 'LR12-447':'trio', 
		'LR13-103':'trio', 'LR13-175':'trio', 'LR13-278':'trio', 'LR14-048':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/id1_mic2_pmg2'
ped_dict = {'LP97-105':'trio', 'LP98-095':'trio', 'LP99-100':'singleton', 'LR00-144':'trio', 'LR01-265':'trio', 'LR01-338':'trio', 'LR02-005':'trio', 
		'LR02-046':'trio', 'LR02-085':'singleton', 'LR03-056':'trio', 'LR03-214':'duo', 'LR04-199':'duo', 'LR04-289':'trio', 'LR06-231':'trio', 
		'LR06-234':'trio', 'LR07-037':'singleton', 'LR07-125':'trio', 'LR07-207':'trio', 'LR08-029':'trio', 'LR08-215':'trio', 'LR08-257':'trio', 
		'LR08-291':'trio', 'LR08-377':'trio', 'LR08-418':'duo', 'LR09-248':'trio', 'LR09-266':'trio', 'LR11-035':'trio', 'LR11-054':'trio', 'LR11-173':'duo', 
		'LR11-189':'trio', 'LR12-053':'trio', 'LR12-204':'duo', 'LR12-214':'duo', 'LR12-224':'duo'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/lis2'
ped_dict = {'LP98-041':'trio', 'LP99-032':'trio', 'LR05-375':'duo', 'LR06-397':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)

work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/meg_gw_0516'
ped_dict = {'LR03-130':'trio', 'LR05-393':'trio', 'LR09-131':'trio', 'LR10-015':'trio', 'LR11-465':'singleton', 'LR12-191_2':'trio', 'LR12-264':'trio', 
		'LR12-291':'trio', 'LR12-316':'trio', 'LR12-339':'trio', 'LR12-444':'trio', 'LR13-068':'trio', 'LR13-069':'trio', 'LR13-180':'trio', 'LR14-190':'trio', 
		'LR14-285':'trio', 'LR16-053':'trio'}
ped_dict = {'LR16-053':'singleton'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/perkin_elmer'
ped_dict = {'LR06-207':'trio', 'LR09-416':'singleton', 'LR10-026':'trio', 'LR10-228':'trio', 'LR10-242':'trio', 'LR11-042':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/rho1'
ped_dict = {'LR11-367':'trio', 'LR12-018':'duo', 'LR12-401':'duo'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/transient_rpt'
ped_dict = {'LP88-024':'trio', 'LR06-254':'trio', 'LR12-165':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/transient1'
ped_dict = {'LP89-036':'trio', 'LR05-007':'trio', 'LR05-060_2':'trio', 'LR06-300':'trio', 'LR07-031':'trio', 'LR08-023':'trio', 'LR09-023':'trio', 'LR10-193':'trio', 'LR10-243':'trio', 
		'LR12-101':'trio', 'LR12-127':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/transient1.3'
ped_dict = {'LR00-220':'trio', 'LR06-355':'trio', 'LR12-343':'trio', 'LR12-379':'trio', 'LR12-436':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/transient2'
ped_dict = {'LP99-206':'trio', 'LR01-224':'trio', 'LR05-035':'trio', 'LR05-404':'trio', 'LR06-029':'trio', 'LR06-099':'trio', 'LR08-049':'trio', 'LR09-105':'trio', 
		'LR09-149':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/uw1'
ped_dict = {'LP92-083':'trio', 'LR00-016':'singleton', 'LR01-279':'singleton', 'LR04-005':'trio', 'LR08-108':'trio', 'LR09-006':'trio', 'LR09-387':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/sherr_exomes_0317'
ped_dict = {'1090':'trio', '1173':'trio', '1175':'trio', '1339':'trio', '1588':'trio', '1629':'trio', '1289':'trio', '1318':'trio', '1324':'trio', '1416':'trio', 
		'1426':'trio', '1512':'trio', '1011':'trio', '1020':'trio', '1200':'trio', '1536':'trio', '1590':'trio', '1201':'trio', '1795':'trio', '2352':'trio', 
		'1036':'trio', '1270':'trio', '1351':'trio', '1250':'trio', '1257':'trio', '1327':'trio', '1510':'trio', '1533':'trio', '1630':'trio', '1810':'trio', 
		'1814':'trio', '1940':'trio', '1978':'trio', '2177':'trio', '2228':'trio', '2301':'trio', '2311':'trio', '2639':'trio', '2644':'trio', '2707':'trio', 
		'2705':'trio', '2756':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
'''
#third batch
'''
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/exomes_after_0317'
ped_dict = {'9267':'trio', 'LR01-194':'trio', 'LR01-381':'trio', 'LR02-291':'trio', 'LR03-395':'trio', 'LR04-239':'trio', 'LR04-290_2':'trio', 'LR05-055':'trio', 
		'LR05-069':'trio', 'LR05-197':'trio', 'LR05-201':'trio', 'LR06-339':'trio', 'LR08-299':'trio', 'LR09-016':'trio', 'LR09-064':'trio', 'LR10-064':'trio', 
		'LR10-084':'trio', 'LR10-246':'trio', 'LR11-003':'trio', 'LR11-025':'trio', 'LR11-124':'trio', 'LR11-151':'singleton', 'LR11-233':'trio', 'LR11-347':'trio', 
		'LR12-172':'trio', 'LR12-188':'trio', 'LR12-225':'trio', 'LR12-304':'trio', 'LR12-308':'trio', 'LR12-346':'duo', 'LR12-354_2':'trio', 'LR12-492':'trio', 
		'LR13-013':'trio', 'LR13-061':'trio', 'LR13-150':'trio', 'LR13-160':'trio', 'LR13-228':'trio', 'LR13-261':'trio', 'LR13-270':'trio', 'LR13-300':'trio', 
		'LR14-031':'trio', 'LR14-071':'trio', 'LR14-211':'trio', 'LR14-221':'trio', 'LR14-291_2':'trio', 'LR15-011':'trio', 'LR15-061':'trio', 'LR15-137':'singleton', 
		'LR15-155':'trio', 'LR15-236':'trio', 'LR16-003':'trio', 'LR16-034':'trio', 'LR16-079':'trio', 'LR16-170':'trio', 'LR16-295':'trio', 'LR16-306':'singleton', 
		'LR16-368':'trio', 'LR16-412':'trio', 'LR16-451':'trio', 'LR16-489':'trio', 'LR17-038':'trio', 'LR17-049':'trio', 'LR17-050':'trio', 'LR17-051':'trio', 
		'LR17-052':'trio', 'LR17-053':'trio', 'LR17-054':'trio', 'LR17-055':'trio', 'LR17-056':'trio', 'LR17-065':'trio', 'LR17-076':'trio', 'LR17-170':'trio', 
		'LR17-260':'duo', 'LR17-332':'trio', 'LR17-361':'trio', 'LR17-432':'trio', 'LR17-438':'trio', 'LR17-439':'duo'}
run_standard_gemini_protocol(work_dir, ped_dict)
'''
#forth batch
'''
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/exomes_after_0616'
ped_dict = {'LR04-316':'trio', 'LR12-305':'trio', 'LR13-356':'trio', 'LR14-226':'trio', 'LR14-252':'trio', 'LR14-255':'trio', 'LR15-272':'trio', 'LR16-029':'duo', 
		'LR16-065':'trio', 'LR16-173':'trio', 'LR16-239':'trio', '160342':'trio', '162107':'trio', 'LR12-459':'trio', 'LR15-121':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/exomes_after_1116'
ped_dict = {'3C-10':'singleton', '3C-11':'singleton', '3C-12':'singleton', '3C-14':'singleton', '3C-2':'trio', '3C-3':'duo', '3C-4':'trio', '3C-6':'trio', 
		'3C-7':'singleton', '3C-8':'trio', 'LR01-312':'trio', 'LR03-264':'trio', 'LR04-036_2':'trio', 'LR04-038':'trio', 'LR05-250':'trio', 'LR06-137':'trio', 
		'LR08-415':'trio', 'LR09-055':'trio', 'LR09-181':'trio', 'LR09-296':'trio', 'LR13-032':'trio', 'LR15-267':'trio', 'LR04-036':'trio', 'LR12-230':'trio', 
		'LR12-476':'trio', 'LR13-321':'trio', 'LR14-092':'trio', 'LR15-003':'trio', 'LR15-004':'trio', 'LR15-081':'trio', 'LR15-097':'trio', 'LR11-383':'trio', 
		'LR15-121_2':'trio', 'LR16-210':'trio', 'LR16-280':'trio', 'LR04-399':'trio', 'LR09-128':'trio', 'LR12-439':'trio', 'LR15-005':'trio', 'LR15-340':'trio', 
		'LR16-228':'trio', 'LR13-084':'singleton', 'LR15-338':'trio', 'LR16-089':'trio', 'LR16-441':'singleton'}
run_standard_gemini_protocol(work_dir, ped_dict)
'''
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/exomes_after_0616'
ped_dict = {'LR14-255':'singleton'}
run_standard_gemini_protocol(work_dir, ped_dict)
##fifth batch
'''
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/exomes_after_1217'
ped_dict = {'LR17-478':'singleton', 'LR02-304':'trio', 'LR04-336':'trio', 'LR08-270':'trio', 'LR13-411':'trio', 'LR14-137_2':'trio', 'LR15-347':'trio', 'LR15-377':'trio', 
		'LR16-087':'trio', 'LR17-075':'trio', 'LR17-106':'trio', 'LR17-172':'trio', 'LR17-433':'trio', 'LR17-490':'trio', 'LR17-505':'trio', 'LR18-002':'trio', 'LR18-014':'trio', 
		'LR18-040':'trio', 'LR18-043':'trio', 'LR18-234':'trio', 'LR18-072':'trio', 'LR18-310':'trio', 'LR18-339':'singleton', 'LR18-340':'singleton', 'LR18-236':'trio', 
		'LR18-030':'duo', 'LR18-031':'trio', 'LP97-070':'trio', 'LR02-434':'trio', 'LR05-393_2':'singleton', 'LR11-216':'trio', 'LR12-213':'trio', 'LR13-241':'trio', 
		'LR14-132':'trio', 'LR14-166':'trio', 'LR14-322':'trio', 'LR14-325':'trio', 'LR15-345':'trio', 'LR16-064':'trio', 'LR16-067':'trio', 'LR16-410':'trio', 
		'LR16-434':'trio', 'LR17-251':'trio', 'LP97-070_2':'trio', 'LR02-434_2':'trio', 'LR04-290':'trio', 'LR06-018':'singleton', 'LR07-201':'trio', 'LR10-016':'trio', 
		'LR11-216_2':'trio', 'LR11-429':'singleton', 'LR12-109':'singleton', 'LR12-191':'trio', 'LR12-241':'trio', 'LR12-265':'trio', 'LR12-266':'trio', 'LR12-317':'trio', 
		'LR12-354':'trio', 'LR12-420':'trio', 'LR12-437':'trio', 'LR13-052':'trio', 'LR13-118':'trio', 'LR13-189':'trio', 'LR13-241_2':'trio', 'LR13-262':'trio', 'LR13-282':'trio', 
		'LR13-364':'trio', 'LR14-026':'trio', 'LR14-041':'trio', 'LR14-066':'trio', 'LR14-093':'trio', 'LR14-097':'trio', 'LR14-102':'trio', 'LR14-132_2':'trio', 'LR14-155':'trio', 
		'LR14-166_2':'trio', 'LR14-223':'trio', 'LR14-266':'trio', 'LR14-291':'trio', 'LR14-322_2':'trio', 'LR14-325_2':'trio', 'LR15-080':'trio', 'LR15-186':'trio', 'LR15-215':'duo', 
		'LR15-225':'trio', 'LR15-284':'trio', 'LR15-345_2':'trio', 'LR16-004':'duo', 'LR16-067_2':'trio', 'LR16-157':'duo', 'LR16-159':'trio', 'LR16-234':'trio', 'LR16-305':'trio', 
		'LR16-360':'trio', 'LR16-407':'trio', 'LR16-432':'trio', 'LR16-434_2':'trio', 'LR17-083':'trio', 'LR17-100':'singleton', 'LR17-101':'singleton', 'LR17-102':'singleton', 
		'LR17-103':'singleton', 'LR17-176':'duo', 'LR17-229':'trio', 'LR17-242':'trio', 'LR17-337':'trio', 'LR15-112':'duo', 'LR17-059':'duo', 'LR13-130':'trio', 'LR13-248':'trio', 
		'LR13-299':'trio', 'LR13-336':'trio', 'LR13-405':'trio', 'LR14-007':'trio', 'LR14-053':'trio', 'LR14-165':'trio', 'LR14-248':'trio', 'LR14-267':'trio', 'LR14-273':'trio', 
		'LR14-292':'trio', 'LR16-227':'trio', 'LR16-293':'trio', 'LR16-300':'trio', 'LR16-421':'trio', 'LR17-217':'trio', '794991':'trio', 'LP99-109':'trio', 'LR01-190':'trio', 
		'LR04-052':'trio', 'LR05-060':'trio', 'LR07-082':'duo', 'LR08-374':'trio', 'LR14-107':'trio', 'LR14-137':'trio', 'LR15-374':'trio', 'LR16-081':'trio', 'LR16-406':'duo', 
		'LR16-030':'trio', 'LR16-031':'trio', 'LR16-172':'trio'}
run_standard_gemini_protocol(work_dir, ped_dict)
'''
'''
##sixth batch
work_dir = '/home/atimms/ngs_data/exomes/working/ghayda_reanalysis_0319/exomes_after_1018'
ped_dict = {'LR18-470':'trio', 'LR18-416':'duo', 'LR18-493':'trio', 'LR18-492':'trio', 'LR18-563':'trio', 'LR05-203':'trio',
		 'LR18-529':'singleton', 'LR12-489':'trio', 'LR19-273':'singleton'}
run_standard_gemini_protocol(work_dir, ped_dict)
'''
