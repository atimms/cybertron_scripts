#!/usr/bin/env python
import os
import subprocess
import glob
import shutil
import sys

'''
requires:
qsub -Iq cdbrmq -l mem=100gb,ncpus=10 -P 2cf4ba3b-f9ef-4cfc-a2c1-be85f5db6730
qsub -Iq cdbrmq -l mem=17gb,ncpus=1 -P 2cf4ba3b-f9ef-4cfc-a2c1-be85f5db6730
'''

##set input variables and parameters
delim = '\t'
threads = '5'


##set input variables and parameters
delim = '\t'
threads = '5'

##programs
bcftools = '/home/atimms/programs/bcftools-1.11/bcftools'
table_annovar = '/home/atimms/programs/annovar_1019/table_annovar.pl'


##annovar parameters
av_genome = 'hg38'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
# av_protocol = ['-protocol', 'refGene,clinvar_20190305,cosmic90_coding,cosmic90_noncoding']
av_protocol = ['-protocol', 'refGene,knownGene,dbnsfp41a,clinvar_20210123,gnomad30_genome']
av_operation = ['-operation', 'g,g,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.', '-vcfinput' ]
extra_header = ['cohort_freq', 'filter', 'info', 'format', 'NWD106517', 'NWD143463', 'NWD154175', 'NWD137485', 'NWD168739', 'NWD166767', 'NWD180579', 'NWD177262', 'NWD195236', 'NWD152229', 'NWD172985', 'NWD194569', 'NWD162784', 'NWD161421', 'NWD215779', 'NWD170420', 'NWD160984', 'NWD128599', 'NWD183506', 'NWD124783', 'NWD215932', 'NWD241600', 'NWD239371', 'NWD259598', 'NWD154254', 'NWD105656', 'NWD195505', 'NWD151355', 'NWD278342', 'NWD194109', 'NWD180893', 'NWD214708', 'NWD244874', 'NWD146928', 'NWD206068', 'NWD232383', 'NWD248077', 'NWD147251', 'NWD197694', 'NWD262109', 'NWD280984', 'NWD118221', 'NWD243367', 'NWD211383', 'NWD113031', 'NWD190729', 'NWD209836', 'NWD247189', 'NWD188965', 'NWD172958', 'NWD287178', 'NWD312195', 'NWD101284', 'NWD140828', 'NWD277162', 'NWD295501', 'NWD330244', 'NWD225096', 'NWD157765', 'NWD124505', 'NWD194601', 'NWD239110', 'NWD265890', 'NWD335421', 'NWD251426', 'NWD326817', 'NWD204367', 'NWD238576', 'NWD199154', 'NWD196048', 'NWD384640', 'NWD268883', 'NWD362624', 'NWD214856', 'NWD163566', 'NWD187086', 'NWD421518', 'NWD360540', 'NWD390204', 'NWD167444', 'NWD429224', 'NWD220071', 'NWD406384', 'NWD108686', 'NWD279693', 'NWD177573', 'NWD355011', 'NWD147406', 'NWD307689', 'NWD353972', 'NWD360848', 'NWD370421', 'NWD257524', 'NWD442543', 'NWD375385', 'NWD348241', 'NWD437304', 'NWD175885', 'NWD142157', 'NWD199671', 'NWD391316', 'NWD308528', 'NWD433762', 'NWD197875', 'NWD105449', 'NWD424245', 'NWD254741', 'NWD259001', 'NWD312042', 'NWD213482', 'NWD406629', 'NWD219802', 'NWD112336', 'NWD350735', 'NWD442666', 'NWD241567', 'NWD182796', 'NWD469914', 'NWD116483', 'NWD402031', 'NWD375215', 'NWD154984', 'NWD470659', 'NWD157378', 'NWD184442', 'NWD325903', 'NWD428903', 'NWD427656', 'NWD397096', 'NWD202851', 'NWD527072', 'NWD132741', 'NWD267017', 'NWD354891', 'NWD406808', 'NWD285938', 'NWD115422', 'NWD205117', 'NWD482024', 'NWD170409', 'NWD201053', 'NWD382751', 'NWD349027', 'NWD473743', 'NWD512549', 'NWD412487', 'NWD241102', 'NWD374227', 'NWD387730', 'NWD130906', 'NWD308907', 'NWD201383', 'NWD550359', 'NWD323492', 'NWD504644', 'NWD272577', 'NWD532882', 'NWD292988', 'NWD359126', 'NWD251301', 'NWD547047', 'NWD486761', 'NWD256214', 'NWD202562', 'NWD411806', 'NWD508063', 'NWD439231', 'NWD496610', 'NWD448474', 'NWD460488', 'NWD201846', 'NWD113907', 'NWD463145', 'NWD107320', 'NWD277715', 'NWD403825', 'NWD475992', 'NWD386557', 'NWD458185', 'NWD117126', 'NWD414754', 'NWD516724', 'NWD483720', 'NWD382420', 'NWD531918', 'NWD126296', 'NWD449106', 'NWD181551', 'NWD559405', 'NWD590418', 'NWD548005', 'NWD615187', 'NWD290367', 'NWD207266', 'NWD129301', 'NWD236227', 'NWD392584', 'NWD320632', 'NWD463667', 'NWD638317', 'NWD444775', 'NWD338702', 'NWD506867', 'NWD189910', 'NWD500227', 'NWD233268', 'NWD347857', 'NWD475652', 'NWD393045', 'NWD496491', 'NWD582968', 'NWD629328', 'NWD561893', 'NWD279860', 'NWD575037', 'NWD558913', 'NWD266069', 'NWD265695', 'NWD550995', 'NWD388992', 'NWD200378', 'NWD262923', 'NWD596337', 'NWD475254', 'NWD472427', 'NWD181479', 'NWD396387', 'NWD462690', 'NWD657292', 'NWD195485', 'NWD589692', 'NWD602941', 'NWD557950', 'NWD227250', 'NWD568604', 'NWD157617', 'NWD506881', 'NWD332379', 'NWD381744', 'NWD599785', 'NWD459192', 'NWD367641', 'NWD185507', 'NWD655267', 'NWD633161', 'NWD336389', 'NWD562686', 'NWD137782', 'NWD403410', 'NWD364824', 'NWD635647', 'NWD362385', 'NWD511935', 'NWD303380', 'NWD354464', 'NWD297116', 'NWD415041', 'NWD656011', 'NWD426024', 'NWD655358', 'NWD639566', 'NWD566518', 'NWD580637', 'NWD111216', 'NWD582010', 'NWD314162', 'NWD728767', 'NWD316160', 'NWD688683', 'NWD711969', 'NWD247578', 'NWD237111', 'NWD696407', 'NWD298155', 'NWD643728', 'NWD435521', 'NWD610591', 'NWD402736', 'NWD302735', 'NWD432045', 'NWD763104', 'NWD688890', 'NWD609009', 'NWD429129', 'NWD764908', 'NWD625449', 'NWD482870', 'NWD573070', 'NWD629676', 'NWD774542', 'NWD593892', 'NWD699349', 'NWD504379', 'NWD749201', 'NWD509210', 'NWD555990', 'NWD351877', 'NWD594846', 'NWD572170', 'NWD427079', 'NWD619433', 'NWD747359', 'NWD515052', 'NWD774221', 'NWD476461', 'NWD596804', 'NWD416851', 'NWD650999', 'NWD677065', 'NWD424326', 'NWD299937', 'NWD785483', 'NWD726903', 'NWD697924', 'NWD347046', 'NWD783464', 'NWD645239', 'NWD290779', 'NWD766429', 'NWD627900', 'NWD618207', 'NWD294004', 'NWD652714', 'NWD659227', 'NWD545229', 'NWD110497', 'NWD437955', 'NWD605469', 'NWD480737', 'NWD673853', 'NWD756099', 'NWD277285', 'NWD415543', 'NWD304766', 'NWD834863', 'NWD725749', 'NWD619790', 'NWD750603', 'NWD619748', 'NWD623175', 'NWD748528', 'NWD417795', 'NWD522992', 'NWD513523', 'NWD769210', 'NWD549871', 'NWD599391', 'NWD813046', 'NWD837230', 'NWD524813', 'NWD774778', 'NWD858031', 'NWD511925', 'NWD753348', 'NWD757501', 'NWD873586', 'NWD766832', 'NWD549095', 'NWD604509', 'NWD155359', 'NWD809752', 'NWD665979', 'NWD490985', 'NWD583889', 'NWD857942', 'NWD887396', 'NWD663932', 'NWD863300', 'NWD869942', 'NWD872277', 'NWD520101', 'NWD898191', 'NWD827563', 'NWD861728', 'NWD461927', 'NWD897542', 'NWD878254', 'NWD776037', 'NWD562995', 'NWD640573', 'NWD549172', 'NWD777457', 'NWD773102', 'NWD495246', 'NWD506080', 'NWD922698', 'NWD646121', 'NWD650067', 'NWD873387', 'NWD689163', 'NWD192630', 'NWD881966', 'NWD846973', 'NWD874216', 'NWD581569', 'NWD531203', 'NWD932551', 'NWD646957', 'NWD587595', 'NWD395441', 'NWD173042', 'NWD857549', 'NWD871339', 'NWD772586', 'NWD623877', 'NWD778133', 'NWD759523', 'NWD905958', 'NWD780209', 'NWD734936', 'NWD717790', 'NWD821993', 'NWD932395', 'NWD797431', 'NWD802887', 'NWD139502', 'NWD593074', 'NWD654247', 'NWD255301', 'NWD586839', 'NWD871893', 'NWD464151', 'NWD602830', 'NWD444909', 'NWD891671', 'NWD680501', 'NWD426987', 'NWD674601', 'NWD540529', 'NWD959224', 'NWD769529', 'NWD607566', 'NWD391260', 'NWD592705', 'NWD964055', 'NWD232569', 'NWD300680', 'NWD333248', 'NWD337116', 'NWD359622', 'NWD393910', 'NWD397737', 'NWD424137', 'NWD427745', 'NWD442840', 'NWD450157', 'NWD456497', 'NWD461148', 'NWD461909', 'NWD464532', 'NWD466937', 'NWD469148', 'NWD470014', 'NWD479072', 'NWD482167', 'NWD485820', 'NWD489577', 'NWD490712', 'NWD494689', 'NWD512837', 'NWD521660', 'NWD523388', 'NWD529793', 'NWD537028', 'NWD539508', 'NWD539918', 'NWD543011', 'NWD545288', 'NWD548898', 'NWD553928', 'NWD556432', 'NWD562339', 'NWD567374', 'NWD568465', 'NWD582270', 'NWD583735', 'NWD612282', 'NWD620072', 'NWD633504', 'NWD640151', 'NWD642366', 'NWD643638', 'NWD649000', 'NWD650829', 'NWD654386', 'NWD655923', 'NWD662677', 'NWD666857', 'NWD669383', 'NWD674152', 'NWD675144', 'NWD679158', 'NWD679813', 'NWD680999', 'NWD683223', 'NWD689669', 'NWD691096', 'NWD699826', 'NWD700258', 'NWD701011', 'NWD705122', 'NWD714615', 'NWD716459', 'NWD719874', 'NWD721157', 'NWD722009', 'NWD727148', 'NWD728463', 'NWD749045', 'NWD750367', 'NWD752573', 'NWD752766', 'NWD754350', 'NWD754553', 'NWD759385', 'NWD763355', 'NWD768186', 'NWD768300', 'NWD775436', 'NWD777673', 'NWD778520', 'NWD779263', 'NWD779965', 'NWD781654', 'NWD783729', 'NWD786446', 'NWD789089', 'NWD791095', 'NWD791684', 'NWD793514', 'NWD795868', 'NWD800626', 'NWD801155', 'NWD801626', 'NWD816080', 'NWD817855', 'NWD819417', 'NWD821282', 'NWD822824', 'NWD823550', 'NWD831552', 'NWD832984', 'NWD835305', 'NWD837577', 'NWD840684', 'NWD843294', 'NWD844412', 'NWD844871', 'NWD845729', 'NWD846259', 'NWD846536', 'NWD846997', 'NWD847072', 'NWD847957', 'NWD849039', 'NWD852195', 'NWD853192', 'NWD854939', 'NWD860129', 'NWD861279', 'NWD861971', 'NWD869379', 'NWD874452', 'NWD875670', 'NWD876199', 'NWD876644', 'NWD879942', 'NWD881030', 'NWD892650', 'NWD892920', 'NWD896397', 'NWD896853', 'NWD899289', 'NWD900321', 'NWD904536', 'NWD906714', 'NWD909195', 'NWD911824', 'NWD912179', 'NWD914437', 'NWD915415', 'NWD915924', 'NWD916171', 'NWD918288', 'NWD922769', 'NWD924608', 'NWD925651', 'NWD928861', 'NWD930344', 'NWD933034', 'NWD934376', 'NWD938588', 'NWD947011', 'NWD947532', 'NWD948417', 'NWD949178', 'NWD950511', 'NWD951431', 'NWD954039', 'NWD954371', 'NWD955193', 'NWD955616', 'NWD957257', 'NWD959156', 'NWD959729', 'NWD961340', 'NWD962202', 'NWD962995', 'NWD963432', 'NWD964195', 'NWD973800', 'NWD974018', 'NWD974183', 'NWD976658', 'NWD976662', 'NWD978533', 'NWD979105', 'NWD981931', 'NWD982856', 'NWD984489', 'NWD985173', 'NWD989436', 'NWD989592', 'NWD991289', 'NWD992137', 'NWD993368', 'NWD993544', 'NWD998627']


##methods
def split_vcf_by_sample(in_vcf, sample_file, out_vcf):
	##get the samples we want, and remove when we don't see a call, using one or multiple threads
	# bcftools_view = subprocess.Popen([bcftools, 'view', '-a', '-Ou', '-S', sample_file, in_vcf], stdout=subprocess.PIPE)
	# bcftools_view2 = subprocess.Popen([bcftools, 'view', '-m', '2', '-O', 'z', '-o', ped_vcf, '-'], stdin=bcftools_view.stdout)
	bcftools_view = subprocess.Popen([bcftools, 'view', '-a', '--threads', '5', '-Ou', '-S', sample_file, in_vcf], stdout=subprocess.PIPE)
	bcftools_view2 = subprocess.Popen([bcftools, 'view', '-m', '2', '--threads', '5', '-O', 'z', '-o', out_vcf, '-'], stdin=bcftools_view.stdout)
	bcftools_view2.wait()


def annotate_vcf_with_annovar(in_vcf, ann_prefix):
	##annotate vcfs with annovar i.e. run_table_annovar
	command = [table_annovar] + av_buildver + [in_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', ann_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()

def multi_to_ann(in_file, out_file):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc == 1:
				header = line[:15] + line[18:20] + [line[40]] + line[54:56] + line[66:84] + extra_header
				out_fh.write(delim.join(header) + '\n')
			else:
				line_out = line[:15] + line[18:20] + [line[40]] + line[54:56] + line[66:85] + line[93:]
				out_fh.write(delim.join(line_out) + '\n')

def filter_exon_rare_vars(in_file, out_file):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc == 1:
				out_fh.write(delim.join(line) + '\n')
			else:
				rg_func = line[5]
				kn_func = line[10]
				rg_exfunc = line[8]
				kn_exfunc = line[13]
				af_info = line[25]

				if af_info == '.':
					af = 0
				else:
					af = float(af_info)
				if rg_func =='exonic' or rg_func == 'splicing' or rg_func == 'exonic;splicing' or kn_func =='exonic' or kn_func == 'splicing' or kn_func == 'exonic;splicing':
					if af < 0.01:
						if rg_exfunc == 'synonymous SNV' and kn_exfunc == 'synonymous SNV':
							pass
						else:
							out_fh.write(delim.join(line) + '\n')


def filter_by_genenames(in_file, genelist, out_file):
	with open(in_file, "r") as in_fh, open(out_file, "w") as out_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc == 1:
				out_fh.write(delim.join(line) + '\n')
			else:
				genes = line[6].split(';') + line[11].split(';')
				if bool(set(genes) & set(genelist)):
					out_fh.write(delim.join(line) + '\n')


def format_vars_by_proband(in_file, out_file):
	proband_dict = {}
	with open(in_file, "r") as in_fh:
		lc = 0
		for line in in_fh:
			line = line.rstrip().split(delim)
			lc += 1
			if lc == 1:
				sample_names = line[42:]
			else:
				gene = line[6]
				var = ':'.join(line[:5])
				genotypes = line[42:]
				gp = -1
				for genotype in genotypes:
					gp += 1
					if genotype == '0/1':
						sample = sample_names[gp]
						if sample in proband_dict:
							if gene not in proband_dict[sample][0]:
								proband_dict[sample][0].append(gene)
							if var not in proband_dict[sample][1]:
								proband_dict[sample][1].append(var)
						else:
							proband_dict[sample] = [[gene], [var]]

	with open(out_file, "w") as out_fh:
		out_fh.write(delim.join(['proband', 'genes', 'vars']) + '\n')
		for p in proband_dict:
			genes = ', '.join(proband_dict[p][0])
			variants = ', '.join(proband_dict[p][1])
			out_fh.write(delim.join([p, genes, variants]) + '\n')

##run methods
working_dir = '/home/atimms/ngs_data/genomes/dave_lisa_chd_genomes_0821'
os.chdir(working_dir)
project_name = 'chd_genomes_0821'
combined_vcf = 'PCGC_CHD_phs001735_TOPMed_WGS_freeze.8.combined.hg38.vcf.gz'
proband_vcf = 'PCGC_CHD_phs001735_TOPMed_WGS_freeze.8.hg38.probands.vcf.gz'
proband_samples = 'proband_ids.txt'
multianno = project_name + '.hg38_multianno.txt'
annotated = project_name + '.hg38_annotated.txt'
exonic_rare_vars = project_name + '.exonic_rare.txt'


##step1. filter vcf for proband samples
# split_vcf_by_sample(combined_vcf, proband_samples, proband_vcf)

##step2. annotate with annovar and sormat
# proband_vcf= 'temp.vcf'
# annotate_vcf_with_annovar(proband_vcf, project_name)
# multi_to_ann(multianno, annotated)

##step3. filter exonic rare vars
# filter_exon_rare_vars(annotated, exonic_rare_vars)

##step4. filter by gene
lisa_genes = ['POMP', 'PSMD6', 'PSMA6', 'PSMD3', 'PSMA7']
chd_genes = ['ABL1', 'ACTC1', 'ACVR1', 'ACVR2B', 'ADAMTS10', 'AFF4', 'ANKRD11', 'ARID1A', 'ARID1B', 'B3GAT3', 'BCOR', 
		'BMPR2', 'BRAF', 'CDK13', 'CFC1', 'CHD4', 'CHD7', 'CHST14', 'CITED2', 'CREBBP', 'CRELD1', 'DLL4', 'DNAH11', 
		'DOCK6', 'EFTUD2', 'EHMT1', 'ELN', 'EP300', 'ESCO2', 'EVC', 'EVC2', 'FBN1', 'FGFR2', 'FLNA', 'FLT4', 'FOXC1', 
		'FOXC2', 'FOXH1', 'FOXP1', 'GATA4', 'GATA5', 'GATA6', 'GDF1', 'GJA1', 'GLI3', 'GPC3', 'HAND1', 'HAND2', 'HDAC8', 
		'HNRNPK', 'HRAS', 'INVS', 'JAG1', 'KANSL1', 'KAT6A', 'KAT6B', 'KDM6A', 'KMT2A', 'KMT2D', 'KRAS', 'KYNU', 'MAP2K1', 
		'MAP2K2', 'MAP3K7', 'MED12', 'MED13L', 'MEIS2', 'MESP1', 'MYBPC3', 'MYH11', 'MYH6', 'MYH7', 'NF1', 'NIPBL', 
		'NKX2-5', 'NKX2-6', 'NODAL', 'NONO', 'NOTCH1', 'NOTCH2', 'NPHP3', 'NPHP4', 'NR2F2', 'NRAS', 'NSD1', 'NUP188', 
		'PBX1', 'PIGL', 'PIGV', 'PITX2', 'PKD1L1', 'PRDM6', 'PRKD1', 'PTPN11', 'RAB23', 'RAD21', 'RAF1', 'RBFOX2', 'RERE', 'RIT1', 'SALL1', 'SALL4', 'SF3B4', 'SHOC2', 'SMAD2', 'SMAD3', 'SMAD4', 'SMAD6', 'SMARCA4', 'SMARCB1', 'SMARCE1', 'SMC1A', 'SMC3', 'SMG9', 'SON', 'SOS1', 'STRA6', 'TAB2', 'TBX1', 'TBX20', 'TBX5', 'TFAP2B', 'TGFBR1', 'TGFBR2', 'TLL1', 'TRAF7', 'TXNL4A', 'UBR1', 'WASHC5', 'ZEB2', 'ZFPM2', 'ZIC3']
carm1 = ['CARM1']
lisa_genes_vars = project_name + '.exonic_rare.lisa_genes.xls'
chd_genes_vars = project_name + '.exonic_rare.chd_genes.xls'
# carm1_vars = project_name + '.exonic_rare.carm1.xls'
# filter_by_genenames(exonic_rare_vars, lisa_genes, lisa_genes_vars)
# filter_by_genenames(exonic_rare_vars, chd_genes, chd_genes_vars)
# filter_by_genenames(exonic_rare_vars, carm1, carm1_vars)

##step5. modify chd vars file to gene/var by sample
chd_genes_vars_formatted = project_name + '.exonic_rare.chd_genes.by_proband.xls'
format_vars_by_proband(chd_genes_vars, chd_genes_vars_formatted)






