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
extra_header = ['cohort_freq', 'filter', 'info', 'format', 'NWD100213', 'NWD101244', 'NWD101284', 'NWD102734', 'NWD105400', 'NWD105449', 'NWD105656', 'NWD106130', 'NWD106517', 'NWD107320', 'NWD108686', 'NWD110497', 'NWD111216', 'NWD112154', 'NWD112336', 'NWD113031', 'NWD113301', 'NWD113907', 'NWD115422', 'NWD116060', 'NWD116483', 'NWD117126', 'NWD117611', 'NWD118221', 'NWD120136', 'NWD120875', 'NWD121407', 'NWD124477', 'NWD124505', 'NWD124783', 'NWD126291', 'NWD126296', 'NWD127271', 'NWD127395', 'NWD128599', 'NWD129301', 'NWD129810', 'NWD130906', 'NWD132083', 'NWD132741', 'NWD135360', 'NWD135366', 'NWD135549', 'NWD135658', 'NWD135825', 'NWD135947', 'NWD136546', 'NWD137485', 'NWD137782', 'NWD137820', 'NWD137997', 'NWD139045', 'NWD139106', 'NWD139502', 'NWD140141', 'NWD140828', 'NWD141296', 'NWD142157', 'NWD143463', 'NWD146928', 'NWD147233', 'NWD147251', 'NWD147406', 'NWD151355', 'NWD152229', 'NWD153211', 'NWD154175', 'NWD154254', 'NWD154984', 'NWD155359', 'NWD157378', 'NWD157617', 'NWD157765', 'NWD160637', 'NWD160941', 'NWD160984', 'NWD161419', 'NWD161421', 'NWD161881', 'NWD162784', 'NWD163226', 'NWD163566', 'NWD164762', 'NWD165662', 'NWD166233', 'NWD166767', 'NWD167119', 'NWD167444', 'NWD168739', 'NWD170409', 'NWD170420', 'NWD171445', 'NWD171826', 'NWD172585', 'NWD172958', 'NWD172985', 'NWD173042', 'NWD173091', 'NWD174133', 'NWD174185', 'NWD174898', 'NWD175885', 'NWD177262', 'NWD177573', 'NWD179493', 'NWD180253', 'NWD180579', 'NWD180893', 'NWD181194', 'NWD181459', 'NWD181479', 'NWD181551', 'NWD181780', 'NWD182796', 'NWD183506', 'NWD184442', 'NWD185507', 'NWD185825', 'NWD185976', 'NWD186407', 'NWD187086', 'NWD188965', 'NWD189910', 'NWD190729', 'NWD191960', 'NWD192630', 'NWD194109', 'NWD194569', 'NWD194601', 'NWD195236', 'NWD195485', 'NWD195505', 'NWD196048', 'NWD196329', 'NWD196925', 'NWD197694', 'NWD197875', 'NWD199154', 'NWD199671', 'NWD200378', 'NWD201053', 'NWD201383', 'NWD201846', 'NWD201904', 'NWD202562', 'NWD202851', 'NWD204367', 'NWD205117', 'NWD205296', 'NWD206068', 'NWD207266', 'NWD207320', 'NWD209240', 'NWD209836', 'NWD210221', 'NWD211347', 'NWD211383', 'NWD212786', 'NWD213482', 'NWD214008', 'NWD214708', 'NWD214767', 'NWD214856', 'NWD215274', 'NWD215779', 'NWD215932', 'NWD219802', 'NWD220071', 'NWD222131', 'NWD223193', 'NWD225096', 'NWD225533', 'NWD226022', 'NWD227250', 'NWD228565', 'NWD229091', 'NWD232383', 'NWD232569', 'NWD233268', 'NWD236227', 'NWD237111', 'NWD237800', 'NWD238576', 'NWD239110', 'NWD239371', 'NWD241102', 'NWD241567', 'NWD241600', 'NWD242722', 'NWD243163', 'NWD243367', 'NWD244865', 'NWD244874', 'NWD247189', 'NWD247578', 'NWD248077', 'NWD249936', 'NWD251301', 'NWD251426', 'NWD254741', 'NWD255052', 'NWD255062', 'NWD255301', 'NWD256214', 'NWD257524', 'NWD258317', 'NWD258902', 'NWD259001', 'NWD259598', 'NWD261223', 'NWD261540', 'NWD262109', 'NWD262923', 'NWD263883', 'NWD264010', 'NWD265695', 'NWD265890', 'NWD266069', 'NWD267017', 'NWD268883', 'NWD269472', 'NWD270867', 'NWD271094', 'NWD272577', 'NWD273675', 'NWD277162', 'NWD277285', 'NWD277405', 'NWD277715', 'NWD278342', 'NWD279693', 'NWD279860', 'NWD280127', 'NWD280265', 'NWD280337', 'NWD280984', 'NWD283465', 'NWD283780', 'NWD285938', 'NWD286240', 'NWD287178', 'NWD288785', 'NWD290367', 'NWD290779', 'NWD290903', 'NWD292988', 'NWD294004', 'NWD294771', 'NWD295501', 'NWD297116', 'NWD298155', 'NWD299143', 'NWD299937', 'NWD300447', 'NWD300680', 'NWD301071', 'NWD302735', 'NWD303380', 'NWD303873', 'NWD304766', 'NWD307461', 'NWD307689', 'NWD308528', 'NWD308907', 'NWD310446', 'NWD312042', 'NWD312195', 'NWD313259', 'NWD314162', 'NWD316160', 'NWD318429', 'NWD318846', 'NWD320632', 'NWD323492', 'NWD323851', 'NWD324386', 'NWD325646', 'NWD325903', 'NWD326817', 'NWD327432', 'NWD328041', 'NWD330244', 'NWD332379', 'NWD332951', 'NWD333248', 'NWD334311', 'NWD334824', 'NWD335421', 'NWD335487', 'NWD335697', 'NWD336389', 'NWD337116', 'NWD337129', 'NWD337496', 'NWD338702', 'NWD341250', 'NWD342236', 'NWD342555', 'NWD342733', 'NWD343018', 'NWD345168', 'NWD346815', 'NWD346819', 'NWD347046', 'NWD347857', 'NWD348241', 'NWD349027', 'NWD349032', 'NWD350735', 'NWD351877', 'NWD351897', 'NWD352068', 'NWD353972', 'NWD354464', 'NWD354866', 'NWD354891', 'NWD355011', 'NWD358869', 'NWD359126', 'NWD359566', 'NWD359622', 'NWD360540', 'NWD360848', 'NWD362114', 'NWD362385', 'NWD362624', 'NWD363255', 'NWD364824', 'NWD364886', 'NWD365846', 'NWD366326', 'NWD366475', 'NWD367641', 'NWD370421', 'NWD371005', 'NWD374048', 'NWD374227', 'NWD374327', 'NWD375215', 'NWD375385', 'NWD376288', 'NWD378087', 'NWD378141', 'NWD378776', 'NWD379975', 'NWD381744', 'NWD382420', 'NWD382751', 'NWD383137', 'NWD384640', 'NWD386306', 'NWD386557', 'NWD387730', 'NWD388628', 'NWD388992', 'NWD389305', 'NWD390204', 'NWD391260', 'NWD391316', 'NWD391750', 'NWD392584', 'NWD393045', 'NWD393910', 'NWD394575', 'NWD395441', 'NWD396387', 'NWD396989', 'NWD397096', 'NWD397737', 'NWD398822', 'NWD402031', 'NWD402082', 'NWD402736', 'NWD403220', 'NWD403410', 'NWD403825', 'NWD404587', 'NWD404643', 'NWD406384', 'NWD406629', 'NWD406808', 'NWD407457', 'NWD409429', 'NWD410040', 'NWD411003', 'NWD411005', 'NWD411806', 'NWD412487', 'NWD414754', 'NWD415041', 'NWD415543', 'NWD416851', 'NWD417795', 'NWD421512', 'NWD421518', 'NWD422181', 'NWD424137', 'NWD424245', 'NWD424326', 'NWD426024', 'NWD426492', 'NWD426987', 'NWD427079', 'NWD427175', 'NWD427184', 'NWD427656', 'NWD427745', 'NWD428903', 'NWD429129', 'NWD429224', 'NWD431586', 'NWD432045', 'NWD433762', 'NWD434856', 'NWD435521', 'NWD437304', 'NWD437955', 'NWD438113', 'NWD439231', 'NWD439646', 'NWD440300', 'NWD442543', 'NWD442666', 'NWD442840', 'NWD444775', 'NWD444909', 'NWD445258', 'NWD446072', 'NWD448474', 'NWD449106', 'NWD450157', 'NWD451158', 'NWD454611', 'NWD456497', 'NWD457761', 'NWD458185', 'NWD459192', 'NWD460488', 'NWD460848', 'NWD461148', 'NWD461909', 'NWD461927', 'NWD462220', 'NWD462690', 'NWD463145', 'NWD463667', 'NWD463986', 'NWD464151', 'NWD464532', 'NWD466937', 'NWD469148', 'NWD469914', 'NWD470014', 'NWD470659', 'NWD472424', 'NWD472427', 'NWD473743', 'NWD474586', 'NWD475254', 'NWD475652', 'NWD475992', 'NWD476461', 'NWD479072', 'NWD480737', 'NWD481226', 'NWD482024', 'NWD482167', 'NWD482316', 'NWD482870', 'NWD483720', 'NWD484711', 'NWD485820', 'NWD486761', 'NWD489577', 'NWD490712', 'NWD490985', 'NWD492559', 'NWD492683', 'NWD494689', 'NWD495246', 'NWD496491', 'NWD496610', 'NWD500227', 'NWD504103', 'NWD504341', 'NWD504379', 'NWD504644', 'NWD506080', 'NWD506867', 'NWD506881', 'NWD507383', 'NWD508063', 'NWD509210', 'NWD509224', 'NWD511289', 'NWD511925', 'NWD511935', 'NWD512549', 'NWD512837', 'NWD513523', 'NWD514581', 'NWD515052', 'NWD516724', 'NWD518738', 'NWD520101', 'NWD520936', 'NWD521660', 'NWD522992', 'NWD523388', 'NWD524813', 'NWD525951', 'NWD527072', 'NWD529159', 'NWD529793', 'NWD530187', 'NWD531203', 'NWD531617', 'NWD531918', 'NWD532882', 'NWD534090', 'NWD535489', 'NWD535738', 'NWD537028', 'NWD537369', 'NWD539172', 'NWD539508', 'NWD539918', 'NWD540529', 'NWD541301', 'NWD541563', 'NWD543011', 'NWD544066', 'NWD544547', 'NWD544957', 'NWD545229', 'NWD545288', 'NWD546216', 'NWD546499', 'NWD547047', 'NWD548005', 'NWD548898', 'NWD549095', 'NWD549172', 'NWD549624', 'NWD549871', 'NWD550359', 'NWD550995', 'NWD553394', 'NWD553648', 'NWD553928', 'NWD554450', 'NWD554891', 'NWD555990', 'NWD556432', 'NWD557950', 'NWD558572', 'NWD558913', 'NWD559405', 'NWD561893', 'NWD562339', 'NWD562588', 'NWD562686', 'NWD562995', 'NWD565491', 'NWD566518', 'NWD566519', 'NWD567374', 'NWD568019', 'NWD568438', 'NWD568465', 'NWD568604', 'NWD569252', 'NWD572170', 'NWD572475', 'NWD573070', 'NWD573847', 'NWD574143', 'NWD574342', 'NWD575037', 'NWD576659', 'NWD576836', 'NWD579904', 'NWD580637', 'NWD581211', 'NWD581569', 'NWD582010', 'NWD582270', 'NWD582968', 'NWD583320', 'NWD583735', 'NWD583889', 'NWD584221', 'NWD585646', 'NWD586839', 'NWD587595', 'NWD588588', 'NWD589692', 'NWD590418', 'NWD591140', 'NWD592705', 'NWD593074', 'NWD593892', 'NWD594512', 'NWD594846', 'NWD595545', 'NWD596337', 'NWD596458', 'NWD596804', 'NWD597844', 'NWD599391', 'NWD599785', 'NWD600333', 'NWD602830', 'NWD602941', 'NWD603112', 'NWD604036', 'NWD604509', 'NWD605469', 'NWD607566', 'NWD608984', 'NWD609009', 'NWD610591', 'NWD611891', 'NWD612282', 'NWD612834', 'NWD614029', 'NWD615187', 'NWD616771', 'NWD618207', 'NWD619433', 'NWD619748', 'NWD619790', 'NWD620072', 'NWD623175', 'NWD623877', 'NWD625449', 'NWD627900', 'NWD629328', 'NWD629448', 'NWD629676', 'NWD631500', 'NWD633161', 'NWD633504', 'NWD635630', 'NWD635647', 'NWD636927', 'NWD637186', 'NWD638317', 'NWD639566', 'NWD639643', 'NWD640151', 'NWD640573', 'NWD641855', 'NWD642366', 'NWD643638', 'NWD643728', 'NWD644210', 'NWD644400', 'NWD645239', 'NWD646121', 'NWD646379', 'NWD646957', 'NWD649000', 'NWD650067', 'NWD650829', 'NWD650999', 'NWD651617', 'NWD652714', 'NWD652869', 'NWD652992', 'NWD653984', 'NWD654247', 'NWD654282', 'NWD654386', 'NWD654806', 'NWD655267', 'NWD655358', 'NWD655923', 'NWD656011', 'NWD657292', 'NWD659227', 'NWD662071', 'NWD662677', 'NWD663932', 'NWD665979', 'NWD666388', 'NWD666424', 'NWD666857', 'NWD669383', 'NWD670719', 'NWD671836', 'NWD672205', 'NWD672462', 'NWD673853', 'NWD674152', 'NWD674601', 'NWD675144', 'NWD675589', 'NWD677065', 'NWD677702', 'NWD679158', 'NWD679813', 'NWD680501', 'NWD680999', 'NWD682543', 'NWD683223', 'NWD687601', 'NWD688683', 'NWD688890', 'NWD689163', 'NWD689669', 'NWD690592', 'NWD691096', 'NWD693760', 'NWD693892', 'NWD696407', 'NWD697768', 'NWD697924', 'NWD699349', 'NWD699826', 'NWD700258', 'NWD701011', 'NWD702133', 'NWD704295', 'NWD705122', 'NWD711818', 'NWD711969', 'NWD713365', 'NWD713897', 'NWD714538', 'NWD714615', 'NWD715170', 'NWD715651', 'NWD716459', 'NWD717242', 'NWD717790', 'NWD719508', 'NWD719848', 'NWD719874', 'NWD720697', 'NWD721157', 'NWD722009', 'NWD723332', 'NWD725032', 'NWD725091', 'NWD725749', 'NWD726903', 'NWD727148', 'NWD728463', 'NWD728733', 'NWD728767', 'NWD733033', 'NWD734936', 'NWD736578', 'NWD738611', 'NWD739039', 'NWD739137', 'NWD742146', 'NWD742510', 'NWD742704', 'NWD744599', 'NWD745306', 'NWD747201', 'NWD747359', 'NWD747936', 'NWD748528', 'NWD749045', 'NWD749201', 'NWD750066', 'NWD750367', 'NWD750603', 'NWD752573', 'NWD752766', 'NWD753348', 'NWD753392', 'NWD754350', 'NWD754553', 'NWD754945', 'NWD755427', 'NWD756099', 'NWD757501', 'NWD758118', 'NWD758422', 'NWD759065', 'NWD759385', 'NWD759523', 'NWD760031', 'NWD760562', 'NWD761615', 'NWD761630', 'NWD763104', 'NWD763273', 'NWD763355', 'NWD764908', 'NWD766429', 'NWD766832', 'NWD768186', 'NWD768300', 'NWD769210', 'NWD769529', 'NWD769895', 'NWD770124', 'NWD772586', 'NWD773102', 'NWD774221', 'NWD774542', 'NWD774778', 'NWD775436', 'NWD776037', 'NWD777167', 'NWD777457', 'NWD777673', 'NWD778133', 'NWD778238', 'NWD778520', 'NWD779263', 'NWD779965', 'NWD780209', 'NWD781654', 'NWD783464', 'NWD783729', 'NWD785483', 'NWD785677', 'NWD786446', 'NWD786492', 'NWD789089', 'NWD790493', 'NWD791095', 'NWD791632', 'NWD791684', 'NWD793514', 'NWD795868', 'NWD797431', 'NWD798390', 'NWD799496', 'NWD799614', 'NWD800626', 'NWD801155', 'NWD801626', 'NWD802887', 'NWD805216', 'NWD807363', 'NWD809752', 'NWD813046', 'NWD814098', 'NWD814668', 'NWD816080', 'NWD816911', 'NWD817855', 'NWD819417', 'NWD820660', 'NWD820851', 'NWD821282', 'NWD821604', 'NWD821739', 'NWD821992', 'NWD821993', 'NWD822824', 'NWD822900', 'NWD823550', 'NWD824056', 'NWD825811', 'NWD827563', 'NWD831095', 'NWD831552', 'NWD832508', 'NWD832984', 'NWD834150', 'NWD834863', 'NWD835305', 'NWD836451', 'NWD837230', 'NWD837577', 'NWD837778', 'NWD837977', 'NWD838859', 'NWD839383', 'NWD839456', 'NWD839918', 'NWD840684', 'NWD842405', 'NWD843294', 'NWD844090', 'NWD844341', 'NWD844412', 'NWD844871', 'NWD845729', 'NWD846259', 'NWD846536', 'NWD846973', 'NWD846997', 'NWD847072', 'NWD847891', 'NWD847957', 'NWD849039', 'NWD850100', 'NWD850358', 'NWD850578', 'NWD850663', 'NWD851078', 'NWD852195', 'NWD853192', 'NWD854795', 'NWD854939', 'NWD857549', 'NWD857942', 'NWD858031', 'NWD859870', 'NWD859907', 'NWD860129', 'NWD861279', 'NWD861728', 'NWD861971', 'NWD862659', 'NWD863300', 'NWD865032', 'NWD869379', 'NWD869942', 'NWD870385', 'NWD871339', 'NWD871893', 'NWD872267', 'NWD872277', 'NWD873387', 'NWD873405', 'NWD873586', 'NWD874216', 'NWD874452', 'NWD875385', 'NWD875670', 'NWD876199', 'NWD876644', 'NWD877738', 'NWD878254', 'NWD878522', 'NWD879942', 'NWD880934', 'NWD881030', 'NWD881232', 'NWD881966', 'NWD883600', 'NWD885054', 'NWD885299', 'NWD887396', 'NWD889562', 'NWD890567', 'NWD891671', 'NWD892650', 'NWD892849', 'NWD892920', 'NWD894144', 'NWD894713', 'NWD896033', 'NWD896397', 'NWD896425', 'NWD896853', 'NWD897542', 'NWD898191', 'NWD899238', 'NWD899289', 'NWD900321', 'NWD902332', 'NWD904155', 'NWD904536', 'NWD905958', 'NWD906714', 'NWD907232', 'NWD907419', 'NWD908438', 'NWD909195', 'NWD910427', 'NWD911068', 'NWD911824', 'NWD912179', 'NWD912399', 'NWD912858', 'NWD913189', 'NWD914437', 'NWD915415', 'NWD915924', 'NWD916037', 'NWD916171', 'NWD918288', 'NWD919632', 'NWD921721', 'NWD921983', 'NWD922698', 'NWD922769', 'NWD924608', 'NWD925651', 'NWD926388', 'NWD928800', 'NWD928861', 'NWD930344', 'NWD930731', 'NWD932395', 'NWD932551', 'NWD932840', 'NWD933034', 'NWD933750', 'NWD934376', 'NWD935750', 'NWD935875', 'NWD936747', 'NWD937099', 'NWD937453', 'NWD938588', 'NWD938973', 'NWD943777', 'NWD944915', 'NWD946382', 'NWD947011', 'NWD947532', 'NWD948110', 'NWD948417', 'NWD949178', 'NWD950511', 'NWD951431', 'NWD953149', 'NWD954039', 'NWD954371', 'NWD955193', 'NWD955616', 'NWD956638', 'NWD956711', 'NWD957257', 'NWD957585', 'NWD958609', 'NWD958879', 'NWD959156', 'NWD959224', 'NWD959688', 'NWD959729', 'NWD961340', 'NWD962202', 'NWD962652', 'NWD962995', 'NWD963432', 'NWD964055', 'NWD964195', 'NWD964770', 'NWD968677', 'NWD969993', 'NWD973067', 'NWD973374', 'NWD973800', 'NWD974018', 'NWD974183', 'NWD976658', 'NWD976662', 'NWD977403', 'NWD978533', 'NWD979105', 'NWD979239', 'NWD980594', 'NWD981120', 'NWD981931', 'NWD982856', 'NWD984489', 'NWD985173', 'NWD987200', 'NWD989436', 'NWD989592', 'NWD991289', 'NWD992137', 'NWD993368', 'NWD993544', 'NWD997376', 'NWD998627'] ##add sample names from vcf

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
working_dir = '/home/atimms/ngs_data/genomes/dave_lisa_chd_genomes_0921'
os.chdir(working_dir)
project_name = 'chd_genomes_0921'
combined_vcf = 'PCGC_CHD_phs001735_TOPMed_WGS_freeze.9b.combined.hg38.c1.vcf.gz'
proband_vcf = 'PCGC_CHD_phs001735_TOPMed_WGS_freeze.9b.probands.hg38.c1.vcf.gz'
proband_samples = 'proband_ids.txt'
multianno = project_name + '.hg38_multianno.txt'
annotated = project_name + '.hg38_annotated.txt'
exonic_rare_vars = project_name + '.exonic_rare.txt'


##step1. filter vcf for proband samples
# split_vcf_by_sample(combined_vcf, proband_samples, proband_vcf)

##step2. annotate with annovar and sormat
# annotate_vcf_with_annovar(proband_vcf, project_name)
multi_to_ann(multianno, annotated)

##step3. filter exonic rare vars
filter_exon_rare_vars(annotated, exonic_rare_vars)

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
carm1_vars = project_name + '.exonic_rare.carm1.xls'
filter_by_genenames(exonic_rare_vars, lisa_genes, lisa_genes_vars)
filter_by_genenames(exonic_rare_vars, chd_genes, chd_genes_vars)
filter_by_genenames(exonic_rare_vars, carm1, carm1_vars)

##step5. modify chd vars file to gene/var by sample
chd_genes_vars_formatted = project_name + '.exonic_rare.chd_genes.by_proband.xls'
format_vars_by_proband(chd_genes_vars, chd_genes_vars_formatted)






