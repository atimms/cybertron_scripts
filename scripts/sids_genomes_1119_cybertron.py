#!/usr/bin/env python
import sys
import subprocess
import os
import glob
import gzip

##modules 
'''
module load biobuilds/2017.11
'''

##parameters
delim = '\t'

##file names etc
fa_file = '/home/atimms/ngs_data/references/hg19/hg19_ucsc.fa'

##programs
gatk = '/cm/shared/apps/GATK/3.7/GenomeAnalysisTK.jar'
table_annovar = '/home/atimms/programs/annovar_0618/table_annovar.pl'
convert_2_annovar = '/home/atimms/programs/annovar/convert2annovar.pl'

#annovar parameters
av_genome = 'hg19'
av_buildver = ['-buildver', av_genome]
av_ref_dir = ['/home/atimms/ngs_data/references/annovar/' + av_genome]
av_protocol = ['-protocol', 'refGene,rmsk,genomicSuperDups,dbnsfp35a,dbnsfp31a_interpro,avsnp147,gnomad211_genome,gnomad211_exome,clinvar_20190305']
av_operation = ['-operation', 'g,r,r,f,f,f,f,f,f']
av_options = ['-otherinfo', '-remove', '-nastring', '.','-vcfinput']


##cases id dict
vg_ids = ['DPXQ372', 'VW8M88B', 'IA4PD32', '3G8D41V', 'SZ8YKSK', 'R9P84BJ', 'U1D61QV', '93JJUT6', 'SS9BTS7', 'YQQN44F', '8G159CU', 'IE986KL', 'UAN8LRY', 
		'HBAWAKX', 'X16XMME', 'B8N26R5', '7VLXKED', 'CYDZKQ1', 'P9Q9KW6', '7RGERVT', 'L6X63ZW', '9CULGUA', '1ZZP8JY', '5YM9XYM', 'BJS8R93', 'Q82C2CL', 
		'NUSZ2WD', 'W3NKD7Z', 'MGDE4YA', '9SLL7T3', 'DEYSE7Y', 'GIX3TM2', 'Z2B53KP', 'U2JDX9Q', '5MI4DGP', 'B1QFERQ', 'NRTN6WP', 'ZT7AD3F', 'GDSJ14I', 
		'IM2L5K1', 'E5NQT7U', '8A7QEVB', 'UTEI9QG', 'RYRAGCH', '2RNQCNR', 'RVRYKCS', 'JDKBM22', 'BDZSWSJ']
lrs = ['LR18-093', 'LR18-098', 'LR18-099', 'LR18-101', 'LR18-102', 'LR18-103', 'LR18-106', 'LR18-122', 'LR18-125', 'LR18-126', 'LR18-166', 'LR18-167', 
		'LR18-168', 'LR18-169', 'LR18-172', 'LR18-173', 'LR18-177', 'LR18-179', 'LR18-180', 'LR18-183', 'LR18-186', 'LR18-187', 'LR18-192', 'LR18-195', 
		'LR18-198', 'LR18-200', 'LR18-244', 'LR18-245', 'LR18-247', 'LR18-250', 'LR18-252', 'LR18-253', 'LR18-255', 'LR18-260', 'LR18-264', 'LR18-272', 
		'LR18-280', 'LR18-285', 'LR18-287', 'LR18-288', 'LR18-294', 'LR18-295', 'LR18-296', 'LR18-297', 'LR18-299', 'LR18-303', 'LR18-304', 'LR18-305']
case_id_dict = dict(zip(vg_ids, lrs))

##case gvcfs to change
case_gvcfs_to_change = ['Sample_DPXQ372-EXT.g.vcf.gz', 'Sample_VW8M88B-EXT.g.vcf.gz', 'Sample_IA4PD32-EXT.g.vcf.gz', 'Sample_3G8D41V-EXT.g.vcf.gz', 'Sample_SZ8YKSK-EXT.g.vcf.gz', 
	'Sample_R9P84BJ-EXT.g.vcf.gz', 'Sample_U1D61QV-EXT.g.vcf.gz', 'Sample_93JJUT6-EXT.g.vcf.gz', 'Sample_SS9BTS7-EXT.g.vcf.gz', 'Sample_YQQN44F-EXT.g.vcf.gz', 'Sample_8G159CU-EXT.g.vcf.gz', 
	'Sample_IE986KL-EXT.g.vcf.gz', 'Sample_UAN8LRY-EXT.g.vcf.gz', 'Sample_HBAWAKX-EXT.g.vcf.gz', 'Sample_X16XMME-EXT.g.vcf.gz', 'Sample_B8N26R5-EXT.g.vcf.gz', 'Sample_7VLXKED-EXT.g.vcf.gz', 
	'Sample_CYDZKQ1-EXT.g.vcf.gz', 'Sample_P9Q9KW6-EXT.g.vcf.gz', 'Sample_7RGERVT-EXT.g.vcf.gz', 'Sample_L6X63ZW-EXT.g.vcf.gz', 'Sample_9CULGUA-EXT.g.vcf.gz', 'Sample_1ZZP8JY-EXT.g.vcf.gz', 
	'Sample_5YM9XYM-EXT.g.vcf.gz', 'Sample_BJS8R93-EXT.g.vcf.gz', 'Sample_Q82C2CL-EXT.g.vcf.gz', 'Sample_NUSZ2WD-EXT.g.vcf.gz', 'Sample_W3NKD7Z-EXT.g.vcf.gz', 'Sample_MGDE4YA-EXT.g.vcf.gz', 
	'Sample_9SLL7T3-EXT.g.vcf.gz', 'Sample_DEYSE7Y-EXT.g.vcf.gz', 'Sample_GIX3TM2-EXT.g.vcf.gz', 'Sample_Z2B53KP-EXT.g.vcf.gz', 'Sample_U2JDX9Q-EXT.g.vcf.gz', 'Sample_5MI4DGP-EXT.g.vcf.gz', 
	'Sample_B1QFERQ-EXT.g.vcf.gz', 'Sample_NRTN6WP-EXT.g.vcf.gz', 'Sample_ZT7AD3F-EXT.g.vcf.gz', 'Sample_GDSJ14I-EXT.g.vcf.gz', 'Sample_IM2L5K1-EXT.g.vcf.gz', 'Sample_E5NQT7U-EXT.g.vcf.gz', 
	'Sample_8A7QEVB-EXT.g.vcf.gz', 'Sample_UTEI9QG-EXT.g.vcf.gz', 'Sample_RYRAGCH-EXT.g.vcf.gz', 'Sample_2RNQCNR-EXT.g.vcf.gz', 'Sample_RVRYKCS-EXT.g.vcf.gz', 'Sample_JDKBM22-EXT.g.vcf.gz', 
	'Sample_BDZSWSJ-EXT.g.vcf.gz']
##control gvcfs to change
control_gvcfs_to_change = ['Sample_BV1XYUL-EXT.g.vcf.gz', 'Sample_M7XITGP-EXT.g.vcf.gz', 'Sample_YX2RKMN-EXT.g.vcf.gz', 'Sample_RJZMFUH-EXT.g.vcf.gz', 'Sample_PBWT8WT-EXT.g.vcf.gz', 
		'Sample_25ZDV2G-EXT.g.vcf.gz', 'Sample_IRY7WMC-EXT.g.vcf.gz', 'Sample_Y19276I-EXT.g.vcf.gz', 'Sample_9KM25YL-EXT.g.vcf.gz', 'Sample_7EXPMGU-EXT.g.vcf.gz', 'Sample_NIWT8YK-EXT.g.vcf.gz', 
		'Sample_YDNWM6Z-EXT.g.vcf.gz', 'Sample_ZY942NL-EXT.g.vcf.gz', 'Sample_TVQ1SV4-EXT.g.vcf.gz', 'Sample_HD1YDNK-EXT.g.vcf.gz', 'Sample_QSZMFW9-EXT.g.vcf.gz', 'Sample_L2F8A3F-EXT.g.vcf.gz', 
		'Sample_HYFIANF-EXT.g.vcf.gz', 'Sample_CWYIVA4-EXT.g.vcf.gz', 'Sample_26UWW3S-EXT.g.vcf.gz', 'Sample_EMV3TSR-EXT.g.vcf.gz', 'Sample_UHXLYCT-EXT.g.vcf.gz', 'Sample_TNT45T4-EXT.g.vcf.gz', 
		'Sample_AKRVUE3-EXT.g.vcf.gz', 'Sample_B7FSEU3-EXT.g.vcf.gz', 'Sample_ATD2WVM-EXT.g.vcf.gz', 'Sample_YWNZURX-EXT.g.vcf.gz', 'Sample_SVIEAEC-EXT.g.vcf.gz', 'Sample_CL1J8B1-EXT.g.vcf.gz', 
		'Sample_YYFY8PN-EXT.g.vcf.gz', 'Sample_576PVIT-EXT.g.vcf.gz', 'Sample_F2GUM9Q-EXT.g.vcf.gz', 'Sample_VKTSS9T-EXT.g.vcf.gz', 'Sample_A8R78VT-EXT.g.vcf.gz', 'Sample_PF44VGN-EXT.g.vcf.gz', 
		'Sample_QW3G4HF-EXT.g.vcf.gz', 'Sample_RZLTHDT-EXT.g.vcf.gz', 'Sample_FPG61PF-EXT.g.vcf.gz', 'Sample_CGAAYCL-EXT.g.vcf.gz', 'Sample_MTADLK7-EXT.g.vcf.gz', 'Sample_FDESA8S-EXT.g.vcf.gz', 
		'Sample_IGFTTKA-EXT.g.vcf.gz', 'Sample_IS99ILP-EXT.g.vcf.gz', 'Sample_HBJNU6D-EXT.g.vcf.gz', 'Sample_T8FXUUX-EXT.g.vcf.gz', 'Sample_4GEZMHM-EXT.g.vcf.gz', 'Sample_8Q1HZWE-EXT.g.vcf.gz', 
		'Sample_XDXP7RF-EXT.g.vcf.gz', 'Sample_1GL82K3-EXT.g.vcf.gz', 'Sample_GICP3NZ-EXT.g.vcf.gz', 'Sample_8F2JCWB-EXT.g.vcf.gz', 'Sample_5Y4M12V-EXT.g.vcf.gz', 'Sample_UHAFECG-EXT.g.vcf.gz', 
		'Sample_RAMIYE7-EXT.g.vcf.gz', 'Sample_R3NW7ET-EXT.g.vcf.gz', 'Sample_B68SMAX-EXT.g.vcf.gz', 'Sample_KDZXW3Z-EXT.g.vcf.gz', 'Sample_IW8KELD-EXT.g.vcf.gz', 'Sample_MFSNMIF-EXT.g.vcf.gz', 
		'Sample_L21L11I-EXT.g.vcf.gz', 'Sample_DB9LUAI-EXT.g.vcf.gz', 'Sample_5J83LZG-EXT.g.vcf.gz', 'Sample_E1R7R8G-EXT.g.vcf.gz', 'Sample_H3D2CMT-EXT.g.vcf.gz', 'Sample_D2WENRB-EXT.g.vcf.gz', 
		'Sample_MG115JD-EXT.g.vcf.gz', 'Sample_LBNI329-EXT.g.vcf.gz', 'Sample_WULHIQ2-EXT.g.vcf.gz', 'Sample_C6MEWBU-EXT.g.vcf.gz', 'Sample_36P2R2Q-EXT.g.vcf.gz', 'Sample_1FSGTK2-EXT.g.vcf.gz', 
		'Sample_RWCSRWK-EXT.g.vcf.gz', 'Sample_9FYGZXJ-EXT.g.vcf.gz', 'Sample_DYMEXA3-EXT.g.vcf.gz', 'Sample_UCS36U9-EXT.g.vcf.gz', 'Sample_4JEBIHC-EXT.g.vcf.gz', 'Sample_F2RJIQ9-EXT.g.vcf.gz', 
		'Sample_WBKH191-EXT.g.vcf.gz', 'Sample_T4AE2BD-EXT.g.vcf.gz', 'Sample_TA9N2TY-EXT.g.vcf.gz', 'Sample_6C4LNY1-EXT.g.vcf.gz', 'Sample_KIABE21-EXT.g.vcf.gz', 'Sample_C5E2FAW-EXT.g.vcf.gz', 
		'Sample_42WCDJS-EXT.g.vcf.gz', 'Sample_WZR1B8L-EXT.g.vcf.gz', 'Sample_SK2S1DK-EXT.g.vcf.gz', 'Sample_HBHW153-EXT.g.vcf.gz', 'Sample_6RCJ3FC-EXT.g.vcf.gz', 'Sample_2H1TBM1-EXT.g.vcf.gz', 
		'Sample_YQA4J7F-EXT.g.vcf.gz', 'Sample_3HGDVJ1-EXT.g.vcf.gz', 'Sample_LYNUGIZ-EXT.g.vcf.gz', 'Sample_UWHESTR-EXT.g.vcf.gz', 'Sample_J6LRK6K-EXT.g.vcf.gz', 'Sample_X8KJVR4-EXT.g.vcf.gz', 
		'Sample_4JRS816-EXT.g.vcf.gz', 'Sample_FBIXHPN-EXT.g.vcf.gz', 'Sample_ZF6BQ6A-EXT.g.vcf.gz', 'Sample_4X2IGIS-EXT.g.vcf.gz', 'Sample_GT6S28L-EXT.g.vcf.gz', 'Sample_JHHXE4D-EXT.g.vcf.gz', 
		'Sample_17ISFMC-EXT.g.vcf.gz', 'Sample_H544FPH-EXT.g.vcf.gz', 'Sample_RN568C1-EXT.g.vcf.gz', 'Sample_ADK9CUJ-EXT.g.vcf.gz', 'Sample_HGZ99N9-EXT.g.vcf.gz', 'Sample_AY4B7U3-EXT.g.vcf.gz', 
		'Sample_AC28ZC2-EXT.g.vcf.gz', 'Sample_9B9UYEE-EXT.g.vcf.gz', 'Sample_63D1I1W-EXT.g.vcf.gz', 'Sample_I41KM4Z-EXT.g.vcf.gz', 'Sample_F9Q6AQM-EXT.g.vcf.gz', 'Sample_2NA5X4X-EXT.g.vcf.gz', 
		'Sample_U66AUUI-EXT.g.vcf.gz', 'Sample_EJVRXS3-EXT.g.vcf.gz', 'Sample_6H9H8YD-EXT.g.vcf.gz', 'Sample_B9JK8SG-EXT.g.vcf.gz', 'Sample_F93MY8G-EXT.g.vcf.gz', 'Sample_6W2W2HN-EXT.g.vcf.gz', 
		'Sample_43H3VZN-EXT.g.vcf.gz', 'Sample_IUSBV56-EXT.g.vcf.gz', 'Sample_CBSXQUE-EXT.g.vcf.gz', 'Sample_B1CNNV1-EXT.g.vcf.gz', 'Sample_TK7LPT2-EXT.g.vcf.gz', 'Sample_QJS4LGM-EXT.g.vcf.gz', 
		'Sample_X5D8979-EXT.g.vcf.gz', 'Sample_RZB4MWA-EXT.g.vcf.gz', 'Sample_ZAD9M6J-EXT.g.vcf.gz', 'Sample_DLVP2AY-EXT.g.vcf.gz', 'Sample_M9TELZU-EXT.g.vcf.gz', 'Sample_MSJCIHJ-EXT.g.vcf.gz', 
		'Sample_ZD5KRL3-EXT.g.vcf.gz', 'Sample_Q1HBJE5-EXT.g.vcf.gz', 'Sample_I29M969-EXT.g.vcf.gz', 'Sample_EAGUMBH-EXT.g.vcf.gz', 'Sample_3H4W627-EXT.g.vcf.gz', 'Sample_JBTD83B-EXT.g.vcf.gz', 
		'Sample_CG8I5BB-EXT.g.vcf.gz', 'Sample_V9PL8QV-EXT.g.vcf.gz', 'Sample_5UUJKGQ-EXT.g.vcf.gz', 'Sample_1LG2YL4-EXT.g.vcf.gz', 'Sample_UI6VQ9U-EXT.g.vcf.gz', 'Sample_NCUHC2T-EXT.g.vcf.gz', 
		'Sample_USMJVSQ-EXT.g.vcf.gz', 'Sample_NFYA6Z7-EXT.g.vcf.gz', 'Sample_4PJUB1W-EXT.g.vcf.gz', 'Sample_H6E52NT-EXT.g.vcf.gz', 'Sample_KYJZR3I-EXT.g.vcf.gz', 'Sample_A7G5LWH-EXT.g.vcf.gz', 
		'Sample_Z3J4T5U-EXT.g.vcf.gz', 'Sample_6YTVGGE-EXT.g.vcf.gz', 'Sample_A7WP6TH-EXT.g.vcf.gz', 'Sample_8ELMJWX-EXT.g.vcf.gz', 'Sample_RBSQWX2-EXT.g.vcf.gz', 'Sample_84CZ2FB-EXT.g.vcf.gz', 
		'Sample_TP2GMU3-EXT.g.vcf.gz', 'Sample_DEN3JQF-EXT.g.vcf.gz', 'Sample_SGQMPD9-EXT.g.vcf.gz', 'Sample_GPNSJ78-EXT.g.vcf.gz', 'Sample_TDXA3C5-EXT.g.vcf.gz', 'Sample_3NRGB39-EXT.g.vcf.gz', 
		'Sample_QVCFZES-EXT.g.vcf.gz', 'Sample_YMBTN7Q-EXT.g.vcf.gz', 'Sample_RUB2SCE-EXT.g.vcf.gz', 'Sample_FC4RNAG-EXT.g.vcf.gz', 'Sample_P5TPIXR-EXT.g.vcf.gz', 'Sample_H6C1G5Q-EXT.g.vcf.gz', 
		'Sample_7EE11F6-EXT.g.vcf.gz', 'Sample_RI22IEY-EXT.g.vcf.gz', 'Sample_UYHCXBN-EXT.g.vcf.gz', 'Sample_5U719ZJ-EXT.g.vcf.gz', 'Sample_6BXERG6-EXT.g.vcf.gz', 'Sample_W88PFRH-EXT.g.vcf.gz', 
		'Sample_AWA88C8-EXT.g.vcf.gz', 'Sample_6R2V7YS-EXT.g.vcf.gz', 'Sample_48MLNI8-EXT.g.vcf.gz', 'Sample_KT2MIJB-EXT.g.vcf.gz', 'Sample_JLFGF3R-EXT.g.vcf.gz', 'Sample_A5H7HDK-EXT.g.vcf.gz', 
		'Sample_IPHXD65-EXT.g.vcf.gz', 'Sample_NQVFZYZ-EXT.g.vcf.gz', 'Sample_DE7S6R4-EXT.g.vcf.gz', 'Sample_CTAWVSY-EXT.g.vcf.gz', 'Sample_I5JMZMH-EXT.g.vcf.gz', 'Sample_E82A8RN-EXT.g.vcf.gz', 
		'Sample_6N8R8GY-EXT.g.vcf.gz', 'Sample_QK21NWV-EXT.g.vcf.gz', 'Sample_5IZ6KJ8-EXT.g.vcf.gz', 'Sample_DNBZKR6-EXT.g.vcf.gz', 'Sample_TNRE1WQ-EXT.g.vcf.gz', 'Sample_M7PLS2G-EXT.g.vcf.gz', 
		'Sample_LKFV22N-EXT.g.vcf.gz', 'Sample_D88H59I-EXT.g.vcf.gz', 'Sample_SNCSRUT-EXT.g.vcf.gz', 'Sample_DJZU9RU-EXT.g.vcf.gz', 'Sample_LD5SKJG-EXT.g.vcf.gz', 'Sample_VA8NK9D-EXT.g.vcf.gz', 
		'Sample_I8CQY4B-EXT.g.vcf.gz', 'Sample_Y9QDK7K-EXT.g.vcf.gz', 'Sample_31I572D-EXT.g.vcf.gz', 'Sample_FP4CKPT-EXT.g.vcf.gz', 'Sample_XSLJD85-EXT.g.vcf.gz', 'Sample_P1V6GZC-EXT.g.vcf.gz', 
		'Sample_XVLV97U-EXT.g.vcf.gz', 'Sample_9YPSMW2-EXT.g.vcf.gz', 'Sample_BR2L3UW-EXT.g.vcf.gz', 'Sample_YI1SXQH-EXT.g.vcf.gz', 'Sample_9FJ8GDF-EXT.g.vcf.gz', 'Sample_RT8ZWX8-EXT.g.vcf.gz', 
		'Sample_9SDN6EU-EXT.g.vcf.gz', 'Sample_95LATDC-EXT.g.vcf.gz', 'Sample_IATZHLI-EXT.g.vcf.gz', 'Sample_54UXA1A-EXT.g.vcf.gz', 'Sample_Q9M68XD-EXT.g.vcf.gz', 'Sample_Z72H2N1-EXT.g.vcf.gz', 
		'Sample_8JWD9XC-EXT.g.vcf.gz', 'Sample_XZK657J-EXT.g.vcf.gz', 'Sample_FFTAZNP-EXT.g.vcf.gz', 'Sample_6XNNIXJ-EXT.g.vcf.gz', 'Sample_KAN5BKG-EXT.g.vcf.gz', 'Sample_QFAQDYF-EXT.g.vcf.gz', 
		'Sample_J61B53L-EXT.g.vcf.gz', 'Sample_A5FFMCA-EXT.g.vcf.gz', 'Sample_2HNCL47-EXT.g.vcf.gz', 'Sample_L23EU3T-EXT.g.vcf.gz', 'Sample_QMU4MDP-EXT.g.vcf.gz', 'Sample_FDC1G7H-EXT.g.vcf.gz', 
		'Sample_7SCWTX5-EXT.g.vcf.gz', 'Sample_UEKTDUA-EXT.g.vcf.gz', 'Sample_97RUHD1-EXT.g.vcf.gz', 'Sample_FKUCU8K-EXT.g.vcf.gz', 'Sample_1PGDULT-EXT.g.vcf.gz', 'Sample_3CCLW4R-EXT.g.vcf.gz', 
		'Sample_L5CELJ2-EXT.g.vcf.gz', 'Sample_6ZEPM17-EXT.g.vcf.gz', 'Sample_1DWNQ7U-EXT.g.vcf.gz', 'Sample_MYVGY1R-EXT.g.vcf.gz', 'Sample_ILZXV5R-EXT.g.vcf.gz', 'Sample_P1BP1WD-EXT.g.vcf.gz', 
		'Sample_FXL1P7P-EXT.g.vcf.gz', 'Sample_YLUWV7C-EXT.g.vcf.gz', 'Sample_S27ZFFE-EXT.g.vcf.gz', 'Sample_95BKYWT-EXT.g.vcf.gz', 'Sample_QGR1WFM-EXT.g.vcf.gz', 'Sample_U7BHRDD-EXT.g.vcf.gz', 
		'Sample_YL5KRNV-EXT.g.vcf.gz', 'Sample_8ZND1X6-EXT.g.vcf.gz', 'Sample_DMY55RJ-EXT.g.vcf.gz', 'Sample_G52YV6D-EXT.g.vcf.gz', 'Sample_Q8YMXF8-EXT.g.vcf.gz', 'Sample_DUZJRS8-EXT.g.vcf.gz', 
		'Sample_N9Y8BGB-EXT.g.vcf.gz', 'Sample_GU4DYN4-EXT.g.vcf.gz', 'Sample_PS3SGZX-EXT.g.vcf.gz', 'Sample_5W8R8IR-EXT.g.vcf.gz', 'Sample_QIA29X5-EXT.g.vcf.gz', 'Sample_DSGIDAQ-EXT.g.vcf.gz', 
		'Sample_XKPPG92-EXT.g.vcf.gz', 'Sample_TRTF1TT-EXT.g.vcf.gz', 'Sample_6HUYN2D-EXT.g.vcf.gz', 'Sample_LIM5U2M-EXT.g.vcf.gz', 'Sample_YRMYZ72-EXT.g.vcf.gz', 'Sample_192US5U-EXT.g.vcf.gz', 
		'Sample_75CKFGE-EXT.g.vcf.gz', 'Sample_QT6VCE4-EXT.g.vcf.gz', 'Sample_BG3NFVU-EXT.g.vcf.gz', 'Sample_LZL5BLL-EXT.g.vcf.gz', 'Sample_1VGEP6P-EXT.g.vcf.gz', 'Sample_SHQZFU2-EXT.g.vcf.gz', 
		'Sample_X913FN5-EXT.g.vcf.gz', 'Sample_W57LQPH-EXT.g.vcf.gz', 'Sample_KLSBV3E-EXT.g.vcf.gz', 'Sample_DCQIMAX-EXT.g.vcf.gz', 'Sample_MS7VTZP-EXT.g.vcf.gz', 'Sample_Y1KJWPC-EXT.g.vcf.gz', 
		'Sample_9J4ZSF4-EXT.g.vcf.gz', 'Sample_3QX2Z3W-EXT.g.vcf.gz', 'Sample_KXFJNLY-EXT.g.vcf.gz', 'Sample_I14ECME-EXT.g.vcf.gz', 'Sample_ITH996U-EXT.g.vcf.gz', 'Sample_JAIBM4Z-EXT.g.vcf.gz', 
		'Sample_24BUJKA-EXT.g.vcf.gz', 'Sample_1BDMDNC-EXT.g.vcf.gz', 'Sample_BXFENAI-EXT.g.vcf.gz', 'Sample_EMXSXP5-EXT.g.vcf.gz', 'Sample_MGM6PJR-EXT.g.vcf.gz', 'Sample_2WGRGNH-EXT.g.vcf.gz', 
		'Sample_I7ZWI4P-EXT.g.vcf.gz', 'Sample_Z68QTNZ-EXT.g.vcf.gz', 'Sample_C9R7RA8-EXT.g.vcf.gz', 'Sample_3LB7SL2-EXT.g.vcf.gz', 'Sample_IKM4F55-EXT.g.vcf.gz', 'Sample_XX5WLQC-EXT.g.vcf.gz', 
		'Sample_PMDFHXK-EXT.g.vcf.gz', 'Sample_TM8CND8-EXT.g.vcf.gz', 'Sample_8ZDN6GL-EXT.g.vcf.gz', 'Sample_WDPAU8D-EXT.g.vcf.gz', 'Sample_4V8EH2X-EXT.g.vcf.gz', 'Sample_K6QL9L2-EXT.g.vcf.gz', 
		'Sample_131SXLY-EXT.g.vcf.gz', 'Sample_3IYE83I-EXT.g.vcf.gz', 'Sample_QNSEHGB-EXT.g.vcf.gz', 'Sample_WRL7MQC-EXT.g.vcf.gz', 'Sample_9WI7YWE-EXT.g.vcf.gz', 'Sample_YTMW5PY-EXT.g.vcf.gz', 
		'Sample_W2MW4RQ-EXT.g.vcf.gz', 'Sample_HA9M882-EXT.g.vcf.gz', 'Sample_39Y2HJY-EXT.g.vcf.gz', 'Sample_EZKYUQK-EXT.g.vcf.gz', 'Sample_DIKPBTQ-EXT.g.vcf.gz', 'Sample_ZQ9RC52-EXT.g.vcf.gz', 
		'Sample_XDVWDQ5-EXT.g.vcf.gz', 'Sample_362IFKJ-EXT.g.vcf.gz', 'Sample_6YI6KZV-EXT.g.vcf.gz', 'Sample_FRXFJ6M-EXT.g.vcf.gz', 'Sample_M2URU1G-EXT.g.vcf.gz', 'Sample_7XUA2FB-EXT.g.vcf.gz', 
		'Sample_5HFCEZF-EXT.g.vcf.gz', 'Sample_VXV88TE-EXT.g.vcf.gz', 'Sample_XE24A81-EXT.g.vcf.gz', 'Sample_Q6T29FJ-EXT.g.vcf.gz', 'Sample_VC7LQRA-EXT.g.vcf.gz', 'Sample_W6FY8RG-EXT.g.vcf.gz', 
		'Sample_IJTD854-EXT.g.vcf.gz', 'Sample_1XT7AP8-EXT.g.vcf.gz', 'Sample_UGRE2UY-EXT.g.vcf.gz', 'Sample_EEVP387-EXT.g.vcf.gz', 'Sample_JI3PVJ9-EXT.g.vcf.gz', 'Sample_26HF8KY-EXT.g.vcf.gz', 
		'Sample_UHZNTRZ-EXT.g.vcf.gz', 'Sample_KFVTPL4-EXT.g.vcf.gz', 'Sample_MBM5UZV-EXT.g.vcf.gz', 'Sample_5CN3421-EXT.g.vcf.gz', 'Sample_GJ5FBN1-EXT.g.vcf.gz', 'Sample_5V2JB1V-EXT.g.vcf.gz', 
		'Sample_NSYYSJG-EXT.g.vcf.gz', 'Sample_SZGYBCR-EXT.g.vcf.gz', 'Sample_S1CGDD3-EXT.g.vcf.gz', 'Sample_V5S36SH-EXT.g.vcf.gz', 'Sample_1B6MM47-EXT.g.vcf.gz', 'Sample_4NWB1IQ-EXT.g.vcf.gz', 
		'Sample_2QRFGK5-EXT.g.vcf.gz', 'Sample_U1ZEAU5-EXT.g.vcf.gz', 'Sample_JKH9A43-EXT.g.vcf.gz', 'Sample_P3TRDFV-EXT.g.vcf.gz', 'Sample_TI92SDJ-EXT.g.vcf.gz', 'Sample_X53JEQQ-EXT.g.vcf.gz', 
		'Sample_CAJXYA9-EXT.g.vcf.gz', 'Sample_JHF6J33-EXT.g.vcf.gz', 'Sample_GVYW1PF-EXT.g.vcf.gz', 'Sample_PP9BAGQ-EXT.g.vcf.gz', 'Sample_SVT45VV-EXT.g.vcf.gz', 'Sample_275XI24-EXT.g.vcf.gz', 
		'Sample_2MEMV3L-EXT.g.vcf.gz', 'Sample_UQ5II99-EXT.g.vcf.gz', 'Sample_M1BQHHY-EXT.g.vcf.gz', 'Sample_LPS1D2Z-EXT.g.vcf.gz', 'Sample_9DH7GFD-EXT.g.vcf.gz', 'Sample_YSGP873-EXT.g.vcf.gz', 
		'Sample_HB775NJ-EXT.g.vcf.gz', 'Sample_VQWLGU1-EXT.g.vcf.gz', 'Sample_VTYLGR3-EXT.g.vcf.gz', 'Sample_D6WQJRZ-EXT.g.vcf.gz', 'Sample_YR6ML8R-EXT.g.vcf.gz', 'Sample_KBINDLS-EXT.g.vcf.gz', 
		'Sample_V5QBCQ7-EXT.g.vcf.gz', 'Sample_G1E4NPI-EXT.g.vcf.gz', 'Sample_1KTIN5X-EXT.g.vcf.gz', 'Sample_RX6IYWL-EXT.g.vcf.gz', 'Sample_IX2BLLE-EXT.g.vcf.gz', 'Sample_VBQC7A3-EXT.g.vcf.gz', 
		'Sample_WZNB7B8-EXT.g.vcf.gz', 'Sample_E16S1AD-EXT.g.vcf.gz', 'Sample_4FL9FHL-EXT.g.vcf.gz', 'Sample_U6GZPB2-EXT.g.vcf.gz', 'Sample_I72NC5Z-EXT.g.vcf.gz', 'Sample_9TJW3WP-EXT.g.vcf.gz', 
		'Sample_3P3IX2K-EXT.g.vcf.gz', 'Sample_XVYCYQN-EXT.g.vcf.gz', 'Sample_VY6AUSQ-EXT.g.vcf.gz', 'Sample_P2TEMX2-EXT.g.vcf.gz', 'Sample_PQ4UCH2-EXT.g.vcf.gz', 'Sample_E61J97H-EXT.g.vcf.gz', 
		'Sample_AJND2TL-EXT.g.vcf.gz', 'Sample_XQ4H1PN-EXT.g.vcf.gz', 'Sample_KTG8SL9-EXT.g.vcf.gz', 'Sample_G4EFJP8-EXT.g.vcf.gz', 'Sample_UHC5IAT-EXT.g.vcf.gz', 'Sample_EQZVMQ5-EXT.g.vcf.gz', 
		'Sample_GFZWI6G-EXT.g.vcf.gz', 'Sample_F2AMQQU-EXT.g.vcf.gz', 'Sample_N35TGYS-EXT.g.vcf.gz', 'Sample_V3LIIST-EXT.g.vcf.gz', 'Sample_XMLHJN9-EXT.g.vcf.gz', 'Sample_7N2W3GW-EXT.g.vcf.gz', 
		'Sample_BIJXYC1-EXT.g.vcf.gz', 'Sample_VAVU59R-EXT.g.vcf.gz', 'Sample_7KQ4HXD-EXT.g.vcf.gz', 'Sample_1X8QUL8-EXT.g.vcf.gz', 'Sample_I27UF4Y-EXT.g.vcf.gz', 'Sample_XTFAL77-EXT.g.vcf.gz', 
		'Sample_WAPYY8P-EXT.g.vcf.gz', 'Sample_EGUN7Q3-EXT.g.vcf.gz', 'Sample_EU5EF9P-EXT.g.vcf.gz', 'Sample_29WZL5S-EXT.g.vcf.gz', 'Sample_BDP41B1-EXT.g.vcf.gz', 'Sample_CN75WBP-EXT.g.vcf.gz', 
		'Sample_L8T1VI2-EXT.g.vcf.gz', 'Sample_NT78JGI-EXT.g.vcf.gz', 'Sample_WPV89SL-EXT.g.vcf.gz', 'Sample_JQABD4S-EXT.g.vcf.gz', 'Sample_4UPQV19-EXT.g.vcf.gz', 'Sample_G13MZ6P-EXT.g.vcf.gz', 
		'Sample_FC2YU86-EXT.g.vcf.gz', 'Sample_QLURVWW-EXT.g.vcf.gz', 'Sample_F47T27G-EXT.g.vcf.gz', 'Sample_NUDGHZD-EXT.g.vcf.gz', 'Sample_IXPH7LS-EXT.g.vcf.gz', 'Sample_CLQVDUH-EXT.g.vcf.gz', 
		'Sample_6QECWGM-EXT.g.vcf.gz', 'Sample_9K2KPVL-EXT.g.vcf.gz', 'Sample_QAFWGWF-EXT.g.vcf.gz', 'Sample_T1GB3UJ-EXT.g.vcf.gz', 'Sample_XUQC86I-EXT.g.vcf.gz', 'Sample_QASD5F9-EXT.g.vcf.gz', 
		'Sample_PIFWFY7-EXT.g.vcf.gz', 'Sample_WCF13BC-EXT.g.vcf.gz', 'Sample_GQL1Q5W-EXT.g.vcf.gz', 'Sample_9E96VE4-EXT.g.vcf.gz', 'Sample_RBF97E8-EXT.g.vcf.gz', 'Sample_XHTFFNJ-EXT.g.vcf.gz', 
		'Sample_WIUTNQX-EXT.g.vcf.gz', 'Sample_7VB8PXT-EXT.g.vcf.gz', 'Sample_FKN4XPQ-EXT.g.vcf.gz', 'Sample_NZMTYXE-EXT.g.vcf.gz', 'Sample_YP5WMNJ-EXT.g.vcf.gz', 'Sample_5KR5YHY-EXT.g.vcf.gz', 
		'Sample_VZ9QX9A-EXT.g.vcf.gz', 'Sample_PKW6YGE-EXT.g.vcf.gz', 'Sample_JFYW2LV-EXT.g.vcf.gz', 'Sample_LR9AWI7-EXT.g.vcf.gz', 'Sample_65UA2H4-EXT.g.vcf.gz', 'Sample_RZ44VC4-EXT.g.vcf.gz', 
		'Sample_VW6X4BX-EXT.g.vcf.gz', 'Sample_9ZUZJFX-EXT.g.vcf.gz', 'Sample_65MAAXX-EXT.g.vcf.gz', 'Sample_2T3724N-EXT.g.vcf.gz', 'Sample_PJ1NWE3-EXT.g.vcf.gz', 'Sample_AGLC2WI-EXT.g.vcf.gz', 
		'Sample_UNN9ZS2-EXT.g.vcf.gz', 'Sample_FRZ8D8X-EXT.g.vcf.gz', 'Sample_Z5VWEND-EXT.g.vcf.gz', 'Sample_32U7T1P-EXT.g.vcf.gz', 'Sample_2WEZLL7-EXT.g.vcf.gz', 'Sample_3MYQ438-EXT.g.vcf.gz', 
		'Sample_RJK4VWG-EXT.g.vcf.gz', 'Sample_46VN9JH-EXT.g.vcf.gz', 'Sample_RF34DV3-EXT.g.vcf.gz', 'Sample_YGP1C7Z-EXT.g.vcf.gz', 'Sample_M7C54IM-EXT.g.vcf.gz', 'Sample_89RTLVW-EXT.g.vcf.gz', 
		'Sample_8ALBNW8-EXT.g.vcf.gz', 'Sample_2DG2ZKC-EXT.g.vcf.gz', 'Sample_VC1ET9E-EXT.g.vcf.gz', 'Sample_MKXJ7IS-EXT.g.vcf.gz', 'Sample_ZLMW5M6-EXT.g.vcf.gz', 'Sample_EYZ8DAQ-EXT.g.vcf.gz', 
		'Sample_MQDSUHV-EXT.g.vcf.gz', 'Sample_1AVK15V-EXT.g.vcf.gz', 'Sample_KSX7F3R-EXT.g.vcf.gz', 'Sample_HYKZ8M4-EXT.g.vcf.gz', 'Sample_P58BSZP-EXT.g.vcf.gz', 'Sample_Q2DULGG-EXT.g.vcf.gz', 
		'Sample_KZCQZ2J-EXT.g.vcf.gz', 'Sample_7163MYU-EXT.g.vcf.gz', 'Sample_Y9DWVPR-EXT.g.vcf.gz', 'Sample_DFJLLRR-EXT.g.vcf.gz', 'Sample_F7LDFR9-EXT.g.vcf.gz', 'Sample_BP8I4D3-EXT.g.vcf.gz', 
		'Sample_TQ8NJDX-EXT.g.vcf.gz', 'Sample_UT6L8B8-EXT.g.vcf.gz', 'Sample_KJSCRKI-EXT.g.vcf.gz', 'Sample_GWKMH5B-EXT.g.vcf.gz', 'Sample_MQQLAHI-EXT.g.vcf.gz', 'Sample_XEEYQ8M-EXT.g.vcf.gz', 
		'Sample_7VNPEFN-EXT.g.vcf.gz', 'Sample_JKTQYMW-EXT.g.vcf.gz', 'Sample_HV88N4K-EXT.g.vcf.gz', 'Sample_YFTHA6M-EXT.g.vcf.gz', 'Sample_RP1PADD-EXT.g.vcf.gz', 'Sample_NAGXXH1-EXT.g.vcf.gz', 
		'Sample_4RAVF3J-EXT.g.vcf.gz', 'Sample_RL2DEEN-EXT.g.vcf.gz', 'Sample_2HY2HKQ-EXT.g.vcf.gz', 'Sample_63NPDHF-EXT.g.vcf.gz', 'Sample_7ZBJKWI-EXT.g.vcf.gz', 'Sample_HQBBUND-EXT.g.vcf.gz', 
		'Sample_LLEGYH6-EXT.g.vcf.gz', 'Sample_ACRI4UI-EXT.g.vcf.gz', 'Sample_K2IBM27-EXT.g.vcf.gz', 'Sample_YFISFP4-EXT.g.vcf.gz', 'Sample_AR7G9VZ-EXT.g.vcf.gz', 'Sample_3M9FYJR-EXT.g.vcf.gz', 
		'Sample_LRP7PIL-EXT.g.vcf.gz', 'Sample_Y8UVI69-EXT.g.vcf.gz', 'Sample_3SEBIJ4-EXT.g.vcf.gz', 'Sample_P8WHCW4-EXT.g.vcf.gz', 'Sample_BCGT9DY-EXT.g.vcf.gz', 'Sample_ZB1SXNQ-EXT.g.vcf.gz', 
		'Sample_BZEDTSE-EXT.g.vcf.gz', 'Sample_TTAPJAZ-EXT.g.vcf.gz', 'Sample_EETX97W-EXT.g.vcf.gz', 'Sample_RTBP2VK-EXT.g.vcf.gz', 'Sample_21Z2Z2R-EXT.g.vcf.gz', 'Sample_UJTE2R1-EXT.g.vcf.gz', 
		'Sample_VYE7V8Y-EXT.g.vcf.gz', 'Sample_4E9EZHZ-EXT.g.vcf.gz', 'Sample_WP45A7U-EXT.g.vcf.gz', 'Sample_V7I2KQ8-EXT.g.vcf.gz', 'Sample_RH7JGCM-EXT.g.vcf.gz', 'Sample_S4QY3TZ-EXT.g.vcf.gz', 
		'Sample_S3KR6B4-EXT.g.vcf.gz', 'Sample_XIL7NNK-EXT.g.vcf.gz', 'Sample_NPCEMFH-EXT.g.vcf.gz', 'Sample_HHL1Q45-EXT.g.vcf.gz', 'Sample_1JK664Z-EXT.g.vcf.gz', 'Sample_KDFGG1Z-EXT.g.vcf.gz', 
		'Sample_6E8EIXE-EXT.g.vcf.gz', 'Sample_JZZWJ2W-EXT.g.vcf.gz', 'Sample_J45FCJG-EXT.g.vcf.gz', 'Sample_1IFY9K4-EXT.g.vcf.gz', 'Sample_Z8M8I4W-EXT.g.vcf.gz', 'Sample_6Z4X2FR-EXT.g.vcf.gz', 
		'Sample_M7FT8F1-EXT.g.vcf.gz', 'Sample_BKBWCAR-EXT.g.vcf.gz', 'Sample_VRFK48K-EXT.g.vcf.gz', 'Sample_KAF6K1A-EXT.g.vcf.gz', 'Sample_JJWF3JA-EXT.g.vcf.gz', 'Sample_MBKC1YK-EXT.g.vcf.gz', 
		'Sample_2DRRU2V-EXT.g.vcf.gz', 'Sample_BNA88AG-EXT.g.vcf.gz', 'Sample_R1XU4C5-EXT.g.vcf.gz', 'Sample_1R7B9KJ-EXT.g.vcf.gz', 'Sample_APGF6TB-EXT.g.vcf.gz', 'Sample_U4PXL9E-EXT.g.vcf.gz', 
		'Sample_PRMTZVL-EXT.g.vcf.gz', 'Sample_VVEVZ8A-EXT.g.vcf.gz', 'Sample_GJXFK4V-EXT.g.vcf.gz', 'Sample_QMHA7D3-EXT.g.vcf.gz', 'Sample_DG5C37L-EXT.g.vcf.gz', 'Sample_6SULFYT-EXT.g.vcf.gz', 
		'Sample_YFZP7NI-EXT.g.vcf.gz', 'Sample_2EM9W37-EXT.g.vcf.gz', 'Sample_9ZS8QDM-EXT.g.vcf.gz', 'Sample_D9V1FQP-EXT.g.vcf.gz', 'Sample_F84FSAR-EXT.g.vcf.gz', 'Sample_LKVELZP-EXT.g.vcf.gz', 
		'Sample_AKH49TM-EXT.g.vcf.gz', 'Sample_QBD5MV4-EXT.g.vcf.gz', 'Sample_CWLNFAH-EXT.g.vcf.gz', 'Sample_AQCY7UM-EXT.g.vcf.gz', 'Sample_AH74ICE-EXT.g.vcf.gz', 'Sample_Z6IFP5J-EXT.g.vcf.gz', 
		'Sample_YMYCYNW-EXT.g.vcf.gz', 'Sample_UXLUVAB-EXT.g.vcf.gz', 'Sample_5YCJ2H4-EXT.g.vcf.gz', 'Sample_5X6C5Y9-EXT.g.vcf.gz', 'Sample_JTE4836-EXT.g.vcf.gz', 'Sample_5M5PD2S-EXT.g.vcf.gz', 
		'Sample_GIJPU85-EXT.g.vcf.gz', 'Sample_5J6DH34-EXT.g.vcf.gz', 'Sample_J93EU5K-EXT.g.vcf.gz', 'Sample_YGKIE8B-EXT.g.vcf.gz', 'Sample_T4HESVJ-EXT.g.vcf.gz', 'Sample_Z89GS87-EXT.g.vcf.gz', 
		'Sample_BHLQSDB-EXT.g.vcf.gz', 'Sample_XQBNN9S-EXT.g.vcf.gz', 'Sample_566WRIT-EXT.g.vcf.gz', 'Sample_XH1M28P-EXT.g.vcf.gz', 'Sample_WYGICR2-EXT.g.vcf.gz', 'Sample_S9419EG-EXT.g.vcf.gz', 
		'Sample_TNPT3VE-EXT.g.vcf.gz', 'Sample_TVNFUUT-EXT.g.vcf.gz', 'Sample_UUH8XCX-EXT.g.vcf.gz', 'Sample_I4VAK6B-EXT.g.vcf.gz', 'Sample_9ALIJX8-EXT.g.vcf.gz', 'Sample_BQPYIV9-EXT.g.vcf.gz', 
		'Sample_WJI9TSA-EXT.g.vcf.gz', 'Sample_HLB7UPN-EXT.g.vcf.gz', 'Sample_V64PVT7-EXT.g.vcf.gz', 'Sample_UA98PBR-EXT.g.vcf.gz', 'Sample_L4P373U-EXT.g.vcf.gz']

##sample ids
case_samples = ['LR18-093', 'LR18-098', 'LR18-099', 'LR18-101', 'LR18-102', 'LR18-103', 'LR18-106', 'LR18-122', 'LR18-125', 'LR18-126', 'LR18-166', 'LR18-167', 'LR18-168', 'LR18-169', 'LR18-172', 
		'LR18-173', 'LR18-177', 'LR18-179', 'LR18-180', 'LR18-183', 'LR18-186', 'LR18-187', 'LR18-192', 'LR18-195', 'LR18-198', 'LR18-200', 'LR18-244', 'LR18-245', 'LR18-247', 'LR18-250', 'LR18-252', 
		'LR18-253', 'LR18-255', 'LR18-260', 'LR18-264', 'LR18-272', 'LR18-280', 'LR18-285', 'LR18-287', 'LR18-288', 'LR18-294', 'LR18-295', 'LR18-296', 'LR18-297', 'LR18-299', 'LR18-303', 'LR18-304', 
		'LR18-305']
control_samples = ['BV1XYUL', 'M7XITGP', 'YX2RKMN', 'RJZMFUH', 'PBWT8WT', '25ZDV2G', 'IRY7WMC', 'Y19276I', '9KM25YL', '7EXPMGU', 'NIWT8YK', 'YDNWM6Z', 'ZY942NL', 'TVQ1SV4', 'HD1YDNK', 'QSZMFW9', 
		'L2F8A3F', 'HYFIANF', 'CWYIVA4', '26UWW3S', 'EMV3TSR', 'UHXLYCT', 'TNT45T4', 'AKRVUE3', 'B7FSEU3', 'ATD2WVM', 'YWNZURX', 'SVIEAEC', 'CL1J8B1', 'YYFY8PN', '576PVIT', 'F2GUM9Q', 'VKTSS9T', 
		'A8R78VT', 'PF44VGN', 'QW3G4HF', 'RZLTHDT', 'FPG61PF', 'CGAAYCL', 'MTADLK7', 'FDESA8S', 'IGFTTKA', 'IS99ILP', 'HBJNU6D', 'T8FXUUX', '4GEZMHM', '8Q1HZWE', 'XDXP7RF', '1GL82K3', 'GICP3NZ', 
		'8F2JCWB', '5Y4M12V', 'UHAFECG', 'RAMIYE7', 'R3NW7ET', 'B68SMAX', 'KDZXW3Z', 'IW8KELD', 'MFSNMIF', 'L21L11I', 'DB9LUAI', '5J83LZG', 'E1R7R8G', 'H3D2CMT', 'D2WENRB', 'MG115JD', 'LBNI329', 
		'WULHIQ2', 'C6MEWBU', '36P2R2Q', '1FSGTK2', 'RWCSRWK', '9FYGZXJ', 'DYMEXA3', 'UCS36U9', '4JEBIHC', 'F2RJIQ9', 'WBKH191', 'T4AE2BD', 'TA9N2TY', '6C4LNY1', 'KIABE21', 'C5E2FAW', '42WCDJS', 
		'WZR1B8L', 'SK2S1DK', 'HBHW153', '6RCJ3FC', '2H1TBM1', 'YQA4J7F', '3HGDVJ1', 'LYNUGIZ', 'UWHESTR', 'J6LRK6K', 'X8KJVR4', '4JRS816', 'FBIXHPN', 'ZF6BQ6A', '4X2IGIS', 'GT6S28L', 'JHHXE4D', 
		'17ISFMC', 'H544FPH', 'RN568C1', 'ADK9CUJ', 'HGZ99N9', 'AY4B7U3', 'AC28ZC2', '9B9UYEE', '63D1I1W', 'I41KM4Z', 'F9Q6AQM', '2NA5X4X', 'U66AUUI', 'EJVRXS3', '6H9H8YD', 'B9JK8SG', 'F93MY8G', 
		'6W2W2HN', '43H3VZN', 'IUSBV56', 'CBSXQUE', 'B1CNNV1', 'TK7LPT2', 'QJS4LGM', 'X5D8979', 'RZB4MWA', 'ZAD9M6J', 'DLVP2AY', 'M9TELZU', 'MSJCIHJ', 'ZD5KRL3', 'Q1HBJE5', 'I29M969', 'EAGUMBH', 
		'3H4W627', 'JBTD83B', 'CG8I5BB', 'V9PL8QV', '5UUJKGQ', '1LG2YL4', 'UI6VQ9U', 'NCUHC2T', 'USMJVSQ', 'NFYA6Z7', '4PJUB1W', 'H6E52NT', 'KYJZR3I', 'A7G5LWH', 'Z3J4T5U', '6YTVGGE', 'A7WP6TH', 
		'8ELMJWX', 'RBSQWX2', '84CZ2FB', 'TP2GMU3', 'DEN3JQF', 'SGQMPD9', 'GPNSJ78', 'TDXA3C5', '3NRGB39', 'QVCFZES', 'YMBTN7Q', 'RUB2SCE', 'FC4RNAG', 'P5TPIXR', 'H6C1G5Q', '7EE11F6', 'RI22IEY', 
		'UYHCXBN', '5U719ZJ', '6BXERG6', 'W88PFRH', 'AWA88C8', '6R2V7YS', '48MLNI8', 'KT2MIJB', 'JLFGF3R', 'A5H7HDK', 'IPHXD65', 'NQVFZYZ', 'DE7S6R4', 'CTAWVSY', 'I5JMZMH', 'E82A8RN', '6N8R8GY', 
		'QK21NWV', '5IZ6KJ8', 'DNBZKR6', 'TNRE1WQ', 'M7PLS2G', 'LKFV22N', 'D88H59I', 'SNCSRUT', 'DJZU9RU', 'LD5SKJG', 'VA8NK9D', 'I8CQY4B', 'Y9QDK7K', '31I572D', 'FP4CKPT', 'XSLJD85', 'P1V6GZC', 
		'XVLV97U', '9YPSMW2', 'BR2L3UW', 'YI1SXQH', '9FJ8GDF', 'RT8ZWX8', '9SDN6EU', '95LATDC', 'IATZHLI', '54UXA1A', 'Q9M68XD', 'Z72H2N1', '8JWD9XC', 'XZK657J', 'FFTAZNP', '6XNNIXJ', 'KAN5BKG', 
		'QFAQDYF', 'J61B53L', 'A5FFMCA', '2HNCL47', 'L23EU3T', 'QMU4MDP', 'FDC1G7H', '7SCWTX5', 'UEKTDUA', '97RUHD1', 'FKUCU8K', '1PGDULT', '3CCLW4R', 'L5CELJ2', '6ZEPM17', '1DWNQ7U', 'MYVGY1R', 
		'ILZXV5R', 'P1BP1WD', 'FXL1P7P', 'YLUWV7C', 'S27ZFFE', '95BKYWT', 'QGR1WFM', 'U7BHRDD', 'YL5KRNV', '8ZND1X6', 'DMY55RJ', 'G52YV6D', 'Q8YMXF8', 'DUZJRS8', 'N9Y8BGB', 'GU4DYN4', 'PS3SGZX', 
		'5W8R8IR', 'QIA29X5', 'DSGIDAQ', 'XKPPG92', 'TRTF1TT', '6HUYN2D', 'LIM5U2M', 'YRMYZ72', '192US5U', '75CKFGE', 'QT6VCE4', 'BG3NFVU', 'LZL5BLL', '1VGEP6P', 'SHQZFU2', 'X913FN5', 'W57LQPH', 
		'KLSBV3E', 'DCQIMAX', 'MS7VTZP', 'Y1KJWPC', '9J4ZSF4', '3QX2Z3W', 'KXFJNLY', 'I14ECME', 'ITH996U', 'JAIBM4Z', '24BUJKA', '1BDMDNC', 'BXFENAI', 'EMXSXP5', 'MGM6PJR', '2WGRGNH', 'I7ZWI4P', 
		'Z68QTNZ', 'C9R7RA8', '3LB7SL2', 'IKM4F55', 'XX5WLQC', 'PMDFHXK', 'TM8CND8', '8ZDN6GL', 'WDPAU8D', '4V8EH2X', 'K6QL9L2', '131SXLY', '3IYE83I', 'QNSEHGB', 'WRL7MQC', '9WI7YWE', 'YTMW5PY', 
		'W2MW4RQ', 'HA9M882', '39Y2HJY', 'EZKYUQK', 'DIKPBTQ', 'ZQ9RC52', 'XDVWDQ5', '362IFKJ', '6YI6KZV', 'FRXFJ6M', 'M2URU1G', '7XUA2FB', '5HFCEZF', 'VXV88TE', 'XE24A81', 'Q6T29FJ', 'VC7LQRA', 
		'W6FY8RG', 'IJTD854', '1XT7AP8', 'UGRE2UY', 'EEVP387', 'JI3PVJ9', '26HF8KY', 'UHZNTRZ', 'KFVTPL4', 'MBM5UZV', '5CN3421', 'GJ5FBN1', '5V2JB1V', 'NSYYSJG', 'SZGYBCR', 'S1CGDD3', 'V5S36SH', 
		'1B6MM47', '4NWB1IQ', '2QRFGK5', 'U1ZEAU5', 'JKH9A43', 'P3TRDFV', 'TI92SDJ', 'X53JEQQ', 'CAJXYA9', 'JHF6J33', 'GVYW1PF', 'PP9BAGQ', 'SVT45VV', '275XI24', '2MEMV3L', 'UQ5II99', 'M1BQHHY', 
		'LPS1D2Z', '9DH7GFD', 'YSGP873', 'HB775NJ', 'VQWLGU1', 'VTYLGR3', 'D6WQJRZ', 'YR6ML8R', 'KBINDLS', 'V5QBCQ7', 'G1E4NPI', '1KTIN5X', 'RX6IYWL', 'IX2BLLE', 'VBQC7A3', 'WZNB7B8', 'E16S1AD', 
		'4FL9FHL', 'U6GZPB2', 'I72NC5Z', '9TJW3WP', '3P3IX2K', 'XVYCYQN', 'VY6AUSQ', 'P2TEMX2', 'PQ4UCH2', 'E61J97H', 'AJND2TL', 'XQ4H1PN', 'KTG8SL9', 'G4EFJP8', 'UHC5IAT', 'EQZVMQ5', 'GFZWI6G', 
		'F2AMQQU', 'N35TGYS', 'V3LIIST', 'XMLHJN9', '7N2W3GW', 'BIJXYC1', 'VAVU59R', '7KQ4HXD', '1X8QUL8', 'I27UF4Y', 'XTFAL77', 'WAPYY8P', 'EGUN7Q3', 'EU5EF9P', '29WZL5S', 'BDP41B1', 'CN75WBP', 
		'L8T1VI2', 'NT78JGI', 'WPV89SL', 'JQABD4S', '4UPQV19', 'G13MZ6P', 'FC2YU86', 'QLURVWW', 'F47T27G', 'NUDGHZD', 'IXPH7LS', 'CLQVDUH', '6QECWGM', '9K2KPVL', 'QAFWGWF', 'T1GB3UJ', 'XUQC86I', 
		'QASD5F9', 'PIFWFY7', 'WCF13BC', 'GQL1Q5W', '9E96VE4', 'RBF97E8', 'XHTFFNJ', 'WIUTNQX', '7VB8PXT', 'FKN4XPQ', 'NZMTYXE', 'YP5WMNJ', '5KR5YHY', 'VZ9QX9A', 'PKW6YGE', 'JFYW2LV', 'LR9AWI7', 
		'65UA2H4', 'RZ44VC4', 'VW6X4BX', '9ZUZJFX', '65MAAXX', '2T3724N', 'PJ1NWE3', 'AGLC2WI', 'UNN9ZS2', 'FRZ8D8X', 'Z5VWEND', '32U7T1P', '2WEZLL7', '3MYQ438', 'RJK4VWG', '46VN9JH', 'RF34DV3', 
		'YGP1C7Z', 'M7C54IM', '89RTLVW', '8ALBNW8', '2DG2ZKC', 'VC1ET9E', 'MKXJ7IS', 'ZLMW5M6', 'EYZ8DAQ', 'MQDSUHV', '1AVK15V', 'KSX7F3R', 'HYKZ8M4', 'P58BSZP', 'Q2DULGG', 'KZCQZ2J', '7163MYU', 
		'Y9DWVPR', 'DFJLLRR', 'F7LDFR9', 'BP8I4D3', 'TQ8NJDX', 'UT6L8B8', 'KJSCRKI', 'GWKMH5B', 'MQQLAHI', 'XEEYQ8M', '7VNPEFN', 'JKTQYMW', 'HV88N4K', 'YFTHA6M', 'RP1PADD', 'NAGXXH1', '4RAVF3J', 
		'RL2DEEN', '2HY2HKQ', '63NPDHF', '7ZBJKWI', 'HQBBUND', 'LLEGYH6', 'ACRI4UI', 'K2IBM27', 'YFISFP4', 'AR7G9VZ', '3M9FYJR', 'LRP7PIL', 'Y8UVI69', '3SEBIJ4', 'P8WHCW4', 'BCGT9DY', 'ZB1SXNQ', 
		'BZEDTSE', 'TTAPJAZ', 'EETX97W', 'RTBP2VK', '21Z2Z2R', 'UJTE2R1', 'VYE7V8Y', '4E9EZHZ', 'WP45A7U', 'V7I2KQ8', 'RH7JGCM', 'S4QY3TZ', 'S3KR6B4', 'XIL7NNK', 'NPCEMFH', 'HHL1Q45', '1JK664Z', 
		'KDFGG1Z', '6E8EIXE', 'JZZWJ2W', 'J45FCJG', '1IFY9K4', 'Z8M8I4W', '6Z4X2FR', 'M7FT8F1', 'BKBWCAR', 'VRFK48K', 'KAF6K1A', 'JJWF3JA', 'MBKC1YK', '2DRRU2V', 'BNA88AG', 'R1XU4C5', '1R7B9KJ', 
		'APGF6TB', 'U4PXL9E', 'PRMTZVL', 'VVEVZ8A', 'GJXFK4V', 'QMHA7D3', 'DG5C37L', '6SULFYT', 'YFZP7NI', '2EM9W37', '9ZS8QDM', 'D9V1FQP', 'F84FSAR', 'LKVELZP', 'AKH49TM', 'QBD5MV4', 'CWLNFAH', 
		'AQCY7UM', 'AH74ICE', 'Z6IFP5J', 'YMYCYNW', 'UXLUVAB', '5YCJ2H4', '5X6C5Y9', 'JTE4836', '5M5PD2S', 'GIJPU85', '5J6DH34', 'J93EU5K', 'YGKIE8B', 'T4HESVJ', 'Z89GS87', 'BHLQSDB', 'XQBNN9S', 
		'566WRIT', 'XH1M28P', 'WYGICR2', 'S9419EG', 'TNPT3VE', 'TVNFUUT', 'UUH8XCX', 'I4VAK6B', '9ALIJX8', 'BQPYIV9', 'WJI9TSA', 'HLB7UPN', 'V64PVT7', 'UA98PBR', 'L4P373U']

##methods
def genotype_gvcfs_and_filter(in_gvcfs, name_prefix):
	vcf_temp0 = name_prefix + 'temp_gatk0.vcf'
	vcf_raw_snps = name_prefix + 'temp_raw_snps.vcf'
	vcf_filtered_snps = name_prefix + 'temp_filtered_snps.vcf'
	vcf_raw_indels = name_prefix + 'temp_raw_indels.vcf'
	vcf_filtered_indels = name_prefix + 'temp_filtered_indels.vcf'
	vcf_temp1 = name_prefix + 'temp_gatk1.vcf'
	vcf_temp2 = name_prefix + 'temp_gatk2.vcf.gz'
	final_vcf = name_prefix + '.gatkHC.vcf.gz'
	gvcf_cmds = []
	for in_gvcf in in_gvcfs:
		gvcf = ['-V', in_gvcf]
		gvcf_cmds.extend(gvcf)
	print gvcf_cmds
	##genotype g.vcfs
	##with nt 
	# command = ['java', '-Xmx100g', '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', fa_file, '-nt', '10'] + gvcf_cmds + ['-o', vcf_temp0]
	##straight up
	'''
	command = ['java', '-Xmx100g', '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', fa_file] + gvcf_cmds + ['-o', vcf_temp0]
	gatk_gg = subprocess.Popen(command)
	gatk_gg.wait()
	'''
	# '''
	##split data in SNPs and indels and apply manual variant filtering (add -L)
	snp_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fa_file, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_snps, '-selectType', 'SNP'])
	snp_cut.wait()
	snp_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fa_file, '-V', vcf_raw_snps, '-o', vcf_filtered_snps, '--filterExpression ', "QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", '--filterName', "snp_filter"])
	snp_vf.wait()
	indel_cut = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'SelectVariants', '-R',fa_file, '-nt', '5', '-V', vcf_temp0, '-o', vcf_raw_indels, '-selectType', 'INDEL'])
	indel_cut.wait()
	indel_vf = subprocess.Popen(['java', '-Xmx100g', '-jar', gatk, '-T', 'VariantFiltration', '-R',fa_file, '-V', vcf_raw_indels, '-o', vcf_filtered_indels, '--filterExpression ', "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0", '--filterName', "indel_filter"])
	indel_vf.wait()
	##combine filtered snps and indels
	combine_var = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineVariants', '-R',fa_file, '-nt', '5', '--variant', vcf_filtered_snps,'--variant', vcf_filtered_indels,'-o', vcf_temp1, '-genotypeMergeOptions', 'UNSORTED'])
	combine_var.wait()
	bgzip_cmd = subprocess.Popen(['bgzip', vcf_temp1])
	bgzip_cmd.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', vcf_temp1 + '.gz'])
	bcf_index.wait()
	#split multi-allelic variants calls in separate lines, and left normalize indels
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'z', '-o', vcf_temp2, vcf_temp1 + '.gz'])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fa_file, '-O', 'z', '-o', final_vcf, vcf_temp2])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', final_vcf])
	bcf_index.wait()
	# '''


def edit_gvcfs_change_name(gvcfs_to_change):
	for gvcf_to_change in gvcfs_to_change:
		vg_id = gvcf_to_change.split('_')[1].split('-')[0]
		if vg_id in case_id_dict:
			sample = case_id_dict[vg_id]
		else:
			sample = vg_id
		print(sample)
		out_gvcf = sample + '.g.vcf'
		##
		with gzip.open(gvcf_to_change, 'rb') as in_fh, open(out_gvcf, "w") as out_fh:
			for line in in_fh:
				if line[0:2] == '#C':
					# print(line)
					line = line.rstrip().split(delim)
					line_out = line[:9] + [sample]
					# print(line_out)
					out_fh.write(delim.join(line_out)+'\n')
				else:
					out_fh.write(line)

def combine_gvcfs(samples, out_vcf):
	gvcf_files = []
	for sample in samples:
		gvcf = ['-V', sample + '.g.vcf']
		gvcf_files.extend(gvcf)
	print gvcf_files
	gatk_hc = subprocess.Popen(['java', '-Xmx20g', '-jar', gatk, '-T', 'CombineGVCFs', '-R', fa_file] + gvcf_files + ['-o', out_vcf])
	gatk_hc.wait()



def split_info_field(info_list):
	indices = [i for i, s in enumerate(info_list) if 'ANNOVAR_DATE' in s]
	# print indices
	i_count = 0
	final_list = []
	for info in info_list:
		# print info
		if i_count > indices[0] and info != 'ALLELE_END':
			info2 = info.split('=')[1]
			#print info2
			final_list.append(info2)
		i_count += 1
	return final_list

def annotate_vars_from_vcf(out_prefix, in_vcf, norm_vcf, samples):
	'''
	vcf_temp1 = 'temp111.vcf'
	#split multi-allelic variants calls in separate lines
	bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-o', vcf_temp1, in_vcf])
	bcf_norm1.wait()
	bcf_norm2 = subprocess.Popen(['bcftools', 'norm', '-f', fa_file, '-O', 'z', '-o', norm_vcf, vcf_temp1])
	bcf_norm2.wait()
	bcf_index = subprocess.Popen(['bcftools', 'index', norm_vcf])
	bcf_index.wait()
	'''
	##annotate vcf file
	command = [table_annovar] + av_buildver + [norm_vcf] + av_ref_dir + av_protocol + av_operation + av_options + ['-out', out_prefix]
	annovar = subprocess.Popen(command)
	annovar.wait()
	##split into individual samples
	post_annovar_vcf = out_prefix + '.' + av_genome + '_multianno.vcf'
	con_ann = subprocess.Popen([convert_2_annovar, '-format', 'vcf4', post_annovar_vcf, '-includeinfo', '-withzyg', '-allsample', '-outfile', out_prefix])
	con_ann.wait()
	
	##format files...
	for sample in samples:
		avinput = out_prefix + '.' + sample + '.avinput'
		outfile = out_prefix + '.' + sample + '.annotated.txt'
		head = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Zygosity', 'Qual', 'Coverage', 'Filter', 'Pos', 'Ref2', 'Alt2', 'Format', 'Info', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'rmsk', 'genomicSuperDups', 'SIFT_score', 'SIFT_converted_rankscore', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_rankscore', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_converted_rankscore', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_converted_rankscore', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_score_rankscore', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_converted_rankscore', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_converted_rankscore', 'PROVEAN_pred', 'VEST3_score', 'VEST3_rankscore', 'MetaSVM_score', 'MetaSVM_rankscore', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_rankscore', 'MetaLR_pred', 'M-CAP_score', 'M-CAP_rankscore', 'M-CAP_pred', 'REVEL_score', 'REVEL_rankscore', 'MutPred_score', 'MutPred_rankscore', 'CADD_raw', 'CADD_raw_rankscore', 'CADD_phred', 'DANN_score', 'DANN_rankscore', 'fathmm-MKL_coding_score', 'fathmm-MKL_coding_rankscore', 'fathmm-MKL_coding_pred', 'Eigen_coding_or_noncoding', 'Eigen-raw', 'Eigen-PC-raw', 'GenoCanyon_score', 'GenoCanyon_score_rankscore', 'integrated_fitCons_score', 'integrated_fitCons_score_rankscore', 'integrated_confidence_value', 'GERP++_RS', 'GERP++_RS_rankscore', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP20way_mammalian', 'phyloP20way_mammalian_rankscore', 'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons20way_mammalian', 'phastCons20way_mammalian_rankscore', 'SiPhy_29way_logOdds', 'SiPhy_29way_logOdds_rankscore', 'Interpro_domain', 'GTEx_V6p_gene', 'GTEx_V6p_tissue', 'Interpro_domain', 'avsnp147', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax', 'cosmic88_coding', 'cosmic88_noncoding']
		head_out = delim.join(head + ['\n'])
		with open(avinput, "r") as av, open(outfile, "w") as final:
			final.write(head_out)
			for line in av:
				line = line.strip('\n').split(delim)
				stuff = line[0:8] + [line[14]] + [line[9]] + line[11:13]
				info = line[15].split(';')
				info_list = split_info_field(info)
				other_stuff = line[16:]
				line_out = delim.join(stuff + other_stuff + info_list +['\n'])
				final.write(line_out)

def filter_combine_format_anntxt_files(out_prefix, cases, controls, all_outfile, filtered_outfile):
	##filter anntxt files
	# '''
	for sample in cases + controls:
		anntxt = out_prefix + '.' + sample + '.annotated.txt'
		filtertxt = out_prefix + '.' + sample + '.exonic_rare.txt'
		with open(anntxt, "r") as at_fh, open(filtertxt, "w") as ft_fh:
			lc, pc = 0, 0
			for line in at_fh:
				line = line.split(delim)
				lc += 1
				if lc == 1:
					ft_fh.write(delim.join(line))
				else:
					gatk_filter = line[8]
					rmsk = line[19]
					segdup = line[20]
					rg_func = line[14]
					exonic_func = line[17]
					genome_af = line[93]
					exome_af = line[110]
					if genome_af == '.':
						genome_af = 0
					if exome_af == '.':
						exome_af = 0
					afs = [float(genome_af), float(exome_af)]
					if gatk_filter == 'PASS' and rmsk == '.' and segdup == '.' and max(afs) <= 0.01:
						if rg_func == 'exonic' or rg_func == 'splicing':
							if exonic_func != 'synonymous_SNV':
								pc +=1
								ft_fh.write(delim.join(line))

		print(sample,lc,pc)
	# '''
	##combine all vars into a dict
	var_dict = {}
	for sample in cases + controls:
		filtertxt = out_prefix + '.' + sample + '.exonic_rare.txt'
		with open(filtertxt, "r") as ft_fh:
			lc = 0
			for line in ft_fh:
				line = line.split(delim)
				lc += 1
				if lc == 1:
					header = ['case_count', 'control_count', 'case_ids', 'control_ids', 'case_quals', 'control_quals', 'case_covs', 'control_covs','case_infos', 'control_infos'] + line[:6] + line[8:13] + line[14:] 
				else:
					var = '_'.join(line[:5])
					qual = line[6]
					cov = line[7]
					info = line[13]
					var_info = line[:6] + line[8:13] + line[14:]
					if var in var_dict:
						if sample in cases:
							var_dict[var][0].append(sample)
							var_dict[var][2].append(qual)
							var_dict[var][4].append(cov)
							var_dict[var][6].append(info)
						elif sample in controls:
							var_dict[var][1].append(sample)
							var_dict[var][3].append(qual)
							var_dict[var][5].append(cov)
							var_dict[var][7].append(info)
					else:
						if sample in cases:
							var_dict[var] = [[sample], [], [qual], [], [cov], [], [info], [], var_info]
						elif sample in controls:
							var_dict[var] = [[], [sample], [], [qual], [], [cov], [], [info], var_info]
	with open(all_outfile, "w") as aout_fh, open(filtered_outfile, "w") as fout_fh:
		aout_fh.write(delim.join(header))
		fout_fh.write(delim.join(header))
		for v in var_dict:
			case_count = len(var_dict[v][0])
			control_count = len(var_dict[v][1])
			case_ids = ','.join(var_dict[v][0])
			control_ids = ','.join(var_dict[v][1])
			case_quals = ','.join(var_dict[v][2])
			control_quals = ','.join(var_dict[v][3])
			case_covs = ','.join(var_dict[v][4])
			control_covs = ','.join(var_dict[v][5])
			case_infos = ','.join(var_dict[v][6])
			control_infos = ','.join(var_dict[v][7])
			line_out = [str(case_count), str(control_count), case_ids, control_ids, case_quals, control_quals, case_covs, control_covs,case_infos, control_infos] + var_dict[v][8]
			aout_fh.write(delim.join(line_out))
			if case_count >= 2 and control_count <= 5:
				fout_fh.write(delim.join(line_out))


##run methods
project = 'sids_genomes_1119'
working_dir = '/home/atimms/ngs_data/genomes/sids_genomes_1019'
os.chdir(working_dir)
case_gvcf = project + '.all_cases.g.vcf'
control_gvcf1 = project + '.controls_1.g.vcf'
control_gvcf2 = project + '.controls_2.g.vcf'
control_gvcf3 = project + '.controls_3.g.vcf'
control_gvcf4 = project + '.controls_4.g.vcf'
control_gvcf5 = project + '.controls_5.g.vcf'
control_gvcf6 = project + '.controls_6.g.vcf'

##change sample name in gvcf and change name
##cases -- done
# edit_gvcfs_change_name(case_gvcfs_to_change[:16])
# edit_gvcfs_change_name(case_gvcfs_to_change[16:32])
# edit_gvcfs_change_name(case_gvcfs_to_change[32:])
##controls -- running
# edit_gvcfs_change_name(control_gvcfs_to_change[:50])
# edit_gvcfs_change_name(control_gvcfs_to_change[50:100])
# edit_gvcfs_change_name(control_gvcfs_to_change[100:150])
# edit_gvcfs_change_name(control_gvcfs_to_change[150:200])
# edit_gvcfs_change_name(control_gvcfs_to_change[200:250])
# edit_gvcfs_change_name(control_gvcfs_to_change[250:300])
# edit_gvcfs_change_name(control_gvcfs_to_change[300:350])
# edit_gvcfs_change_name(control_gvcfs_to_change[350:400])
# edit_gvcfs_change_name(control_gvcfs_to_change[400:450])
# edit_gvcfs_change_name(control_gvcfs_to_change[450:500])
# edit_gvcfs_change_name(control_gvcfs_to_change[500:550])
# edit_gvcfs_change_name(control_gvcfs_to_change[550:])

##combine gvcfs in batches
##cases -- going
# combine_gvcfs(case_samples, case_gvcf)
##controls -- run when ready, check all gvcfs 
# combine_gvcfs(control_samples[:100], control_gvcf1)
# combine_gvcfs(control_samples[100:200], control_gvcf2)
# combine_gvcfs(control_samples[200:300], control_gvcf3)
# combine_gvcfs(control_samples[300:400], control_gvcf4)
# combine_gvcfs(control_samples[400:500], control_gvcf5)
# combine_gvcfs(control_samples[500:], control_gvcf6)

##genotype gvcf
gvcfs_to_genotype = [case_gvcf, control_gvcf1, control_gvcf2, control_gvcf3, control_gvcf4, control_gvcf5, control_gvcf6]
# genotype_gvcfs_and_filter(gvcfs_to_genotype, project)

##annotatate variants
variant_vcf = project + '.gatkHC.vcf.gz'
normalized_vcf = project + '.gatkHC_normalized.vcf.gz'
all_samples = case_samples + control_samples
all_exonic_rare_vars = project + '.combined.exonic_rare.xls'
filtered_exonic_rare_vars = project + '.filtered.exonic_rare.xls'

##test
# project = 'test'

# variant_vcf = project + '.gatkHC.vcf.gz'
# normalized_vcf = project + '.gatkHC_normalized.vcf.gz'
# sample_dict = ['163499', '163500']
annotate_vars_from_vcf(project, variant_vcf, normalized_vcf, all_samples)

# case_samples = ['LR18-093', 'LR18-098', 'LR18-099']
# control_samples = ['BV1XYUL', 'M7XITGP', 'YX2RKMN']
filter_combine_format_anntxt_files(project, case_samples, control_samples, all_exonic_rare_vars, filtered_exonic_rare_vars)










