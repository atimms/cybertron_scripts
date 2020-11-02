#!/usr/bin/env python
import sys
import subprocess
import os


##parms
delim = '\t'



##methods

def get_dbsnp_locations(dbsnp_file, infile, outfile, reg_file):
    ##make dbSNP dict
    dbsnp_dict = {}
    with open(dbsnp_file, "r") as in_fh:
         for line in in_fh:
              line = line.rstrip().split(delim)
              rs = line[5]
              pos = line[:3]
              dbsnp_dict[rs] = pos

    print('dict made')

    ##query dict
    with open(infile, "r") as in_fh, open(outfile, "w") as out_fh, open(reg_file, "w") as reg_fh:
         for line in in_fh:
              marker = line.rstrip()
              if marker in dbsnp_dict:
                    out_fh.write(delim.join([marker] + dbsnp_dict[marker]) + '\n')
                    reg_fh.write(delim.join(dbsnp_dict[marker][:2]) + '\n')
              else:
                    out_fh.write(marker + '\n')

    print('dict queried')

def subset_vcf_bcftools(in_vcf, out_vcf, reg_file):
    bcftools_view = subprocess.Popen(['bcftools', 'view', '-O', 'z', '-o', 'temp.vcf.gz', '-R', reg_file, in_vcf])
    bcftools_view.wait()
    bcf_norm1 = subprocess.Popen(['bcftools', 'norm', '-m-both', '-O', 'v', '-o', out_vcf, 'temp.vcf.gz'])
    bcf_norm1.wait()


def format_vcf(in_vcf, outfile, marker_file):
    markers_list = []
    with open(marker_file, "r") as mark_fh:
        for line in mark_fh:
            marker = line.rstrip()
            markers_list.append(marker)
    with open(in_vcf, "r") as in_fh, open(outfile, "w") as out_fh:
        for line in in_fh:
            if line[:2] != '##':
                if line[0] == '#':
                    line = line.rstrip().split(delim)
                    header = line[:5] + line[9:]
                    out_fh.write(delim.join(header) + '\n')
                else:
                    line = line.rstrip().split(delim)
                    snp_info = line[:5]
                    rs_number = line[2]
                    ref = line[3]
                    alt = line[4]
                    vcf_genotypes = line[9:]
                    genotypes = []
                    if rs_number in markers_list:
                        for geno1 in vcf_genotypes:
                            geno = geno1.split(':')[0]
                            if geno == '0/0':
                                genotype = ref + '/' + ref
                                genotypes.append(genotype)
                            elif geno == '0/1':
                                genotype = ref + '/' + alt
                                genotypes.append(genotype)
                            elif geno == '1/1':
                                genotype = alt + '/' + alt
                                genotypes.append(genotype)
                            elif geno == './.' or '1/0':
                                genotype = ''
                                genotypes.append(genotype)
                            else:
                                print(geno, 'is a weird one', line[:5], geno1)
                        line_out = snp_info + genotypes
                        out_fh.write(delim.join(line_out) + '\n')

##run methods
working_dir = "/home/atimms/ngs_data/misc/daniela_genotype_snps_0920"
os.chdir(working_dir)

##add genomic location to rs numbers
avsnp147_file = '/home/atimms/ngs_data/references/annovar/hg19/hg19_avsnp147.txt'
markers = 'markers_wanted.txt'
markers_with_rs_locations = 'markers_with_dbsnp_info.txt'
rs_regions_file = 'dbsnp_regions.txt'
exome_vcf = 'microtia_exomes_recall_0720.gatkHC.vcf.gz'
exome_rs_vcf = 'microtia_exomes_recall_0720.gatkHC.markers.vcf'
genome_vcf = 'luquetti_grc_wgs_combined.HF.final.vcf.gz'
genome_rs_vcf = 'luquetti_grc_wgs_combined.markers.vcf'
genome_genotypes = 'genome_genotypes_markers_0920.txt'

##get marker locations from annovar dbSNP file
# get_dbsnp_locations(avsnp147_file, markers, markers_with_rs_locations, rs_regions_file)
##make vcf with just the markers (don't worry about the exomes)
# subset_vcf_bcftools(exome_vcf, exome_rs_vcf, rs_regions_file)
# subset_vcf_bcftools(genome_vcf, genome_rs_vcf, rs_regions_file)
##format vcf into genotyes
format_vcf(genome_rs_vcf, genome_genotypes, markers)


