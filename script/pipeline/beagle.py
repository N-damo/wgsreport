#!/usr/bin/env python3
# coding:utf8

import os
import sys
import pysam
import logging
import gzip
import re
import subprocess

base_dir = os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))  # dir/to/script/pipeline
main_root = os.path.dirname(base_dir)  # dir/to/wgsreport

"""
基因组填补后的数据过滤没有同一标准，一般用MAF，AR2或DR2来判断，AR2和DR2高度相关，选择其中一个即可。MAF大于0.05的基因型填补效果最好，也可根据DR2值来判断，一般在0.3-0.8属于合适阈值。注意的是，beagle5.0没有AR2统计值。
现有流程不对imputation结果做qc
"""
class Beagle(object):

    beagle = '/share/data1/local/bin/beagle.27Jan18.7e1.jar'
    bref = '/share/data3/lianlin/soft/bin/wgs/eas/eas.allgrch38.dropdup.genotype.bref'
    genetic_map = '/share/data1/PublicProject/BEAGLE'

    def __init__(self, vcf):
        self.vcf = vcf

    def parse_vcf(self):
        vcf = pysam.VariantFile(self.vcf)
        for rec in vcf.fetch():
            chrom = rec.chrom.replace('chr', '')
            pos = str(rec.pos)
            rsid = rec.id
            if rsid == None:
                rsid = '.'
            ref = rec.ref
            alt = ','.join(rec.alts)
            qual = str(rec.qual)
            filter_ = '.'
            info = '.'
            format_ = 'GT:PL'
            record = ''
            for sampleRecord in rec.samples.values():
                try:
                    GT = sampleRecord['GT']
                    PL = sampleRecord['PL']
                    # print(PL)
                    if GT == (None, None):
                        GT = './.'
                    else:
                        GT = '/'.join([str(i) for i in GT])
                    if PL == (None,):
                        PL = '.'
                    else:
                        PL = ','.join([str(i) for i in PL])
                except KeyError:
                    GT = './.'
                    PL = '.'

                record += GT+':'+PL+'\t'
            line = '\t'.join([chrom, pos, rsid, ref, alt,
                              qual, filter_, info, format_, record])
            yield line

            # break
        vcf.close()

    def vcf_head(self):
        vcf = pysam.VariantFile(self.vcf)
        head = list(vcf.header.samples)
        vcf.close()
        return '\t'.join(head)

    def vcf_out(self):
        head = self.vcf_head()
        with gzip.open('multiplex1.vcf.gz', 'wt') as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write(
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
            f.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">\n')
            f.write(
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(head))
            for line in self.parse_vcf():
                f.write(line+'\n')
        subprocess.call('tabix multiplex1.vcf.gz', shell=True)

    def sequence_split(self):
        vcfheader = os.path.join(main_root, 'file/vcfheader.vcf')
        splitjar = os.path.join(base_dir, 'pipeline/splitvcf.jar')
        allsplit = open(os.path.join(base_dir, 'pipeline/allsplit.bat')).read()
        allsplit = re.sub('vcfheader', vcfheader, allsplit)
        allsplit = re.sub('splitvcf.jar', splitjar, allsplit)
        open('allsplit.bat', 'wt').write(allsplit)
        subprocess.call('sh allsplit.bat multiplex1.vcf.gz hmmgl', shell=True)

    def pbs(self):
        hmmgl = os.listdir('./')
        hmmgl = [file for file in hmmgl if file.endswith('vcf.gz')]
        for i in range(0, 25):
            vcf = [file for file in hmmgl if file.endswith(
                'vcf.gz') and file.startswith('hmmgl'+str(i)+'.')]
            if len(vcf) != 0:
                if i < 22:
                    plink = 'chr' + 'i'
                else:
                    plink = 'chrX'
                for j in vcf:
                    section = j.split('.')[1]
                    cmd1 = "#!/bin/bash\n#$ -N imputation\n#$ -pe smp 8\n#$ -q all.q\n#$ -cwd\n"
                    cmd2 = "java -Djava.io.tmpdir=./ -Xss5m -Xmx128g -jar {beagle} gtgl={vcf} ref={bref} map={map}/plink.{plink}.GRCh38.map out=mid_{vcf} lowmem=true ;\nwait ;\n".format(
                        beagle=self.beagle, bref=self.bref, vcf=j, map=self.genetic_map, plink=plink)
                    print(cmd2)
                    open('impjob{chrom}.{section}.pbs'.format(
                        chrom=i, section=section), 'wt').write(''.join([cmd1, cmd2]))


if __name__ == "__main__":
    vcf = sys.argv[1]
    Beagle(vcf).vcf_out()
    Beagle(vcf).sequence_split()
    Beagle(vcf).pbs()
