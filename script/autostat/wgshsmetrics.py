#!/usr/bin/env python3
#coding:utf-8

import os
import sys
from .pyfastp import Basic_QC
from .pyCollectMetrics import cov_flag_merge
from .pyVcfStat import VcfStat
from .pySV import SV
from .pyCNV import CNV
from .pyCircos import Circos


base_dir=os.path.dirname(os.path.abspath(__file__))


def stat(sample_list=None):
    if sample_list == None:
        raise SystemError('CNV stat need sample list,cnv file looks like sample_cnvnator.vcf.gz and tabix file also need.\n')
    genecode='/Users/linlian/genome/gencode.v31.annotation.bed'
    makewindow='/Users/linlian/Desktop/reports/200kb.bed'
    vcf = 'annovar.hg38_multianno.vcf'
    sv='delly.vcf.gz'
    Basic_QC('1.qc')
    cov_flag_merge('2.mapping')
    VcfStat(vcf, 'SNP', '3.SNP')
    VcfStat(vcf, 'INDEL', '4.INDEL')
    SV(sv,'5.SV')
    for sample in sample_list:
        cnv=sample+'_cnvnator.vcf.gz'
        CNV(cnv,genecode,'6.CNV')
    circos_path=os.path.join(os.path.dirname(base_dir),'circos')
    Circos(vcf=vcf,window=makewindow,variant_type='SNP',directory='3.SNP',conf=circos_path)
    Circos(vcf=vcf,window=makewindow,variant_type='INDEL',directory='4.INDEL',conf=circos_path)