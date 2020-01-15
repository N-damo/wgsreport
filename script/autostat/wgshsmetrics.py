#!/usr/bin/env python3
# coding:utf-8

import os
import sys
from pyfastp import Basic_QC
from pyCollectMetrics import cov_flag_merge
from pyVcfStat import VcfStat
from pySV import SV
from pyCNV import CNV
from pyCircos import Circos
import logging
import subprocess
import pandas as pd 


base_dir = os.path.dirname(os.path.abspath(__file__))


def stat(sample_list=None):
    if sample_list == None:
        raise SystemError(
            'CNV stat need sample list,cnv file looks like sample_cnvnator.vcf.gz and tabix file also need.\n')
    genecode = '/share/data3/lianlin/1kg/gene/gencode.v31.annotation.bed'
    makewindow = '/share/data3/lianlin/soft/bin/wgs/200kb.bed'
    vcf = 'annovar.hg38_multianno.vcf'
    sv = 'delly.vcf.gz'
    # Basic_QC('qc')
    # cov_flag_merge('mapping')
    # VcfStat(vcf, 'SNP', 'SNP')
    # VcfStat(vcf, 'INDEL', 'INDEL')
    # SV(sv, 'SV')

    for sample in sample_list:
        cnv = sample+'_cnvnator.vcf.gz'
        CNV(cnv, genecode, 'CNV')
        
    if len(sample_list) == 1:
        pass
    else:
        df1=pd.read_csv(os.path.join('CNV',sample_list[0]+'cnv_stat.csv'))
        for sample in sample_list:
            df2=pd.read_csv(os.path.join('CNV',sample+'cnv_stat.csv'))
            df1=pd.concat([df1,df2])
        df=df1.drop_duplicates()
        df.to_csv('CNV/cnv_stat.csv',index=False,header=True)
    # circos_path = os.path.join(os.path.dirname(base_dir), 'circos')
    # Circos(vcf=vcf, window=makewindow, variant_type='SNP',
    #        directory='SNP', conf=circos_path)
    # Circos(vcf=vcf, window=makewindow, variant_type='INDEL',
    #        directory='INDEL', conf=circos_path)


def sample_get(file):
    sample_list = []
    with open(file, 'rt') as f:
        for k, v in enumerate(f):
            if k == 0:
                working_space = v.strip()
            else:
                sample = v.strip().split()[-1]
                if sample not in sample_list:
                    sample_list.append(sample)
                else:
                    pass
    return sample_list, os.path.abspath(working_space)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
                        format='%(levelname)s:%(asctime)s:%(message)s')
    input = sys.argv[1]
    sample_list, working_space = sample_get(input)
    os.chdir(working_space)
    stat(sample_list=sample_list)
