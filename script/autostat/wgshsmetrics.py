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


def stat():
    genecode = '/share/data3/lianlin/1kg/gene/gencode.v31.annotation.bed'
    vcf = 'annovar.hg38_multianno.vcf'
    sv = 'delly.vcf.gz'
    cnv='cnvnator.vcf.gz'
    # Basic_QC('qc')
    # cov_flag_merge('mapping')
    # VcfStat(vcf, 'SNP', 'SNP')
    # VcfStat(vcf, 'INDEL', 'INDEL')
    # SV(sv, 'SV')
    CNV(cnv, genecode, 'CNV')
 


if __name__ == '__main__':
    stat()
