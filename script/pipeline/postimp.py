#!/usr/bin/env python3
# coding:utf-8

import sys
import os
import subprocess


def merge(prefix):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    mergejar = os.path.join(base_dir, 'mergevcf.jar')
    for line in sys.stdin:
        line = line.strip()
        chrom = line.split()[0].split('.')[0].replace('mid_hmmgl', '')
        if chrom not in ['22', '23', '24']:
            chrom = 1+int(chrom)
        else:
            chrom = 'X'
        cmd = 'java -jar {mergejar} {chrom} {line}|bgzip >chr{chrom}.vcf.gz ;'.format(
            mergejar=mergejar, chrom=chrom, line=line)
        print(cmd)
        subprocess.call(cmd, shell=True)
        subprocess.call('tabix chr{chrom}.vcf.gz ;'.format(
            chrom=chrom), shell=True)
    subprocess.call('ls chr*.vcf.gz >vcf.list ;', shell=True)
    subprocess.call(
        'vcf-concat -f vcf.list |vcf-sort -t ./ -c|bgzip >{prefix}_all.vcf.gz ;'.format(prefix=prefix), shell=True)
    subprocess.call('tabix {prefix}_all.vcf.gz ;'.format(
        prefix=prefix), shell=True)


if __name__ == '__main__':
    prefix = sys.argv[1]
    merge(prefix)
