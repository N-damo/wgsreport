#!/usr/bin/env python3
# coding:utf-8

import os
import pandas as pd
import sys
from linecache import getline
import subprocess
from collections import defaultdict


def nesteddict():
    return defaultdict(nesteddict)


def lb_sample_match():
    sample_dict = {}
    files = [f for f in os.listdir('./') if f.endswith('aln.bam')]
    for f in files:
        library, sampleID, *_ = f.split('.')
        if library not in sample_dict:
            sample_dict[library] = sampleID
        else:
            raise SystemError(
                'can not accept multi libraryID during test time\n')
    return sample_dict


def get_line(file, pattern):
    with open(file, 'rt') as f:
        for index, value in enumerate(f):
            if value.startswith(pattern):
                return index + 1
            else:
                continue


def dedup_log():
    files = [f for f in os.listdir('./') if f.endswith('dedup.log')]
    sample = nesteddict()
    for log in files:
        line = get_line(log, 'LIBRARY') + 1
        rec = getline(log, line).strip().split()
        library, UNPAIRED_READS_EXAMINED, READ_PAIRS_EXAMINED, PERCENT_DUPLICATION = rec[
            0], rec[1], rec[2], rec[-2]
        sample[library]['unmerged_reads'] = UNPAIRED_READS_EXAMINED
        sample[library]['merged_reads'] = READ_PAIRS_EXAMINED
        sample[library]['pcr_duplication'] = PERCENT_DUPLICATION
    df = pd.DataFrame.from_dict(sample).T
    df = df.reset_index()
    df = df.rename(columns={'index': 'sample'})
    df.to_csv('de.csv', header=True, index=False)
    return df


def bwa_mem():
    bwamem = [f for f in os.listdir('./') if f.endswith('bwamem')]
    sample_dict = {}
    for sample in bwamem:
        sampleID = sample.split('.')[0]
        line = get_line(
            sample, '[M::mem_pestat] analyzing insert size distribution for orientation FR...') + 3
        line = getline(sample, line)
        len_std = line.split(':')[-1].strip().strip('(').strip(')')
        tlen, std = len_std.split(',')
        sample_dict[sampleID] = [tlen, std]
    df = pd.DataFrame.from_dict(sample_dict).T
    df = df.reset_index()
    df.columns = ['sample', 'mean_tlen_unmerged', 'sd_tlen_unmerged']
    df.to_csv('bwa.csv', header=True, index=False)
    return df


def flagstatx():
    df = pd.read_table('combined_bqsr.flagstatx', header=24)
    df = df[df['2'] != 'fail']
    df['# 1'] = df['# 1'].map(lambda x: x.split('_')[-1])
    df = df[['# 1', '19']]
    df.columns = ['sample', 'percentage_proper_mapped']
    df.to_csv('flag.csv', header=True, index=False)
    return df


def cov_stat():
    df = pd.read_table('combined_bqsr.covstat', header=0)
    df = df[df['QC'] != 'fail']
    chrs = []
    for i in range(1, 23):
        chrs.append('chr'+str(i))
    chrs.append('chrX')
    chrs.append('chrY')
    df = df[df['# chr'].isin(chrs)]
    df = df.reset_index()
    df = df.drop('index', axis=1)
    df = df.drop('QC', axis=1)
    df_mean = df.loc[0:22].mean()
    df_mean = df_mean.reset_index()
    df_mean.columns = ['sample', 'coverage_avg']
    x = df.loc[22, ].reset_index()
    x = x.loc[1:len(x) - 2]
    y = df.loc[23, ].reset_index()
    y = y.loc[1:len(y) - 2]
    xToy = pd.merge(x, y)
    xToy['coverage_X2Y'] = xToy[22]/xToy[23]
    xToy.columns = ['sample', '22', '23', 'coverage_X2Y']
    xToy = xToy[['sample', 'coverage_X2Y']]
    df = pd.merge(df_mean, xToy)
    df['sample'] = df['sample'].map(lambda x: x.split('_')[-1])
    df.to_csv('cov.csv', header=True, index=False)
    return df


def predict_gender(gender):
    if gender > 1:
        return '女'
    else:
        return '男'


if __name__ == '__main__':
    dirpath = sys.argv[1]
    os.chdir(dirpath)
    batch_no = os.path.basename(os.path.abspath('./'))
    columns = ['样本号(订单号)', 'batch_no', 'deposit_cram', 'deposit_vcf', 'unmerged_reads', 'merged_reads', 'mean_tlen_unmerged', 'sd_tlen_unmerged', 'pcr_duplication', 'percentage_proper_mapped',
               'coverage_X2Y', 'sex_from_data', 'contamination_MT', 'coverage_avg', 'coverage_80', 'imputation_error']
    sample_dict = lb_sample_match()
    dedup = dedup_log()
    insert_size = bwa_mem()
    proper_mapped = flagstatx()
    coverage = cov_stat()
    l = [insert_size, proper_mapped, coverage]
    for i in l:
        dedup = pd.merge(dedup, i)
    dedup['样本号(订单号)'] = dedup['sample'].map(lambda x: sample_dict[x])
    dedup['sex_from_data'] = dedup['coverage_X2Y'].map(predict_gender)
    dedup['coverage_80'] = ''
    dedup['contamination_MT'] = ''
    dedup['imputation_error'] = ''
    dedup['batch_no'] = batch_no
    dedup['deposit_cram'] = '/share/data4/deposit/CRAM/{}_combined_bqsr.cram'.format(batch_no)
    dedup['deposit_vcf'] = '/share/data4/deposit/VCF/{}.vcf.gz'.format(batch_no)
    dedup = dedup[columns]
    dedup = dedup.reindex(columns=columns)
    dedup.to_excel('sample_info_merge.xlsx', header=True,
                   index=False, encoding='utf_8_sig')
