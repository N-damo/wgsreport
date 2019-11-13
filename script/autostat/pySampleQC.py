#!/usr/bin/env python3
#coding:utf-8

import os
import pandas as pd

"""
mergy pyCollectMetrics.py and pyfastp.py result to qc samples.if the result of qc is not good,maybe need re-library construct and re-sequencing.
the suitable sample qc still not right now, maybe a few weeks later.
"""

class SampleQC(object):
    columns=['Samples','InsertSize_mean(bp)','Insert_std(bp)','Mapped_rate(%)','PCR_duplication(%)','MEAN_COVERAGE(X)','MEDIAN_COVERAGE(X)','Coverage_at_least_1X(%)','Coverage_at_least_5X(%)']
    def __init__(self):
        if os.path.exists('1.qc/after_filtering.csv') and os.path.exists('2.mapping/mapping_stat.csv'):
            pass
        else:
            raise OSError('no 1.qc/after_filtering.csv or 2.mapping/mapping_stat.csv')
        self.qc=pd.read_csv('1.qc/after_filtering.csv',header=0)
        self.mapping=pd.read_csv('2.mapping/mapping_stat.csv',header=0)
        self.qc_samples=self.qc['Samples']
        self.mapping_samples=self.mapping['Samples']
        if len(set(self.qc_samples) - set(self.mapping_samples)) != 0:
            raise SystemError('sample does not match perfectly.\n')
        self.df=pd.merge(self.qc,self.mapping,on='Samples')
        self.df=self.df[self.columns]
        self.df=self.df.reindex(columns=self.columns)
        self.df.to_csv('SampleQC.csv',index=False,header=True)


if __name__ == '__main__':
    SampleQC()