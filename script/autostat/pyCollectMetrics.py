#!/usr/bin/env python3
#!coding:utf-8

import pandas as pd 
from .coverage import CollectWgsMetrics
from .collect_insert_size import CollectInsertSizeMetrics
from .markdup import MarkDuplicate


"""
gatk CollectInsertSizeMetrics(insert_size_metrics.txt),MarkDuplicate(dedup.log),CollectWgsMetrics(wgs_metrics.txt)
samtools bedcov 100kb.bed B1810012.aln.bam|awk '{print $1,$2,$3,$4/($3-$2)}'|tr ' ' '\t' >bedgraph.bed
"""


def cov_flag_merge(module):
    insert=CollectInsertSizeMetrics(module).insert_size()
    dedup=MarkDuplicate(module).dedup_log()
    metrics=CollectWgsMetrics(module).parse_metric()
    df=pd.merge(dedup,metrics,on='Samples')
    df=pd.merge(df,insert,on='Samples')
    expect_columns_order=['Samples','InsertSize_mean(bp)','Insert_std(bp)','Mapped_rate(%)','PCR_duplication(%)','MEAN_COVERAGE(X)','MEDIAN_COVERAGE(X)','Coverage_at_least_1X(%)', 'Coverage_at_least_5X(%)', 'Coverage_at_least_10X(%)', 'Coverage_at_least_20X(%)','Coverage_at_least_30X(%)']
    df=df.reindex(columns=expect_columns_order)
    df.to_csv('{}/mapping_stat.csv'.format(module),index=False,header=True)


if __name__ == '__main__':
    cov_flag_merge('2.mapping')



