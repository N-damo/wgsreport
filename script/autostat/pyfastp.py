#!/usr/bin/env python3
# coding:utf-8

import json
import os
import sys
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from metrics_search import Search, nesteddict

"""
fastp log result,specially end with json and starts with sample's name.
"""


class Basic_QC(object):

    def __init__(self, module):
        self.module = module
        self.collect_file = Search('fastp').match()
        self.create_log()
        self.out_dataframe()

    def create_log(self):
        if os.path.exists(self.module):
            pass
        else:
            os.mkdir(self.module)
        for sample in self.collect_file:
            direc = os.path.join(self.module, sample)
            if os.path.exists(direc):
                pass
            else:
                os.mkdir(direc)

    def json_in(self, file):
        with open(file, 'rt') as fin:
            qc = json.load(fin)
        return qc

    def read_json(self):
        sequencing_qc = nesteddict()
        for sample in self.collect_file:
            qc = self.json_in(self.collect_file[sample])
            sequencing_qc[sample] = qc
        return sequencing_qc

    def quality_png(self, df, read, sample):
        df = pd.DataFrame.from_dict(df, orient='index').T
        df['pos'] = [i + 1 for i in df.index]
        df = df.reindex(columns=['pos', 'A', 'T', 'C', 'G', 'mean'])
        df.to_csv('{module}/{sample}/{read}.base_quality.csv'.format(
            module=self.module, sample=sample, read=read), index=False, header=True)
        plt.clf()
        plt.plot('pos', 'A', data=df, color='skyblue',
                 linewidth=2, linestyle='dashed')
        plt.plot('pos', 'T', data=df, color='red', linewidth=2)
        plt.plot('pos', 'C', data=df, color='black',
                 linewidth=2, linestyle='dashed')
        plt.plot('pos', 'G', data=df, color='green',
                 linewidth=2, linestyle='dashed')
        plt.plot('pos', 'mean', data=df, color='yellow',
                 linewidth=2, linestyle='dashed')
        plt.title('base quality')
        plt.legend(bbox_to_anchor=(1.05, 1),
                   loc="upper left", ncol=1, borderaxespad=0)
        plt.tight_layout()
        plt.savefig('{module}/{sample}/{read}.base_quality.png'.format(
            module=self.module, sample=sample, read=read), format='png', dpi=200)
        plt.close()

    def content_png(self, df, read, sample):
        df = pd.DataFrame.from_dict(df, orient='index').T
        df['pos'] = [i + 1 for i in df.index]
        df = df.reindex(columns=['pos', 'A', 'T', 'C', 'G', 'N', 'GC'])
        df.to_csv('{module}/{sample}/{read}.base_content.csv'.format(
            module=self.module, sample=sample, read=read), index=False, header=True)
        plt.clf()
        plt.plot('pos', 'A', data=df, color='skyblue',
                 linewidth=2, linestyle='dashed')
        plt.plot('pos', 'T', data=df, color='red', linewidth=2)
        plt.plot('pos', 'C', data=df, color='black',
                 linewidth=2, linestyle='dashed')
        plt.plot('pos', 'G', data=df, color='green',
                 linewidth=2, linestyle='dashed')
        plt.plot('pos', 'N', data=df, color='yellow',
                 linewidth=2, linestyle='dashed')
        plt.plot('pos', 'GC', data=df, color='blue',
                 linewidth=2, linestyle='solid')
        plt.title('base content & GC content')
        plt.legend(bbox_to_anchor=(1.05, 1),
                   loc="upper left", ncol=1, borderaxespad=0)
        plt.tight_layout()
        plt.savefig('{module}/{sample}/{read}.base_content.png'.format(
            module=self.module, sample=sample, read=read), format='png', dpi=200)
        plt.close()

    def out_dataframe(self):
        before_filtering_dataframe, after_filtering_dataframe, filtering_result_dataframe = self.raw_reads_content()
        before_filtering_pd = pd.DataFrame.from_dict(
            before_filtering_dataframe, orient='index')
        before_filtering_pd = before_filtering_pd.reset_index()
        after_filtering_pd = pd.DataFrame.from_dict(
            after_filtering_dataframe, orient='index')
        after_filtering_pd = after_filtering_pd.reset_index()
        filtering_result_pd = pd.DataFrame.from_dict(
            filtering_result_dataframe, orient='index')
        filtering_result_pd = filtering_result_pd.reset_index()
        before_filtering_pd = before_filtering_pd.rename(
            columns={'index': 'Samples'})
        columns = ['Samples', 'Raw_total_reads', 'Raw_total_bases', 'Raw_q20_bases', 'Raw_q30_bases',
                   'Raw_q20_rate', 'Raw_q30_rate', 'Raw_read1_mean_length', 'Raw_read2_mean_length', 'Raw_gc_content']
        before_filtering_pd = before_filtering_pd.reindex(columns=columns)
        before_filtering_pd.to_csv(
            '{}/before_filtering.csv'.format(self.module), index=False, header=True)
        after_filtering_pd = after_filtering_pd.rename(
            columns={'index': 'Samples'})
        columns = ['Samples', 'Clean_total_reads', 'Clean_total_bases', 'Clean_q20_bases', 'Clean_q30_bases',
                   'Clean_q20_rate', 'Clean_q30_rate', 'Clean_read1_mean_length', 'Clean_read2_mean_length', 'Clean_gc_content']
        after_filtering_pd = after_filtering_pd.reindex(columns=columns)
        after_filtering_pd.to_csv(
            '{}/after_filtering.csv'.format(self.module), index=False, header=True)
        filtering_result_pd = filtering_result_pd.rename(
            columns={'index': 'Samples'})
        columns = ['Samples', 'passed_filter_reads', 'low_quality_reads',
                   'too_many_N_reads', 'too_short_reads', 'too_long_reads']
        filtering_result_pd = filtering_result_pd.reindex(columns=columns)
        filtering_result_pd.to_csv(
            '{}/filtering_result.csv'.format(self.module), index=False, header=True)

    def raw_reads_content(self):
        before_filtering_dataframe = nesteddict()
        after_filtering_dataframe = nesteddict()
        filtering_result_dataframe = nesteddict()
        sequencing_qc = self.read_json()
        for sample in sequencing_qc:
            try:
                before_filtering = sequencing_qc[sample]['summary']['before_filtering']
                before_filtering = {
                    'Raw_'+key: before_filtering[key] for key in before_filtering}
                after_filtering = sequencing_qc[sample]['summary']['after_filtering']
                after_filtering = {
                    'Clean_'+key: after_filtering[key] for key in after_filtering}
                filtering_result = sequencing_qc[sample]['filtering_result']
                before_filtering_dataframe[sample] = before_filtering
                after_filtering_dataframe[sample] = after_filtering
                filtering_result_dataframe[sample] = filtering_result
                # quality and content png
                read1_after_filtering_quality_curves = sequencing_qc[
                    sample]['read1_after_filtering']['quality_curves']
                read1_after_filtering_content_curves = sequencing_qc[
                    sample]['read1_after_filtering']['content_curves']
                read2_after_filtering_quality_curves = sequencing_qc[
                    sample]['read2_after_filtering']['quality_curves']
                read2_after_filtering_content_curves = sequencing_qc[
                    sample]['read2_after_filtering']['content_curves']
                self.quality_png(
                    read1_after_filtering_quality_curves, 'read1', sample)
                self.content_png(
                    read1_after_filtering_content_curves, 'read1', sample)
                self.quality_png(
                    read2_after_filtering_quality_curves, 'read2', sample)
                self.content_png(
                    read2_after_filtering_content_curves, 'read2', sample)
            except KeyError:
                continue

        return before_filtering_dataframe, after_filtering_dataframe, filtering_result_dataframe


if __name__ == '__main__':
    Basic_QC('1.qc')
