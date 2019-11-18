#!/usr/bin/env python3
# coding:utf-8

import os
import pandas as pd
from collections import defaultdict
import pysam
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
from pybedtools import BedTool
import json


"""
delly del dup inv bnd ins
delly call -g hunman.fa -x exclude.bed -o delly.bcf in1.bam
delly filer -p -o delly.pass.bcf delly.bcf
bcftools view -O z -o delly.vcf.gz delly.pass.bcf
"""


class SV(object):
    sv_col = ['Samples', 'DEL', 'DUP', 'INS', 'INV', 'BND']

    def __init__(self, vcf, module):
        self.vcf = vcf
        self.module = module
        if os.path.exists(self.module):
            pass
        else:
            os.mkdir(self.module)
        self.sv_format()

    def source_dict(self, sample_list):
        dict_init = defaultdict(dict)
        for sample in sample_list:
            dict_init[sample]['DEL'] = 0
            dict_init[sample]['DUP'] = 0
            dict_init[sample]['BND'] = 0
            dict_init[sample]['INS'] = 0
            dict_init[sample]['INV'] = 0
            dict_init[sample]['DEL_length'] = []
            dict_init[sample]['DUP_length'] = []
            dict_init[sample]['BND_length'] = []
            dict_init[sample]['INS_length'] = []
            dict_init[sample]['INV_length'] = []
        return dict_init

    def parse_vcf(self):
        vcf = pysam.VariantFile(self.vcf)
        sample_list = list(vcf.header.samples)
        sample_dict = self.source_dict(sample_list)
        for rec in vcf.fetch():
            self.chrom = rec.chrom
            start = rec.start
            end = rec.stop
            sv_length = end-start
            filter_ = list(rec.filter)[0]
            if filter_ == 'PASS':
                svtype = rec.info['SVTYPE']
                for sampleRecord in rec.samples.values():
                    GT = sampleRecord['GT']
                    if GT == (0, 0) or GT == (None, None):
                        continue
                    else:
                        if svtype == 'DEL':
                            sample_dict[sampleRecord.name]['DEL'] += 1
                            sample_dict[sampleRecord.name]['DEL_length'].append(
                                sv_length)
                        if svtype == 'DUP':
                            sample_dict[sampleRecord.name]['DUP'] += 1
                            sample_dict[sampleRecord.name]['DUP_length'].append(
                                sv_length)
                        if svtype == 'INS':
                            sample_dict[sampleRecord.name]['INS'] += 1
                            sample_dict[sampleRecord.name]['INS_length'].append(
                                sv_length)
                        if svtype == 'INV':
                            sample_dict[sampleRecord.name]['INV'] += 1
                            sample_dict[sampleRecord.name]['INV_length'].append(
                                sv_length)
                        if svtype == 'BND':
                            sample_dict[sampleRecord.name]['BND'] += 1
                            sample_dict[sampleRecord.name]['BND_length'].append(
                                sv_length)
        vcf.close()
        return sample_dict

    def sv_format(self):
        sample_dict = self.parse_vcf()
        df = pd.DataFrame.from_dict(sample_dict, orient='index')
        df = df.reset_index()
        df = df.rename(columns={'index': 'Samples'})
        with open('test.json', 'wt') as f:
            json.dump(sample_dict, f)
        df = df.reindex(columns=self.sv_col)
        df.to_csv('{module}/sv_stat.csv'.format(module=self.module),
                  header=True, index=False)
        self.sv_plot(sample_dict)

    def sv_plot(self, sample_dict):
        d = ['DEL_length', 'DUP_length', 'INS_length', 'INV_length', 'BND_length']
        for sample in sample_dict:
            if os.path.exists(os.path.join(self.module, sample)):
                pass
            else:
                os.mkdir(os.path.join(self.module, sample))
            sv = sample_dict[sample]
            plt.clf()
            n = 0
            fig = plt.figure(figsize=(7, 7), tight_layout=True)
            palate = plt.get_cmap('Set1')
            for i in d:
                n += 1
                plt.subplot(3, 2, n)
                a = self.group(sv[i])
                df = pd.DataFrame.from_dict(a, orient='index')
                df = df.reset_index()
                df.columns = ['group', 'count']
                plt.bar(df['group'], df['count'], label=i, color=palate(n))
                plt.axis(ymin=0)
                plt.title('{} distribution'.format(i))
                plt.savefig('{module}/{sample}/sv_length_distribution.png'.format(
                    module=self.module, sample=sample), dpi=200, format='png')
            plt.close()

    def group(self, d):
        a = {}
        a['0-1k'] = 0
        a['1-5k'] = 0
        a['5-10k'] = 0
        a['over 10k'] = 0
        if len(d) == 0:
            return a
        else:
            for i in d:
                if 0 < i <= 1e3:
                    a['0-1k'] += 1
                if 1e3 < i <= 5e3:
                    a['1-5k'] += 1
                if 5e3 < i <= 10e3:
                    a['5-10k'] += 1
                if i > 10e3:
                    a['over 10k'] += 1
            return a


if __name__ == '__main__':
    vcf = sys.argv[1]
    SV(vcf, '5.SV')
