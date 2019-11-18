#!/usr/bin/env python3
# coding:utf-8

import pysam
import os
import sys
import pandas as pd
from collections import defaultdict
import json
from pybedtools import BedTool
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import re


class Circos():
    def __init__(self, vcf=None, window=None, variant_type=None, directory=None, conf=None):
        self.vcf = vcf
        self.window = window
        self.variant_type = variant_type
        self.dir = directory
        self.conf = conf
        self.format_stat()
        #subprocess.call('circos -conf /Users/linlian/Documents/GitHub/wgs_stat_html/scripts/circos/circos.conf -outputdir {} -outputfile {}.png'.format(self.dir,'marker'),shell=True)

    def chr_list_init(self):
        chr_list = []
        for i in range(1, 23):
            chr_list.append('chr'+str(i))
        chr_list.append('chrX')
        chr_list.append('chrY')
        if os.path.exists(self.dir):
            pass
        else:
            os.mkdir(self.dir)
        return chr_list

    def source_dict(self, sample_list):
        dict_init = defaultdict(dict)
        chr_list = self.chr_list_init()
        sample_list.append('all')
        for sample in sample_list:
            for chrom in chr_list:
                dict_init[sample][chrom] = []
        return dict_init

    def snp_or_indel(self, ref: str, alt: tuple) -> str:
        alt = sorted(alt, key=lambda x: len(x), reverse=False)[-1]
        if len(ref) > 1 or len(alt) > 1:
            if len(ref) == len(alt):
                if len(set(list(ref)) - set(list(alt))) >= 1:
                    print('unlegal variant record in {} {}.Maybe this is a snp or mnp,not indel record.if has multi-allelic snp,the stat will be not accuracy.\n'.format(self.chrom, self.position))
            else:
                return 'INDEL'
        else:
            return 'SNP'

    def parse_vcf(self):
        vcf = pysam.VariantFile(self.vcf)
        sample_list = list(vcf.header.samples)
        sample_dict = self.source_dict(sample_list)
        for rec in vcf.fetch():
            self.chrom = rec.chrom
            self.position = rec.pos
            ref = rec.ref
            alt = rec.alts
            if self.snp_or_indel(ref, alt) == self.variant_type:
                sample_dict['all'][self.chrom].append(self.position)
                for sampleRecorde in rec.samples.values():
                    GT = sampleRecorde['GT']
                    if GT == (0, 0) or GT == (None, None):
                        continue
                    else:
                        sample_dict[sampleRecorde.name][self.chrom].append(
                            self.position)
        return sample_dict

    def format_stat(self):
        sample_dict = self.parse_vcf()
        self.save_dict(sample_dict)

    def save_dict(self, all_dict):
        chr_list = self.chr_list_init()
        with open('test.json', 'wt') as f:
            json.dump(all_dict, f)
        for sample in all_dict:
            if sample == 'all':
                pass
            else:
                if os.path.exists(os.path.join(self.dir, sample)):
                    pass
                else:
                    os.mkdir(os.path.join(self.dir, sample))
            partial = all_dict[sample]
            df1 = pd.DataFrame({'chrom': ['chr1' for i in range(
                len(partial['chr1']))], 'end': partial['chr1']})
            for chrom in chr_list[1:]:
                length = len(partial[chrom])
                df2 = pd.DataFrame(
                    {'chrom': [chrom for i in range(length)], 'end': partial[chrom]})
                df1 = pd.concat([df1, df2], axis=0)
            df1['start'] = df1['end']-1
            df1 = df1.reindex(columns=['chrom', 'start', 'end'])
            df1['start'] = df1['start'].astype(int)
            df1['end'] = df1['end'].astype(int)
            if sample != 'all':
                df1.to_csv('{module}/{sample}/{variant}_pos.bed'.format(module=self.dir,
                                                                        sample=sample, variant=self.variant_type), index=False, header=False, sep='\t')
            else:
                df1.to_csv('{module}/all_{variant}_pos.bed'.format(module=self.dir,
                                                                   variant=self.variant_type), index=False, header=False, sep='\t')
            self.bed_intersect(sample)

    def bed_intersect(self, sample):
        if sample == 'all':
            bed_root = '{module}/all_{variant}_pos.bed'.format(
                module=self.dir, variant=self.variant_type)
        else:
            bed_root = '{module}/{sample}/{variant}_pos.bed'.format(
                module=self.dir, sample=sample, variant=self.variant_type)
        bed = BedTool(bed_root)
        window = BedTool(self.window)
        # calculate the number of snp/indel in 200kb window
        intersect = window.intersect(bed, c=True)
        #intersect=intersect.filter(lambda x:x[3] !='0')
        intersect.moveto(bed_root)
        outputdir = os.path.dirname(bed_root)
        self.circos()
        self.plot_conf(sample)
        cmd = 'circos -conf circos2.conf -outputdir {outputdir} -outputfile {sample}_{variant}.circos.png'.format(
            outputdir=outputdir, sample=sample, variant=self.variant_type)
        subprocess.call(cmd, shell=True)
        self.heatmap_plot(bed_root)

    def circos(self):
        circos_conf = os.path.join(self.conf, 'circos.conf')
        line = self.pattern_search(circos_conf, 'karyotype')
        karyotype = os.path.join(self.conf, 'karyotype.human.hg38.txt')
        route = 'karyotype = ' + karyotype
        self.edit_conf(circos_conf, 'circos.conf', line, route)
        if os.path.exists('ticks.conf'):
            pass
        else:
            subprocess.call(
                'cp {}/ticks.conf ticks.conf'.format(self.conf), shell=True)
            subprocess.call(
                'cp {}/ideogram.conf ideogram.conf'.format(self.conf), shell=True)

    def edit_conf(self, conf, output, line, route):
        with open(output, 'wt') as fout:
            with open(conf, 'rt') as fin:
                for k, v in enumerate(fin):
                    rec = v.strip()
                    if k == line:
                        fout.write(route+'\n')
                    else:
                        fout.write(rec+'\n')
            fin.close()

    def plot_conf(self, sample):
        conf = os.path.join(
            self.conf, 'plots{}_histogram.conf'.format(self.variant_type))
        line = self.pattern_search(conf, 'file')
        if sample == 'all':
            conf_edit = 'file = {module}/all_{m}_pos.bed'.format(
                m=self.variant_type, module=self.dir)
        else:
            conf_edit = 'file = {module}/{sample}/{m}_pos.bed'.format(
                m=self.variant_type, module=self.dir, sample=sample)
        self.edit_conf(conf, 'plots{}_histogram.conf'.format(
            self.variant_type), line, conf_edit)

        snp_line = self.pattern_search(
            'circos.conf', '<<include plotsSNP_histogram.conf>>')
        indel_line = self.pattern_search(
            'circos.conf', '<<include plotsINDEL_histogram.conf>>')
        if self.variant_type == 'SNP':
            self.edit_conf('circos.conf', 'circos2.conf', indel_line, '')
        elif self.variant_type == 'INDEL':
            self.edit_conf('circos.conf', 'circos2.conf', snp_line, '')

    def pattern_search(self, conf, pattern):
        with open(conf, 'rt') as f:
            for k, v in enumerate(f):
                p = re.search(pattern, v)
                if p == None:
                    continue
                else:
                    return k

    def heatmap_plot(self, bed):
        plt.clf()
        df = pd.read_table(bed, sep='\t', header=None)
        df.columns = ['chrom', 'start', 'end', 'number']
        df2 = pd.DataFrame()
        for chr in df['chrom'].drop_duplicates():
            data = df[df['chrom'] == chr]
            data = data.reset_index()
            data = data.drop('index', axis=1)
            df2[chr] = data[data['chrom'] == chr]['number']
        df2 = df2.fillna(0)
        df2_norm_col = (df2-df2.mean())/df2.std()
        fig = plt.figure(figsize=(20, 15), tight_layout=True)
        sns.heatmap(df2_norm_col, robust=True, cmap='viridis')
        plt.title('Genome Distribution Heatmap of Molecular Markers')
        dirname = os.path.dirname(bed)
        plt.savefig(os.path.join(dirname, '{}_marker_density.png'.format(
            self.variant_type)), format='png', bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    vcf = sys.argv[1]
    window = sys.argv[2]
    Circos(vcf, window, 'SNP', '3.SNP')
    Circos(vcf, window, 'INDEL', '4.INDEL')
