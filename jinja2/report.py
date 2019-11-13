#!/usr/bin/env python3
# coding:utf-8

from jinja2 import Environment, FileSystemLoader
import os
import sys
import json
import pandas as pd
from collections import defaultdict
import re


class HTML(object):
    def __init__(self, template, prefix, dict_sub=None,template_dir='templates'):
        self.template_dir=os.path.abspath(template_dir)
        self.template_file = template
        self.dict_sub = dict_sub
        self.dir_test()
        self.prefix = prefix
        self.out_dir = os.path.abspath(os.path.dirname(prefix))
        if os.path.exists(self.out_dir):
            pass
        else:
            os.mkdir(self.out_dir)

    def render__(self):
        env = Environment(loader=FileSystemLoader(self.template_dir))
        env = env.get_template(self.template_file)
        print(self.template_file)
        if self.dict_sub != None:
            content = env.render(self.dict_sub)
        else:
            content = env.render()
        static_file = open(self.prefix, 'wt',encoding='utf-8').write(content)

    def dir_test(self):
        if os.path.exists('page'):
            pass
        else:
            os.mkdir('page')


class Index(object):
    def __init__(self):
        self.render()

    def render(self, template='index.html'):
        return HTML(template, './index.html').render__()


class Introduction(object):
    def __init__(self):
        self.background()
        self.general_stat()
        self.sample_info()
        self.software_used()

    def background(self, template='项目背景/background.html'):
        return HTML(template, 'page/项目背景/background.html').render__()

    def general_stat(self, template='项目背景/general_stat.html'):
        fastp = pd.read_csv('1.qc/after_filtering.csv')
        clean_data = fastp['Clean_total_bases'].sum()
        avg_clean_data = fastp['Clean_total_bases'].mean()
        q30_percent = fastp['Clean_q30_rate'].mean()
        bwa = pd.read_csv('2.mapping/mapping_stat.csv')
        avg_mapping_depth = bwa['MEAN_COVERAGE(X)'].mean()
        avg_mapping_rate = bwa['Mapped_rate(%)'].mean()
        coverage = bwa['Coverage_at_least_1X(%)'].mean()
        snp = pd.read_csv('3.SNP/SNP.stat.csv')
        avg_snp_number = snp['Total_SNPs'].mean()
        indel = pd.read_csv('4.INDEL/INDEL.stat.csv')
        avg_indel_number = indel['Total_INDELs'].mean()
        cnv=pd.read_csv('5.CNV/cnv_stat.csv')
        total_cnv=cnv['TotalDEL'] + cnv['TotalDUP']
        avg_cnv_number=total_cnv.mean()
        dict_sub = {'clean_data': clean_data, 'avg_clean_data': avg_clean_data, 'q30_percent': q30_percent,
                    'avg_mapping_depth': avg_mapping_depth, 'avg_mapping_rate': avg_mapping_rate, 'coverage': coverage, 'avg_snp_number': avg_snp_number, 'avg_indel_number': avg_indel_number, 'avg_cnv_number': avg_cnv_number}
        return HTML(template, 'page/项目背景/general_stat.html', dict_sub=dict_sub).render__()

    def sample_info(self, template='项目背景/sample_info.html'):
        return HTML(template, 'page/项目背景/sample_info.html').render__()

    def software_used(self, template='项目背景/software_used.html'):
        dict_sub = {}
        soft = pd.read_csv('software.csv', encoding='gbk')
        for i in range(len(soft)):
            dict_sub[i] = soft.loc[i, ]
        return HTML(template, 'page/项目背景/software_used.html', dict_sub={'dict_sub':dict_sub}).render__()


class Data_qc(object):
    def __init__(self):
        self.background()
        self.qc()


    def background(self,template='测序数据质控/background.html'):
        return HTML(template,'page/测序数据质控/background.html').render__()

    def qc_tuple(self,sub,df):
        for sample in df['Samples']:
            sample_df=df[df['Samples']==sample]
            for col in df.columns[1:]:
                value=sample_df[col].values[0]
                sub[sample].append((re.sub(r'\s+','_',col),value))
        return sub

    def qc(self,template='测序数据质控/qc.html'):
        dict_sub=defaultdict(list)
        fastp = pd.read_csv('1.qc/before_filtering.csv')
        sample_list=fastp['Samples'].to_list()
        dict_sub=self.qc_tuple(dict_sub,fastp)
        fastp = pd.read_csv('1.qc/after_filtering.csv')
        dict_sub=self.qc_tuple(dict_sub,fastp)

        return HTML(template,'page/测序数据质控/qc.html',dict_sub={'sample_list':sample_list,'dict_sub':dict_sub}).render__()


class BWA_mem(object):
    def __init__(self):
        self.background()
        self.mapping_stat()

    def bwa_tuple(self,sub,df):
        for sample in df['Samples']:
            sample_df=df[df['Samples']==sample]
            for col in df.columns[1:]:
                value=sample_df[col].values[0]
                sub[sample].append((re.sub(r'\s+','_',col),value))
        return sub


    def background(self,template='基因组比对/background.html'):
        return HTML(template,'page/基因组比对/background.html').render__()


    def mapping_stat(self,template='基因组比对/mapping_stat.html'):
        dict_sub=defaultdict(list)
        bwa = pd.read_csv('2.mapping/mapping_stat.csv')
        sample_list=bwa['Samples'].to_list()
        dict_sub=self.bwa_tuple(dict_sub,bwa)
        return HTML(template,'page/基因组比对/mapping_stat.html',dict_sub={'sample_list':sample_list,'dict_sub':dict_sub}).render__()


class SNV(object):
    def __init__(self):
        self.background()
        self.variant_stat()

    def snv_tuple(self,sub,df):
        for sample in df['Samples']:
            sample_df=df[df['Samples']==sample]
            for col in df.columns[1:]:
                value=sample_df[col].values[0]
                sub[sample].append((re.sub(r'\s+','_',col),value))
        return sub

    def background(self,template='短变异位点检测/background.html'):
        return HTML(template,'page/短变异位点检测/background.html').render__()

    def variant_stat(self,template='短变异位点检测/variant_stat.html'):
        snp_dict_sub=defaultdict(list)
        snp = pd.read_csv('3.SNP/SNP.stat.csv')
        sample_list=snp['Samples'].to_list()
        snp_dict_sub=self.snv_tuple(snp_dict_sub,snp) 

        indel_dict_sub=defaultdict(list)
        indel = pd.read_csv('4.INDEL/INDEL.stat.csv')
        indel_dict_sub=self.snv_tuple(indel_dict_sub,indel)

        snp_annotation_func=defaultdict(list)
        snp_annotation=pd.read_csv('3.SNP/SNP.func.csv') 
        snp_annotation_func=self.snv_tuple(snp_annotation_func,snp_annotation)

        indel_annotation_func=defaultdict(list)
        indel_annotation=pd.read_csv('4.INDEL/INDEL.func.csv')
        indel_annotation_func=self.snv_tuple(indel_annotation_func,indel_annotation)

        return HTML(template,'page/短变异位点检测/variant_stat.html',dict_sub={'sample_list':sample_list,'snp_dict_sub':snp_dict_sub,'indel_dict_sub':indel_dict_sub,'snp_annotation_func':snp_annotation_func,'indel_annotation_func':indel_annotation_func}).render__()


class SV(object):
    def __init__(self):
        self.background()

    def background(self,template='SV变异检测/background.html'):
        return HTML(template,'page/SV变异检测/background.html').render__()

    def sv_stat(self,template='SV变异检测/sv_stat.html'):
        return HTML(template).render__()


class CNV(object):
    def __init__(self):
        self.background()
        self.cnv_stat()

    def cnv_tuple(self,sub,df):
        for sample in df['Samples']:
            sample_df=df[df['Samples']==sample]
            for col in df.columns[1:]:
                value=sample_df[col].values[0]
                sub[sample].append((re.sub(r'\s+','_',col),value))
        return sub

    def background(self,template='CNV变异检测/background.html'):
        return HTML(template,'page/CNV变异检测/background.html').render__()

    def cnv_stat(self,template='CNV变异检测/cnv_stat.html'):
        cnv=pd.read_csv('5.CNV/cnv_stat.csv')
        sample_list=cnv['Samples'].to_list()
        dict_sub=defaultdict(list)
        dict_sub=self.cnv_tuple(dict_sub,cnv)
        return HTML(template,'page/CNV变异检测/cnv_stat.html',dict_sub={'dict_sub':dict_sub,'sample_list':sample_list}).render__()


class SSR(object):
    def background(self,template='SSR标记分析/background.html'):
        return HTML(template).render__()

    def ssr_stat(self,template='SSR标记分析/ssr_stat.html'):
        return HTML(template).render__()

class Marker(object):
    def __init__(self):
        self.background()
        self.density()

    def background(self,template='标记分布可视化/background.html'):
        return HTML(template,'page/标记分布可视化/background.html').render__()

    def density(self,template='标记分布可视化/marker_stat.html'):
        dict_sub=defaultdict(list)
        snp = pd.read_csv('3.SNP/SNP.stat.csv')
        dict_sub=snp['Samples'].to_list()
        return HTML(template,'page/标记分布可视化/marker_stat.html',dict_sub={'dict_sub':dict_sub}).render__()

if __name__ == '__main__':
    Index()
    Introduction()
    Data_qc()
    BWA_mem()
    SNV()
    SV()
    CNV()
    Marker()