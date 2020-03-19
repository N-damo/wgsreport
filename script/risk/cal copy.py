#!/usr/bin/env python3
#coding:utf-8


import json
import os
import sys
import re
import numpy as np
import pysam
import subprocess
from collections import defaultdict
import math
import mysql_db
from scipy.stats import norm
import pandas as pd
import shutil


def nesteddict():
    return defaultdict(nesteddict)



def read_json(file):
    with open(file,'rt') as f:
        for i in f:
            d=json.loads(i.strip())
    return d


def out_json(save_name,json_file):
    with open(save_name,'wt') as f:
        json.dump(json_file,f,ensure_ascii=False)


def trim_vcf(vcf,locus):
    rs=[]
    with open('gwas_id.txt','wt') as f:
        for value in locus.values():
            for rsnumber in value.keys():
                if rsnumber not in rs:
                    rs.append(rsnumber) 
                    f.write(rsnumber+'\n')
    cmd='vcftools --gzvcf {} --recode --snps {} -c |bgzip > sub.vcf.gz && tabix -f sub.vcf.gz'.format(vcf,'gwas_id.txt')
    subprocess.call(cmd,shell=True)



class Beta(object):
    def __init__(self, vcf, all_summary_data):
        self.ori_vcf=vcf
        self.vcf = 'sub.vcf.gz'
        self.all_summary_data = all_summary_data


    def source_dict(self):
        sample_beta = defaultdict(dict)
        vcf = pysam.VariantFile(self.vcf)
        for trait in self.all_summary_data:
            for sample in vcf.header.samples:
                sample_beta[trait][sample] = 0
        vcf.close()
        return sample_beta        


    def parse_vcf(self):
        vcf_in=pysam.VariantFile(self.vcf)
        sample_beta=self.source_dict()
        for rec in vcf_in.fetch():
            rsid = rec.id
            sample_alleles = {s.name: list(set(s.alleles)) for s in rec.samples.values()}
            for trait in self.all_summary_data:
                summary_data = self.all_summary_data[trait]
                if not rsid in summary_data:
                    continue
                else:
                    #gwas[rsid] = [riskallele, beta, pvalue, freq]
                    riskallele = summary_data[rsid][0]
                    beta = float(summary_data[rsid][1])
                    freq = float(summary_data[rsid][3])
                    for sample in sample_alleles:
                        allele_dup = sample_alleles[sample]
                        if len(allele_dup) == 2:  # double == het genotype
                            x = 1
                        elif allele_dup[0] == riskallele:  # single == ref or alt,if equal to riskallele,x=2
                            x = 2
                        else:
                            x = 0
                        # beta*(x-2f)/sqrt(2f(1-f))
                        glm_beta = beta*(x-(2*freq)) / \
                            math.sqrt(2*freq*(1.0-freq))
                        sample_beta[trait][sample] += glm_beta
        vcf_in.close()
        return sample_beta


    def infer_gender(self):
        sample_x_hom = defaultdict(int)
        sample_gender = {}
        vcf = pysam.VariantFile(self.ori_vcf)
        total = 0
        for rec in vcf.fetch('X', 2781478, 155701383):
            total += 1
            sample_alleles = {s.name: len(set(s.alleles))
                              for s in rec.samples.values()}
            for sample in sample_alleles:
                #print(sample_alleles[sample])
                if sample_alleles[sample] == 1:
                    sample_x_hom[sample] += 1
        for sample in vcf.header.samples:
            x_hom = sample_x_hom[sample]/total
            if x_hom > 0.98:
                sample_gender[sample] = 0  # male
            else:
                sample_gender[sample] = 1  # female
        vcf.close()
        return sample_gender


class Diseaese(object):
    
    def __init__(self,sample_beta,sample_gender,gwas_beta_avg_std):
        self.sample_beta=sample_beta
        self.sample_gender=sample_gender
        self.gwas_beta_avg_std=gwas_beta_avg_std
        if os.path.exists('risk.csv'):
            os.remove('risk.csv')
        self.save_file=open('risk.csv','ta')
        self.risk=nesteddict()


    def dataset2chinese(self):
        sample_disease_beta={}
        sample_beta_avg_std={}
        for basic_data_name in self.sample_beta:
            table=mysql_db.Available.select().where(mysql_db.Available.basic_data_name == basic_data_name)
            table=list(table.dicts())
            if len(table) == 0:
                continue
            else:
                for record in table:
                    chinese_name=record['chinese']
                    sample_disease_beta[chinese_name]=self.sample_beta[basic_data_name]
                    sample_beta_avg_std[chinese_name]=self.gwas_beta_avg_std[basic_data_name]
        return sample_disease_beta,sample_beta_avg_std


    def get_dict_list(self,chinese_name,*table):
        context=[]
        for tb in table:
            t=tb.select().where(tb.chinese_name == chinese_name)
            dict_list=list(t.dicts())
            context += dict_list
            if len(dict_list) != 0:
                major_field=tb._meta.table_name
                return major_field,context
            else:
                major_field='未知'
        return major_field,context



    def body_character(self):
        sample_disease_beta,sample_beta_avg_std=self.dataset2chinese()
        for chinese_name in sample_disease_beta:
            beta_avg=sample_beta_avg_std[chinese_name][0]
            beta_std=sample_beta_avg_std[chinese_name][1]
            major_field,table=self.get_dict_list(chinese_name,mysql_db.BodyChracteristic,mysql_db.Metabolism,mysql_db.DrugReaction,mysql_db.IndividualChracteristic)
            if len(table) == 0:
                continue
            else:
                record=table[0]
                sex=record.get('sex','2')
                classification=record.get('classification','其他')
                for sample in sample_disease_beta[chinese_name]:
                    beta=float(sample_disease_beta[chinese_name][sample])
                    gender=self.sample_gender[sample]
                    if sex == '2':
                        #超过2倍标准差为高，-2-2倍就是正常，反之就是低
                        suggest=self.z_score(beta,beta_avg,beta_std)
                    elif sex == '0':
                        if gender == '1':
                            continue
                        else:
                            suggest=self.z_score(beta,beta_avg,beta_std)
                    elif sex == '1':
                        if gender == '0':
                            continue
                        else:
                            suggest=self.z_score(beta,beta_avg,beta_std)
                    cnmorbidity=''
                    individual_morbidity='' 
                    risk_value=''
                    self.risk[sample][chinese_name]['major_field']=major_field
                    self.risk[sample][chinese_name]['classification']=classification
                    self.risk[sample][chinese_name]['cnmorbidity']=cnmorbidity
                    self.risk[sample][chinese_name]['individual_morbidity']=individual_morbidity
                    self.risk[sample][chinese_name]['risk_value']=risk_value
                    self.risk[sample][chinese_name]['suggest']=suggest
                    print(sample,gender, chinese_name,major_field,classification,cnmorbidity,individual_morbidity, risk_value, suggest, file=self.save_file,sep=',')


    def z_score(self,beta,beta_avg,beta_std):
        z=(beta-beta_avg)/beta_std
        if z >=2:
            level='偏高'
        elif -2 < z < 2:
            level='正常'
        else:
            level='偏低'
        return level




    # male=morbidity1,morbidity3
    # female=morbidity2,morbidity4
    # sex=0,only male
    # sex=1,only female
    #sex=2,male and female4
    def disease(self):
        sample_disease_beta,sample_beta_avg_std=self.dataset2chinese()
        for chinese_name in sample_disease_beta:
            beta_avg=sample_beta_avg_std[chinese_name][0]
            beta_std=sample_beta_avg_std[chinese_name][1]
            major_field,table=self.get_dict_list(chinese_name,mysql_db.CommonDisease,mysql_db.InfectiousDiseases,mysql_db.CancerDisease)
            if len(table) == 0:
                continue
            else:
                record=table[0]
                sex=record.get('sex','2')
                classification=record.get('classification','其他')
                gene_factor=record.get('gene_factor','')
                morbidity1=record.get('morbidity1','')
                morbidity2=record.get('morbidity2','')
                morbidity3=record.get('morbidity3','')
                morbidity4=record.get('morbidity4','')
                if gene_factor != '':
                    gene_factor=float(gene_factor)
                else:
                    gene_factor=1
                for sample in sample_disease_beta[chinese_name]:
                    beta=float(sample_disease_beta[chinese_name][sample])
                    gender=self.sample_gender[sample]
                    cnmorbidity=self.common_morbidity(sex,gender,morbidity1,morbidity2,morbidity3,morbidity4)
                    if cnmorbidity:
                        t = norm.ppf(1-cnmorbidity)  # 人群发病率对应beta值
                        new_beta = (beta - beta_avg)/beta_std
                        x = new_beta
                        _ = 1-norm.cdf(t-x)
                        individual_morbidity = _*gene_factor  # average h = 20%
                        times = individual_morbidity/cnmorbidity
                        if times > 1:
                            original_risk_value = 1+math.log(times, math.e**2)
                        else:
                            original_risk_value = times
                        # if times >1,risk_value may still lower than 1
                        # (1+math.log(math.e, math.e**2))=1.5
                        risk_value = original_risk_value / \
                            (1+math.log(math.e, math.e**2))
                        suggest = self.report_value(
                            risk_value, cnmorbidity, individual_morbidity)
                        risk_value=round(risk_value,2)
                        self.risk[sample][chinese_name]['major_field']=major_field
                        self.risk[sample][chinese_name]['classification']=classification
                        self.risk[sample][chinese_name]['cnmorbidity']=cnmorbidity
                        self.risk[sample][chinese_name]['individual_morbidity']=individual_morbidity
                        self.risk[sample][chinese_name]['risk_value']=risk_value
                        self.risk[sample][chinese_name]['suggest']=suggest
                        print(sample,gender, chinese_name,major_field,classification,cnmorbidity,individual_morbidity, risk_value, suggest, file=self.save_file,sep=',')


    def common_morbidity(self,sex,gender,morbidity1,morbidity2,morbidity3,morbidity4):
        if sex == '2':
            if gender == 0:
                if morbidity1 != '':
                    cnmorbidity=float(morbidity1)
                else:
                    cnmorbidity=float(morbidity3)
            else:
                if morbidity2 != '':
                    cnmorbidity=float(morbidity2)
                else:
                    cnmorbidity=float(morbidity4)
        elif sex == '0':
            if gender == 0:
                if morbidity1 != '':
                    cnmorbidity=float(morbidity1)
                else:
                    cnmorbidity=float(morbidity3)
            else:
                cnmorbidity=False
            
        elif sex == '1':
            if gender == 1:
                if morbidity2 != '':
                    cnmorbidity=float(morbidity2)
                else:
                    cnmorbidity=float(morbidity4)
            else:
                cnmorbidity=False
        return cnmorbidity

    def report_value(self, risk_value, cnmorbidity, individual_morbidity):
        cal_cnmorbidity_level = self.morbidity_level(cnmorbidity)
        cal_individual_morbidity_level = self.morbidity_level(individual_morbidity)
        if risk_value == 1:
            if cal_individual_morbidity_level > cal_cnmorbidity_level:
                cal_individual_morbidity_level = cal_cnmorbidity_level
        individual_risk_value_level = self.individual_risk_value(risk_value)
        if cal_individual_morbidity_level == 3 and individual_risk_value_level == 3:
            return '重点关注'
        elif cal_individual_morbidity_level == 3 or individual_risk_value_level == 3:
            return '适度关注'
        elif individual_risk_value_level == 2 or cal_individual_morbidity_level == 2:
            return '略微关注'
        else:
            return '正常对待'

    def morbidity_level(self, morbidity):
        if morbidity >= 1e-1:
            level = 3
        elif 1e-3 <= morbidity < 1e-1:
            level = 2
        else:
            level = 1
        return level

    def individual_risk_value(self, risk_value):
        # (1+math.log(50*math.e, math.e**2))/(1+math.log(math.e, math.e**2)) ≈ 2.3
        if risk_value < 1:
            level = 1
        elif 1 <= risk_value < (1+math.log(50*math.e, math.e**2))/(1+math.log(math.e, math.e**2)):
            level = 2
        else:
            level = 3
        return level

    

if __name__ == '__main__':
    if os.path.exists('risk'):
        pass
    else:
        os.mkdir('risk')
    trait_locus='/Users/linlian/Documents/GitHub/wgsreport/script/risk/data/trait_locus.json'
    beta_avg_std='/Users/linlian/Documents/GitHub/wgsreport/script/risk/data/beta_avg_std.json'
    test_vcf='gmb26.vcf.gz'
    gwas_locus=read_json(trait_locus)
    gwas_beta_avg_std=read_json(beta_avg_std)
    # trim_vcf(test_vcf,gwas_locus)
    # sample_beta=Beta(test_vcf,gwas_locus).parse_vcf()
    # sample_gender=Beta(test_vcf,gwas_locus).infer_gender()
    # out_json('sample_beta.json',sample_beta)
    # out_json('sample_gender.json',sample_gender)
    sample_beta=read_json('/Users/linlian/Documents/GitHub/wgsreport/script/risk/sample_beta.json')
    sample_gender=read_json('/Users/linlian/Documents/GitHub/wgsreport/script/risk/sample_gender.json')
    disease=Diseaese(sample_beta,sample_gender,gwas_beta_avg_std)
    disease.disease()
    disease.body_character()
    out_json('risk.json',disease.risk)
    df=pd.read_csv('risk.csv')
    df.columns=['Samples','gender','chinese_name','major_field','classification','cnmorbidity','individual','risk_value','suggest']
    for i in df.groupby('Samples'):
        sample=i[0]
        dataframe=i[1]
        if os.path.exists(os.path.join('risk',sample)):
            pass
        else:
            os.mkdir(os.path.join('risk',sample))
            dataframe.to_csv(os.path.join('risk',sample,'risk.csv'),index=False,encoding='utf_8_sig')
    shutil.copy('risk.json','risk')
    df.to_csv(os.path.join('risk','risk.csv'),index=False,encoding='utf_8_sig')
