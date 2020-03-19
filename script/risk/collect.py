import pysam
import json
from collections import defaultdict
import subprocess
import numpy as np
import os
import mysql_db
import pandas as pd 
import shutil


def read_json(f):
    with open(f, 'rt') as f:
        for i in f:
            l = json.loads(i.strip())
    return l


def out_json(save_name, json_file):
    with open(save_name, 'wt') as f:
        json.dump(json_file, f, ensure_ascii=False)


def trim_vcf(vcf, locus):
    rs = []
    with open('gwas_id.txt', 'wt') as f:
        for value in locus.values():
            rs += value.keys()
        rs = list(set(rs))
        f.write('\n'.join(rs))
    print('ok')
    # cmd = 'plink2 --vcf {} --recode vcf --extract {} --out sub && bgzip -f sub.vcf && tabix -f sub.vcf.gz'.format(
    #     vcf, 'gwas_id.txt')
    cmd='vcftools --gzvcf {} --recode --snps {} -c |bgzip > sub.vcf.gz && tabix -f sub.vcf.gz'.format(vcf,'gwas_id.txt')
    print(cmd)
    subprocess.call(cmd, shell=True)


def nesteddict():
    return defaultdict(nesteddict)


class Beta(object):
    def __init__(self, vcf, all_summary_data):
        self.ori_vcf = vcf
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
        vcf_in = pysam.VariantFile(self.vcf)
        sample_beta = self.source_dict()
        for rec in vcf_in.fetch():
            rsid = rec.id
            sample_alleles = {s.name: list(set(s.alleles))
                              for s in rec.samples.values()}
            for trait in self.all_summary_data:
                summary_data = self.all_summary_data[trait]
                if not rsid in summary_data:
                    continue
                else:
                    #gwas[rsid] = [riskallele, beta, pvalue, freq]
                    riskallele = summary_data[rsid][0]
                    beta = float(summary_data[rsid][1])
                    for sample in sample_alleles:
                        allele_dup = sample_alleles[sample]
                        if len(allele_dup) == 2:  # double == het genotype
                            x = 1
                        # single == ref or alt,if equal to riskallele,x=2
                        elif allele_dup[0] == riskallele:
                            x = 2
                        else:
                            x = 0
                        glm_beta = x*beta
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
                # print(sample_alleles[sample])
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

    def __init__(self, sample_beta, sample_gender, gwas_beta_avg_std,module='risk'):
        self.sample_beta = sample_beta
        self.sample_gender = sample_gender
        self.gwas_beta_avg_std = gwas_beta_avg_std
        self.module=module
        if os.path.exists('risk'):
            pass
        else:
            os.mkdir('risk')
        self.risk = nesteddict()


    def dataset2chinese(self):
        sample_disease_beta = {}
        sample_beta_avg_std = {}
        for basic_data_name in self.sample_beta:
            table = mysql_db.Available.select().where(
                mysql_db.Available.basic_data_name == basic_data_name)
            table = list(table.dicts())
            if len(table) == 0:
                continue
            else:
                for record in table:
                    chinese_name = record['chinese']
                    sample_disease_beta[chinese_name] = self.sample_beta[basic_data_name]
                    sample_beta_avg_std[chinese_name] = self.gwas_beta_avg_std[basic_data_name]
        return sample_disease_beta, sample_beta_avg_std

    def get_dict_list(self, chinese_name, *table):
        context = []
        for tb in table:
            t = tb.select().where(tb.chinese_name == chinese_name)
            dict_list = list(t.dicts())
            context += dict_list
            if len(dict_list) != 0:
                major_field = tb._meta.table_name
                return major_field, context
            else:
                major_field = '未知'
        return major_field, context

    def disease_risk(self):
        sample_disease_beta, sample_beta_avg_std = self.dataset2chinese()
        with open('risk.csv','wt') as fout:
            fout.write(','.join(['Samples','gender','chinese_name','major_field','classification','suggest'])+'\n')
            for chinese_name in sample_disease_beta:
                beta_avg = sample_beta_avg_std[chinese_name][0]
                beta_std = sample_beta_avg_std[chinese_name][1]
                major_field, table = self.get_dict_list(chinese_name, mysql_db.CommonDisease, mysql_db.InfectiousDiseases, mysql_db.CancerDisease,
                                                        mysql_db.BodyChracteristic, mysql_db.Metabolism, mysql_db.DrugReaction, mysql_db.IndividualChracteristic)
                if len(table) == 0:
                    continue
                else:
                    record = table[0]
                    sex = record.get('sex', '2')
                    classification = record.get('classification', '其他')
                    for sample in sample_disease_beta[chinese_name]:
                        beta = float(sample_disease_beta[chinese_name][sample])
                        gender = self.sample_gender[sample]
                        if sex == '2':
                            # 超过1.5倍标准差为高，-1.5-1.5倍就是正常，反之就是低
                            suggest = self.z_score(beta, beta_avg, beta_std)
                        elif sex == '0':
                            if gender == '1':
                                continue
                            else:
                                suggest = self.z_score(beta, beta_avg, beta_std)
                        elif sex == '1':
                            if gender == '0':
                                continue
                            else:
                                suggest = self.z_score(beta, beta_avg, beta_std)
                        #cnmorbidity = ''
                        #individual_morbidity = ''
                        #risk_value = ''
                        #self.risk[sample][chinese_name]['major_field'] = major_field
                        #self.risk[sample][chinese_name]['classification'] = classification
                        #self.risk[sample][chinese_name]['cnmorbidity'] = cnmorbidity
                        #self.risk[sample][chinese_name]['individual_morbidity'] = individual_morbidity
                        #self.risk[sample][chinese_name]['risk_value'] = risk_value
                        #self.risk[sample][chinese_name]['suggest'] = suggest
                        out=[sample, gender, chinese_name, major_field, classification, suggest]
                        fout.write(','.join(map(str,out))+'\n')


    def z_score(self, beta, beta_avg, beta_std):
        z = (beta-beta_avg)/beta_std
        if z >= 1.5:
            level = '偏高'
        elif -1.5 < z < 1.5:
            level = '正常'
        else:
            level = '偏低'
        return level


def avg_std(gwas):
    cal_avg_std = {}
    for t in gwas:
        values = gwas[t]
        beta = list(values.values())
        avg = np.mean(beta)
        std = np.std(beta, ddof=1)
        cal_avg_std[t] = [avg, std]
    return cal_avg_std


if __name__ == '__main__':
    vcf='gmb26.vcf.gz'
    #vcf = '/share/data1/PublicProject/BEAGLE/eas.allgrch38.genotype.vcf.gz'
    #gwas_locus = read_json('gwas_summary_stat.json')
    # trim_vcf(vcf,gwas_locus)
    # sample_beta=Beta(vcf,gwas_locus).parse_vcf()
    # sample_gender=Beta(vcf,gwas_locus).infer_gender()
    # out_json('sample_beta.json',sample_beta)
    # out_json('sample_gender.json',sample_gender)
    # sample_beta = read_json('sample_beta.json')
    # cal_avg_std = avg_std(sample_beta)
    # out_json('beta_avg_std.json', cal_avg_std)
    gwas_locus = read_json('gwas_summary_stat.json')
    trim_vcf(vcf,gwas_locus)
    sample_beta=Beta(vcf,gwas_locus).parse_vcf()
    sample_gender=Beta(vcf,gwas_locus).infer_gender()
    gwas_beta_avg_std=read_json('beta_avg_std.json')
    disease=Diseaese(sample_beta,sample_gender,gwas_beta_avg_std)
    disease.disease_risk()
    df=pd.read_csv('risk.csv')
    df.to_csv(os.path.join('risk','risk.csv'),index=False,encoding='utf_8_sig')
    for i in df.groupby('Samples'):
        sample=i[0]
        dataframe=i[1]
        if os.path.exists(os.path.join('risk',sample)):
            pass
        else:
            os.mkdir(os.path.join('risk',sample))
            dataframe.to_csv(os.path.join('risk',sample,'risk.csv'),index=False,encoding='utf_8_sig')
    