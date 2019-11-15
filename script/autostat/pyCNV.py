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


"""
install AnnotSV.tcl,bedtools,delly
annotation: AnnotSV.tcl
~/soft/AnnotSV_2.2/sv/bin/AnnotSV/AnnotSV.tcl -SVinputFile -outputFile ./sv.csv -genomeBuild GRCh38 
structure variant calling: delly 
"""

def nesteddict():
    return defaultdict(nesteddict)

class CNV(object):

    #stat_col =['Samples','TotalDEL','TotalDUP','Total_DEL_length(bp)','Total_DUP_length(bp)','Overlapped gene']
    stat_col =['Samples','TotalDEL','TotalDUP','Overlapped gene','over 1KB DEL','over 5KB DEL','over 10KB DEL','over 100KB DEL','over 1000KB DEL','over 1KB DUP','over 5KB DUP','over 10KB DUP','over 100KB DUP','over 1000KB DUP']
    def __init__(self,cnv,genecode,module):
        self.vcf = os.path.abspath(cnv)
        self.genecode = genecode
        self.module = module
        self.sample_stat = defaultdict(dict)   
        if os.path.exists(self.module):
            pass
        else:
            os.mkdir(self.module)

        self.cnv_format()
        self.sample_stat_out()



    def sample_stat_out(self):
        cnv_stat = pd.DataFrame.from_dict(self.sample_stat,orient='index')
        cnv_stat = cnv_stat.reset_index()
        cnv_stat = cnv_stat.rename(columns={'index':'Samples'})
        cnv_stat = cnv_stat.reindex(columns=self.stat_col)
        cnv_stat.to_csv('{module}/cnv_stat.csv'.format(module=self.module),index=False,header=True)

        
    def source_dict(self,sample_list):
        samples_sv = nesteddict()
        for sample in sample_list:
            samples_sv[sample]['Chrom'] = []
            samples_sv[sample]['Start']= []
            samples_sv[sample]['End']= []
            samples_sv[sample]['Length(bp)']= []
            samples_sv[sample]['Type'] = []
            samples_sv[sample]['CN'] = []
        return samples_sv


    def parse_vcf(self):
        vcf = pysam.VariantFile(self.vcf)
        sample_list=list(vcf.header.samples)
        samples_sv =self.source_dict(sample_list)
        for rec in vcf.fetch():
            sv_type = rec.info['SVTYPE']
            if sv_type == 'DEL' or sv_type == 'DUP':
                chrom =rec.chrom
                start = rec.pos
                end = rec.stop
                cnv_length = end - start
                filter__=list(rec.filter)[0]
                if filter__ == 'PASS':
                    if cnv_length > 1000:
                        for sampleRecord in rec.samples.values():
                            GT=sampleRecord['GT']
                            CN=sampleRecord['CN']
                            if CN == 2:
                                continue
                            else:
                                if GT == (0,0) or GT == (None,None):
                                    continue
                                else:
                                    samples_sv[sampleRecord.name]['Chrom'].append(chrom)
                                    samples_sv[sampleRecord.name]['Start'].append(start)
                                    samples_sv[sampleRecord.name]['End'].append(end)
                                    samples_sv[sampleRecord.name]['Length(bp)'].append(cnv_length)
                                    samples_sv[sampleRecord.name]['Type'].append(sv_type)
                                    samples_sv[sampleRecord.name]['CN'].append(CN)
                    
        vcf.close()
        return samples_sv


    def cnv_format(self):
        samples_sv = self.parse_vcf()
        for sample in samples_sv:
            if os.path.exists(os.path.join(self.module,sample)):
                pass
            else:
                os.mkdir(os.path.join(self.module,sample))
            df = pd.DataFrame()
            df['Samples'] = [sample for i in range(len(samples_sv[sample]['Chrom']))] 
            df['Chrom'] = samples_sv[sample]['Chrom']
            df['Start'] = samples_sv[sample]['Start']
            df['End'] =samples_sv[sample]['End']
            df['Length(bp)'] = samples_sv[sample]['Length(bp)']
            df['Type'] = samples_sv[sample]['Type']
            df['CN'] = samples_sv[sample]['CN']
            df.to_csv('{module}/{sample}/cnv_igv.seg'.format(module=self.module,sample=sample),index=False,header=True,sep='\t')
            df=df.drop('Samples',axis=1)
            #df=df.drop('CN',axis=1)
            df.to_csv('{module}/{sample}/cnv.bed'.format(module=self.module,sample=sample),index=False,header=False,sep='\t')
            self.cnv_plot(df,sample)
            cnv = BedTool('{module}/{sample}/cnv.bed'.format(module=self.module,sample=sample))
            gene = BedTool(self.genecode)
            intersect = cnv.intersect(gene,wb=True,f=0.5)
            intersect_gene_count = self.gene_count(intersect)
            #e.g chr1	1406909	1406998	18358	DUP	chr1	1406909	1406998	MRPL20
            intersect.moveto('{module}/{sample}/cnv.annotatedGenecodeV31.bed'.format(module=self.module,sample=sample))
            self.cnv_stat(df,sample,intersect_gene_count)
   
        #print(self.sample_stat)
       


    def gene_count(self,intersect):
        gene_count = []
        for feature in intersect:
            gene=feature[-1]
            if gene not in gene_count:
                gene_count.append(gene)
        return len(gene_count)


    def stat_length_distribution(self,length_dis):
        a={}
        a['>1k']=0
        a['>5k']=0
        a['>10k']=0
        a['>100k']=0
        a['>1000k']=0
        if len(length_dis) == 0:
            return a['>1k'],a['>5k'],a['>10k'],a['>100k'],a['>1000k']
        else:
            for i in length_dis:
                kb = i/1000
                if kb > 1:
                    a['>1k'] += 1
                if kb > 5:
                    a['>5k'] += 1
                if kb > 10:
                    a['>10k'] += 1
                if kb > 100:
                    a['>100k'] += 1
                if kb > 1000:
                    a['>1000k'] += 1
            return a['>1k'],a['>5k'],a['>10k'],a['>100k'],a['>1000k']



    def cnv_stat(self,df,sample,gene_count):
        del_length = df[df['Type'] == 'DEL']['Length(bp)']
        delk1,delk5,delk10,delk100,delk1000 = self.stat_length_distribution(del_length)
        dup_length = df[df['Type'] == 'DUP']['Length(bp)']
        dupk1,dupk5,dupk10,dupk100,dupk1000 = self.stat_length_distribution(dup_length)
        self.sample_stat[sample]['TotalDEL'] = len(df[df['Type'] == 'DEL'])
        self.sample_stat[sample]['over 1KB DEL'] = delk1
        self.sample_stat[sample]['over 5KB DEL'] = delk5
        self.sample_stat[sample]['over 10KB DEL'] = delk10
        self.sample_stat[sample]['over 100KB DEL'] = delk100
        self.sample_stat[sample]['over 1000KB DEL'] = delk1000
        #self.sample_stat[sample]['Total_DEL_length(bp)'] = sum(df[df['Type'] == 'DEL']['Length(bp)'])
        self.sample_stat[sample]['TotalDUP'] = len(df[df['Type'] == 'DUP'])
        self.sample_stat[sample]['over 1KB DUP'] = dupk1
        self.sample_stat[sample]['over 5KB DUP'] = dupk5
        self.sample_stat[sample]['over 10KB DUP'] = dupk10
        self.sample_stat[sample]['over 100KB DUP'] = dupk100
        self.sample_stat[sample]['over 1000KB DUP'] = dupk1000
        #self.sample_stat[sample]['Total_DUP_length(bp)'] = sum(df[df['Type'] == 'DUP']['Length(bp)'])
        self.sample_stat[sample]['Overlapped gene'] = gene_count
        

    def cnv_plot(self,df,sample):
        plt.clf()
        sns.violinplot(y="Length(bp)", x="Type", hue='Type', data=df, palette="Pastel1")
        plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", ncol=1, borderaxespad=0)
        plt.title('CNV distribution')
        plt.tight_layout()
        plt.savefig('{module}/{sample}/CNV_distribution.png'.format(module=self.module,sample=sample),format='png')
        plt.close()



if __name__ == '__main__':
    cnv = sys.argv[1]
    #genecode = sys.argv[2]
    genecode='/Users/linlian/genome/gencode.v31.annotation.bed'
    CNV(cnv,genecode,'6.CNV')






