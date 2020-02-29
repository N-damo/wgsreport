#!/usr/bin/env python3
# coding:utf-8

import os
import pandas as pd
from collections import defaultdict
import pysam
import sys
import matplotlib.pyplot as plt
import shutil
import json
from gprofiler import GProfiler
import numpy as np 

"""
in-house script to stat snp and indel distribution.Typically, the vcf should annotated by annovar,it should annotated with 1000genome database and refGene database.
if the vcf contatin multi-allelic snp,the statistic result will be not very accuracy. Maybe you can compare it with RTG.jar VCFSTATS tools to check your result if you needed.
snp and indel should coexist in the same vcf,if not you could mergy them with gatk,bcftools or other available tools.
annovar command:perl /share/data1/src/annovar/table_annovar.pl han.vcf.gz /share/data2/leon/annovar/humandb -buildver hg38 -out test -otherinfo -remove -protocol refGene,cytoBand,exac03,gnomad_exome,gnomad_genome,ALL.sites.2015_08,clinvar_20170905 -operation g,r,f,f,f,f,f --vcfinput -nastring . -polish -thread 5
"""


def nesteddict():
    return defaultdict(nesteddict)


class VcfStat(object):

    #snp_col = ['Samples','Total_variants','Passing filter','Unpassing filter','Total_SNPs','Missing genotype','Same as reference','Homozygous','Heterozygous','Ti','Tv','Ti/Tv','Fraction of 1000genome(%)','Fraction of dbsnp(%)','Novel','Intergenic','Intronic','Exonic','Splicing','ncRNA','Downstream','Upstream','UTR3','UTR5','Unknown']
    snp_col = ['Samples', 'Total_SNPs', 'Missing_genotype', 'Homozygous', 'Heterozygous', 'Ti', 'Tv', 'TivsTv',
               'Intergenic', 'Intronic', 'Exonic', 'Splicing', 'ncRNA', 'Downstream', 'Upstream', 'UTR3', 'UTR5', 'Unknown']
    indel_col = ['Samples', 'Total_INDELs', 'Missing_genotype', 'TotalInsertion', 'TotalDeletion', 'Homozygous', 'Heterozygous',
                 'Intergenic', 'Intronic', 'Exonic', 'Splicing', 'ncRNA', 'Downstream', 'Upstream', 'UTR3', 'UTR5', 'Unknown']
    func_col = ['Samples', 'synonymous', 'nonsynonymous', 'stopgain', 'stoploss', 'frameshift_insertion', 'frameshift_deletion',
                'frameshift_block_substitution', 'nonframeshift_insertion', 'nonframeshift_deletion', 'nonframeshift_block_substitution', 'ExonicFunc_Unknown']
    #total_col = list(set(snp_col + indel_col + func_col))

    def __init__(self, vcf, variant_type, module):
        self.vcf = os.path.abspath(vcf)
        self.variant_type = variant_type
        self.module = module
        if os.path.exists(self.module):
            pass
        else:
            os.mkdir(self.module)
        self.format_stat()

    def snp_or_indel(self, ref: str, alt: tuple) -> str:
        alt = sorted(alt, key=lambda x: len(x), reverse=False)[-1]
        if len(ref) > 1 or len(alt) > 1:
            if len(ref) == len(alt):
                if len(set(list(ref)) - set(list(alt))) >= 1:
                    print('unlegal variant record in {} {}.Maybe this is a snp or mnp,not indel record.if has multi-allelic snp,the stat will be not accuracy.\n'.format(self.chr, self.pos))
            else:
                return 'INDEL'
        else:
            return 'SNP'

    def source_dict(self, sample_list):
        source_dict = defaultdict(dict)
        for sample in sample_list:
            source_dict[sample]['Passing filter'] = 0
            source_dict[sample]['Unpassing filter'] = 0
            source_dict[sample]['Homozygous'] = 0
            source_dict[sample]['Heterozygous'] = 0
            source_dict[sample]['Total_{}s'.format(self.variant_type)] = 0
            source_dict[sample]['Deletion'] = []
            source_dict[sample]['Insertion'] = []
            if self.variant_type == 'SNP':
                source_dict[sample]['Ti'] = 0
                source_dict[sample]['Tv'] = 0
            source_dict[sample]['Missing_genotype'] = 0
            source_dict[sample]['Same as reference'] = 0
            source_dict[sample]['Intergenic'] = 0
            source_dict[sample]['Intronic'] = 0
            source_dict[sample]['Exonic'] = 0
            source_dict[sample]['Downstream'] = 0
            source_dict[sample]['Upstream'] = 0
            source_dict[sample]['UTR3'] = 0
            source_dict[sample]['UTR5'] = 0
            source_dict[sample]['ncRNA'] = 0
            source_dict[sample]['Splicing'] = 0
            source_dict[sample]['Unknown'] = 0
            #source_dict[sample]['1000genome'] = 0
            #source_dict[sample]['Novel'] = 0
            #source_dict[sample]['dbsnp'] = 0
            source_dict[sample]['synonymous'] = 0
            source_dict[sample]['nonsynonymous'] = 0
            source_dict[sample]['stopgain'] = 0
            source_dict[sample]['stoploss'] = 0
            source_dict[sample]['frameshift_insertion'] = 0
            source_dict[sample]['frameshift_deletion'] = 0
            source_dict[sample]['frameshift_block_substitution'] = 0
            source_dict[sample]['nonframeshift_insertion'] = 0
            source_dict[sample]['nonframeshift_deletion'] = 0
            source_dict[sample]['nonframeshift_block_substitution'] = 0
            source_dict[sample]['ExonicFunc_Unknown'] = 0
            source_dict[sample]['depth'] = []
            source_dict[sample]['quality'] = []
            source_dict[sample]['gene'] = []
        return source_dict

    # def snp_mnp_deletion_insertion_indel(self,ref:str,alt:str):
    #     ref_len = len(ref)
    #     alt_len = len(alt)
    #     if ref_len == alt_len:
    #         ref = list(ref)
    #         alt = list(alt)
    #         if len(set(ref) - set(alt)) == 1:
    #             return 'SNP'
    #         else:
    #             return 'MNP'
    #     elif ref_len > alt_len :
    #         head = ref[0:alt_len]
    #         if head == alt:
    #             return 'Deletion'
    #         else:
    #             return 'INDEL'
    #     else:
    #         head = alt[0:ref_len]
    #         if head == ref:
    #             return 'Insertion'
    #         else:
    #             return 'INDEL'

    def indel_distribution(self, alleles):
        ref = alleles[0]
        alt = alleles[1]
        ref_len = len(ref)
        alt_len = len(alt)
        if ref_len > alt_len:
            return 'Deletion', ref_len - alt_len
        elif ref_len < alt_len:
            return 'Insertion', alt_len - ref_len
        else:
            return 'Insertion', alt_len

    def ti_tv(self, alleles, ref):
        ti = {'A': 'G', 'G': 'A', 'T': 'C', 'C': 'T'}
        #tv={'A':'T','A':'C','G':'T','G':'C','T':'G','T':'A','C':' A','C':'G'}
        allele1 = alleles[0]
        allele2 = alleles[1]
        try:
            if allele1 == allele2:
                if ti[ref] == allele1:
                    return 'ti', 2
                else:
                    return 'tv', 2
            else:
                if ti[ref] == allele1 or ti[ref] == allele2:
                    return 'ti', 1
                else:
                    return 'tv', 1
        except KeyError:
            raise SystemError(
                'unlegal variant: {} {}\n'.format(self.chr, self.pos))

    def parse_vcf(self):
        vcf = pysam.VariantFile(self.vcf)
        self.sample_list = list(vcf.header.samples)
        samples_stat = self.source_dict(self.sample_list)
        indel_length = nesteddict()
        self.snp_or_indel_record = 0
        for rec in vcf.fetch():
            self.chr = rec.chrom
            self.pos = rec.pos
            ref = rec.ref
            alt = rec.alts
            self.rsid = rec.id
            try:
                filter_ = list(rec.filter)[0]
            except IndexError:
                #print('filter stat unknow,set it to be PASS stat.\n')
                filter_ = 'PASS'
            try:
                #kilogenome_ALL = rec.info['ALL.sites.2015_08']
                Func_refGene = rec.info['Func.refGene'][0]
                ExonicFunc_refGene = rec.info['ExonicFunc.refGene'][0]
                gene=rec.info['Gene.refGene'][0]
                
            except KeyError:
                #kilogenome_ALL = 'error'
                Func_refGene = 'error'
                ExonicFunc_refGene = 'error'
                gene='error'
                # print(
                #     'can not find ALL.sites.2015_08 or Func.refGene or ExonicFunc.refGene.The relative record will be 0.\n')
            #print(gene)
            if self.snp_or_indel(ref, alt) == self.variant_type:
                self.snp_or_indel_record += 1
                if filter_ == 'PASS':
                    for sampleRecord in rec.samples.values():
                        # EQUAL TO TOTAL_SNP(INDEL)S + SAME AS REFERENCE + MISSING GENOTYPE
                        samples_stat[sampleRecord.name]['Passing filter'] += 1
                        alleles = sampleRecord.alleles
                        try:
                            GT = sampleRecord['GT']
                            # print(GT,sampleRecord.name,self.pos)
                            DP = sampleRecord['DP']
                            GQ = sampleRecord['GQ']
                        except KeyError:
                            continue
                        if GT.__eq__((0, 0)):
                            samples_stat[sampleRecord.name]['Same as reference'] += 1
                            continue
                        elif alleles == (None, None):
                            samples_stat[sampleRecord.name]['Missing_genotype'] += 1
                            continue
                        else:
                            if gene not in samples_stat[sampleRecord.name]['gene'] and gene.find('x3b') == -1:
                                samples_stat[sampleRecord.name]['gene'].append(gene)
                            samples_stat[sampleRecord.name]['depth'].append(DP)
                            samples_stat[sampleRecord.name]['quality'].append(GQ)
                            samples_stat[sampleRecord.name]['Total_{}s'.format(
                                self.variant_type)] += 1  # EUQAL TO TOTALINSERT + TOTALDELETE
                            if self.variant_type == 'INDEL':
                                indel_stat, indel_length = self.indel_distribution(
                                    alleles)
                                samples_stat[sampleRecord.name][indel_stat].append(
                                    indel_length)
                            #######################################
                            if len(set(alleles)) == 1:
                                samples_stat[sampleRecord.name]['Homozygous'] += 1
                            else:
                                samples_stat[sampleRecord.name]['Heterozygous'] += 1
                            #################################
                            if self.variant_type == 'SNP':
                                tiOrtv, number = self.ti_tv(alleles, ref)
                                if tiOrtv == 'ti':
                                    samples_stat[sampleRecord.name]['Ti'] += number
                                else:
                                    samples_stat[sampleRecord.name]['Tv'] += number
                            ################################
                            if Func_refGene == 'intergenic':
                                samples_stat[sampleRecord.name]['Intergenic'] += 1
                            elif Func_refGene == 'intronic':
                                samples_stat[sampleRecord.name]['Intronic'] += 1
                            elif Func_refGene == 'exonic':
                                samples_stat[sampleRecord.name]['Exonic'] += 1
                            elif Func_refGene == 'downstream':
                                samples_stat[sampleRecord.name]['Downstream'] += 1
                            elif Func_refGene == 'upstream':
                                samples_stat[sampleRecord.name]['Upstream'] += 1
                            elif Func_refGene == 'UTR3':
                                samples_stat[sampleRecord.name]['UTR3'] += 1
                            elif Func_refGene == 'UTR5':
                                samples_stat[sampleRecord.name]['UTR5'] += 1
                            elif Func_refGene.startswith('ncRNA'):
                                samples_stat[sampleRecord.name]['ncRNA'] += 1
                            elif Func_refGene == 'splicing':
                                samples_stat[sampleRecord.name]['Splicing'] += 1
                            else:
                                samples_stat[sampleRecord.name]['Unknown'] += 1
                            #########################
                            # if kilogenome_ALL != '.':
                            #     samples_stat[sampleRecord.name]['1000genome'] += 1
                            #######################
                            # if self.rsid == '.':
                            #     samples_stat[sampleRecord.name]['Novel'] += 1
                            # else:
                            #     samples_stat[sampleRecord.name]['dbsnp'] += 1
                            ##############################
                            if ExonicFunc_refGene == 'synonymous_SNV':
                                samples_stat[sampleRecord.name]['synonymous'] += 1
                            elif ExonicFunc_refGene == 'nonsynonymous_SNV':
                                samples_stat[sampleRecord.name]['nonsynonymous'] += 1
                            elif ExonicFunc_refGene == 'stopgain':
                                samples_stat[sampleRecord.name]['stopgain'] += 1
                            elif ExonicFunc_refGene == 'stoploss':
                                samples_stat[sampleRecord.name]['stoploss'] += 1
                            elif ExonicFunc_refGene == 'frameshift_insertion':
                                samples_stat[sampleRecord.name]['frameshift_insertion'] += 1
                            elif ExonicFunc_refGene == 'frameshift_deletion':
                                samples_stat[sampleRecord.name]['frameshift_deletion'] += 1
                            elif ExonicFunc_refGene == 'frameshift_block_substitution':
                                samples_stat[sampleRecord.name]['frameshift_block_substitution'] += 1
                            elif ExonicFunc_refGene == 'nonframeshift_insertion':
                                samples_stat[sampleRecord.name]['nonframeshift_insertion'] += 1
                            elif ExonicFunc_refGene == 'nonframeshift_deletion':
                                samples_stat[sampleRecord.name]['nonframeshift_deletion'] += 1
                            elif ExonicFunc_refGene == 'nonframeshift_block_substitution':
                                samples_stat[sampleRecord.name]['nonframeshift_block_substitution'] += 1
                            else:
                                samples_stat[sampleRecord.name]['ExonicFunc_Unknown'] += 1
        vcf.close()
        samples_stat = self.samples_stat_multi(samples_stat)
        return samples_stat

    def samples_stat_multi(self, samples_stat):
        for sample in samples_stat:
            if self.variant_type == 'SNP':
                samples_stat[sample]['TivsTv'] = samples_stat[sample]['Ti'] / \
                    samples_stat[sample]['Tv']
            else:
                samples_stat[sample]['TotalDeletion'] = len(
                    samples_stat[sample]['Deletion'])
                samples_stat[sample]['TotalInsertion'] = len(
                    samples_stat[sample]['Insertion'])
            # samples_stat[sample]['Fraction of 1000genome(%)'] = samples_stat[sample]['1000genome'] / \
            #     samples_stat[sample]['Total_{}s'.format(self.variant_type)]
            # samples_stat[sample]['Fraction of dbsnp(%)'] = samples_stat[sample]['dbsnp'] / \
            #     samples_stat[sample]['Total_{}s'.format(self.variant_type)]
            samples_stat[sample]['Unpassing filter'] = self.snp_or_indel_record - \
                samples_stat[sample]['Passing filter']
        return samples_stat

    def get_gene_list(self,samples_stat):
        for sample in samples_stat:
            gene=samples_stat[sample]['gene']
            if len(gene) == 0:
                continue
            else:
                gp = GProfiler(user_agent='ExampleTool', return_dataframe=True)
                df=gp.profile(organism='hsapiens',query=gene)
                go=df[df['native'].str.contains('GO')]
                go.to_csv('{module}/{sample}/GO_FuncTerm.csv'.format(module=self.module,sample=sample),header=True,index=False,sep=',')
                self.plot_go(go,sample,'GO')
                kegg=df[df['native'].str.contains('KEGG')]
                kegg.to_csv('{module}/{sample}/KEGG_FuncTerm.csv'.format(module=self.module,sample=sample),header=True,index=False,sep=',')
                self.plot_go(kegg,sample,'KEGG')
                df=gp.convert(organism='hsapiens',query=gene,target_namespace='ENTREZGENE_ACC')
                df.to_csv('{module}/{sample}/Entrez_Gene_converted.csv'.format(module=self.module,sample=sample),header=True,index=False,sep=',')
                with open('{module}/{sample}/gene_list.txt'.format(module=self.module,sample=sample),'wt') as f:
                    f.write('\n'.join(gene))


    def plot_go(self,df,sample,db):
        df=df.reset_index()
        df['log_pvalue']=-np.log10(df['p_value'])
        df=df.sort_values('log_pvalue',ascending=False)
        df['label']=df['name']+'  '+df['native']
        plt.clf()
        plt.figure(figsize=(10,10))
        ax=plt.subplot(111)
        ax.barh('label','log_pvalue',data=df.loc[0:20,],color=['skyblue', 'red', 'orange', 'purple', 'cyan'])
        ax.set_xlabel('-log10(P_value)',fontsize=15)
        ax.set_ylabel('Function',labelpad=50,fontsize=15)
        ax.set_title('{} Enrichment Function'.format(db),fontsize=20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.axis(ymin=0.5)
        plt.savefig('{module}/{sample}/{db}_FuncTerm.png'.format(module=self.module,sample=sample,db=db),format='png',dpi=200,bbox_inches='tight')
        plt.close()


    def format_stat(self):
        samples_stat = self.parse_vcf()
        df = pd.DataFrame.from_dict(samples_stat, orient='index')
        df = df.reset_index()
        df = df.rename(columns={'index': 'Samples'})
        if self.variant_type == 'SNP':
            distribution = df[self.snp_col]
        else:
            distribution = df[self.indel_col]
        function = df[self.func_col]
        distribution.to_csv('{}/{}.stat.csv'.format(self.module,
                                                    self.variant_type), index=False, header=True)  # sample level
        function.to_csv('{}/{}.func.csv'.format(self.module, self.variant_type),
                        index=False, header=True)  # sample level
        self.indel_length_format(samples_stat, 'Deletion')
        self.indel_length_format(samples_stat, 'Insertion')
        self.site_depth_acumulative_plot(samples_stat, 'depth')
        self.site_depth_acumulative_plot(samples_stat, 'quality')
        self.get_gene_list(samples_stat)

    def site_depth_acumulative_plot(self, samples_stat, stat):
        # with open('sample_stats.json','wt') as f:
        #     json.dump(samples_stat,f)
        for sample in samples_stat:
            try:
                os.mkdir(os.path.join(self.module, sample))
            except OSError:
                pass
            depth = samples_stat[sample][stat]
            df = pd.DataFrame({stat: depth})
            df2 = {}
            for group in df.groupby(stat):
                df2[group[0]] = len(group[1])
            df2 = pd.DataFrame.from_dict(df2, orient='index')
            df2 = df2.reset_index()
            df2.columns = [stat, 'dense']
            df2['acumulative'] = 0
            df2 = df2.sort_values(stat)
            for i in range(len(df2)):
                df2.loc[i, 'acumulative'] = sum(df2.loc[:i, 'dense'])
            df2['acumulative'] = df2['acumulative']/len(df)
            plt.clf()
            plt.plot(stat, 'acumulative', data=df2, color='skyblue')
            plt.fill_between(stat, 'acumulative', data=df2, color='skyblue')
            plt.xlabel('{}'.format(stat))
            plt.ylabel('acumulative {} proportion(%)'.format(
                self.variant_type))
            plt.title('{} {} acumulative'.format(self.variant_type, stat))
            plt.savefig("{module}/{sample}/{stat}_acumulative.png".format(
                module=self.module, sample=sample, stat=stat), format='png', bbox_inches='tight', dpi=200)
            plt.close()
            df2 = df2[[stat, 'acumulative']]
            df2.to_csv('{module}/{sample}/{prefix}_acumulative.csv'.format(module=self.module,
                                                                           sample=sample, prefix=stat), index=False, header=True)

    def indel_length_format(self, samples_stat, stat):
        if self.variant_type == 'INDEL':
            indel_collect = nesteddict()
            for sample in samples_stat:
                try:
                    os.mkdir(os.path.join(self.module, sample))
                except OSError:
                    pass

                indel_stat = samples_stat[sample][stat]
                keys = set(indel_stat)
                for k in keys:
                    count = indel_stat.count(k)
                    indel_collect[stat][k] = count
                indel_length_distribution = pd.DataFrame.from_dict(
                    indel_collect, orient='index').T
                indel_length_distribution = indel_length_distribution.reset_index()
                indel_length_distribution = indel_length_distribution.rename(
                    columns={'index': stat, stat: 'count'})
                indel_length_distribution.to_csv('{module}/{sample}/{stat}_distribution.csv'.format(
                    module=self.module, sample=sample, stat=stat), index=False, header=True)
                self.indel_length_plot(indel_length_distribution, stat, sample)
                self.indel_length_cumulative_plot(
                    indel_length_distribution, stat, sample)
                # os.chdir(sample)

    def indel_length_plot(self, df, prefix, sample):
        plt.clf()
        plt.plot(prefix, 'count', data=df)
        plt.xlabel('{} length(bp)'.format(prefix))
        plt.ylabel('count')
        plt.title('{} length distribution'.format(prefix))
        plt.savefig("{module}/{sample}/{prefix}.length_distribution.png".format(
            module=self.module, sample=sample, prefix=prefix), format='png')
        plt.close()

    def indel_length_cumulative_plot(self, df, prefix, sample):
        df['acumulative'] = 0
        #df['count'] = df['count'].astype(int)
        df = df.sort_values(prefix)
        df2 = {}
        for group in df.groupby(prefix):
            df2[group[0]] = len(group[1])
        df2 = pd.DataFrame.from_dict(df2, orient='index')
        df2 = df2.reset_index()
        df2.columns = [prefix, 'dense']
        df2['acumulative'] = 0
        df2 = df2.sort_values(prefix)
        for i in range(len(df2)):
            df2.loc[i, 'acumulative'] = sum(df2.loc[:i, 'dense'])
        df2['acumulative'] = df2['acumulative']/len(df)
        plt.clf()
        plt.plot(prefix, 'acumulative', data=df2, color='skyblue')
        plt.fill_between(prefix, 'acumulative', data=df2, color='skyblue')
        plt.xlabel('{} length(bp)'.format(prefix))
        plt.ylabel('acumulative proportion(%)')
        plt.title('{} length acumulative'.format(prefix))
        plt.savefig("{module}/{sample}/{prefix}.length_acumulative.png".format(
            module=self.module, sample=sample, prefix=prefix), bbox_inches='tight', dpi=200, format='png')
        plt.close()
        df = df[[prefix, 'acumulative']]
        df.to_csv('{module}/{sample}/{prefix}_acumulative.csv'.format(module=self.module,
                                                                      sample=sample, prefix=prefix), index=False, header=True)


if __name__ == '__main__':
    vcf = sys.argv[1]
    VcfStat(vcf, 'SNP', '3.SNP')
    VcfStat(vcf, 'INDEL', '4.INDEL')
