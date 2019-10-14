#!/usr/bin/env python2
#coding:utf-8
import gzip
import os
import multiprocessing
from multiprocessing import Pool
import sys
import subprocess



class Postimputation(object):
    def __init__(self,chromosome):
        self.chromosome=chromosome
        self.list_file()
        self.vcf_merge()

    def list_file(self):
        #wgs.chr14.0.vcf.gz
        files=os.listdir('./')
        self.vcf_list=[]
        for f in files:
            if f.startswith('mid_wgs.{}.'.format(self.chromosome)) and f.endswith('vcf.gz'):
                self.vcf_list.append(f)
        if len(self.vcf_list) == 0:
            print 'there is no mid_wgs.{}.vcf.gz files'.format(self.chromosome)
            os._exit()
        self.vcf_list=sorted(self.vcf_list,key=lambda x:int(x.split('.')[2]),reverse=False)
        #print self.vcf_list

    def first_vcf(self):
        self.result=[]
        with gzip.open(self.vcf_list[0],'r') as f:
            for i in f:
                self.result.append(i)

    def vcf_append(self):
        self.first_vcf()
        if len(self.vcf_list) >= 2:
            for vcf in self.vcf_list[1:]:
                with gzip.open(vcf,'r') as f:
                    n=1
                    for i in f:
                        if i.startswith('#'):
                            continue
                        else:
                            if n > 3000:
                                #print i
                                self.result.append(i)
                                n += 1
                            else:
                                n += 1
                    
                            
                       
    def vcf_merge(self):
        self.vcf_append()
        with gzip.open("final.{}.vcf.gz".format(self.chromosome),'w') as f:
            for i in self.result:
                f.write(i)

def all_merge(batch):
    #cmd='bcftools merge --threads 32 -o {}.vcf.gz -O z -l vcf_list.txt'.format(batch)
    cmd="(less final.chr1.vcf.gz|head -30;less final.chr*vcf.gz|grep -v '^#')|vcf-sort -c|bgzip >{}.vcf.gz".format(batch)
    subprocess.call(cmd,shell=True)
    

def chromosomes():
    chrs = []
    for i in range(1, 23):
        chr = 'chr' + str(i)
        chrs.append(chr)
    chrs.append('chrX_part1')
    chrs.append('chrX')
    chrs.append('chrX_part2')
    return chrs

if __name__ == '__main__':  
    batch=sys.argv[1]
    pool = Pool()
    chrs=chromosomes()
    for chromosome in chrs:
        pool.apply_async(Postimputation, args=(chromosome,))
    pool.close()
    pool.join()
    all_merge(batch)
    #Postimputation('chrX_part1')
