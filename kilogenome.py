#!/usr/bin/env python2
#coding:utf-8
import gzip
import os
import sys
import json
from collections import defaultdict

class Freq(object):
    def __init__(self,vcf):
        self.vcf=vcf
        self._open()
        self.get_head()
		
    def _open(self):
        if os.path.splitext(self.vcf)[-1] == '.gz':
            self.input=gzip.open(self.vcf)
        else:
            self.input=open(self.vcf)

    def get_head(self):
        self.vcfhead=[]
        for i in self.input:
            if i.startswith('##'):
                self.vcfhead.append(i)
            elif i.startswith('#CHROM'):
                self.samples=i.strip().split('\t')
                self.vcfhead.append(i)
            else:
                break

    def pop(self):
        with open('/share/data3/lianlin/soft/bin/wgs/1kgenome.json','r') as f:
            for i in f:
                pops=json.loads(i)
        return pops

    def f(self):
        #print self.samples
        pops=self.pop()
        population=defaultdict(list)
        for p in pops:
            samples=pops[p]
            for s in samples:
                #print s
                pos=[i for i,v in enumerate(self.samples) if v.find(s) != -1][0]
                #print pos
                population[p].append(pos)
        #print population['EAS']
        return population

    def parse_vcf(self):
        for i in self.input:
            if i.startswith('#'):
                continue
            else:
                line=i.strip().split()
                rsid=line[2]
                if rsid == '.':
                    continue
                else:
                    ref=line[3]
                    alt=line[4]
                    if ref not in ['A','C','G','T'] or alt not in ['A','C','G','T']:
                        continue
                    else:
                        if len(ref) > 1 or len(alt) > 1:
                            continue
                        else:
                            yield line

    def cal_freq(self,kilo):
        pf={}
        for p in kilo:
            AA=kilo[p].count('0|0')
            Aa=kilo[p].count('0|1')
            aa=kilo[p].count('1|1')
            AC=2.0*aa + Aa
            AN=2*(AA+Aa+aa)
            #print p,AA,Aa,aa
            f=AC/AN
            pf[p]=round(f,4)
        return pf

    def std_freq(self,pf):
        info=[]
        for p in pf:
            #print p
            _=p+'='+str(pf[p])
            info.append(_)
        
        return info

    def filter(self):
        population=self.f()
        for line in self.parse_vcf():
            part1='\t'.join(line[0:7])
            _format=line[8]
            kilo=defaultdict(list)
            for p in population: 
                samples_pos=population[p]
                for sp in samples_pos:
                    genotype=line[sp]
                    kilo[p].append(genotype)
            pf=self.cal_freq(kilo)
            info=';'.join(self.std_freq(pf))
            if pf['EAS'] >= 0.001:
                _=[line[i] for i in population['EAS']]
                part2='\t'.join(_)
                all=part1+'\t'+info+'\t'+_format+'\t'+part2
                yield all


    def write_out(self):
        prefix=os.path.basename(self.vcf).split('.')[0]
        with open('/share/data3/lianlin/soft/bin/wgs/eas/eas.{}.kilogenome.phase3.genotype.vcf'.format(prefix),'w') as f:
            self.vcfhead.insert(5,'##INFO=<ID=EAS,Number=A,Type=Float,Description="altelative Allele frequency in the EAS populations in the range (0,1)">\n')
            self.vcfhead.insert(6,'##INFO=<ID=EUR,Number=A,Type=Float,Description="altelative Allele frequency in the EAS populations in the range (0,1)">\n')
            self.vcfhead.insert(7,'##INFO=<ID=AFR,Number=A,Type=Float,Description="altelative Allele frequency in the EAS populations in the range (0,1)">\n')
            self.vcfhead.insert(8,'##INFO=<ID=AMR,Number=A,Type=Float,Description="altelative Allele frequency in the EAS populations in the range (0,1)">\n')
            self.vcfhead.insert(9,'##INFO=<ID=SAS,Number=A,Type=Float,Description="altelative Allele frequency in the EAS populations in the range (0,1)">\n')
            population=self.f()
            head=[self.samples[i] for i in population['EAS']]
            self.vcfhead[-1]='\t'.join(self.samples[0:9]) + '\t' + '\t'.join(head) +'\n'
            print len(self.vcfhead[-1].split())
            for h in self.vcfhead:
                f.write(h)
            for i in self.filter():
                #print i
                f.write(i+'\n')
if __name__ == '__main__':
    vcf=sys.argv[1]
    Freq(vcf).write_out()
            






        
