#!/usr/bin/env python3
#coding:utf-8
import os
import subprocess



class SV(object):
    reference = "/share/data1/genome/hs38DH.fa"
    delly='/share/data3/lianlin/soft/delly/delly'
    excludeTemplates='/share/data3/lianlin/soft/delly/excludeTemplates/human.hg38.excl.tsv'

    def __init__(self,groups:dict,working_space:str):
        self.groups=groups
        self.working_space=os.path.abspath(working_space)
 
    
    def dellySV(self):
        bams=[]
        if len(self.groups) >1:
            for sample in self.groups:
                bam=sample+'.dedup.bam'
                bams.append(bam)
            step1='{tool} call -x {excludeTemplates} -o sv.bcf -g {reference} {bams} ;'.format(tool=self.delly,reference=self.reference,bams=' '.join(bams),excludeTemplates=self.excludeTemplates)
            step2='{tool} filter -f germline -p -o delly.bcf sv.bcf'.format(tool=self.delly)
        else:
            step1='{tool} call -x {excludeTemplates} -o sv.bcf -g {reference} {sample}.dedup.bam ;'.format(tool=self.delly,reference=self.reference,sample=list(self.groups.keys())[0],excludeTemplates=self.excludeTemplates)
            step2='{tool} filter -p -o delly.bcf sv.bcf ;'.format(tool=self.delly)
        step3='bcftools view -O z -o delly.vcf.gz delly.bcf ;'
        step4='tabix delly.vcf.gz ;'
        sge = "#$ -N delly\n#$ -pe smp 5\n#$ -q all.q\n#$ -cwd\nset -e\ncd {working_space}\nsource ~/.bash_profile\n".format(working_space=self.working_space)
        with open('{}/sv.bat'.format(self.working_space),'wt') as f:
            f.write(sge)
            f.write('\n'.join([step1,step2,step3,step4]))
            f.write('wait ;')
        
