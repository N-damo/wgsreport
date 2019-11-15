#!/usr/bin/env python3
#coding:utf-8
import os
import subprocess



class CNV(object):
    reference = "/share/data1/genome/hs38DH.fa"
    chrom="/share/data3/lianlin/genome/hg38"
    cnvnator='/share/data3/lianlin/soft/CNVnator-master/cnvnator'
    cnvnator2vcf='/share/data3/lianlin/soft/CNVnator-master/cnvnator2VCF.pl'
    bin_size=1000
    
    def __init__(self,groups:dict,working_space:str):
        self.groups=groups
        self.working_space=os.path.abspath(working_space)
          
    def CnvnatorCNV(self):
        for sample in self.groups:
            root=sample+'.root'
            if os.path.exists(os.path.join(self.working_space,root)):
                subprocess.call('rm root',shell=True)
            step1="{tool} -root {sample}.root -tree {sample}.dedup.bam -chrom $(seq -f 'chr%g' 1 22) chrX chrY ;".format(tool=self.cnvnator,sample=sample)
            step2='{tool} -root {sample}.root -his {bin_size} -fasta {reference} ;'.format(tool=self.cnvnator,sample=sample,reference=self.reference,bin_size=self.bin_size)
            step3='{tool} -root {sample}.root -stat {bin_size} ;'.format(tool=self.cnvnator,sample=sample,bin_size=self.bin_size)
            step4='{tool} -root {sample}.root -partition {bin_size} ;'.format(tool=self.cnvnator,sample=sample,bin_size=self.bin_size)
            step5='{tool} -root {sample}.root -call {bin_size} >{sample}.out ;'.format(tool=self.cnvnator,sample=sample,bin_size=self.bin_size)
            step6='{tool} -prefix {sample} -reference hg38 {sample}.out {chrom} |bgzip >{sample}.vcf.gz ;'.format(tool=self.cnvnator2vcf,sample=sample,chrom=self.chrom)
            step7='tabix {sample}.vcf.gz ;'.format(sample=sample)
            sge = "#$ -N {sample}.cnvnator\n#$ -pe smp 5\n#$ -q all.q\n#$ -cwd\nset -e\ncd {working_space}\nsource ~/.bash_profile\n".format(sample=sample,working_space=self.working_space)
            with open('{}/{}_cnv.bat'.format(self.working_space,sample),'wt') as f:
                f.write(sge)
                f.write('\n'.join([step1,step2,step3,step4,step5,step6,step7]))
                f.write('wait ;')