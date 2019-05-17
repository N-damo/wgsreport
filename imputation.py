#!/use/bin/env python2
# coding:utf-8
import gzip
import os
import multiprocessing
from multiprocessing import Pool
import sys
import subprocess


class Imputation(object):
    vcfhead = ['##fileformat=VCFv4.2',
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    X = {"chrX:1-2781479": "chrX_part1", "chrX:2781480-155701382": "chrX",
         "chrX:155701383-156040895": "chrX_part2"}
    # x染色体需要分三块进行

    def __init__(self, vcf, chromosome):
        self.vcf = vcf
        self.chromosome = chromosome
        self._open()
        self.get_head()
        self._block()

    def _open(self):
        if os.path.splitext(self.vcf)[-1] == '.gz':
            self.input = gzip.open(self.vcf, 'r')
        else:
            self.input = open(self.vcf, 'r')

    def get_head(self):
        for i in self.input:
            if i.startswith('#CHROM'):
                self.vcfhead.append(i.strip())
                break

    def std(self, lines):
        result = []
        for line in lines:
            part = line.split('\t')
            part[0] = part[0].replace('chr', '')
            part[6] = '.'
            part[7] = '.'
            part[8] = 'GT:PL'
            for i in range(len(part[9:])):
                gt = part[9+i].split(':')[0]
                pl = part[9+i].split(':')[-1]
                if gt == './.':
                    pl = '.'
                _ = ':'.join([gt, pl])
                part[9+i] = _
            result.append('\t'.join(part))
        return result

    def parse_vcf(self):
        chr_list = []
        part = []
        if self.chromosome != 'chrX':
            cmd = 'tabix {} {}'.format(self.vcf, self.chromosome)
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            stdout = p.communicate()[0]
            result = self.std(stdout.strip().split('\n'))
            chr_list.append(result)
            part.append(self.chromosome)
        elif self.chromosome == 'chrX':
            for i in self.X:
                cmd = 'tabix {} {}'.format(self.vcf, i)
                p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                stdout = p.communicate()[0]
                if stdout != '':
                    result = self.std(stdout.strip().split('\n'))
                    chr_list.append(result)
                    part.append(self.X[i])
                else:
                    print '{} {} tabix result is bulk'.format(i, self.X[i])
        return chr_list, part

    def _block(self):
        l, part = self.parse_vcf()
        for p in range(len(l)):
            chr_list = l[p]
            n = len(chr_list)
            block_nm = 50000  # 每一个block 50000位点数
            block_cross = 3000  # 相邻block之间交叉位点数3000
            start = 0
            block_it = 0
            print self.chromosome, part[p], block_nm, block_cross, n
            while (block_nm + start) <= n:  # 右半部分下标有可能溢出
                block = chr_list[start:block_nm+start]
                block_name = 'wgs.{}.{}.vcf.gz'.format(
                    part[p], block_it)
                self.out(block, block_name)
                self.pbs(block_name, part[p], block_it)
                start = start + block_nm-block_cross
                block_it += 1
            else:
                block = chr_list[start:]
                block_name = 'wgs.{}.{}.vcf.gz'.format(
                    part[p], block_it)
                self.out(block, block_name)
                self.pbs(block_name, part[p], block_it)

    def out(self, block, block_name):
        with gzip.open(block_name, 'w') as f:
            merge = self.vcfhead+block
            for i in merge:
                f.write(i+'\n')

    def pbs(self, block_name, part, block_it):
        self.workspace = os.path.abspath('./')
        cmd1 = "#!/bin/bash\n#$ -N wgs.{}.{}.pbs\n#$ -pe smp 8\n#$ -q all.q\n#$ -cwd\n".format(
            part, block_it)
        cmd2 = "set -e\ncd {}\n".format(self.workspace)
        cmd3 = "java -Djava.io.tmpdir=/share/data1/tmp/ -Xss5m -Xmx128g -jar /share/data1/local/bin/beagle.27Jan18.7e1.jar gtgl={vcf}.vcf.gz ref=/share/data1/PublicProject/BEAGLE/eas.allgrch38.genotype.bref map=/share/data1/PublicProject/BEAGLE/plink.{plink}.GRCh38.map out=mid_{vcf} lowmem=true ;\nwait ;".format(
            vcf=block_name.replace('.vcf.gz', ''), plink=part)
        with open('wgs.{}.{}.pbs'.format(part, block_it), 'w') as f:
            f.write(cmd1+cmd2+cmd3)


def shuf():
    cmd = 'for i in $(ls *.pbs);do echo qsub $i \;;done >alljob.bat;chmod +x alljob.bat'
    subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    vcf = sys.argv[1]
    pool = Pool()
    chrs = []
    for i in range(1, 23):
        chr = 'chr' + str(i)
        chrs.append(chr)
    chrs.append('chrX')
    for chromosome in chrs:
        pool.apply_async(Imputation, args=(vcf, chromosome,))
    pool.close()
    pool.join()
    # Imputation(sys.argv[1],'chrX',sys.argv[2])
    shuf()
