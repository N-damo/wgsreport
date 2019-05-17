#!/usr/bin/env python
import pandas as pd
import subprocess
import argparse
import os
import sys
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(
        description='Germline short variant discovery (SNPs + Indels)')
    req = parser.add_argument_group('required arguments')
    opt = parser.add_argument_group('optional arguments')
    req.add_argument('-i', '--input', help='input your file', required=True)
    req.add_argument('-o', '--output', help='working directory', required=True)
    req.add_argument('-f', '--final_name',
                     help='final result name', required=True)
    opt.add_argument('-t', '--thread',
                     help='choose thread mode,the default is 10', default='10')
    opt.add_argument('-m', '--method',
                     help='wgs or wes,the default is wes', default='wes')
    opt.add_argument('-ip', '--padding',
                     help='Amount of padding (in bp) to add to each interval you are including, default=0.', default='0')
    opt.add_argument('-v', '--variant_filter',
                     help='hardfilter or vqsr,the default is hardfilter', default='hardfilter')
    opt.add_argument(
        '-q', '--queue', help='choose compution node,the default is all.q', default='all.q')
    args = parser.parse_args()
    return args


class sequencing(object):

    gatk = '/share/data1/src/gatk/gatk'
    dbsnp = '/share/data1/PublicProject/GATK_bundle/dbsnp150_chr.vcf.gz'
    indel = '/share/data1/PublicProject/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
    snp = '/share/data1/PublicProject/GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz'
    bed = '/share/data2/leon/S07084713_Regions.bed'
    gatkbed = '/share/data1/PublicProject/GATK_bundle/wgs_calling_regions.hg38.bed '
    interval = '/share/data2/leon/list.interval_list'
    reference = '/share/data1/genome/hs38DH.fa'
    Nextera = '-f CTGTCTCTTATACACATCTCCGAGCCCACGAGACIIIIIIIIATCTCGTATGCCGTCTTCTGCTTG -s CTGTCTCTTATACACATCTGACGCTGCCGACGAIIIIIIIIGTGTAGATCTCGGTGGTCGCCGTATCATT -c TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
    Truseq = '-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIIATCTCGTATGCCGTCTTCTGCTTG -s AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTIIIIIIIIGTGTAGATCTCGGTGGTCGCCGTATCATT -c ACACTCTTTCCCTACACGACGCTCTTCCGATCT'

    def __init__(self):
        self.df = self.input()
        self.group = self.getgroup()
        self.first_step, self.second_step = self.combine()
        self._write_first()
        self._write_second()

    def input(self):
        df = pd.read_table(args.input, header=None)
        return df

    def getgroup(self):
        l = []
        group = self.df.groupby(5)
        for groups in group:
            l.append(groups[1])
        return l

    def number(self):
        if len(self.df) == 1:
            return True
        else:
            return False

    def trim(self):
        l = []
        for i in self.df.index.values:
            fq1 = self.df.iloc[i, 0]
            fq2 = self.df.iloc[i, 1]
            run = self.df.iloc[i, 3]
            sample = self.df.iloc[i, 5]
            if self.df.iloc[i, 2] == 'Truseq':
                cmd = 'java -jar /share/data1/local/bin/trimmomatic-0.38.jar PE -threads {thread} -phred33 {fq1} {fq2} {run}.{sample}_r1.fq.gz {run}.{sample}_r1_unpaired.fq.gz {run}.{sample}_r2.fq.gz {run}.{sample}_r2_unpaired.fq.gz ILLUMINACLIP:/share/data1/local/bin/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;'.format(
                    thread=args.thread, fq1=fq1, fq2=fq2, run=run, sample=sample)
                l.append(cmd)
            elif self.df.iloc[i, 2] == 'Nextera':
                cmd = 'java -jar /share/data1/local/bin/trimmomatic-0.38.jar PE -threads {thread} -phred33 {fq1} {fq2} {run}.{sample}_r1.fq.gz {run}.{sample}_r1_unpaired.fq.gz {run}.{sample}_r2.fq.gz {run}.{sample}_r2_unpaired.fq.gz ILLUMINACLIP:/share/data1/local/bin/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;'.format(
                    thread=args.thread, fq1=fq1, fq2=fq2, run=run, sample=sample)
                l.append(cmd)
            else:
                print '{} :\nadaptor error, Nextera or Truseq accept,not nextera or truseq or others'.format(
                    self.df.iloc[i, 2])
                sys.exit()
        return l

    def bwa(self):
        l = []
        for i in self.df.index.values:
            run = self.df.iloc[i, 3]
            sample = self.df.iloc[i, 5]
            library = self.df.iloc[i, 4]
            cmd = '/share/data1/src/bwa/bwakit/seqtk mergepe {run}.{sample}_r1.fq.gz {run}.{sample}_r2.fq.gz |bwa mem -p -t {thread} -R"@RG\\tID:{run}\\tLB:{library}\\tSM:{sample}\\tPU:flowcell\\tPL:Illumina" /share/data1/genome/hs38DH.fa - 2>{run}.log.bwamem |k8 /share/data1/src/bwa/bwakit/bwa-postalt.js -p {run}.hla /share/data1/genome/hs38DH.fa.alt | samtools sort -@ 2 -m8g - -o {run}.{sample}.aln.bam ;'.format(
                thread=args.thread, run=run, library=library, sample=sample)
            l.append(cmd)
        return l

    def MarkDuplicates(self):
        l = []
        for i in range(len(self.group)):
            bams = ' -I '.join(list(self.group[i]
                                    [3]+'.'+self.group[i][5]+'.aln.bam'))
            sample = self.group[i][5].unique()[0]
            cmd = '/share/data1/src/gatk/gatk MarkDuplicates -I {bams} -O {sample}.dedup.bam -M {sample}.dedup.log -PG null --TMP_DIR ~/tmp/{sample} ;'.format(
                bams=bams, sample=sample)
            l.append(cmd)
        return l

    def index(self):
        l = []
        for i in range(len(self.group)):
            sample = self.group[i][5].unique()[0]
            cmd = 'samtools index -@ 4 {sample}.dedup.bam ;'.format(
                sample=sample)
            l.append(cmd)
        return l

    def Collect(self):
        l = []
        for i in range(len(self.group)):
            sample = self.group[i][5].unique()[0]
            if args.method == 'wes':
                cmd = '/share/data1/src/gatk/gatk CollectHsMetrics -I {sample}.dedup.bam -O {sample}_hs_metrics.txt -BI {interval} -TI {interval} -R /share/data1/genome/hs38DH.fa &'.format(
                    sample=sample, interval=self.interval)
                l.append(cmd)
            elif args.method == 'wgs':
                cmd = '/share/data1/src/gatk/gatk CollectWgsMetrics -I {sample}.dedup.bam -O {sample}_wgs_metrics.txt -R /share/data1/genome/hs38DH.fa &'.format(
                    sample=sample)
                l.append(cmd)
        return l

    def BaseRecalibrator(self):
        l = []
        for i in range(len(self.group)):
            sample = self.group[i][5].unique()[0]
            if args.method == 'wes':
                cmd = '/share/data1/src/gatk/gatk BaseRecalibrator -I {sample}.dedup.bam -O {sample}.recal.table --known-sites {snp} --known-sites {indel} -L {bed} -ip {ip} -R /share/data1/genome/hs38DH.fa ;'.format(
                    sample=sample, bed=self.bed, ip=args.padding, snp=self.snp, indel=self.indel)
                l.append(cmd)
            elif args.method == 'wgs':
                cmd = '/share/data1/src/gatk/gatk BaseRecalibrator -I {sample}.dedup.bam -O {sample}.recal.table --known-sites {snp} --known-sites {indel} -R /share/data1/genome/hs38DH.fa ;'.format(
                    sample=sample, snp=self.snp, indel=self.indel)
                l.append(cmd)
        return l

    def ApplyBQSR(self):
        l = []
        for i in range(len(self.group)):
            sample = self.group[i][5].unique()[0]
            if args.method == 'wes':
                cmd = '/share/data1/src/gatk/gatk ApplyBQSR -I {sample}.dedup.bam -O {sample}.dedup.bqsr.bam -bqsr {sample}.recal.table -L {bed} -ip {ip} -R /share/data1/genome/hs38DH.fa ;'.format(
                    sample=sample, bed=self.bed, ip=args.padding)
                l.append(cmd)
            elif args.method == 'wgs':
                cmd = '/share/data1/src/gatk/gatk ApplyBQSR -I {sample}.dedup.bam -O {sample}.dedup.bqsr.bam -bqsr {sample}.recal.table -R /share/data1/genome/hs38DH.fa ;'.format(
                    sample=sample)
                l.append(cmd)
        return l

    def HaplotypeCaller(self):
        l = []
        for i in range(len(self.group)):
            sample = self.group[i][5].unique()[0]
            if args.method == 'wes':
                cmd = '/share/data1/src/gatk/gatk HaplotypeCaller -I {sample}.dedup.bqsr.bam -O {sample}.HC.g.vcf.gz --emit-ref-confidence GVCF --dbsnp {dbsnp} -L {bed} -ip {ip} -R /share/data1/genome/hs38DH.fa ;'.format(
                    sample=sample, bed=self.bed, ip=args.padding, dbsnp=self.dbsnp)
                l.append(cmd)
            elif args.method == 'wgs':
                chr = []
                for i in range(1, 23, 1):
                    chr.append('chr' + str(i))
                chr.append('chrX')
                chr.append('chrY')
                chr.append('chrM')
                for i in chr:
                    cmd = '/share/data1/src/gatk/gatk HaplotypeCaller -I {sample}.dedup.bqsr.bam -O {sample}.{chromosome}.HC.g.vcf.gz --emit-ref-confidence GVCF --dbsnp {dbsnp} -R /share/data1/genome/hs38DH.fa &'.format(
                        sample=sample, chromosome=i, dbsnp=self.dbsnp)
                    l.append(cmd)
                m = []
                for i in chr:
                    g = '{sample}.{chromosome}.HC.g.vcf.gz'.format(
                        sample=sample, chromosome=i)
                    m.append(g)
                cmd = '/share/data1/src/gatk/gatk CombineGVCFs -V {vcfs} -O {sample}.HC.g.vcf.gz -R /share/data1/genome/hs38DH.fa ;'.format(
                    vcfs=' -V '.join(m), sample=sample)
                l.append(cmd)
        return l

    def covstat(self):
        l = []
        if args.method == 'wgs':
            for i in range(len(self.group)):
                sample = self.group[i][5].unique()[0]
                cmd = 'sam covstat {sample}.dedup.bqsr.bam > {sample}.covstat &'.format(
                    sample=sample)
                l.append(cmd)
        elif args.method == 'wes':
            l = []
        return l

    def cram(self):
        l = []
        for i in range(len(self.group)):
            sample = self.group[i][5].unique()[0]
            cmd = 'samtools view -C -T /share/data1/genome/hs38DH.fa -@ 4 -o {sample}.dedup.bqsr.cram {sample}.dedup.bqsr.bam &'.format(
                sample=sample)
            l.append(cmd)
        return l

    def GenomicsDBImport(self):
        l = []
        inputfile = ' -V '.join(self.df[5].unique()+'.HC.g.vcf.gz')
        if args.method == 'wgs':
            for i in range(1, 23):
                cmd = '/share/data1/src/gatk/gatk GenomicsDBImport -V {vcfs} --genomicsdb-workspace-path chr{n}.database --intervals chr{n} --reader-threads 6 &'.format(
                    vcfs=inputfile, n=i)
                l.append(cmd)
            for i in range(5, 25, 5):
                l.insert(i, 'wait ;')
            l.append('wait ;')
            chr = ['chrX', 'chrY', 'chrM']
            for i in chr:
                cmd = '/share/data1/src/gatk/gatk GenomicsDBImport -V {vcfs} --genomicsdb-workspace-path {chromosome}.database --intervals {chromosome} --reader-threads 6 &'.format(
                    vcfs=inputfile, chromosome=i)
                l.append(cmd)
            l.append('wait ;')
        elif args.method == 'wes':
            if len(self.df[5].unique()) == 1:
                l = []
            else:
                cmd = '/share/data1/src/gatk/gatk CombineGVCFs -V {vcfs} -O cohort.g.vcf.gz -R /share/data1/genome/hs38DH.fa ;'.format(
                    vcfs=inputfile)
                l.append(cmd)
        return l

    def GenotypeGVCFs(self):
        l = []
        if args.method == 'wes':
            if len(self.df[5].unique()) == 1:
                cmd = '/share/data1/src/gatk/gatk GenotypeGVCFs -V {vcfs} -O han.vcf.gz -R /share/data1/genome/hs38DH.fa --dbsnp {dbsnp} -L {bed} -ip {ip} ;'.format(
                    vcfs=' -V '.join(self.df[5].unique()+'.HC.g.vcf.gz'), dbsnp=self.dbsnp, bed=self.bed, ip=args.padding)
                l.append(cmd)
            else:
                cmd = '/share/data1/src/gatk/gatk GenotypeGVCFs -V cohort.g.vcf.gz -O han.vcf.gz -R /share/data1/genome/hs38DH.fa --dbsnp {dbsnp} -L {bed} -ip {ip} ;'.format(
                    dbsnp=self.dbsnp, bed=self.bed, ip=args.padding)
                l.append(cmd)
        elif args.method == 'wgs':
            for i in range(1, 23):
                cmd = '/share/data1/src/gatk/gatk GenotypeGVCFs -V gendb://chr{n}.database -O han.chr{n}.vcf.gz -R /share/data1/genome/hs38DH.fa --dbsnp {dbsnp} &'.format(
                    n=i, dbsnp=self.dbsnp)
                l.append(cmd)
            for i in range(5, 25, 5):
                l.insert(i, 'wait ;')
            l.append('wait ;')
            chr = ['chrX', 'chrY', 'chrM']
            for i in chr:
                cmd = '/share/data1/src/gatk/gatk GenotypeGVCFs -V gendb://{chromosome}.database -O han.{chromosome}.vcf.gz -R /share/data1/genome/hs38DH.fa --dbsnp {dbsnp} &'.format(
                    chromosome=i, dbsnp=self.dbsnp)
                l.append(cmd)
            l.append('wait ;')
            c = []
            for i in range(1, 23):
                c.append('han.chr{}.vcf.gz'.format(i))
            cmd = "/share/data1/src/gatk/gatk MergeVcfs -I {vcfs} -I han.chrX.vcf.gz -I han.chrY.vcf.gz -I han.chrM.vcf.gz -O han.vcf.gz ;".format(
                vcfs=' -I '.join(c))
            l.append(cmd)
        return l

    def variantfilter(self):
        if args.variant_filter == 'hardfilter':
            l = []
            cmd = '/share/data1/src/gatk/gatk SelectVariants -select-type SNP -V han.vcf.gz -O han.snp.vcf.gz ;'
            l.append(cmd)
            cmd = '/share/data1/src/gatk/gatk SelectVariants -select-type INDEL -V han.vcf.gz -O han.indel.vcf.gz ;'
            l.append(cmd)
            cmd = '/share/data1/src/gatk/gatk VariantFiltration -V han.snp.vcf.gz --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "my_snp_Filter" -O han.snp.hardfilter.vcf.gz ;'
            l.append(cmd)
            cmd = '/share/data1/src/gatk/gatk VariantFiltration -V han.indel.vcf.gz --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "my_indel_Filter" -O han.indel.hardfilter.vcf.gz ;'
            l.append(cmd)
            cmd = "/share/data1/src/gatk/gatk MergeVcfs -R /share/data1/genome/hs38DH.fa -I han.snp.hardfilter.vcf.gz -I han.indel.hardfilter.vcf.gz -O {finalname}.snp.indel.hardfilter.vcf.gz ;".format(
                finalname=args.final_name)
            l.append(cmd)
        elif args.variant_filter == 'vqsr':
            l = []
            cmd = '/share/data1/src/gatk/gatk VariantRecalibrator -R /share/data1/genome/hs38DH.fa -V han.vcf.gz -resource hapmap,known=false,training=true,truth=true,prior=15.0:/share/data1/PublicProject/GATK_bundle/hapmap_3.3.hg38.vcf.gz -resource omni,known=false,training=true,truth=false,prior=12.0:/share/data1/PublicProject/GATK_bundle/1000G_omni2.5.hg38.vcf.gz -resource 1000G,known=false,training=true,truth=false,prior=10.0:/share/data1/PublicProject/GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz -resource dbsnp,known=true,training=false,truth=false,prior=2.0:/share/data1/PublicProject/GATK_bundle/dbsnp150_chr.vcf.gz -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -tranche 99.0 -mode SNP -O output.snp.recal -tranches-file output.snp.tranches -rscript-file output.snp.plots.R ;'
            l.append(cmd)
            cmd = '/share/data1/src/gatk/gatk ApplyVQSR -R /share/data1/genome/hs38DH.fa -V han.vcf.gz -O han.snp.vqsr.vcf.gz --ts-filter-level 99.0 -tranches-file output.snp.tranches -recal-file output.snp.recal -mode SNP ;'
            l.append(cmd)
            cmd = '/share/data1/src/gatk/gatk VariantRecalibrator -R /share/data1/genome/hs38DH.fa -V han.snp.vqsr.vcf.gz  -mode INDEL -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -resource dbsnp,known=true,training=false,truth=false,prior=2.0:/share/data1/PublicProject/GATK_bundle/dbsnp150_chr.vcf.gz -resource mills,known=true,training=true,truth=true,prior=12.0:/share/data1/PublicProject/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --max-gaussians 6 -O output.indel.recal -tranches-file output.indel.tranches -rscript-file output.indel.plots.R ;'
            l.append(cmd)
            cmd = '/share/data1/src/gatk/gatk ApplyVQSR -R /share/data1/genome/hs38DH.fa -V han.snp.vqsr.vcf.gz -O {finalname}.snp.indel.vqsr.vcf.gz -tranches-file output.indel.tranches -recal-file output.indel.recal -mode INDEL ;'.format(
                finalname=args.final_name)
            l.append(cmd)
        return l

    def annotation(self):
        l = []
        if args.method == 'wes':
            cmd = '/share/data1/src/annovar/table_annovar.pl {finalname}.snp.indel.hardfilter.vcf.gz /share/data1/src/annovar/humandb/ -buildver hg38 -out myanno -otherinfo -remove -protocol refGene,cytoBand,exac03,gnomad_exome,clinvar_20170905,dbnsfp33a, -operation g,r,f,f,f,f --vcfinput -nastring . -polish -thread 30'.format(
                finalname=args.final_name)
            l.append(cmd)
        elif args.method == 'wgs':
            l = []

        return l

    def combine(self):
        trim_adap = self.trim()
        bwa = self.bwa()
        dedup = self.MarkDuplicates()
        index = self.index()
        collect = self.Collect()
        bqsr = self.BaseRecalibrator()
        applybqsr = self.ApplyBQSR()
        hp = self.HaplotypeCaller()
        GenomicsDBImport = self.GenomicsDBImport()
        genotype = self.GenotypeGVCFs()
        VariantFiltration = self.variantfilter()
        cram = self.cram()
        covstat = self.covstat()
        annotation = self.annotation()
        first_step = pd.DataFrame(trim_adap + bwa + dedup + index +
                                  collect + bqsr + applybqsr + cram + covstat + hp)
        second_step = GenomicsDBImport + genotype + VariantFiltration + annotation
        return first_step, second_step

    def _write_first(self):
        for i in self.df[5].unique():
            with open(i, 'w') as f:
                f.write('#$ -N ' + i + '\n')
                f.write('#$ -pe smp ' + args.thread + '\n')
                f.write('#$ -q ' + args.queue + '\n')
                f.write('#$ -cwd\n')
                f.write('set -e\n')
                f.write('cd ' + args.output + '\n')
                f.write('source ~/.bash_profile\n')
                for j in self.first_step[self.first_step[0].str.contains(i)].values:
                    for line in j:
                        if '.covstat' in line:
                            f.write(line + '\n')
                            f.write('wait ;\n')
                        elif '-L chr10' in line:
                            f.write(line + '\n')
                            f.write('wait ;\n')
                        elif '-L chr20' in line:
                            f.write(line + '\n')
                            f.write('wait ;\n')
                        elif '-L chrM' in line:
                            f.write(line + '\n')
                            f.write('wait ;\n')
                        else:
                            f.write(line + '\n')
                f.write('wait ;')

        with open('job.bat', 'w') as f:
            for i in self.df[5].unique():
                f.write('qsub {}\n'.format(i))
        cmd = 'chmod +x job.bat'
        subprocess.call(cmd, shell=True)

    def _write_second(self):
        with open('genotype', 'w') as f:
            f.write('#$ -N genotype\n')
            f.write('#$ -pe smp 32\n')
            f.write('#$ -q ' + args.queue + '\n')
            f.write('#$ -cwd\n')
            f.write('set -e\n')
            f.write('cd ' + args.output + '\n')
            f.write('source ~/.bash_profile\n')
            for i in self.second_step:
                f.write(i + '\n')
            f.write('wait ;')


if __name__ == '__main__':
    args = get_args()
    sequencing()
