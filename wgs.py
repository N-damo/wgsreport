#!/usr/bin/env python3
# coding:utf-8

import os
import time
import sys
import subprocess
from collections import defaultdict
import argparse
from argparse import ArgumentParser
import shutil


def get_args():
    parser = ArgumentParser(
        description='wgs analysis workflow made by guomai biotechonology')
    parser.add_argument('--input', '-i', help='provide sample info txt',
                        required=True, type=extant_file, metavar='FILE')
    parser.add_argument(
        '--run', '-r', choices={'on','off'},help='if you choose run on,please add -r argument.the default is off', required=False, default='off')
    args = parser.parse_args()
    return args


def extant_file(file):
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError(
            '{file} dose not exist'.format(file=file))
    return file


class Wgs(object):

    trimmomatic = "/share/data1/local/bin/trimmomatic-0.38.jar"
    gatk = "/share/data1/src/gatk/gatk"
    gatk3_8 = "/share/data1/local/bin/GenomeAnalysisTK.jar"
    seqtk = "/share/data1/src/bwa/bwakit/seqtk"
    alt = "/share/data1/genome/hs38DH.fa.alt"
    reference = "/share/data1/genome/hs38DH.fa"
    bed = "/share/data1/PublicProject/GATK_bundle/wgs_calling_regions.hg38.bed"
    dbsnp = "/share/data1/PublicProject/GATK_bundle/dbsnp150_chr.vcf.gz"
    indel = "/share/data1/PublicProject/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    indel2 = "/share/data1/PublicProject/GATK_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz"
    snp = "/share/data1/PublicProject/GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    alleles = "/share/data1/PublicProject/BEAGLE/eas.allgrch38.sites_chr.vcf.gz"

    def __init__(self, input):
        self.input = input
        self.group()

    def create_working_space(self):
        try:
            if not os.path.exists(self.working_space):
                os.mkdir(self.working_space)
                sys.stdout.write('mkdir {}  directory\n'.format(
                    self.working_space))
            else:
                os.removedirs(self.working_space)
                os.mkdir(self.working_space)
                sys.stdout.write('delete this empty directory {}\n'.format(
                    self.working_space))
        except OSError:
            sys.stdout.write('mkdir {}  directory.\n'.format(
                    self.working_space))
            sys.stdout.write('directory is not empty,make sure first:{}.\nexit!\n'.format(self.working_space))
            sys.exit()
            # n=1
            # while n <= 3:
            #     command = input(
            #     'delete this none-empty directory?: {}\nyou can choose yes/no:'.format(self.working_space))
            #     command=command.lower()
            #     if command == 'yes' or command == 'y':
            #         sys.stdout.write(
            #             'delete this none-empty directory: {}'.format(self.working_space))
            #         shutil.rmtree(self.working_space)
            #         os.mkdir(self.working_space)
            #         sys.stdout.write('mkdir {}  directory'.format(
            #             self.working_space))
            #         break
            #     elif command == 'no' or command == 'n':
            #         sys.stdout.write(
            #             'remain this none-empty directory.All bat job will be running in this directory:{}'.format(self.working_space))
            #         break
            #     else:
            #         sys.stdout.write('bad input,plese retry yes/no')
            #         n += 1
            # else:
            #     sys.stdout.write('please try later')
            #     sys.exit()


    def group(self):
        self.groups = defaultdict(list)
        self.working_space = ''
        row_line = 0
        with open(self.input, 'r') as f:
            for i in f:
                if row_line > 0:
                    row_line += 1
                    line = i.strip().split()
                    if len(line) == 5:
                        sample = line[4]
                        self.groups[sample].append(line)
                    else:
                        sys.stdout.write('it seem has less or more than five items in the {} row,it may be a blank row,so pass it and continue\n'.format(
                            row_line))
                else:
                    self.working_space = i.strip().split()
                    assert len(
                        self.working_space) == 1, 'please use batch number,e.g. gmb1,gmb2 ...'
                    self.working_space = self.working_space[0]
                    # if not os.path.exists(self.working_space):
                    #     os.mkdir(self.working_space)
                    self.create_working_space()
                    self.working_space = os.path.abspath(self.working_space)
                    row_line += 1

    def info_list_by_sample(self):
        for sample in self.groups:
            info_list = self.groups[sample]
            # each sample may contain more than one library id,these sample will be merge in markduplicate step.
            yield info_list

    def trim(self):
        for info_list in self.info_list_by_sample():
            part = []
            for parallel_sample in info_list:
                fq1, fq2, adapter, library, sample = parallel_sample
                # only accept nextera or truseq adapter
                if adapter.lower() == 'truseq':
                    trim = 'java -jar {tool} PE -threads 5 -phred33 {fq1} {fq2} {library}.{sample}_r1.fq.gz {library}.{sample}_r1_unpaired.fq.gz {library}.{sample}_r2.fq.gz {library}.{sample}_r2_unpaired.fq.gz ILLUMINACLIP:{trimmomatic_dir}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &'.format(
                        tool=self.trimmomatic, trimmomatic_dir=os.path.dirname(self.trimmomatic), fq1=fq1, fq2=fq2, library=library, sample=sample)
                elif adapter.lower() == 'nextera':
                    trim = 'java -jar {tool} PE -threads 5 -phred33 {fq1} {fq2} {library}.{sample}_r1.fq.gz {library}.{sample}_r1_unpaired.fq.gz {library}.{sample}_r2.fq.gz {library}.{sample}_r2_unpaired.fq.gz ILLUMINACLIP:{trimmomatic_dir}/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &'.format(
                        tool=self.trimmomatic, trimmomatic_dir=os.path.dirname(self.trimmomatic), fq1=fq1, fq2=fq2, library=library, sample=sample)
                else:
                    sys.stdout.write('{} :\nadaptor error, only Nextera or Truseq accept,lower or upper dose not matter\n'.format(
                        adapter))
                    sys.exit()
                part.append(trim)
            part.append('wait ;')
            yield part, info_list

    def bwa(self):
        for part, info_list in self.trim():
            for parallel_sample in info_list:
                *_, library, sample = parallel_sample
                bwa = 'mergeps.o {library}.{sample}_r1.fq.gz {library}.{sample}_r2.fq.gz |bwa mem -p -t 5 -R"@RG\\tID:{batch_library}\\tLB:{library}\\tSM:{sample}\\tPU:flowcell\\tPL:Illumina" {reference} - 2>{library}.{sample}.log.bwamem |k8 {bwa_dir}/bwa-postalt.js -p {library}.hla {alt} | samtools sort -@ 2 -m8g - -o {library}.{sample}.aln.bam ;'.format(
                    reference=self.reference, bwa_dir=os.path.dirname(self.seqtk), alt=self.alt, library=library, batch_library=os.path.basename(self.working_space)+'_'+library, sample=sample)
                part.append(bwa)
            yield part, info_list

    def MarkDup(self):
        for part, info_list in self.bwa():
            bams = []
            for parallel_sample in info_list:
                *_, library, sample = parallel_sample
                bams.append("{library}.{sample}.aln.bam".format(
                    library=library, sample=sample))
            MarkDuplicates = '{tool} MarkDuplicates -I {bams}  -O {sample}.dedup.bam -M {sample}.dedup.log -PG null --TMP_DIR ~/tmp/{sample} ;'.format(
                tool=self.gatk, bams=' -I '.join(bams), sample=bams[0].split('.')[1])
            part.append(MarkDuplicates)
            yield part, info_list

    def trim_bwa_dedup_indel_realign(self):
        for part, info_list in self.MarkDup():
            *_, sample = info_list[0]
            index = 'samtools index -@ 4 {sample}.dedup.bam ;'.format(
                sample=sample)
            RealignerTargetCreator = "java -Xmx32g -jar {tool} -T RealignerTargetCreator -R {reference} -I {sample}.dedup.bam --known {indel} --known {indel2} -o {sample}.forIndelRealigner.intervals -nt 5 ;".format(
                tool=self.gatk3_8, sample=sample, reference=self.reference, indel=self.indel, indel2=self.indel2)
            IndelRealigner = "java -Xmx32g -jar {tool} -T IndelRealigner -R {reference} -I {sample}.dedup.bam -known {indel} -known {indel2} -targetIntervals {sample}.forIndelRealigner.intervals -o {sample}._realign.bam ;".format(
                tool=self.gatk3_8, sample=sample, reference=self.reference, indel=self.indel, indel2=self.indel2)
            index2 = 'samtools index -@ 4 {sample}._realign.bam ;'.format(
                sample=sample)
            part.append(index)
            part.append(RealignerTargetCreator)
            part.append(IndelRealigner)
            part.append(index2)
            sge = "#$ -N {sample}\n#$ -pe smp 6\n#$ -q all.q\n#$ -cwd\nset -e\ncd {working_space}\nsource ~/.bash_profile\n".format(
                sample=sample, working_space=self.working_space)
            with open(os.path.join(self.working_space, sample+'.bat'), 'w') as f:
                f.write(sge)
                f.write('\n'.join(part)+'\n')
                f.write('wait ;')

    def genotype(self):
        bams = ['{sample}._realign.bam'.format(
            sample=sample) for sample in self.groups]
        BaseRecalibrator = '{tool} BaseRecalibrator -I {bams} -O recal.table --known-sites {snp} --known-sites {indel} -L {bed}  -R {reference} ;'.format(
            tool=self.gatk, bed=self.bed, snp=self.snp, indel=self.indel, reference=self.reference, bams=' -I '.join(bams))
        ApplyBQSR = '{tool} ApplyBQSR -I {bams} -O combined_bqsr.bam -bqsr recal.table -R {reference} ;'.format(
            tool=self.gatk, reference=self.reference, bams=' -I '.join(bams))
        # DepthOfCoverage = "java -Xmx32g -jar {tool} -T DepthOfCoverage -R {reference} -o coverage -I combined_bqsr.bam -L {bed} -ct 1 -ct 3 -ct 10 -ct 20 -omitBaseOutput &".format(
        #     tool=self.gatk3_8, bed=self.bed, reference=self.reference)
        flagstatx = "sam flagstatx combined_bqsr.bam > combined_bqsr.flagstatx &"
        covstat = "sam covstat combined_bqsr.bam > combined_bqsr.covstat &"
        cram = "samtools view -C -T {reference} -@ 4 -o {batch}_combined_bqsr.cram combined_bqsr.bam &".format(
            reference=self.reference, batch=os.path.basename(self.working_space))
        UnifiedGenotyper = "java -Xmx32g -jar {tool} -T UnifiedGenotyper -R {reference} -I combined_bqsr.bam -o multiplex.vcf.gz --dbsnp {dbsnp} --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES --alleles {allels} ;".format(
            tool=self.gatk3_8, reference=self.reference, dbsnp=self.dbsnp, allels=self.alleles)
        preimput = "zcat multiplex.vcf.gz | preimput.o | bgzip > multiplex1.vcf.gz ;"
        tabix = "tabix multiplex1.vcf.gz ;"
        allsplit = "mkdir workspace ;cd {workspace}; allsplit.bat ../multiplex1.vcf.gz hmmgl".format(
            workspace=os.path.join(self.working_space, 'workspace'))
        alljobs = "ls hmmgl*.vcf.gz | shuf | pbsgen.o $PWD impjob > alljobs.bat"
        sge = "#$ -N {batch}\n#$ -pe smp 32\n#$ -q all.q\n#$ -cwd\nset -e\ncd {working_space}\nsource ~/.bash_profile\n".format(
            working_space=self.working_space, batch=os.path.basename(self.working_space))
        # part = [BaseRecalibrator, ApplyBQSR, DepthOfCoverage, flagstatx,
        #         covstat, cram, UnifiedGenotyper, preimput, tabix, allsplit, alljobs]
        part = [BaseRecalibrator, ApplyBQSR, flagstatx,
                covstat, cram, UnifiedGenotyper, preimput, tabix, allsplit, alljobs]
        with open(os.path.join(self.working_space, 'genotype.bat'), 'w') as f:
            f.write(sge)
            f.write('\n'.join(part)+'\n')
            f.write('wait ;')
        #subprocess.call('chmod +x genotype.bat', shell=True)

    def ssh_check(self):
        p = subprocess.Popen('hostname', shell=True, stdout=subprocess.PIPE)
        hostname = p.communicate()[0].strip().decode()
        return hostname

    def run_qsub_step1_step2(self):
        assert self.ssh_check() == 'guomai1', 'please login to guomai1 firstly'
        job_ID = []
        os.chdir(self.working_space)
        sys.stdout.write('changing working directory to your provide {}\n'.format(
            self.working_space))
        if os.path.exists('gm.log'):
            os.remove('gm.log')
        subprocess.call(
            'date +"%Y-%m-%d %H:%M:%S" >> gm.log', shell=True)
        subprocess.call('echo begin step1>> gm.log', shell=True)
        for sample in self.groups:
            cmd = 'qsub {sample}.bat ;'.format(
                sample=sample)
            print(cmd)
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            step1_ID = p.communicate()[0].split()[2]
            job_ID.append(step1_ID)
        while 1:
            qstat_ID = self.run_check_by_time()
            _ = set(job_ID) & set(qstat_ID)
            if len(_) == 0:
                subprocess.call(
                    'date +"%Y-%m-%d %H:%M:%S" >> gm.log', shell=True)
                subprocess.call('echo end step1>> gm.log', shell=True)
                cmd = 'qsub genotype.bat'
                sys.stdout.write(cmd)
                p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                step2_ID = p.communicate()[0].split()[2]
                subprocess.call(
                    'date +"%Y-%m-%d %H:%M:%S" >> gm.log', shell=True)
                subprocess.call('echo begin genotyping>> gm.log', shell=True)
                break
            else:
                time.sleep(60)
        return step2_ID

    def run_pbs_step3(self, step2_ID):
        while 1:
            qstat_ID = self.run_check_by_time()
            if step2_ID not in qstat_ID:
                subprocess.call(
                    'date +"%Y-%m-%d %H:%M:%S" >> gm.log', shell=True)
                subprocess.call('echo end genotyping>> gm.log', shell=True)
                files = os.listdir('workspace')
                pbs = [i for i in files if i.endswith('.pbs')]
                pbs_job_ID = []
                for sample in pbs:
                    cmd = 'qsub {sample} ;'.format(sample=os.path.join('workspace',sample))
                    p = subprocess.Popen(
                        cmd, shell=True, stdout=subprocess.PIPE)
                    ID = p.communicate()[0].split()[2]
                    pbs_job_ID.append(ID)
                    print(cmd)
                subprocess.call(
                    'date +"%Y-%m-%d %H:%M:%S" >> gm.log', shell=True)
                subprocess.call('echo begin imputation>> gm.log', shell=True)
                break
            else:
                time.sleep(60)
        return pbs_job_ID

    def run_merge_step4(self, pbs_job_ID):
        while 1:
            qstat_ID = self.run_check_by_time()
            _ = set(pbs_job_ID) & set(qstat_ID)
            if len(_) == 0:
                subprocess.call(
                    'date +"%Y-%m-%d %H:%M:%S" >> gm.log', shell=True)
                subprocess.call(
                    'echo completed imputation and merge vcf,genelize final.vcf.gz locally>> gm.log', shell=True)
                subprocess.call('./asnoerror.bat', shell=True)
                # add mtDNA haplogroup prediction
                break
            else:
                time.sleep(60)

    def run_check_by_time(self):
        cmd = "qstat|sed -n '3,$p'|cut -d ' ' -f3"
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        qstat_ID = p.communicate()[0].split(b'\n')
        return qstat_ID


def asnoerror(working_space):
    working_space = os.path.basename(working_space)
    with open(os.path.join(working_space, 'asnoerror.bat'), 'w') as f:
        f.write("#!/bin/bash"+'\n')
        f.write('cd workspace'+'\n')
        f.write('gatherpostname.bat mid_hmmgl | postimpbatch.o final | sh'+'\n')
        #f.write('vcftools --gzvcf final_all.vcf.gz --recode --snps /share/data1/tmp/gwas-rsid.txt -c |bgzip > {batch}_report.vcf.gz'.format(batch=working_space)+'\n')
        f.write('mv final_all.vcf.gz {batch}.vcf.gz'.format(
            batch=working_space)+'\n')


if __name__ == "__main__":
    args = get_args()
    job = args.input
    result = Wgs(job)
    asnoerror(result.working_space)
    result.trim_bwa_dedup_indel_realign()
    result.genotype()
    if args.run == 'wet_run':
        sys.stdout.write('auto qsub will be running in guomai1 right now.\n')
        step2_ID = result.run_qsub_step1_step2()
        pbs_job_ID = result.run_pbs_step3(step2_ID)
        result.run_merge_step4(pbs_job_ID)
    elif args.run == 'dry_run':
        sys.stdout.write('cause you choose dry_run,all job bat will be in {}.\n'.format(result.working_space))
    else:
        sys.stdout.write('only support dry_run or wet_run.The defatult is dry_run.\n')
        sys.exit()
