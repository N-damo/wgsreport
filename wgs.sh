#!/bin/bash
#coding:utf-8
set -e

trimmomatic="/share/data1/local/bin/trimmomatic-0.38.jar"
gatk="/share/data1/src/gatk/gatk"
seqtk="/share/data1/src/bwa/bwakit/seqtk"
annovar="/share/data1/src/annovar/"
dbsnp="/share/data1/PublicProject/GATK_bundle/dbsnp150_chr.vcf.gz"
indel="/share/data1/PublicProject/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
snp="/share/data1/PublicProject/GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
reference="/share/data1/genome/hs38DH.fa"
eas_allele="/share/data1/PublicProject/BEAGLE/eas.allgrch38.sites_chr.vcf.gz"
imputation="/share/data3/lianlin/soft/bin/wgs/imputation.py"
postimputation="/share/data3/lianlin/soft/bin/wgs/postimputation.py"

if [ $1 == -h ]
then
echo 'usage:'
echo -e 'prepare input file:\nfq1\tfq2\tadapter\tlibraryID\tsampleID\nfastq file must contain absolute path.\nAdapter only accept Truseq or Nextera.The file should be splited by tab,not space.'
echo 'step1: cat input|exom.sh $PWD output_prefix'
echo 'step2: ./qsub.bat,the job will be qsubed into sge cluster.You need to ssh guomai1 before running step2.'
echo 'step3: when the end of the genotyping jobs, run "nohup ./genotype.bat &"'
exit
elif [ $# != 2 ]
then
echo -e 'argument number error,please check it\nplease type exom.sh -h'
exit
fi

final=$1
pwd=$2


i=0
while read line
do
    f1=$(echo ${line}|awk '{print $1}')
    f2=$(echo ${line}|awk '{print $2}')
    adapter=$(echo ${line}|awk '{print $3}')
    library=$(echo ${line}|awk '{print $4}')
    run=$library
    sample=$(echo ${line}|awk '{print $5}')

    if [ $adapter = 'Truseq' ]
    then
        trim='java -jar '$trimmomatic' PE -threads 4 -phred33 '$f1' '$f2' '${run}'.'${sample}'_r1.fq.gz '${run}'.'${sample}'_r1_unpaired.fq.gz '${run}'.'${sample}'_r2.fq.gz '${run}'.'${sample}'_r2_unpaired.fq.gz ILLUMINACLIP:/share/data1/local/bin/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;'
    elif [ $adapter = 'Nextera' ]
    then
        trim='java -jar '$trimmomatic' PE -threads 4 -phred33 '$f1' '$f2' '${run}'.'${sample}'_r1.fq.gz '${run}'.'${sample}'_r1_unpaired.fq.gz '${run}'.'${sample}'_r2.fq.gz '${run}'.'${sample}'_r2_unpaired.fq.gz ILLUMINACLIP:/share/data1/local/bin/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;'
    else
        echo $adapter
        echo ''${adapter}' error,adapter is not Truseq or Nextera'
        exit
    fi
    samples[$i]=${sample}
    i=`expr $i + 1`
    bwa=''${seqtk}' mergepe '${run}.${sample}'_r1.fq.gz '${run}.${sample}'_r2.fq.gz |bwa mem -p -t 4 -R "@RG\\tID:'${run}'\\tLB:'${library}'\\tSM:'${sample}'\\tPU:flowcell\\tPL:Illumina" '${reference}' - 2>'${run}'.log.bwamem |k8 /share/data1/src/bwa/bwakit/bwa-postalt.js -p '${run}'.hla '${reference}'.alt | samtools sort -@ 2 -m8g - -o '${run}.${sample}'.aln.bam ;'
    markduplicate=''${gatk}' MarkDuplicates -I '${run}'.'${sample}'.aln.bam -O '${sample}'.dedup.bam -M '${sample}'.dedup.log -PG null --TMP_DIR ~/tmp/'${sample}' ;'
    index='samtools index -@ 4 '${sample}'.dedup.bam ;'
    CollectWgsMetrics=''${gatk}' CollectWgsMetrics -I '${sample}'.dedup.bam -O '${sample}'_hs_metrics.txt -R '${reference}' &'
    BaseRecalibrator=''${gatk}' BaseRecalibrator -I '${sample}'.dedup.bam -O '${sample}'.recal.table --known-sites '$indel' --known-sites '$snp'  -R '${reference}' ;'
    ApplyBQSR=''${gatk}' ApplyBQSR -I '${sample}'.dedup.bam -O '${sample}'.dedup.bqsr.bam -bqsr '${sample}'.recal.table -R '${reference}' ;'
    HaplotypeCaller=''${gatk}' HaplotypeCaller  --genotyping-mode GENOTYPE_GIVEN_ALLELES --output-mode EMIT_ALL_SITES -I '${sample}'.dedup.bqsr.bam -O '${sample}'.HC.vcf.gz --dbsnp '$dbsnp'  -R '${reference}' --alleles '$eas_allele' ;'
    sge="#$ -N ${sample}\n#$ -pe smp 32\n#$ -q all.q\n#$ -cwd\nset -e\ncd $pwd\nsource ~/.bash_profile\n"
    echo -e "$sge\n${trim}\n$bwa\n$markduplicate\n$index\n$CollectWgsMetrics\n$BaseRecalibrator\n$ApplyBQSR\n$HaplotypeCaller\nwait ;" > ${sample}.bat

done

if [ -f 'genotype.bat' ]
then
    rm 'genotype.bat'
fi

if [ $i == 1 ]
then
    mergebam=''
    cram='samtools view -C -T '${reference}' -@ 4 -o '${final}'_bqsr.cram '${sample}'.dedup.bqsr.bam &' 

else
    mergebam='samtools merge -f -@ 32 -O bam '${final}_bqsr.bam' '$(echo ${samples[@]}|sed 's/\s/.dedup.bqsr.bam /g')'.dedup.bqsr.bam ;'
    cram='samtools view -C -T '${reference}' -@ 4 -o '${final}'_bqsr.cram '${final}'_bqsr.bam &' 
fi
sge="#$ -N imputation\n#$ -pe smp 32\n#$ -q all.q\n#$ -cwd\nset -e\ncd $pwd\nsource ~/.bash_profile\n"
MergeVcfs=''${gatk}' MergeVcfs -R '${reference}' -I '$(echo ${samples[@]}|sed 's/\s/.HC.vcf.gz /g'|sed 's/\s/ -V /g')'.HC.vcf.gz -O '${final}'.vcf.gz ;'

beagle_work="$pwd/beagle"

if [ -d $beagle_work ]
then
    rm -rf $beagle_work
fi
    


beagle='mkdir '$beagle_work' ;\ncd '$beagle_work' ;\n'$imputation' ../'${final}'.vcf.gz ;'
echo -e "$sge\n$mergebam\n$cram\n$MergeVcfs\n$beagle\nwait ;"|sed '/^$/d' > genotype.bat
echo -e 'cd '$beagle_work' ;\n'$postimputation' '$final' ;\nwait ;' >asnoerror.bat
chmod +x asnoerror.bat
chmod +x genotype.bat


if [ -f 'qsub.bat' ]
then
    rm qsub.bat
fi

echo ${samples[@]}|tr ' ' '\n'|sed 's/^/qsub /'|sed 's/$/.bat/' > qsub.bat
chmod +x qsub.bat






