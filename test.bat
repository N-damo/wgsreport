#!/bin/bash
#$ -N NGSjob
#$ -pe smp 32
#$ -q all.q
#$ -cwd
cd /share/data2/leon/run12
source ~/.bash_profile
leeHomMulti -t 3 -fq1 /share/data2/gene/GMB12/S12011_Test_S108_L006_R1_001.fastq.gz -fq2 /share/data2/gene/GMB12/S12011_Test_S108_L006_R2_001.fastq.gz -fqo S12011_Test -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIIATCTCGTATGCCGTCTTCTGCTTG -s AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTIIIIIIIIGTGTAGATCTCGGTGGTCGCCGTATCATT -c ACACTCTTTCCCTACACGACGCTCTTCCGATCT &
wait;
mergeps.o S12011_Test_r1.fq.gz S12011_Test_r2.fq.gz S12011_Test.fq.gz | bwa  mem -p -t24 -R'@RG\tID:GMB12_S12011\tLB:S12011\tSM:S12011_Test\tPU:flowcell\tPL:Illumina' /share/data1/genome/hs38DH.fa - 2>S12011_Test.log.bwamem |  k8 /share/data1/src/bwa/bwakit/bwa-postalt.js -p S12011_Test.hla /share/data1/genome/hs38DH.fa.alt | samtools sort -@ 2 -m8g - -o S12011_Test.aln.bam;
/share/data1/src/gatk/gatk MarkDuplicates  -I S12011_Test.aln.bam -O combined_dedup.bam -M dedup.log -PG null ;
samtools index -@ 4 combined_dedup.bam ;
java -Xmx32g -jar /share/data1/local/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /share/data1/genome/hs38DH.fa -I combined_dedup.bam --known /share/data1/PublicProject/GATK_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz --known /share/data1/PublicProject/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -o forIndelRealigner.intervals -nt 8 ;
java -Xmx32g -jar /share/data1/local/bin/GenomeAnalysisTK.jar -T IndelRealigner -R /share/data1/genome/hs38DH.fa -I combined_dedup.bam -known /share/data1/PublicProject/GATK_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz -known /share/data1/PublicProject/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -targetIntervals forIndelRealigner.intervals -o combined_realign.bam ;
java -Xmx32g -jar /share/data1/local/bin/GenomeAnalysisTK.jar -T DepthOfCoverage -R /share/data1/genome/hs38DH.fa -o coverage -I combined_realign.bam -L /share/data1/PublicProject/GATK_bundle/wgs_calling_regions.hg38.bed -ct 1 -ct 3 -ct 10 -ct 20 -omitBaseOutput &
/share/data1/src/gatk/gatk BaseRecalibrator -I combined_realign.bam -O recal.table -R /share/data1/genome/hs38DH.fa --known-sites /share/data1/PublicProject/GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites /share/data1/PublicProject/GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz -L /share/data1/PublicProject/GATK_bundle/wgs_calling_regions.hg38.bed;
/share/data1/src/gatk/gatk ApplyBQSR -I combined_realign.bam -O combined_bqsr.bam -bqsr recal.table -R /share/data1/genome/hs38DH.fa;
sam flagstatx combined_bqsr.bam > combined_bqsr.flagstatx &
sam covstat combined_bqsr.bam > combined_bqsr.covstat &
samtools view -C -T /share/data1/genome/hs38DH.fa -@ 4 -o GMB12_bqsr.cram combined_bqsr.bam &
java -Xmx32g -jar /share/data1/local/bin/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /share/data1/genome/hs38DH.fa -I combined_bqsr.bam -o multiplex.vcf.gz --dbsnp /share/data1/PublicProject/GATK_bundle/dbsnp150_chr.vcf.gz --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES --alleles /share/data1/PublicProject/BEAGLE/eas.allgrch38.sites_chr.vcf.gz ;
zcat multiplex.vcf.gz | preimput.o | bgzip > multiplex1.vcf.gz ;
tabix multiplex1.vcf.gz ;
mkdir workspace; cd workspace; allsplit.bat ../multiplex1.vcf.gz hmmgl
ls hmmgl*.vcf.gz | shuf | pbsgen.o $PWD impjob > alljobs.bat 
chmod +x alljobs.bat
chmod +x ../asnoerror.bat
wait;
