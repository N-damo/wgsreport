#!/usr/bin/bash
#plink2 --vcf ~/soft/bin/wgs/b38/1kg.phase3.v5a.norm_snps.vcf.gz --biallelic-only --make-bed --out 1kg --set-missing-var-ids @:#:\$1:\$2
#plink2 --vcf /share/data4/deposit/VCF/gmb25.vcf.gz --biallelic-only --make-bed --out gmb25 --maf 0.05 --set-missing-var-ids @:#:\$1:\$2

if [ $# != 1 ]
then
echo the only argument is vcf file route
echo sh prepare.sh /path/to/gmb1.vcf.gz
exit
elif [ $1 == -h ]
then
echo usage:sh prepare.sh /path/to/something.vcf.gz
exit
elif [ $1 == --help ]
then
exit
fi


batch=$(basename $1|cut -d '.' -f1)
plink2 --vcf $1 --biallelic-only --make-bed --out ${batch} --set-missing-var-ids @:#:\$1:\$2
plink2 --bfile /share/data3/lianlin/admixture/1kg --bmerge ${batch} --make-bed --out merge_test
plink2 --bfile /share/data3/lianlin/admixture/1kg --exclude merge_test-merge.missnp --make-bed  --out 1kg_subset 
plink2 --bfile ${batch} --exclude merge_test-merge.missnp --make-bed  --out ${batch}_subset 
plink2 --bfile 1kg_subset --bmerge ${batch}_subset --make-bed --out panel.1kg
plink2 --bfile panel.1kg  --indep-pairwise 50 10 0.1 
plink2 --bfile panel.1kg --extract plink.prune.in --make-bed --out panel_prune.1kg
python3 /share/data3/lianlin/admixture/admix.py /share/data3/lianlin/admixture/1kgenome.txt panel_prune.1kg.fam panel_prune.1kg
cmd="#$ -N genotye\n#$ -pe smp 20\n#$ -q all.q\n#$ -cwd\nset -e\ncd $PWD\nsource ~/.bash_profile\n" 
echo -e $cmd >sge
for i in {1..22};do plink2 --bfile panel_prune.1kg --chr $i --make-bed --out panel_prune.1kg.chr$i && cp panel_prune.1kg.pop panel_prune.1kg.chr$i.pop;done ;
for i in {1..22};do (cat sge; echo admixture  --supervised -j20 panel_prune.1kg.chr$i.bed 26) >chr$i.pbs;done ;
ls *.pbs|while read i;do qsub $i;done ;

