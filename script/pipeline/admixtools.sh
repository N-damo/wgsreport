#!/usr/bin/bash

source ~/.bash_profile
source activate gatk



if [ $# != 1 ]
then
echo the only argument is vcf file route
echo sh admixtools.sh /path/to/gmb1.vcf.gz
exit
elif [ $1 == -h ]
then
echo usage:sh admixtools.sh /path/to/something.vcf.gz
exit
elif [ $1 == --help ]
then
exit
fi


batch=$(basename $1|cut -d '.' -f1)
plink2 --vcf /share/data3/lianlin/chimp/v37.2.1240K/panel_filter.vcf.gz --make-bed --out panel --set-missing-var-ids @:#:\$1:\$2
cut -f2 panel.bim >panel.snp
plink2 --vcf $1 --extract panel.snp --make-bed --out $batch --set-missing-var-ids @:#:\$1:\$2
plink2 --bfile panel --bmerge $batch --make-bed --out merge_test
plink2 --bfile panel --exclude merge_test-merge.missnp --make-bed --out panel_drop
plink2 --bfile $batch --exclude merge_test-merge.missnp --make-bed --out ${batch}_drop
plink2 --bfile panel_drop --bmerge ${batch}_drop --make-bed --recode --out test
/share/data1/src/EIG-6.1.4/bin/convertf -p /share/data3/lianlin/chimp/v37.2.1240K/par.txt
awk '{print $1,"U",$1}' test.fam >test.ind
R CMD BATCH /share/data3/lianlin/chimp/v37.2.1240K/qpF4ratio.R
