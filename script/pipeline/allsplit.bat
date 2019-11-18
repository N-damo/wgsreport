#!/bin/bash
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 1 ) | java -jar splitvcf.jar 1 50000 3000 $2"0" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 2 ) | java -jar splitvcf.jar 2 50000 3000 $2"1" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 3 ) | java -jar splitvcf.jar 3 50000 3000 $2"2" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 4 ) | java -jar splitvcf.jar 4 50000 3000 $2"3" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 5 ) | java -jar splitvcf.jar 5 50000 3000 $2"4" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 6 ) | java -jar splitvcf.jar 6 50000 3000 $2"5" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 7 ) | java -jar splitvcf.jar 7 50000 3000 $2"6" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 8 ) | java -jar splitvcf.jar 8 50000 3000 $2"7" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 9 ) | java -jar splitvcf.jar 9 50000 3000 $2"8" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 10 ) | java -jar splitvcf.jar 10 50000 3000 $2"9" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 11 ) | java -jar splitvcf.jar 11 50000 3000 $2"10" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 12 ) | java -jar splitvcf.jar 12 50000 3000 $2"11" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 13 ) | java -jar splitvcf.jar 13 50000 3000 $2"12" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 14 ) | java -jar splitvcf.jar 14 50000 3000 $2"13" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 15 ) | java -jar splitvcf.jar 15 50000 3000 $2"14" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 16 ) | java -jar splitvcf.jar 16 50000 3000 $2"15" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 17 ) | java -jar splitvcf.jar 17 50000 3000 $2"16" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 18 ) | java -jar splitvcf.jar 18 50000 3000 $2"17" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 19 ) | java -jar splitvcf.jar 19 50000 3000 $2"18" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 20 ) | java -jar splitvcf.jar 20 50000 3000 $2"19" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 21 ) | java -jar splitvcf.jar 21 50000 3000 $2"20" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 22 ) | java -jar splitvcf.jar 22 50000 3000 $2"21" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 X:1-2781479 ) | java -jar splitvcf.jar X:1-2781479 50000 3000 $2"22" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 X:2781480-155701382 ) | java -jar splitvcf.jar X:2781480-155701382 50000 3000 $2"23" &
( cat vcfheader ; zcat $1 | head -5000|grep '^#CHROM'; tabix $1 X:155701383-156040895 ) | java -jar splitvcf.jar X:155701383-156040895 50000 3000 $2"24" &
wait