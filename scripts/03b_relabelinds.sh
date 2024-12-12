#!/bin/bash

sed -i 's/\t/:/g' badgoodsamplist.txt

while read NAME
do
  bcftools view -h $NAME.vcf | grep 'QUE' | sed 's/\t/\n/g' | sed '/QUE/!d' > "$NAME-vcf-samps.txt"
done < filteredvcflist.txt

while read NAME
do
  while read LINE
  do
    grep ^$LINE badgoodsamplist.txt | sed 's/.*://' >> "new-$NAME-vcf-samps.txt"
  done < "$NAME-vcf-samps.txt"
done < filteredvcflist.txt


# bcftools query -l indsrem_DP8.vcf | while read LINE
# do
#   grep ^$LINE ~/qmac/analyses/calling/badgoodsamplist.txt | sed 's/.*://g' >> new-DP8-vcf-samps.txt
# done
#
# bcftools reheader --threads 96 -s new-DP8-vcf-samps.txt -o reheaded_indsrem_DP8.vcf indsrem_DP8.vcf

parallel 'bcftools reheader --threads 5 -s new-{}-vcf-samps.txt {}.vcf > {}-reheaded.vcf' :::: filteredvcflist.txt

while read NAME
do
  awk '{print $1}' ../$NAME.fam | sed 's/QUE//g' | sed 's/_/\t(0)\t/g' | sed 's/mac/1/g;s/alb/2/g;s/bic/3/g;s/ste/4/g;s/mue/5/g' | nl > "new-$NAME-prelimclumpp.txt"
done < ../totalvcflist.txt

while read NAME
do
  sed 's/.*:/:/g' "${NAME}_clumpp.Q" >> "${NAME}_frac_only_clumpp.Q"
done < clumppQlist.txt

while read NAME
do
  for i in {2..10}
  do
    paste "new-$NAME-prelimclumpp.txt" $(echo `echo $NAME | sed 's/filtered-genome-//' | sed 's/-/_/g'`"_K${i}_frac_only_clumpp.Q") > "$(echo `echo $NAME | sed 's/filtered-genome-//' | sed 's/-/_/g'`"_K${i}_updatedinds_clumpp.Q")"
  done
done < ../totalvcflist.txt

while read NAME
do
  zip -r "${NAME}_updatedinds_clumpp.zip" . -i "${NAME}_K"\*"_updatedinds_clumpp.Q"
done < clumppQlistnoK.txt

#
# while read NAME
# do
#   echo `echo $NAME | sed 's/filtered-genome-//' | sed 's/-/_/'`"_clumpp.Q"
# done < ../totalvcflist.txt
#
#
# while read NAME
# do
#   while read $a $b $c $d $e $f
#     do
#       echo "`grep ^$a goodbadsamplist.txt | sed 's/.*://'` `grep ^$a goodbadsamplist.txt | sed 's/.*://'` $c $d $e $f" >> new-$NAME.fam
#   done < $NAME.fam
# done < totalvcflist.txt
#
# while read $a $b $c $d $e $f
#   do
#     echo "`grep ^$a goodbadsamplist.txt | sed 's/.*://'` `grep ^$a goodbadsamplist.txt | sed 's/.*://'` $c $d $e $f" >> new-pruned-filtered-genome-DP8.fam
# done < pruned-filtered-genome-DP8.fam
#
# while read $a
# do
#   echo $a
# done < pruned-filtered-genome-DP8.fam
#
#
#
#
#
# while read NAME
# do
#   a=`awk '{print $1}'`
#   cat $a | while read $samp
#     do
#       echo "`grep ^$samp goodbadsamplist.txt | sed 's/.*://'` `grep ^$samp goodbadsamplist.txt | sed 's/.*://'` `grep $samp $NAME | sed `
#   done
# done < totalvcflist.txt
