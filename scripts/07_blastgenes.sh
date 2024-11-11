#!/bin/bash

# make directory for

mkdir ~/qmac/analyses/blast
cd ~/qmac/analyses/blast

# make new fasta file for each target sequence

csplit -kfs genesequences.fa '/^>/' '{*}'
rm s00

for i in [^gene]*
do
  mv $i "$(head -1 $i | sed 's/>//g')_target.fa"
done

# make namelist of target sequence fasta files

for i in *_target.fa
do
  echo ${i%.fa} >> targetlist.txt
done

# make BLAST databases for Q. lobata and Q. robur genomes - make sure genomes are in directory

makeblastdb -dbtype nucl -in QlobataReference.fna -title "Q. lobata Reference Genome" -logfile Qlobatadatabase.log -taxid 97700
makeblastdb -dbtype nucl -in QroburReference.fa -title "Q. robur Reference Genome" -logfile Qroburdatabase.log -taxid 38942

# BLAST target queries against Q. lobata genome and retain best hit

while read NAME
do
  blastn -query "$NAME-target.fa" -db QlobataReference.fna -num_alignments 1 -outfmt "10 qacc length sseqid evalue bitscore score sstart send sstrand" | head -1 >> prelimbed.txt
done < targetlist.txt

# parallelized version:

# parallel 'blastn -query {} -db QlobataReference.fna -outfmt "10 qacc length sseqid evalue bitscore score sstart send sstrand" | head -1 >> prelimbed.txt' :::: targetsequencelist.txt

# get positions/chr of each gene for grabbing sequence from Q. lobata genome

while IFS="," read a b c d e f g h i
do
  if [[ $g -gt $h ]]
  then
    echo "$c $h $g $a . $i" >> Qlobatagenes.bed
  else
    echo "$c $g $h $a . $i" >> Qlobatagenes.bed
  fi
done < prelimbed.txt
sed -i 's/ /\t/g;s/plus/+/g;s/minus/-/g' Qlobatagenes.bed

# get gene sequences from Q. lobata reference

bedtools getfasta -fi QlobataReference.fna -bed Qlobatagenes.bed -name -fo | csplit -f recip -ks - '/^>/' '{*}'

for i in recip*
do
  mv $i "$(head -1 $i | sed 's/>//;s/\.2.*/.2/')_recip.fa"
done

# reciprocal best BLAST of best Q. lobata result against Q. robur genome

while read NAME
do
  sed -i '1s/::.*//' "${NAME}_recip.fa"
  blastn -query "${NAME}_recip.fa" -db QroburReference.fa -outfmt "10 qacc length sseqid evalue bitscore score sstart send sstrand" | head -1 >> prelimrecip.txt
done < targetlist.txt

# parallelized version:

# parallel 'blastn -query {}_recip.fa -db QroburReference.fa -outfmt "10 qacc length sseqid evalue bitscore score sstart send sstrand" | head -1 >> prelimrecip.txt' :::: targetlist.txt

# check loci that failed RBB

blastn -query T0299670.2-target.fa -db QlobataReference.fna
blastn -query T0200130.2-target.fa -db QlobataReference.fna
blastn -query T0160680.2-target.fa -db QlobataReference.fna

blastn -query T0299670.2-target.fa -db QroburReference.fa
blastn -query T0200130.2-target.fa -db QroburReference.fa
blastn -query T0160680.2-target.fa -db QroburReference.fa

sed 's/^Qrob_//g;1d;s/Chr/Qrob_Chr/g' quercus_hybseq_genes.txt | awk '{print $3,$4,$5,$1}' - > Qroburreglist.txt

while read a b c d
do
  if [[ $b -gt $c ]]
  then
    if [[ $(expr length "$a") == 9 ]]
    then
      echo `echo $a | sed "s/$a/${a%[0-9]}0${a#*Chr}/"`$'\t'$c$'\t'$b$'\t'$d >> Qroburreglist.bed
    else
      echo $a$'\t'$c$'\t'$b$'\t'$d >> Qroburreglist.bed
    fi
    echo
  else
    if [[ $(expr length "$a") == 9 ]]
    then
      echo `echo $a | sed "s/$a/${a%[0-9]}0${a#*Chr}/"`$'\t'$b$'\t'$c$'\t'$d >> Qroburreglist.bed
    else
      echo $a$'\t'$b$'\t'$c$'\t'$d >> Qroburreglist.bed
    fi
  fi
done < Qroburreglist.txt

bedtools getfasta -fi QroburReference.fa -bed Qroburreglist.bed -name | csplit -f "fullgene" -ks - '/^>/' '{*}'

while read NAME
do
  if [[ $(blastn -query "$NAME-target.fa" -db "$NAME-full.fa") == 1 ]]
  then

  else

  fi
done < targetlist.txt

# v3.2 from NCBI

grep '>' QlobataReference.fna | while read LINE
do
  echo $LINE | sed 's/>//;s/Que.*SW786//;s/un.*Assembly//;s/,.*//;s/\s\+/:/;s/chromosome/chr/;s/:.*chr2_un/:chr2_un/;s/\s//' >> Qlobataregconvlist.txt
done

cp Qlobata.v3.0.PCG.COPY.bed Qlobata.v3.0.PCG.bed
while read LINE
do
  sed -i "s/${LINE#*:}\s/${LINE%:*}\t/g" Qlobata.v3.0.PCG.bed
done < Qlobataregconvlist.txt

awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6}' Qlobata.v3.0.PCG.bed > QlobataPCGreduced.bed

while read LINE
do
  sed "s/${LINE#*:}$/${LINE%:*}/g" QlobataRef3.0.fa
done < Qlobataregconvlist.txt

bedtools getfasta -s -name -fi QlobataRef3.0.fa -bed QlobataPCGreduced.bed > QlobataPCG.fa

# Scq3eQI_473, Scq3eQI_539, Scq3eQI_790, Scq3eQI_853, Scq3eQI_943, Scq3eQI_1014 not found in v3.2 ref genome

makeblastdb -dbtype nucl -in QlobataPCG.fa -title "Q. lobata v3.0 PCG" \
-logfile QlobataPCGdatabase.log -taxid 97700

while read NAME
do
  blastn -query "$NAME-target.fa" -db QlobataPCG.fa -outfmt "10 qacc length sseqid evalue bitscore score sstart send sstrand" | head -1 >> prelimbedPCG.txt
done < targetlist.txt

awk -F',' '{print $1}' prelimbedPCG.txt > targetlist-PCGblasted.txt

grep -v -f targetlist-PCGblasted.txt targetlist.txt | while read NAME
do
  blastn -query "$NAME-target.fa" -db QlobataRef3.0.fa -outfmt "10 qacc length sseqid evalue bitscore score sstart send sstrand" | head -1 >> missing-genes.txt
done

# T0078670.2,1386,NC_044905.1,0.0,2499,1353,46875137,46873752,minus
# T0236590.2,343,NC_044906.1,6.66e-59,233,126,47834876,47834571,minus
# T0654080.2,366,NC_044907.1,0.0,660,357,19284257,19283892,minus
# T0716070.2,312,NC_044911.1,4.46e-159,566,306,62598905,62599216,plus

while read NAME
do
  grep "$(echo "$NAME" | awk -F',' '{print $3}')" QlobataPCGreduced.bed | grep "\s$(echo "$NAME" | awk -F',' '{print $7}' | cut -c1-2 )" - >> potential-missing-genes.txt
done < missing-genes.txt
awk -F',' '{OFS="\t";if($7 < $8){print $3,$7,$8,$1}else{print $3,$8,$7,$1}}' missing-genes.txt | cat potential-missing-genes.txt - | sort -k1,1 -k2,2 -

# manually examine unmapped targets and loci

awk -F',' '{OFS=","; print $1,$3}' prelimbedPCG.txt | sed 's/(.*//g;s/:\+/,/g;s/-/,/g' - | awk -F',' '{OFS="\t";print $3,$4,$5,$1"_"$2}' - > PCG-3.0.bed
awk '{print $4}' PCG-3.0.bed > QrobtoQlobconvlist.txt

bedtools getfasta -fi ../QlobataRef3.0.fa -bed ../PCG-3.0.bed -nameOnly | csplit -z -k -s -n 3 -b "%02d.fa" -f "PCG-3.0-" - '/^>/' '{*}'
for i in PCG-3.0-*.fa
do
  mv "$i" PCG-3.0-"$(sed -n '1p' "$i" | sed 's/^>//')".fa
done

makeblastdb -dbtype nucl -in ~/qmac/analyses/align/index/QlobataReference.fna -title "Q. lobata v3.2 Reference Genome" \
-logfile Qlobata3.2database.log -taxid 97700

while read NAME
do
  blastn -query "PCG-3.0-$NAME.fa" -db ~/qmac/analyses/align/index/QlobataReference.fna -outfmt "10 qacc length sseqid evalue bitscore score sstart send sstrand" | head -1 >> prelimbedPCG3.2.txt
done < QrobtoQlobconvlist.txt

awk -F',' '{OFS="\t"; print $3,$7,$8,$1,".",":"}' prelimbedPCG3.2.txt | while read LINE
do
  echo "$LINE" | sed "s/:/$(grep "$(echo "$LINE" | cut -f4 | sed 's/.*\.2_//')" Qlobata.v3.0.PCG.bed | cut -f6 -)/" >> sorted-Qlobatagenes-PCG.bed
done

awk '{print $2}' ../subsample/DP8.bim | awk -F':' '{print $1,$2,$2}' > DP8-variants.bed

for i in {1..20}
do
  sort -R QrobtoQlobconvlist.txt | head -15 | while read LINE
  do
    printf "$(grep $(echo "$LINE" | sed 's/\.2_.*//') Qroburgenenames.txt)\n$(grep $(echo "$LINE" | sed 's/.*\.2_//') Qlobatagenenames.txt)\n\n"
  done
done

# manually check if random subset of gene names match

while read NAME
do
  grep "$NAME" sorted-Qlobatagenes-PCG.bed | sort -n -k2 - > "$NAME-PCG.bed"
  bedtools closest -d -N -io -a "$NAME-PCG.bed" -b "$NAME-PCG.bed" | awk '{print $4,$10,$13}' - >> genedistprelim.txt
done < chrlist.txt

sort -n -k3 genedistprelim.txt > sorted-genedist.txt

# manually inspect via IGV to determine appropriate borders
