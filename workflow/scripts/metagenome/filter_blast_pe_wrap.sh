#!/usr/bin/env bash

R1=$1
R2=$2
threads=$3
out=$4
taxid=$5
SGREP=$6

rm -f $R1.*
rm -f $R2.*

LC_ALL=C
grep -v "GI_NOT_FOUND" $R1 | sort > $R1.sorted
grep -v "GI_NOT_FOUND" $R2 | sort > $R2.sorted

cut -f1 $R1.sorted | cut -f1 -d'_' | sort -u > $R1.readsID
cut -f1 $R2.sorted | cut -f1 -d'_' | sort -u > $R2.readsID

count=1

comm -12 $R1.readsID $R2.readsID > $R1.mutual.reads
touch $R1.count.$count

while read RID
do
  i=$(wc -l $R1.count.$count | cut -f1 -d' ')
  if (( $i > 5000 ))
  then
    Rscript scripts/metagenome/filter_paired_end_blast.R $R1.count.$count $R2.count.$count $out.$count $taxid
    count=$((count+1))
  fi
  $SGREP $RID $R1.sorted >> $R1.count.$count
  $SGREP $RID $R2.sorted >> $R2.count.$count
  #echo "LC_ALL=C grep $RID $R1 > $R1.$RID; LC_ALL=C grep $RID $R2 > $R2.$RID; Rscript scripts/metagenome/filter_paired_end_blast.R $R1.$RID $R2.$RID $out.$RID $taxid" >> $R1.prompt

done < $R1.mutual.reads

Rscript scripts/metagenome/filter_paired_end_blast.R $R1.count.$count $R2.count.$count $out.$count $taxid

cat $out.* > $out

rm $R1.* $R2.* $out.*


#for ((j=1; j<=count; j++))
#do
#  echo "Rscript scripts/metagenome/filter_paired_end_blast.R $R1.count.$j $R2.count.$j $out.$j $taxid" >> $R1.prompt
#done
#perl_fork_univ.pl $R1.prompt $threads
