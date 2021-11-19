#!/usr/bin/env bash

R1=$1
R2=$2
basename=$3
taxid_lengths=$4
threads=$5
out=$6

random-string() {
        cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w ${1:-10} | head -n 1
}

dir=$(random-string)"$basename""tmp"

mkdir -p $dir

t=$(( threads / 2 ))

echo "pigz -dc -p $t $R1 | sort -S 2G --parallel=$t -T $dir > $dir/$basename.R1.outfmt6.sorted; cut -d'_' -f1 $dir/$basename.R1.outfmt6.sorted | cut -f1 | sort -S 2G --parallel=$t -T $dir -u > $dir/$basename.R1.readids" >> $dir/prompt
echo "pigz -dc -p $t $R2 | sort -S 2G --parallel=$t -T $dir > $dir/$basename.R2.outfmt6.sorted; cut -d'_' -f1 $dir/$basename.R2.outfmt6.sorted | cut -f1 | sort -S 2G --parallel=$t -T $dir -u > $dir/$basename.R2.readids" >> $dir/prompt

parallel -j $threads < $dir/prompt

comm -12 <(sort -S 3G --parallel=$t -T $dir $dir/$basename.R1.readids) <(sort -S 3G --parallel=$t -T $dir $dir/$basename.R2.readids) > $dir/READIDS

echo "we sorting now homies"

Rscript ./scripts/metagenome/filter_paired_end_blast5.R $taxid_lengths $dir/READIDS $dir $threads $basename

echo -e "taxid\tqseqid\tstrand\tsseqid\tpident_R1\tlength_R1\tmismatch_R1\tgapopen_R1\tqstart_R1\tqend_R1\tsstart_R1\tsend_R1\tevalue_R1\tbitscore_R1\tqlen_R1\tpident_R2\tlength_R2\tmismatch_R2\tgapopen_R2\tqstart_R2\tqend_R2\tsstart_R2\tsend_R2\tevalue_R2\tbitscore_R2\tqlen_R2\tgenome_len\teffective_length" > $out

find $dir -name "*final" -type f -exec cat {} + > $out
####
rm -r $dir
