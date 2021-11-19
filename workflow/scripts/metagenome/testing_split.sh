start_time="$(date -u +%s)"

R2='sample_output/blast/GA_genome/R2/urine001.R2.outfmt6.gz'
R1='sample_output/blast/GA_genome/R1/urine001.R1.outfmt6.gz'
taxid_lengths='databases/GenomeDB/taxids_lengths.txt'
threads=10
mkdir -p 'GA'
mkdir -p 'GA/TEMP_urine001/'

pigz -dc -p $threads $R1 | while read line
do
  id=$(echo $line | cut -f1 | cut -f1 -d'_')
  echo $line >> GA/TEMP_urine001/R1.$id
  echo $line >> GA/TEMP_urine001/R1_IDS
done"

end_time="$(date -u +%s)"

elapsed="$(($end_time-$start_time))"
echo "Total of $elapsed seconds elapsed R1"

pigz -dc -p $threads | while read line
do
  id=$(echo $line | cut -f1 | cut -f1 -d'_')
  echo $line >> GA/TEMP_urine001/R2.$id
  echo $line >> GA/TEMP_urine001/R2_IDS
done

end_time="$(date -u +%s)"

elapsed="$(($end_time-$start_time))"
echo "Total of $elapsed seconds elapsed for R2"

comm -12 <(sort -u GA/TEMP_urine001/R1_IDS) <(sort -u GA/TEMP_urine001/R2_IDS) > GA/TEMP_urine001_COMMON



cat TEMP_urine001_COMMON | while read line
do
  Rscript scripts/metagenome/filter_paired_end_blast2.R GA/TEMP_urine001/R1.$line GA/TEMP_urine001/R2.$line test $taxid_lengths
done

end_time="$(date -u +%s)"

elapsed="$(($end_time-$start_time))"
echo "Total of $elapsed seconds elapsed for process"
