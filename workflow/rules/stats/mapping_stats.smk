rule mapping_stats:
    input:
        original_r1 = CLUSTER+config['DATA']+'samples/{sample}_R1.fastq.gz',
        trimmed_r1 = 'sample_output/trim/{sample}_R2_trim.fastq',
        raw_bam = 'sample_output/aligned/raw_aligned/{sample}.bam',
        deduped_bam = 'sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam',
        deduped_bai = 'sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam.bai',
    output:
        stats = 'sample_output/stats/{sample}_mapping_stats.txt'
    params:
        mappable_hg19='2948611470',
        mappable_chr21 = '40088623'
    shell:
        """
        original_reads=$(pigz -dc {input.original_r1} | wc -l)
        original_reads=$((original_reads / 4))
        trimmed_reads=$(wc -l {input.trimmed_r1} | cut -d' ' -f1)
        trimmed_reads=$((trimmed_reads / 4))
        mapped_reads=$(samtools view -c {input.raw_bam})
        mapped_reads=$((mapped_reads / 2))
        deduped_reads=$(samtools view -c {input.deduped_bam})
        deduped_reads=$((deduped_reads / 2))
        mapping_eff=$(echo "scale=2;$mapped_reads/$trimmed_reads" | bc)
        deduped_frac=$(echo "scale=2;$deduped_reads/$mapped_reads" | bc)
        depth_var=$(samtools depth {input.deduped_bam} | awk '{{sum+=$3}} END {{print sum/{params.mappable_hg19}" "NR/{params.mappable_hg19}}}')
        depth=$(echo $depth_var | cut -f1 -d' ')
        bp_frac=$(echo $depth_var | cut -f2 -d' ')

        echo -e "SAMPLE\tNUM_READS\tREADS_AFTER_TRIM\tALIGNED\tMAPPING_EFFICIENCY\tDEDUPED_READS\tDEDUP_FRAC\tDEPTH\tFRACTION_PB_BP" > {output.stats}
        echo -e "{wildcards.sample}\t$original_reads\t$trimmed_reads\t$mapped_reads\t$mapping_eff\t$deduped_reads\t$deduped_frac\t$depth\t$bp_frac" >> {output.stats}
        """
