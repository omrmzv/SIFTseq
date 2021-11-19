rule phix_alignment:
    input:
        R1 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz',
        R2 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz'
    output:
        bam = temp('sample_output/decontaminate_unmapped/phix_aligned/bam/{sample}.bam'),
        unmapped_R1 = 'sample_output/decontaminate_unmapped/phix_aligned/unmapped/{sample}_unmapped_R1.fastq',
        unmapped_R2 = 'sample_output/decontaminate_unmapped/phix_aligned/unmapped/{sample}_unmapped_R2.fastq'
    params:
        genome_prefix = 'references/phix/phix',
        unmapped_pe = 'sample_output/decontaminate_unmapped/phix_aligned/unmapped/{sample}_unmapped_R%.fastq',
        prep_and_seq_type = get_sample_info
    shell:
        """
        bowtie2 -x {params.genome_prefix} \
                -1 {input.R1} -2 {input.R2} \
                --local --very-sensitive-local \
                --un-conc {params.unmapped_pe} | samtools view -bh - > {output.bam}
        """


rule subset_FQ:
    input:
        unmapped_R1 = 'sample_output/decontaminate_unmapped/phix_aligned/unmapped/{sample}_unmapped_R1.fastq',
        unmapped_R2 = 'sample_output/decontaminate_unmapped/phix_aligned/unmapped/{sample}_unmapped_R2.fastq',
        R1 = '/workdir/apc88/COFEESEQ/Data/samples/{sample}_R1.fastq.gz',
        R2 = '/workdir/apc88/COFEESEQ/Data/samples/{sample}_R2.fastq.gz'
    output:
        tempR1 = temp('{sample}_R1.fq'),
        tempR2 = temp('{sample}_R2.fq'),
        R1='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}_R1.fastq',
        R2='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}_R2.fastq'
    shell:
        """
        zcat {input.R1} | sed 's/ 1/_1/g' > {output.tempR1}
        zcat {input.R2} | sed 's/ 2/_2/g' > {output.tempR2}
        grep "@" {input.unmapped_R1} | sed 's/^@//g' | ./software/seqtk/seqtk subseq {output.tempR1} - > {output.R1}
        grep "@" {input.unmapped_R2} | sed 's/^@//g' | ./software/seqtk/seqtk subseq {output.tempR2} - > {output.R2}
        """

rule adapter_and_quality_filter:
    input:
        'sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}_R1.fastq',
        'sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}_R2.fastq',
        adapter_file = get_adapter_file
    output:
        R1='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.R1.fastq',
        R2='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.R2.fastq'
    shell:
        """
        ./software/bbmap/bbduk.sh in1={input[0]} in2={input[1]} \
                out1={output.R1} out2={output.R2} -Xmx1g \
                ref={input[2]} \
                threads=1 \
                maq=32
        """

rule dedupe:
    input:
        R1='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.R1.fastq',
        R2='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.R2.fastq'
    output:
        R1='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.dedupe.R1.fastq',
        R2='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.dedupe.R2.fastq'
    resources:
        mem_mb = 20000
    log: 'logs/dedupe/{sample}.dedupe'
    shell:
        """
        ./software/bbmap/clumpify.sh -Xmx65000m in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} dedupe &> {log}
        """

rule FLAHS:
    input:
        R1='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.dedupe.R1.fastq',
        R2='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.dedupe.R2.fastq'
    output:
        combined ='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.dedupe.extendedFrags.fastq',
        R1_not_combined = 'sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.dedupe.notCombined_1.fastq'
    resources:
        mem_mb = 20000
    log: 'logs/FLASH/{sample}.flash'
    shell:
        """
        /programs/FLASH2/flash2 -q -M75 -O -o {wildcards.sample}.trim.dedupe -d ./sample_output/decontaminate_unmapped/original_fq_unmapped/ {input.R1} {input.R2} &> {log}
        """

rule combine:
    input:
        combined ='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.dedupe.extendedFrags.fastq',
        R1_not_combined = 'sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.dedupe.notCombined_1.fastq'
    output:
        catted = 'sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.dedupe.merged.cat.fastq'
    shell:
        """
        sed '1~4s/@.*/&-COMBINED/' {input.combined} > {output.catted}
        sed '1~4s/@.*/&-R1SINGLEEND/' {input.R1_not_combined} >> {output.catted}
        """

rule decontaminate:
    input:
        combined = 'sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.trim.dedupe.merged.cat.fastq',
        CT = '/workdir/apc88/WGBS_pipeline/references/for_decontamination/CT.fa',
        GA = '/workdir/apc88/WGBS_pipeline/references/for_decontamination/GA.fa'
    output:
        insilico_CT = 'sample_output/decontaminate_unmapped/kmer_decontaminate/{sample}.CT.fastq',
        pass1 = 'sample_output/decontaminate_unmapped/kmer_decontaminate/{sample}.CT.pass1.fastq',
        pass2 = 'sample_output/decontaminate_unmapped/kmer_decontaminate/{sample}.CT.pass2.fastq',
        decon = 'sample_output/decontaminate_unmapped/kmer_decontaminate/{sample}.decon.fa',
        deconCT = 'sample_output/decontaminate_unmapped/kmer_decontaminate/{sample}.CT.decon.fa'
    resources:
        mem_mb=80000
    log: 'logs/decontaminate/{sample}.decon'
    shell:
        """
        python workflow/scripts/metagenome/insilico_conversion2.py -i {input.combined} -o {output.insilico_CT} -r R1
        {BBDUK} in={output.insilico_CT} \
                out={output.pass1} \
                -Xmx{resources.mem_mb}m \
                prealloc=t rcomp=f \
                ref={input.GA} k=50 &> {log}
        {BBDUK} in={output.pass1} \
                out={output.pass2} \
                -Xmx{resources.mem_mb}m \
                prealloc=t rcomp=f \
                ref={input.CT} k=50 &>> {log}
        grep "@" {output.pass2} | sed 's/^@//g' | ./software/seqtk/seqtk subseq {input.combined} - | \
            fastq_to_fasta -Q33 -o {output.decon}

        fastq_to_fasta -Q33 -i {output.pass2} -o {output.deconCT}
        """
