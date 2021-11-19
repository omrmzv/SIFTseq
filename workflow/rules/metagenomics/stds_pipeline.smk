rule subset_FQ_stds:
    input:
        unmapped_R1 = 'sample_output/decontaminate_unmapped/phix_aligned/unmapped/{sample}_unmapped_R1.fastq',
        unmapped_R2 = 'sample_output/decontaminate_unmapped/phix_aligned/unmapped/{sample}_unmapped_R2.fastq',
        R1 = 'Data/samples/{sample}_R1.fastq.gz',
        R2 = 'Data/samples/{sample}_R2.fastq.gz'
    output:
        tempR1 = temp('{sample}.stds_R1.fq'),
        tempR2 = temp('{sample}.stds_R2.fq'),
        R1='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.stds_R1.fastq',
        R2='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.stds_R2.fastq'
    shell:
        """
        zcat {input.R1} | sed 's/ 1/_1/g' > {output.tempR1}
        zcat {input.R2} | sed 's/ 2/_2/g' > {output.tempR2}
        grep "@" {input.unmapped_R1} | sed 's/^@//g' | sed 's/ 1/_1/g' | ./software/seqtk/seqtk subseq {output.tempR1} - > {output.R1}
        grep "@" {input.unmapped_R2} | sed 's/^@//g' | sed 's/ 2/_2/g' | ./software/seqtk/seqtk subseq {output.tempR2} - > {output.R2}
        """

rule adapter_and_quality_filter_stds:
    input:
        'sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.stds_R1.fastq',
        'sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.stds_R2.fastq',
        adapter_file = get_adapter_file
    output:
        R1='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.stds.trim.R1.fastq',
        R2='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.stds.trim.R2.fastq'
    shell:
        """
        ./software/bbmap/bbduk.sh in1={input[0]} in2={input[1]} \
                out1={output.R1} out2={output.R2} -Xmx1g \
                ref={input[2]} \
                threads=1 \
                maq=32
        """
rule dedupe_stds:
    input:
        R1='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.stds.trim.R1.fastq',
        R2='sample_output/decontaminate_unmapped/original_fq_unmapped/{sample}.stds.trim.R2.fastq'
    output:
        R1='sample_output/decontaminate_unmapped/stds/original_fq_unmapped/{sample}.stds.trim.dedupe.R1.fastq',
        R2='sample_output/decontaminate_unmapped/stds/original_fq_unmapped/{sample}.stds.trim.dedupe.R2.fastq'
    resources:
        mem_mb = 20000
    log: 'logs/dedupe/{sample}.dedupe'
    shell:
        """
        ./software/bbmap/clumpify.sh -Xmx65000m in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} dedupe &> {log}
        """

rule FLAHS_stds:
    input:
        R1='sample_output/decontaminate_unmapped/stds/original_fq_unmapped/{sample}.stds.trim.dedupe.R1.fastq',
        R2='sample_output/decontaminate_unmapped/stds/original_fq_unmapped/{sample}.stds.trim.dedupe.R2.fastq'
    output:
        combined ='sample_output/decontaminate_unmapped/stds/original_fq_unmapped/{sample}.stds.trim.dedupe.extendedFrags.fastq',
        R1_not_combined = 'sample_output/decontaminate_unmapped/stds/original_fq_unmapped/{sample}.stds.trim.dedupe.notCombined_1.fastq'
    resources:
        mem_mb = 20000
    log: 'logs/FLASH/{sample}.flash'
    shell:
        """
        /programs/FLASH2/flash2 -q -M75 -O -o {wildcards.sample}.stds.trim.dedupe -d ./sample_output/decontaminate_unmapped/stds/original_fq_unmapped/ {input.R1} {input.R2} &> {log}
        """

rule combine_stds:
    input:
        combined ='sample_output/decontaminate_unmapped/stds/original_fq_unmapped/{sample}.stds.trim.dedupe.extendedFrags.fastq',
        R1_not_combined = 'sample_output/decontaminate_unmapped/stds/original_fq_unmapped/{sample}.stds.trim.dedupe.notCombined_1.fastq'
    output:
        catted = 'sample_output/decontaminate_unmapped/stds/original_fq_unmapped/{sample}.stds.trim.dedupe.merged.cat.fastq'
    shell:
        """
        sed '1~4s/@.*/&-COMBINED/' {input.combined} > {output.catted}
        sed '1~4s/@.*/&-R1SINGLEEND/' {input.R1_not_combined} >> {output.catted}
        """

rule decontaminate_stds:
    input:
        combined = 'sample_output/decontaminate_unmapped/stds/original_fq_unmapped/{sample}.stds.trim.dedupe.merged.cat.fastq',
        REF = '/workdir/apc88/WGBS_pipeline/references/hg19_Bowtie2/hg19.fa'
    output:
        pass0 = 'sample_output/decontaminate_unmapped/stds/kmer_decontaminate/{sample}.stds.pass.fastq',
        decon = 'sample_output/decontaminate_unmapped/stds/kmer_decontaminate/{sample}.stds.decon.fa'
    resources:
        mem_mb=80000
    log: 'logs/decontaminate/{sample}.decon'
    shell:
        """
        {BBDUK} in={input.combined} \
                out={output.pass0} \
                -Xmx{resources.mem_mb}m \
                prealloc=t \
                ref={input.REF} k=50 &> {log}

        grep "@" {output.pass0} | sed 's/^@//g' | ./software/seqtk/seqtk subseq {input.combined} - | \
            fastq_to_fasta -Q33 -o {output.decon}
        """


rule hs_blastn_stds:
    input:
        fa = 'sample_output/decontaminate_unmapped/stds/kmer_decontaminate/{sample}.stds.decon.fa',
        db = 'databases/blast/stds/NCBIGenomes06.fna',
        gi_to_taxid = 'databases/blast/NCBIGenomes06.gis.taxids',
        grammy_db_proof = 'logs/grammy/grammy_prompts_CT'
    output:
        outfmt6 = 'sample_output/stds/blast/{sample}.outfmt6'
    threads: 8
    resources:
        mem_mb=80000
    shell:
        """
		{HSBLASTN} align -query {input.fa} \
                    	-db {input.db} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 > {output.outfmt6}
		"""


rule filter_blast_stds: #remove hits from human, and messy GIs and loow quality BLAST
    input:
        blast_outfmt6 = 'sample_output/stds/blast/{sample}.outfmt6',
        human='/workdir/apc88/COFEESEQ/databases/GenomeDB/human.gis',
        gis = '/workdir/apc88/COFEESEQ/databases/blast/notkept.gis2',
        clean_fasta = 'sample_output/decontaminate_unmapped/stds/kmer_decontaminate/{sample}.stds.decon.fa'
    output:
        human_blasted = 'sample_output/stds/blast_filtered/{sample}.human',
        blast_filtered = 'sample_output/stds/blast_filtered/{sample}.sorted.outfmt6'
    run:
        import os
        import dnaio
        os.system("mkdir -p sample_output/stds/blast_filtered")
        print("1")
        os.system("mkdir -p sample_output/stds/blast_filtered/stds")
        print("2")

        cmd = f"grep -f {input.human} {input.blast_outfmt6} | cut -f1 | sort -u > {output.human_blasted}"
        os.system(cmd)
        human_readIDs = [line.rstrip('\n') for line in open(output.human_blasted)]
        bad_GIs = [line.rstrip('\n') for line in open(input.gis)]
        with open(input.blast_outfmt6) as f, open(output.blast_filtered, 'w') as w:
            for hit in f:
                #if wildcards.conversion == "GA_genome":
                #    print(hit)
                read_id = hit.strip().split('\t')[0]
                gi = hit.strip().split('\t')[1].split('|')[1]
                gi = 'gi|'+gi+'|'
                mol_len = float(hit.strip().split('\t')[3])
                mapped_len = float(hit.strip().split('\t')[12])
                ratio = mol_len/mapped_len
                if (mol_len >= 50) and (ratio > 0.9):
                    if gi not in bad_GIs:
                        if read_id not in human_readIDs:
                            w.write(hit)

rule get_relevant_stds:
    input:
        outfmt6 = 'sample_output/stds/blast_filtered/{sample}.sorted.outfmt6',
        clean_fasta = 'sample_output/decontaminate_unmapped/stds/kmer_decontaminate/{sample}.stds.decon.fa',
        human='/workdir/apc88/COFEESEQ/databases/GenomeDB/human.gis',
        gis = '/workdir/apc88/COFEESEQ/databases/blast/notkept.gis2'
    output:
        cleantblat1 = 'sample_output/stds/grammy/{sample}/{sample}.tblat.1'
    shell:
        """
        grep "Plus/Plus" {input.outfmt6} > {output.cleantblat1} || true
        grep "Plus/Minus" {input.outfmt6} >> {output.cleantblat1} || true
        """

rule grammy_clean_stds:
    input:
        clean_fasta = 'sample_output/decontaminate_unmapped/stds/kmer_decontaminate/{sample}.stds.decon.fa',
        cleantblat1 = 'sample_output/stds/grammy/{sample}/{sample}.tblat.1'
    output:
        nonhumanfa_gz = 'sample_output/stds/grammy/{sample}/{sample}.fa.gz',
        nonhumanfasta_gz = temp('sample_output/stds/grammy/{sample}/{sample}.fasta.gz'),
        rdt = 'sample_output/stds/grammy/{sample}/{sample}.rdt',
        mtx = 'sample_output/stds/grammy/{sample}/{sample}.mtx',
        lld = 'sample_output/stds/grammy/{sample}/{sample}.lld',
        btp = 'sample_output/stds/grammy/{sample}/{sample}.btp',
        est = 'sample_output/stds/grammy/{sample}/{sample}.est',
        gra = 'sample_output/stds/grammy/{sample}/{sample}.gra',
        avl = 'sample_output/stds/grammy/{sample}/{sample}.avl'
    resources: mem_mb=1
    shell: #used to be cat R1 R2 | sed ....
        """
        if [ $(wc -l {input.cleantblat1} | cut -d' ' -f1) -gt 1 ]
        then
            cut -f1 {input.cleantblat1} | sort -u | ./software/seqtk/seqtk subseq {input.clean_fasta} - | gzip -1 > {output.nonhumanfa_gz}
            cd sample_output/stds/grammy/{wildcards.sample}
            python2.7 {GRAMMY_RDT} -t illumina . .
            python2.7 {GRAMMY_PRE} -q "40,75,-5" {wildcards.sample} {GRAMMY_REF_FASTA}
            python2.7 {GRAMMY_EM} -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx
            python2.7 {GRAMMY_POST} {wildcards.sample}.est {GRAMMY_REF_FASTA} {wildcards.sample}.btp
            cd ../../../
        else
            touch {output.nonhumanfa_gz} {output.nonhumanfasta_gz}
            touch {output.rdt} {output.mtx} {output.lld} {output.btp} {output.est} {output.gra} {output.avl}
        fi
		"""

rule annotate_grammy_stds:
	input:
		rdt = 'sample_output/stds/grammy/{sample}/{sample}.rdt',
		mtx = 'sample_output/stds/grammy/{sample}/{sample}.mtx',
		lld = 'sample_output/stds/grammy/{sample}/{sample}.lld',
		btp = 'sample_output/stds/grammy/{sample}/{sample}.btp',
		est = 'sample_output/stds/grammy/{sample}/{sample}.est',
		gra = 'sample_output/stds/grammy/{sample}/{sample}.gra',
		avl = 'sample_output/stds/grammy/{sample}/{sample}.avl',
		tblat1 = 'sample_output/stds/grammy/{sample}/{sample}.tblat.1',
		LUT = 'LUTGrammy/taxids_names_lengths_tax.tab'
	output:
		tab='sample_output/stds/grammy/{sample}/{sample}.tab',
		anno = 'sample_output/stds/grammy/{sample}/{sample}.grammy.tab'
	params:
		DIR='sample_output/',
		DB='stds/grammy/{sample}/',
        stats = 'sample_output/stats/{sample}_mapping_stats.txt',
		LUT='LUTGrammy/taxids_names_lengths_tax.tab'
	shell:
		"""
        if [ $(cut -f1 {input.tblat1} | sort -u | wc -l) -gt 2 ]
        then
            Rscript workflow/scripts/metagenome/filter_gra_file.R {input.gra} {output.tab} {wildcards.sample}
            Rscript workflow/scripts/metagenome/annotate_grammy_apc.R {output.tab} {input.tblat1} {params.stats} {input.LUT} {output.anno}
        else
            touch {output.tab}
            head -n1 {input.LUT} > {output.anno}
        fi
		"""
