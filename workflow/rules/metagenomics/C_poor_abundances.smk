rule get_C_poor:
    input:
        fasta = 'sample_output/decontaminate_unmapped/kmer_decontaminate/{sample}.decon.fa'
    output:
        clean_fasta = 'sample_output/C_poor/fasta/{sample}.cpoor.fa'
    run:
        import os
        import dnaio
        os.system("mkdir -p sample_output/C_poor")
        os.system("mkdir -p sample_output/C_poor/fasta")

        with dnaio.open(input.fasta) as f, dnaio.open(output.clean_fasta, mode = 'w') as w:
            for record in f:
                seq = str(record.sequence)
                if (seq.count('CG')==0) and (seq.count('C')<=4):
                    w.write(record)

rule hs_blastn_wgbs_cpoor:
    input:
        clean_fasta = 'sample_output/C_poor/fasta/{sample}.cpoor.fa',
        db_CT = 'databases/blast/CT_conversion/NCBIGenomes06_CT.fna',
        db_GA = 'databases/blast/GA_conversion/NCBIGenomes06_GA.fna',
        gi_to_taxid = 'databases/blast/NCBIGenomes06.gis.taxids',
        grammy_db_proof = 'logs/grammy/grammy_prompts_CT'
    output:
        outfmt6_CT = 'sample_output/C_poor/blast/CT_genome/{sample}.outfmt6',
        outfmt6_GA = 'sample_output/C_poor/blast/GA_genome/{sample}.outfmt6'
    threads: 8
    resources:
        mem_mb=80000
    shell:
        """
		{HSBLASTN} align -query {input.clean_fasta} \
                    	-db {input.db_CT} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 > {output.outfmt6_CT}
        {HSBLASTN} align -query {input.clean_fasta} \
                    	-db {input.db_GA} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 > {output.outfmt6_GA}
		"""

rule filter_blast_cpoor: #remove hits from human, and messy GIs and loow quality BLAST
    input:
        blast_outfmt6 = 'sample_output/C_poor/blast/{conversion}/{sample}.outfmt6',
        human='/workdir/apc88/COFEESEQ/databases/GenomeDB/human.gis',
        gis = '/workdir/apc88/COFEESEQ/databases/blast/notkept.gis2',
        clean_fasta = 'sample_output/C_poor/fasta/{sample}.cpoor.fa'
    output:
        human_blasted = 'sample_output/C_poor/blast_filtered/{conversion}/{sample}.human',
        blast_filtered = 'sample_output/C_poor/blast_filtered/{conversion}/{sample}.sorted.outfmt6'
    run:
        import os
        import dnaio
        os.system("mkdir -p sample_output/C_poor/blast_filtered")
        os.system("mkdir -p sample_output/C_poor/blast_filtered/CT_genome")
        os.system("mkdir -p sample_output/C_poor/blast_filtered/GA_genome")

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

rule get_relevant_cpoor:
    input:
        outfmt6_CT = 'sample_output/C_poor/blast_filtered/CT_genome/{sample}.sorted.outfmt6',
        outfmt6_GA = 'sample_output/C_poor/blast_filtered/GA_genome/{sample}.sorted.outfmt6',
        clean_fasta = 'sample_output/C_poor/fasta/{sample}.cpoor.fa',
        human='/workdir/apc88/COFEESEQ/databases/GenomeDB/human.gis',
        gis = '/workdir/apc88/COFEESEQ/databases/blast/notkept.gis2'
    output:
        cleantblat1 = 'sample_output/C_poor/grammy/{sample}/{sample}.tblat.1'
    shell:
        """
        grep "Plus/Plus" {input.outfmt6_CT} > {output.cleantblat1} || true
        grep "Plus/Minus" {input.outfmt6_GA} >> {output.cleantblat1} || true
        """


rule grammy_clean_cpoor:
    input:
        clean_fasta = 'sample_output/C_poor/fasta/{sample}.cpoor.fa',
        cleantblat1 = 'sample_output/C_poor/grammy/{sample}/{sample}.tblat.1'
    output:
        nonhumanfa_gz = 'sample_output/C_poor/grammy/{sample}/{sample}.fa.gz',
        nonhumanfasta_gz = temp('sample_output/C_poor/grammy/{sample}/{sample}.fasta.gz'),
        rdt = 'sample_output/C_poor/grammy/{sample}/{sample}.rdt',
        mtx = 'sample_output/C_poor/grammy/{sample}/{sample}.mtx',
        lld = 'sample_output/C_poor/grammy/{sample}/{sample}.lld',
        btp = 'sample_output/C_poor/grammy/{sample}/{sample}.btp',
        est = 'sample_output/C_poor/grammy/{sample}/{sample}.est',
        gra = 'sample_output/C_poor/grammy/{sample}/{sample}.gra',
        avl = 'sample_output/C_poor/grammy/{sample}/{sample}.avl'
    resources: mem_mb=1
    shell: #used to be cat R1 R2 | sed ....
        """
        if [ $(wc -l {input.cleantblat1} | cut -d' ' -f1) -gt 1 ]
        then
            cut -f1 {input.cleantblat1} | sort -u | ./software/seqtk/seqtk subseq {input.clean_fasta} - | gzip -1 > {output.nonhumanfa_gz}
            cd sample_output/C_poor/grammy/{wildcards.sample}
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

rule annotate_grammy_cpoor:
	input:
		rdt = 'sample_output/C_poor/grammy/{sample}/{sample}.rdt',
		mtx = 'sample_output/C_poor/grammy/{sample}/{sample}.mtx',
		lld = 'sample_output/C_poor/grammy/{sample}/{sample}.lld',
		btp = 'sample_output/C_poor/grammy/{sample}/{sample}.btp',
		est = 'sample_output/C_poor/grammy/{sample}/{sample}.est',
		gra = 'sample_output/C_poor/grammy/{sample}/{sample}.gra',
		avl = 'sample_output/C_poor/grammy/{sample}/{sample}.avl',
		tblat1 = 'sample_output/C_poor/grammy/{sample}/{sample}.tblat.1',
		LUT = 'LUTGrammy/taxids_names_lengths_tax.tab'
	output:
		tab='sample_output/C_poor/grammy/{sample}/{sample}.tab',
		anno = 'sample_output/C_poor/grammy/{sample}/{sample}.grammy.tab'
	params:
		DIR='sample_output/',
		DB='grammy/{sample}/',
        stats = 'sample_output/stats/{sample}_mapping_stats.txt',
		LUT='LUTGrammy/taxids_names_lengths_tax.tab'
	shell:
		"""
        if [ $(wc -l {input.tblat1} | cut -d' ' -f1) -gt 2 ]
        then
            Rscript workflow/scripts/metagenome/filter_gra_file.R {input.gra} {output.tab} {wildcards.sample}
            echo "aa"
            Rscript workflow/scripts/metagenome/annotate_grammy_apc.R {output.tab} {input.tblat1} {params.stats} {input.LUT} {output.anno}
        else
            touch {output.tab}
            head -n1 {input.LUT} > {output.anno}
        fi
		"""
