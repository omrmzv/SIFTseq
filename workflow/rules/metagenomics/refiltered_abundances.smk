rule refilter:
    input:
        genome_fasta = 'NCBIGenomes06.fna.1',
        unfiltered_grammy = 'sample_output/unfiltered/grammy/{sample}/{sample}.grammy.tab',
        unfiltered_tblat1_file = 'sample_output/unfiltered/grammy/{sample}/{sample}.tblat.1',
        filtered_grammy = 'sample_output/C_poor/grammy/{sample}/{sample}.grammy.tab',
        filtered_tblat1_file = 'sample_output/C_poor/grammy/{sample}/{sample}.tblat.1',
        gi_to_tax = 'databases/blast/NCBIGenomes06.gis.taxids'
    output:
        filt_comp = 'sample_output/refiltered/{sample}/figures/unfilt_vs_filt_AdjBlast.png',
        refiltered_tblat1 = 'sample_output/refiltered/{sample}/{sample}.tblat.1'
    params:
        outdir = './sample_output/refiltered/{sample}/'
    shell:
        """
        if [ $(wc -l {input.filtered_tblat1_file} | cut -d' ' -f1) -gt 2 ]
        then
            python ./workflow/scripts/metagenome/notebook.py \
                {input.genome_fasta} \
                {input.unfiltered_grammy} \
                {input.unfiltered_tblat1_file} \
                {input.filtered_grammy} \
                {input.filtered_tblat1_file} \
                {input.gi_to_tax} \
                {params.outdir} \
                {output.refiltered_tblat1}
        else
            touch {output.filt_comp}
            touch {output.refiltered_tblat1}
        fi
        """

rule grammy_clean_refiltered:
    input:
        clean_fasta = 'sample_output/C_poor/fasta/{sample}.cpoor.fa',
        cleantblat1 = 'sample_output/refiltered/{sample}/{sample}.tblat.1'
    output:
        nonhumanfa_gz = 'sample_output/refiltered/{sample}/{sample}.fa.gz',
        nonhumanfasta_gz = temp('sample_output/refiltered/{sample}/{sample}.fasta.gz'),
        rdt = 'sample_output/refiltered/{sample}/{sample}.rdt',
        mtx = 'sample_output/refiltered/{sample}/{sample}.mtx',
        lld = 'sample_output/refiltered/{sample}/{sample}.lld',
        btp = 'sample_output/refiltered/{sample}/{sample}.btp',
        est = 'sample_output/refiltered/{sample}/{sample}.est',
        gra = 'sample_output/refiltered/{sample}/{sample}.gra',
        avl = 'sample_output/refiltered/{sample}/{sample}.avl'
    resources: mem_mb=1
    shell: #used to be cat R1 R2 | sed ....
        """
        if [ $(wc -l {input.cleantblat1} | cut -d' ' -f1) -gt 1 ]
        then
            cut -f1 {input.cleantblat1} | sort -u | ./software/seqtk/seqtk subseq {input.clean_fasta} - | gzip -1 > {output.nonhumanfa_gz}
            cd sample_output/refiltered/{wildcards.sample}
            python2.7 {GRAMMY_RDT} -t illumina . .
            python2.7 {GRAMMY_PRE} -q "40,75,-5" {wildcards.sample} {GRAMMY_REF_FASTA}
            python2.7 {GRAMMY_EM} -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx
            python2.7 {GRAMMY_POST} {wildcards.sample}.est {GRAMMY_REF_FASTA} {wildcards.sample}.btp
            cd ../../
        else
            touch {output.nonhumanfa_gz} {output.nonhumanfasta_gz}
            touch {output.rdt} {output.mtx} {output.lld} {output.btp} {output.est} {output.gra} {output.avl}
        fi
		"""

rule annotate_grammy_refiltered:
	input:
		rdt = 'sample_output/refiltered/{sample}/{sample}.rdt',
		mtx = 'sample_output/refiltered/{sample}/{sample}.mtx',
		lld = 'sample_output/refiltered/{sample}/{sample}.lld',
		btp = 'sample_output/refiltered/{sample}/{sample}.btp',
		est = 'sample_output/refiltered/{sample}/{sample}.est',
		gra = 'sample_output/refiltered/{sample}/{sample}.gra',
		avl = 'sample_output/refiltered/{sample}/{sample}.avl',
		tblat1 = 'sample_output/refiltered/{sample}/{sample}.tblat.1',
		LUT = 'LUTGrammy/taxids_names_lengths_tax.tab'
	output:
		tab='sample_output/refiltered/{sample}/{sample}.tab',
		anno = 'sample_output/refiltered/{sample}/{sample}.grammy.tab'
	params:
		DIR='sample_output/',
		DB='grammy/{sample}/',
        stats = 'sample_output/stats/{sample}_mapping_stats.txt',
		LUT='LUTGrammy/taxids_names_lengths_tax.tab'
	shell:
		"""
        if [ $(wc -l {input.tblat1} | cut -d' ' -f1) -gt 1 ]
        then
            Rscript workflow/scripts/metagenome/filter_gra_file.R {input.gra} {output.tab} {wildcards.sample}
            Rscript workflow/scripts/metagenome/annotate_grammy_apc.R {output.tab} {input.tblat1} {params.stats} {input.LUT} {output.anno}
        else
            touch {output.tab}
            head -n1 {input.LUT} > {output.anno}
        fi
		"""
