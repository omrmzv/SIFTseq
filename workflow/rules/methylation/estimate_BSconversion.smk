rule estimate_BSconversion:
	input:
		mapped = 'sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam',
		genome = get_reference_genome_fasta
	output:
		mr = temp('sample_output/{sample}.mr'),
        bsrate = 'sample_output/conversion_rates/{sample}.bsrate.txt',
		sample_rate = 'sample_output/conversion_rates/{sample}.bsconversion'
	shell:
		"""
        {METHPIPETOMR} -m general -o {output.mr} -L 500 {input.mapped}
		{METHPIPEBSRATE} -c {input.genome} -o {output.bsrate} {output.mr}
        X=$(head -1 {output.bsrate})
        echo -e "{wildcards.sample}\t$X" > {output.sample_rate}
        """

rule estimate_BSconversion2:
	input:
		mapped = 'sample_output/aligned/all_chr/{sample}_all_chr.nonCG_filtered.bam',
		genome = get_reference_genome_fasta
	output:
		mr = temp('sample_output/{sample}.CHfilt.mr'),
        bsrate = 'sample_output/conversion_rates/{sample}.CHfilt.bsrate.txt',
		sample_rate = 'sample_output/conversion_rates/{sample}.CHfilt.bsconversion'
	shell:
		"""
        {METHPIPETOMR} -m general -o {output.mr} -L 500 {input.mapped}
		LC_ALL=C sort -T ./ -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o {output.mr} {output.mr}
		{METHPIPEBSRATE} -c {input.genome} -o {output.bsrate} {output.mr}
        X=$(head -1 {output.bsrate})
        echo -e "{wildcards.sample}\t$X" > {output.sample_rate}
        """
