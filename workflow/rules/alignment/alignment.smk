rule pe_bismark:
	input:
		genome = get_reference_genome,
		r1p = 'sample_output/trim/{sample}_R1_trim.fastq',
		r2p = 'sample_output/trim/{sample}_R2_trim.fastq'
	output:
		bam = temp('sample_output/pe_bisulfite_aligned/raw_aligned/{sample}.bam'),
		unmapped_R1 = temp('sample_output/pe_bisulfite_aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz'),
		unmapped_R2 = temp('sample_output/pe_bisulfite_aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz')
	log: 'logs/alignment/{sample}.alignment.log'
	threads: alignment_threads
	params:
		outdir = 'sample_output/pe_bisulfite_aligned/raw_aligned/'
	shell:
		"""
		threads=$(( {threads} / 2 ))
		{BISMARK} --genome {input.genome} \
					--parallel $threads \
					--quiet \
					--unmapped \
					-o {params.outdir} \
					-1 {input.r1p} \
					-2 {input.r2p}
		mv {params.outdir}{wildcards.sample}_R1_trim_bismark_bt2_pe.bam {output.bam}
		mv {params.outdir}{wildcards.sample}_R1_trim_bismark_bt2_PE_report.txt {log}

		mv {params.outdir}{wildcards.sample}_R1_trim.fastq_unmapped_reads_1.fq.gz {output.unmapped_R1}
		mv {params.outdir}{wildcards.sample}_R2_trim.fastq_unmapped_reads_2.fq.gz {output.unmapped_R2}
		"""

rule filter_bismark:
	input:
		bam = 'sample_output/pe_bisulfite_aligned/raw_aligned/{sample}.bam'
	output:
		bismark_dup = temp('sample_output/pe_bisulfite_aligned/raw_aligned/{sample}.deduplicated.bam'),
		mapped_all_chr=temp('sample_output/pe_bisulfite_aligned/all_chr/{sample}_mapped_all_chr.bam'),
		mapped_all_chr_bai = temp('sample_output/pe_bisulfite_aligned/all_chr/{sample}_mapped_all_chr.bam.bai'),
		dedup_file = temp('sample_output/pe_bisulfite_aligned/raw_aligned/{sample}.deduplication_report.txt')
	params:
		mapQ = '15',
		outdir = 'sample_output/pe_bisulfite_aligned/raw_aligned/'
	log: 'logs/deduplication/{sample}.log'
	threads: filter_bam_threads
	shell:
		"""
		{RMDUPS} -o {wildcards.sample} -p --output_dir {params.outdir} --bam {input.bam} &> {log}
		samtools sort -@ {threads} {output.bismark_dup} | samtools view -h -q {params.mapQ} -o {output.mapped_all_chr}
		samtools index {output.mapped_all_chr}
		"""

rule pe_bowtie2:
	input:
		genome = get_reference_genome_fasta,
		r1p = 'sample_output/trim/{sample}_R1_trim.fastq',
		r2p = 'sample_output/trim/{sample}_R2_trim.fastq',
	output:
		bam = 'sample_output/pe_stdseq_aligned/raw_aligned/{sample}.bam',
		unmapped_R1 = 'sample_output/pe_stdseq_aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz',
		unmapped_R2 = 'sample_output/pe_stdseq_aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz'
	log: 'logs/alignment/{sample}.alignment.log'
	params:
		unmapped_pe = 'sample_output/pe_stdseq_aligned/unmapped/{sample}_pe_unmapped_R%.fastq.gz'
	threads: alignment_threads
	shell:
		"""
		genome_prefix=$(echo {input.genome} | cut -f1 -d'.')
		bowtie2 -x $genome_prefix \
				-1 {input.r1p} -2 {input.r2p} \
				-p {threads} -q --score-min L,0,-0.2 --ignore-quals --no-mixed --no-discordant --dovetail --maxins 500 \
				--un-conc-gz {params.unmapped_pe} | samtools view -f3 -bh - > {output.bam}
		"""

rule filter_bowtie2:
	input:
		bam = 'sample_output/pe_stdseq_aligned/raw_aligned/{sample}.bam'
	output:
		mapped_all_chr=temp('sample_output/pe_stdseq_aligned/all_chr/{sample}_mapped_all_chr.bam'),
		mapped_all_chr_bai = temp('sample_output/pe_stdseq_aligned/all_chr/{sample}_mapped_all_chr.bam.bai')
	params:
		mapQ = '15'
	threads: filter_bam_threads
	shell:
		"""
		samtools sort -n -@ {threads} {input.bam} | \
			samtools fixmate -m -@ {threads} - - | \
			samtools sort -@ {threads} | \
			samtools markdup -r -@ {threads} - - | samtools sort -@ {threads} -o {output.mapped_all_chr}
		samtools index {output.mapped_all_chr}
		"""

def aggregate_bam_files(wildcards):
	"""
	we have this function to figure out which inputs are hg19, bisulfite, etc.
	"""
	sample_name, prep_type, seq_mode, abundance_control, sample_type, seq_type = get_sample_info(wildcards)

	if "2x" in seq_mode and seq_type == "bisulfite":
		bam = 'sample_output/pe_bisulfite_aligned/raw_aligned/{sample}.bam'
		mapped_all_chr = 'sample_output/pe_bisulfite_aligned/all_chr/{sample}_mapped_all_chr.bam'
		mapped_all_chr_bai = 'sample_output/pe_bisulfite_aligned/all_chr/{sample}_mapped_all_chr.bam.bai'
		unmapped_R1 = 'sample_output/pe_bisulfite_aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz'
		unmapped_R2 = 'sample_output/pe_bisulfite_aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz'

	if "2x" in seq_mode and seq_type == "standard":
		bam = 'sample_output/pe_stdseq_aligned/raw_aligned/{sample}.bam'
		mapped_all_chr = 'sample_output/pe_stdseq_aligned/all_chr/{sample}_mapped_all_chr.bam'
		mapped_all_chr_bai = 'sample_output/pe_stdseq_aligned/all_chr/{sample}_mapped_all_chr.bam.bai'
		unmapped_R1 = 'sample_output/pe_stdseq_aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz'
		unmapped_R2 = 'sample_output/pe_stdseq_aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz'

	return[bam, mapped_all_chr, mapped_all_chr_bai, unmapped_R1, unmapped_R2]


rule aggregate_bams:
	input:
		aggregate_bam_files
	output:
		bam = 'sample_output/aligned/raw_aligned/{sample}.bam',
		mapped_all_chr = 'sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam',
		mapped_all_chr_bai ='sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam.bai',
		unmapped_R1 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz',
		unmapped_R2 ='sample_output/aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz'
	shell:
		"""
		cp {input[0]} {output.bam}
		cp {input[1]} {output.mapped_all_chr}
		cp {input[2]} {output.mapped_all_chr_bai}
		cp {input[3]} {output.unmapped_R1}
		cp {input[4]} {output.unmapped_R2}
		"""

rule filter_nonconversion:
	input:
		mapped_all_chr = 'sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam'
	output:
		namesorted = temp('sample_output/aligned/all_chr/{sample}_all_chr.bam'),
		filtered = 'sample_output/aligned/all_chr/{sample}_all_chr.nonCG_filtered.bam',
		removed = 'sample_output/aligned/all_chr/{sample}_all_chr.nonCG_removed_seqs.bam',
		txt = temp('sample_output/aligned/all_chr/{sample}_all_chr.non-conversion_filtering.txt')
	log: 'logs/filterCH/{sample}.txt'
	shell:
		"""
		samtools sort -n {input.mapped_all_chr} -o {output.namesorted}
		{FILTERNONCONV} -p --threshold 3 --minimum_count 3 {output.namesorted}
		cp {output.txt} {log}
		"""
