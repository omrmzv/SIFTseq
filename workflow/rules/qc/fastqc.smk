rule fqc:
	input:
		get_fastq_reads
	output:
		'sample_output/fastqc/{sample}_R1_fastqc.html',
		'sample_output/fastqc/{sample}_R2_fastqc.html'
	threads: fastqc_threads
	params:
		outdir = 'sample_output/fastqc/',
		prep_and_seq_type = get_sample_info
	log: 'logs/fastqc/{sample}.fqc.log'
	shell:
		"""
		fastqc {input[0]} {input[1]} -t {threads} --outdir {params.outdir} &>{log}
		"""
