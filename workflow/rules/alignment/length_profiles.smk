rule length_profiles:
	input:
		bam = 'sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam',
		bai = 'sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam.bai'
	output:
		lengths='sample_output/Lengths/{sample}_aligned.lengths.gz',
		length_counts = 'sample_output/Lengths/{sample}_aligned.lengths.counts',
		pdf='sample_output/Lengths/{sample}_FragsHistogram.pdf'
	shell:
		"""
		samtools view {input.bam} | awk -v OFS='\t' '{{if ($9>0 && $9<1000) print $1,$9}}' | gzip -9 > {output.lengths}
        zcat {output.lengths} | cut -f2 | sort | uniq -c > {output.length_counts}
		Rscript workflow/scripts/length_profiles/HistogramFragmentLengths.R {output.length_counts} {output.pdf}
		"""
