################################################################################
# FUNCTIONS
################################################################################
def string_split(string, delimiter, number_of_splits):
    s=string.split(delimiter, number_of_splits)
    return(s)

def get_fastq_reads(wcrds):
    """
    We get the FQs based on if we got SE or PE
    """
    with open(SEQUENCING_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_name, prep_type, seq_mode, abundance_control, sample_type, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_name:
                seq = seq_mode.split('x')[0]
                prep = prep_type
    if seq == '1': #single-end
        return[CLUSTER+config['DATA']+'samples/{sample}_R1.fastq.gz', CLUSTER+config['DATA']+'samples/{sample}_R1.fastq.gz'] #second R1 is a dumb placeholder
    if seq == '2': #paired-end
        return[CLUSTER+config['DATA']+'samples/{sample}_R1.fastq.gz', CLUSTER+config['DATA']+'samples/{sample}_R2.fastq.gz']

    raise(ValueError('Could not figure paired end or single end'))

def get_sample_info(wcrds):
    """
    We get the FQs based on if we got SE or PE
    """
    with open(SEQUENCING_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_name, prep_type, seq_mode, abundance_control, sample_type, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_name:
                return[sample_name, prep_type, seq_mode, abundance_control, sample_type, seq_type]

def get_reference_genome(wcrds):
    sample_name, prep_type, seq_mode, abundance_control, sample_type, seq_type = get_sample_info(wcrds)
    if sample_type == "hg19":
        if seq_type == "bisulfite":
            dir = HG19METH
        if seq_type == "standard":
            dir = HG19STD
    if sample_type == "microbial_control" and seq_type == "bisulfite":
        dir = CTLMETH
    if sample_type == "phix":
        if seq_type == "bisulfite":
            dir = PHIXMETH
        if seq_type == "standard":
            dir = PHIXSTD
    if sample_type == "mm10":
        if seq_type == "bisulfite":
            dir = MM10METH
        if seq_type == "standard":
            dir = MM10STD
    return(dir)

def get_reference_genome_fasta(wcrds):
    sample_name, prep_type, seq_mode, abundance_control, sample_type, seq_type = get_sample_info(wcrds)
    if sample_type == "hg19":
        if seq_type == "bisulfite":
            dir = HG19METH + 'hg19.fa'
        if seq_type == "standard":
            dir = HG19STD + 'hg19.fa'
    if sample_type == "microbial_control" and seq_type == "bisulfite":
        dir = CTLMETH + 'mic.fa'
    if sample_type == "phix":
        if seq_type == "bisulfite":
            dir = PHIXMETH + 'phix.fa'
        if seq_type == "standard":
            dir = PHIXSTD + 'phix.fa'
    if sample_type == "mm10":
        if seq_type == "bisulfite":
            dir = MM10METH + 'mm10.fa'
        if seq_type == "standard":
            dir = MM10STD + 'mm10.fa'
    return(dir)

def get_adapter_file(wcrds):
    """
    For trimming
    """
    with open(SEQUENCING_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_name, prep_type, seq_mode, abundance_control, sample_type, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_name:
                prep = prep_type
    if prep == "MEYER_SSLP":
        adapt = MEYER_ADAPTORS
    if prep == "SRSLY_SSLP":
        adapt = SRSLY_ADAPTORS
    if prep == "SWIFT_ACCEL":
        adapt = SA_ADAPTORS
    if prep == "MEYER_SRSLY":
        adapt = MEYER_SRSLY_ADAPTORS

    return(adapt)

def aget_training_samples():
    """
    Get the samples from CFS/training_samples.tsv
    """
    sams2 = []
    print("hi")
    with open('CFS/training_samples.tsv') as f:
        next(f)
        for line in f:
            sams2.append(line.split('\t')[0])
    print(sams2)
    print(len(sams2))
    return(sams2)

def is_spring(wcrds):
    """
    Determine if the input is FQ or SPRING
    """
    if os.path.exists(CLUSTER + config['DATA']+'samples/{sample}.spring') and not os.path.exists(CLUSTER+config['DATA']+'samples/{sample}_R1.fastq.gz'):
        return(True)
    else:
        return(False)
