import argparse
from pyfaidx import Fasta
import dnaio
import random

def rev_complement(string):
    d = {'A': 'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    reverse = string[::-1]
    rev_comp = ''
    for i in reverse:
        rev_comp += d.get(i, 'N')

    return(rev_comp)

def get_cytosines(r1_start, r1_end, r2_start, r2_end, record1, record2, reference_sequence, strand):
    # use length and positions to determine what regions to get in the case of circular genomes
    R1 = record1.sequence[(r1_start-1):(r1_end)]
    R2 = rev_complement(record2.sequence)[(r2_start-1):(r2_end)]
    qlen1 = len(R1)
    qlen2 = len(R2)
    qual1 = record1.qualities[(r1_start-1):(r1_end)]
    qual2 = record2.qualities[::-1][(r2_start-1):(r2_end)]
    if qlen1 + qlen2 < len(reference_sequence):
        #no overlap
        molecule = R1 + (len(reference_sequence) - qlen1 - qlen2)*'N' + R2
    else:
        #with overlap ughhhhhh
        # any overlap would have their start and stops ordered as A/B/A/B
        #assume R1 is at 5' and R2 is at 3'
        R1_r15p = R1 + (len(reference_sequence) - len(R1)) * 'N'
        qual1_r15p = qual1 + (len(reference_sequence) - len(R1)) * '!'
        R2_r15p = (len(reference_sequence) - len(R2)) * 'N' + R2
        qual2_r15p = (len(reference_sequence) - len(R2)) * '!' + qual2

        R1_5p_mol = ''
        R1_5p_score = 0

        for a,b, q_a, q_b, REF in zip(R1_r15p, R2_r15p, qual1_r15p, qual2_r15p, reference_sequence):
            if a == b:
                R1_5p_mol += a
            elif a=='N' and b != 'N':
                R1_5p_mol += b
            elif a != 'N' and b =='N':
                R1_5p_mol +=a
            else:
                if ord(q_a) > ord(q_b):
                    R1_5p_mol +=a
                elif ord(q_b) > ord(q_a):
                    R1_5p_mol +=b
                else:
                    R1_5p_mol += random.choice([a, b])
            if R1_5p_mol[-1].replace('C', 'T') == REF:
                R1_5p_score+=1

        R1_r25p = R1 + (len(reference_sequence) - len(R1)) * 'N'
        qual1_r25p = qual1 + (len(reference_sequence) - len(R1)) * '!'
        R2_r25p = (len(reference_sequence) - len(R2)) * 'N' + R2
        qual2_r25p = (len(reference_sequence) - len(R2)) * '!' + qual2
        R2_5p_mol = ''
        R2_5p_score = 0
        for a,b, q_a, q_b, REF in zip(R1_r25p, R2_r25p, qual1_r25p, qual2_r25p, reference_sequence):
            if a == b:
                R2_5p_mol += a
            elif a=='N' and b != 'N':
                R2_5p_mol += b
            elif a != 'N' and b =='N':
                R2_5p_mol +=a
            else:
                if ord(q_a) > ord(q_b):
                    R2_5p_mol +=a
                elif ord(q_b) > ord(q_a):
                    R2_5p_mol +=b
                else:
                    R2_5p_mol += random.choice([a, b])
            if R2_5p_mol[-1].replace('C', 'T') == REF:
                R2_5p_score+=1

        if R1_5p_score >= R2_5p_score:
            molecule = R1_5p_mol
        else:
            molecule = R2_5p_mol

    mismatch_string = ''
    for bp, ref in zip(molecule, reference_sequence):
        if bp == 'N':
            mismatch_string+= 'N'
        elif bp == 'C' and ref == 'C':
            mismatch_string+='C'
        elif bp == 'T' and ref == 'C':
            mismatch_string += 'T'
        else:
            mismatch_string += '.'

    #print('reference sequence:')
    #print(reference_sequence)
    #print(f'length reference sequence: {len(reference_sequence)}')
    #print('mismatch string: ')
    #print(mismatch_string)
    #print(f'length reference sequence: {len(mismatch_string)}')
    #print(f'read positions: {r1_start}, {r1_end}, {r2_start}, {r2_end}')

    return(mismatch_string, molecule, reference_sequence)

def get_reference_seq(sseqid, ref_r1_start, ref_r1_end, ref_r2_start, ref_r2_end, frag_length, fasta, strand, genome_length, record1):
    positions = [ref_r1_start, ref_r1_end, ref_r2_start, ref_r2_end]
    positions.sort()
    if positions[3]-positions[0]+1 == frag_length:
        #reads do not land on end/start of circular genome
        if strand == 'Plus/Plus':
            reference_sequence = fasta[sseqid][(positions[0]-1):(positions[3])].seq
        if strand == 'Plus/Minus':
            reference_sequence = fasta[sseqid][(positions[0]-1):(positions[3])].seq
            reference_sequence = rev_complement(reference_sequence)
    else:
        if strand == 'Plus/Plus':
            reference_sequence = fasta[sseqid][(positions[2]-1):(genome_length)].seq + fasta[sseqid][(0):(positions[1])].seq
        if strand == 'Plus/Minus':
            reference_sequence = fasta[sseqid][(positions[2]-1):(genome_length)].seq + fasta[sseqid][(0):(positions[1])].seq
            reference_sequence = rev_complement(reference_sequence)

    return(reference_sequence)


def mismatch(entry, record1, record2, fasta):
    read_id = entry.split('\t')[1]
    strand = entry.split('\t')[2]
    sseqid = entry.split('\t')[3]

    r1_start = entry.split('\t')[8]
    r1_end = entry.split('\t')[9]
    ref_r1_start = entry.split('\t')[10]
    ref_r1_end = entry.split('\t')[11]

    r2_start = entry.split('\t')[19]
    r2_end = entry.split('\t')[20]
    ref_r2_start = entry.split('\t')[21]
    ref_r2_end = entry.split('\t')[22]
    frag_length = entry.strip().split('\t')[27]
    genome_length = entry.split('\t')[26]



    reference_sequence = get_reference_seq(sseqid, int(ref_r1_start), int(ref_r1_end), int(ref_r2_start), int(ref_r2_end), int(frag_length), fasta, strand, int(genome_length), record1)
    mismatch_string, molecule, ref_seq = get_cytosines(int(r1_start), int(r1_end), int(r2_start), int(r2_end), record1, record2, reference_sequence, strand)
    #print(f'ref postiions: {ref_r1_start}, {ref_r1_end}, {ref_r2_start}, {ref_r2_end}')

    #print(f'fragment length noted in file is {frag_length}')
    #print(f'length of mismatch string is {len(mismatch_string)}')
    #print(f'genome length is {genome_length}')
    #print('+'*50)
    #new_entry = entry.strip('\n') + '\t' + mismatch_string + '\n'
    new_entry = entry.strip('\n') + '\t' + mismatch_string + '\t' + molecule + '\t' + ref_seq + '\n'
    return(new_entry)


def main(r1, r2, tblatpe, output, ref):
    fasta = Fasta(ref)
    with open(tblatpe) as f, open(output, 'w') as w, dnaio.open(r1) as fqr1, dnaio.open(r2) as fqr2:
        w.write('\t'.join(['taxid', 'qseqid', 'strand', 'sseqid', \
                            'pident_R1', 'length_R1', 'mismatch_R1', 'gapopen_R1', 'qstart_R1', \
                            'qend_R1', 'sstart_R1', 'send_R1', 'evalue_R1', 'bitscore_R1', 'qlen_R1', \
                            'pident_R2', 'length_R2', 'mismatch_R2', 'gapopen_R2', 'qstart_R2', \
                            'qend_R2', 'sstart_R2', 'send_R2', 'evalue_R2', 'bitscore_R2', 'qlen_R2', \
                            'genome_len', 'effective_length', 'mismatch_string', 'molecule', 'reference_sequence'])+ '\n')
        for entry, record1, record2 in zip(f, fqr1, fqr2):
            w.write(mismatch(entry, record1, record2, fasta))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", help="r1 r2 paired", nargs='+')
    parser.add_argument("--tblat",help="tblat paired end file")
    parser.add_argument('-o',"--out", help="tblat paired end with CT tags")
    parser.add_argument('--ref', help = 'reference fasta')

    args = parser.parse_args()
    reads = args.reads
    if len(reads) ==2:
        r1 = reads[0]
        r2 = reads[1]
    else:
        r1 = reads
    tblatpe = args.tblat
    output = args.out
    ref = args.ref
    main(r1, r2, tblatpe, output, ref)
