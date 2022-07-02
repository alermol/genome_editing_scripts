"""
Script takes baf file containing mutated reads, classified reads by mutation
type (SNP, insertion or deletion) considering mutation only in defined positions
range and write a bam file for each type of mutation.
"""

import os
import re
from pathlib import Path

import pysam
from Bio import SeqIO


def parser_resolve_path(path):
    return Path(path).resolve()


def get_ref_name(ref_fasta: Path) -> str:
    ref = SeqIO.read(ref_fasta, 'fasta')
    return ref.id


def get_ref_len(ref_fasta: Path) -> int:
    ref = SeqIO.read(ref_fasta, 'fasta')
    return len(ref.seq)


def split_reads(bam_file: Path,
                ref_fasta: Path,
                start_pos: int,
                end_pos: int,
                threads: int):
    """Function exctract unique reads with mutation in specific positions"""
    total_reads, snp_reads, ins_reads, del_reads = 0, 0, 0, 0
    with pysam.AlignmentFile(bam_file, 'rb', threads=threads) as bam, \
        pysam.AlignmentFile(f'{bam_file}.SNP', 'wb',
                            reference_lengths=[get_ref_len(ref_fasta)],
                            reference_names=[get_ref_name(ref_fasta)]) as snp_bam, \
        pysam.AlignmentFile(f'{bam_file}.ins', 'wb',
                            reference_lengths=[get_ref_len(ref_fasta)],
                            reference_names=[get_ref_name(ref_fasta)]) as ins_bam, \
        pysam.AlignmentFile(f'{bam_file}.del', 'wb',
                            reference_lengths=[get_ref_len(ref_fasta)],
                            reference_names=[get_ref_name(ref_fasta)]) as del_bam:
        for column in bam.pileup(contig=get_ref_name(ref_fasta=ref_fasta),
                                 start=start_pos,
                                 stop=end_pos,
                                 fastafile=pysam.FastaFile(ref_fasta),
                                 max_depth=100000000,
                                 truncate=True):
            reads = column.pileups
            mutations = column.get_query_sequences(mark_matches=True,
                                                   add_indels=True)
            for read, mut in zip(reads, mutations):
                if re.search(r'^[ACGTNacgtn]$', str(mut)):
                    snp_bam.write(read.alignment)
                    snp_reads += 1
                    total_reads += 1
                elif re.search(r'\+[0-9]+[ACGTNacgtn]+', str(mut)):
                    ins_bam.write(read.alignment)
                    ins_reads += 1
                    total_reads += 1
                elif re.search(r'-[0-9]+[ACGTNacgtn]+', str(mut)):
                    del_bam.write(read.alignment)
                    del_reads += 1
                    total_reads += 1
    print(f'{total_reads} were recorded of these:')
    print(f'  {snp_reads} containing SNPs')
    print(f'  {ins_reads} containing insertions')
    print(f'  {del_reads} containing deletions')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('bam_file',
                        help='path to sorted bam file',
                        type=parser_resolve_path)
    parser.add_argument('ref_fasta',
                        help='path to reference fasta file (single record)',
                        type=parser_resolve_path)
    parser.add_argument('start_pos',
                        help='start position for analisys (0-based)',
                        type=int)
    parser.add_argument('end_pos',
                        help='end position for analisys (0-based)',
                        type=int)

    parser.add_argument('--threads', metavar='',
                        help='number of threads (default: 1)',
                        type=int, default=1)

    args = parser.parse_args()

    pysam.index(f'{args.bam_file}')
    split_reads(bam_file=args.bam_file,
                ref_fasta=args.ref_fasta,
                start_pos=args.start_pos,
                end_pos=args.end_pos,
                threads=args.threads)
    for file in [f'{args.bam_file}.SNP',
                 f'{args.bam_file}.ins',
                 f'{args.bam_file}.del']:
        pysam.sort('-o', f'{file}.sorted', f'{file}')
        os.remove(f'{file}')
        os.rename(f'{file}.sorted', f'{file}.bam'.replace('.bam', '', 1))
