"""
This script extracts unique mutated reads from NGS deep amplicon sequencing
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


def extract_mutated_reads(bam_file: Path,
                          ref_fasta: Path,
                          output_bam: Path,
                          start_pos: int,
                          end_pos: int,
                          min_mapping_qual: int,
                          max_indel_length: int,
                          threads: int):
    """Function exctract unique reads with mutation in specific positions"""
    indel_pattern = re.compile(
        f'[.,][+-][0-9]+[atgcnATGCN]{{{max_indel_length},}}')
    unique_sequences = list()
    total_reads = 0
    with pysam.AlignmentFile(bam_file, 'rb', threads=threads) as bam, \
        pysam.AlignmentFile(output_bam, 'wb',
                            reference_lengths=[get_ref_len(ref_fasta)],
                            reference_names=[get_ref_name(ref_fasta)]) as mutated_bam:
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
                if ((mut != '.') &
                        (mut != ',') &
                        (re.search(indel_pattern, str(mut)) is None) &
                        (not read.alignment.is_duplicate) &
                        (read.alignment.mapping_quality > min_mapping_qual) &
                        (not read.alignment.is_supplementary) &
                        (read.alignment.is_mapped)):
                    if read.alignment.query_sequence in unique_sequences:
                        continue
                    else:
                        mutated_bam.write(read.alignment)
                        total_reads += 1
                        unique_sequences.append(read.alignment.query_sequence)
    print(f'{total_reads} reads were recorded in {output_bam}')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('bam_file',
                        help='path to sorted bam file',
                        type=parser_resolve_path)
    parser.add_argument('ref_fasta',
                        help='path to reference fasta file (single record)',
                        type=parser_resolve_path)
    parser.add_argument('output_bam',
                        help='path to output bam file for mutated reads',
                        type=parser_resolve_path)
    parser.add_argument('start_pos',
                        help='start position for analisys (0-based)',
                        type=int)
    parser.add_argument('end_pos',
                        help='end position for analisys (0-based)',
                        type=int)

    parser.add_argument('--min_map_qual', metavar='',
                        help='minimal mapping quality for read to consider (default: 20)',
                        type=int, default=20)
    parser.add_argument('--max_indel_len', metavar='',
                        help='maximum length of indel to filter out (default: 40)',
                        type=int, default=40)
    parser.add_argument('--threads', metavar='',
                        help='number of threads (default: 1)',
                        type=int, default=1)

    args = parser.parse_args()

    pysam.index(f'{args.bam_file}')
    extract_mutated_reads(bam_file=args.bam_file,
                          ref_fasta=args.ref_fasta,
                          output_bam=args.output_bam,
                          start_pos=args.start_pos,
                          end_pos=args.end_pos,
                          min_mapping_qual=args.min_map_qual,
                          max_indel_length=args.max_indel_len,
                          threads=args.threads)
    pysam.sort('-o', f'{args.output_bam}.sorted', f'{args.output_bam}')
    os.remove(f'{args.output_bam}')
    os.rename(f'{args.output_bam}.sorted', f'{args.output_bam}')
