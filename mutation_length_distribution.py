"""
Script takes bam file containing mutated reads and count number of mutation of
each length from 1 to >=10 regardless to mutation type.
"""

import re
from pathlib import Path

import pysam
from Bio import SeqIO


def parser_resolve_path(path):
    return Path(path).resolve()


def get_ref_name(ref_fasta: Path) -> str:
    ref = SeqIO.read(ref_fasta, 'fasta')
    return ref.id


def extract_mutations_info(bam_file: Path,
                           ref_fasta: Path,
                           start_pos: int,
                           end_pos: int,
                           threads: int,
                           header: bool):
    lengths = []
    with pysam.AlignmentFile(bam_file, 'rb', threads=threads) as bam:
        for column in bam.pileup(contig=get_ref_name(ref_fasta=ref_fasta),
                                 start=start_pos,
                                 stop=end_pos,
                                 fastafile=pysam.FastaFile(ref_fasta),
                                 max_depth=100000000,
                                 truncate=True):
            for pos in column.get_query_sequences(mark_matches=True,
                                                  add_indels=True):
                if (pos == '.') | (pos == ','):
                    continue
                elif re.search(r'^[ACGTNacgtn]$', str(pos)):
                    lengths.append(1)
                elif re.search(r'[.,][+-][0-9]+[atgcnATGCN]+', str(pos)):
                    lengths.append(
                        int(re.search(r'[0-9]+', str(pos)).group(0)))
    strict_lengths = [lengths.count(i) for i in range(1, 10)]
    length_gt10 = str(len(lengths) - sum(strict_lengths))
    if not header:
        header = '\t'.join([str(i) for i in range(1, 10)])
        print(f"{header}\t>=10")
    print('\t'.join([str(i) for i in strict_lengths] + [length_gt10]))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('bam_file',
                        help='path to sorted bam file containing unique mutated reads (result of extract_mutated_reads.py)',
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
    parser.add_argument('--header',
                        help='add header to results output',
                        action='store_false')

    args = parser.parse_args()

    pysam.index(f'{args.bam_file}')
    extract_mutations_info(bam_file=args.bam_file,
                           ref_fasta=args.ref_fasta,
                           start_pos=args.start_pos,
                           end_pos=args.end_pos,
                           threads=args.threads,
                           header=args.header)
