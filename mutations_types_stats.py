"""
Script extracts information about types of mutation in unique mutated reads
that cover specific position
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
                           threads: int):
    insertion, deletion, snp = 0, 0, 0
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
                elif re.search(r'^[ACGTNacgtn]$', str(pos)):  # SNP
                    snp += 1
                elif re.search(r'\+[0-9]+[ACGTNacgtn]+', str(pos)):
                    insertion += 1
                elif re.search(r'-[0-9]+[ACGTNacgtn]+', str(pos)):
                    deletion += 1
    print(f's\ti\td\n{snp}\t{insertion}\t{deletion}')


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

    args = parser.parse_args()

    pysam.index(f'{args.bam_file}')
    extract_mutations_info(bam_file=args.bam_file,
                           ref_fasta=args.ref_fasta,
                           start_pos=args.start_pos,
                           end_pos=args.end_pos,
                           threads=args.threads)
