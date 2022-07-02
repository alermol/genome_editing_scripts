# Genome editing

Repository contains scripts for analysis of NGS data acquired from samples after genome editing

### ***extract_mutated_reads.py***
Script takes bam file and create another bam file containing reads that 
covers define positions of reference and somehow differ from reference 
(mutated reads) because SNP, insertion or deletion

```
usage: extract_mutated_reads.py [-h] [--min_map_qual] [--max_indel_len] 
[--threads] [--unique] bam_file ref_fasta output_bam start_pos end_pos

positional arguments:
  bam_file          path to sorted bam file
  ref_fasta         path to reference fasta file (single record)
  output_bam        path to output bam file for mutated reads
  start_pos         start position for analisys (0-based)
  end_pos           end position for analisys (0-based)

optional arguments:
  -h, --help        show this help message and exit
  --min_map_qual    minimal mapping quality for read to consider (default: 20)
  --max_indel_len   maximum length of indel to filter out (default: 40)
  --threads         number of threads (default: 1)
  --unique          extract only reads with unique sequence (default: False)

```
