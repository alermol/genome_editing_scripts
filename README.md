# Genome editing

Repository contains scripts for analysis of NGS data acquired from samples after genome editing.

### ***extract_mutated_reads.py***
Script takes bam file and create another bam file containing reads that covers define positions of reference and somehow differ from reference in this positions (mutated reads) because of SNP, insertion or deletion. Optionally only unique by sequence reads could be remained in resulted bam file.

```None
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
  --unique          extract only reads with unique sequences
```

**Following scripts should be used for bam file containing mutated reads only which is output of this script**

### ***split_bam_by_mutation_type.py***

Script takes baf file containing mutated reads, classified reads by mutation type (SNP, insertion or deletion) considering mutation only in defined positions range and write a bam file for each type of mutation.

```None
usage: split_bam_by_mutation_type.py [-h] [--threads] bam_file ref_fasta 
       start_pos end_pos

positional arguments:
  bam_file    path to sorted bam file
  ref_fasta   path to reference fasta file (single record)
  start_pos   start position for analisys (0-based)
  end_pos     end position for analisys (0-based)

optional arguments:
  -h, --help  show this help message and exit
  --threads   number of threads (default: 1)

```
### ***mutations_types_stats.py***

Script takes bam file containing mutated reads and count number of mutations for each type (SNP, insertion and deletion) inside defined positions range.

```None
usage: mutations_types_stats.py [-h] [--threads] bam_file ref_fasta start_pos 
       end_pos

positional arguments:
  bam_file    path to sorted bam file containing unique mutated reads (result of extract_mutated_reads.py)
  ref_fasta   path to reference fasta file (single record)
  start_pos   start position for analisys (0-based)
  end_pos     end position for analisys (0-based)

optional arguments:
  -h, --help  show this help message and exit
  --threads   number of threads (default: 1)
```

### ***mutation_length_distribution.py***

Script takes bam file containing mutated reads and count number of mutation of each length from 1 to >=10 regardless to mutation type.

```None
usage: mutation_length_distribution.py [-h] [--threads] bam_file ref_fasta 
       start_pos end_pos

positional arguments:
  bam_file    path to sorted bam file containing unique mutated reads (result of extract_mutated_reads.py)
  ref_fasta   path to reference fasta file (single record)
  start_pos   start position for analisys (0-based)
  end_pos     end position for analisys (0-based)

optional arguments:
  -h, --help  show this help message and exit
  --threads   number of threads (default: 1)
```
