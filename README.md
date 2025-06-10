# K-Hill

This repository is a Rust implementation of the k-mer Hill statistic proposed in [Narechania et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11529847). It is based on the [K-Hill Python](https://github.com/deanbobo/khill) implementation. This code takes inspiration from the k-mer sketching approaches in [skani](https://github.com/bluenote-1577/skani) and [Finch](https://github.com/onecodex/finch-rs).

K-Hill is a method and software package that can quantify molecular diversity in pangenomes. The method draws on information theory (Shannon Diversity) to quantify richness and evenness of K-mers between genomes in groups of samples. The approach is computationally efficient - not relying on databases, alignments, or genome graphs.

# Running K-Hill

K-Hill can be run two ways:
1. with --directory which will measure K-Hill across all genomic FASTA files in the specified directory
2. with --file_table which will measure K-Hill across specified groups of genomes

The `--file_table` input file should be a tab separated values (TSV) file with two columns indicating the group of each genome and the path to a genomic FASTA file. For example:

```
# group_id  fasta_file_path
groupA  /path/to/genome1.fna
groupA  /path/to/genome2.fna
groupB  /path/to/genome3.fna
groupB  /path/to/genome4.fna
```
