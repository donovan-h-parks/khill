# K-Hill

This repository is a Rust implementation of the k-mer Hill statistic proposed in [Narechania et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11529847). It is based on the [K-Hill Python](https://github.com/deanbobo/khill) implementation though does not produce identical results as a difference hashing function
is used. This code takes inspiration from the k-mer sketching approaches in [skani](https://github.com/bluenote-1577/skani) and [Finch](https://github.com/onecodex/finch-rs).

K-Hill is a method and software package that can quantify molecular diversity in pangenomes. The method draws on information theory (Shannon Diversity) to quantify richness and evenness of K-mers between genomes in groups of samples. The approach is computationally efficient - not relying on databases, alignments, or genome graphs.

# Running K-Hill

K-Hill can be run two ways:
1. with --input-dir which will measure K-Hill across all genomic FASTA files (*.fa, *.fna, *.fasta) in the specified directory
2. with --genome-group-table which will measure K-Hill across specified groups of genomes

The `genome-group-table` input file should be a tab separated values (TSV) file with two columns indicating the group of each genome and the path to a genomic FASTA file. For example:

```
# group_id  fasta_file_path
groupA  /path/to/genome1.fna
groupA  /path/to/genome2.fna
groupB  /path/to/genome3.fna
groupB  /path/to/genome4.fna
```

K-Hill benefits substantially from using multiple threads which can be specified with the `--threads` flag. By default, k-hill runs with a k-mer length (--kmer_length) of 19 and a scaling factor (--scale) of 100 (i.e. k-hill is applied to sketches containing ~1% of all k-mers). 

# Install

## Building K-Hill from Source

Ensure you have [Rust](https://www.rust-lang.org/tools/install) installed. Then, build the project with:

```!sh
git clone https://github.com/your-username/khill.git
cd khill
cargo build --release
```

The compiled binary will be located at `target/release/khill`.

## Pre-build Executable

A pre-built executable is provided for x86-64 Linux systems which can be obtained using:

```!sh
wget https://github.com/donovan-h-parks/khill/releases/download/<version>/khill
chmod +x khill
./khill -h
```

__Important__: the pre-built executable can run substantially slower (2 to 3x) than building a native executable
