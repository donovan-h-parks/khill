use std::path::Path;

/// Extracts genome identifier from a given sequence file path by removing common file extensions.
pub fn genome_id_from_filename(seq_file: &Path) -> String {
    let mut genome_id = seq_file.file_name().unwrap().to_string_lossy().to_string();

    if genome_id.ends_with(".gz") {
        genome_id = genome_id.replace(".gz", "");
    }

    if genome_id.ends_with(".fq") {
        genome_id = genome_id.replace(".fq", "");
    } else if genome_id.ends_with(".fna") {
        genome_id = genome_id.replace(".fna", "");
    } else if genome_id.ends_with(".fa") {
        genome_id = genome_id.replace(".fa", "");
    } else if genome_id.ends_with(".fasta") {
        genome_id = genome_id.replace(".fasta", "");
    } else if genome_id.ends_with(".fastq") {
        genome_id = genome_id.replace(".fastq", "");
    }

    genome_id
}