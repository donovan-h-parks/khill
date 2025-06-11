
//! Main entry point for the khill application.
//!
//! This file handles command-line parsing, logging setup, input validation, and orchestrates
//! the computation of k-hill statistics and genome entropy for groups of genomic FASTA files.
//! It supports input via a directory of FASTA files or a TSV file specifying genome groups.
//! Results are written to output files in the specified directory.

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use ahash::{HashMap, HashMapExt};
use anyhow::{Result};
use clap::Parser;
use log::info;

use crate::cli::Cli;
use crate::logging::setup_logger;
use crate::khill::khill;
use crate::progress::progress_bar;
use crate::sketch_params::SketchParams;

mod cli;
pub mod logging;
pub mod progress;
pub mod khill;
pub mod sketch_params;
pub mod frac_min_hash;
pub mod hashing;
pub mod io_utils;

/// Common initialization required by all commands.
fn init(threads: usize) -> Result<()> {
    const VERSION: &str = env!("CARGO_PKG_VERSION");
    info!("{} v{}", env!("CARGO_PKG_NAME"), VERSION);
    info!("{}", env::args().collect::<Vec<String>>().join(" "));

    info!("Using {} threads.", threads);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;

    Ok(())
}

/// Parse a TSV file containing genome groups and the path to their genomic FASTA files.
fn parse_genome_groups_file(file_path: &PathBuf) -> Result<HashMap<String, Vec<PathBuf>>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    
    // process each lines
    let mut groups: HashMap<String, Vec<PathBuf>> = HashMap::new();
    for line in reader.lines() {
        let line = line?;

        // skip comment lines starting with #
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 2 {
            return Err(anyhow::anyhow!("Invalid TSV format: each line must have exactly 2 columns (group_id, path)"));
        }
        
        let group = fields[0].to_string();
        let path = PathBuf::from(fields[1]);

        groups.entry(group)
            .or_default()
            .push(path);
    }
    
    Ok(groups)
}

fn main() -> Result<()> {
    let start = Instant::now();

    let args = Cli::parse();

    setup_logger(&args.out_dir)?;

    init(args.threads)?;

    // determine if input is being specified via a directory or a file table

    let groups = if let Some(genome_group_table) = args.genome_group_table {
        info!("Using genome group file: {}", genome_group_table.display());
        parse_genome_groups_file(&genome_group_table)?
    } else if let Some(input_dir) = args.input_dir {
        info!("Using input directory: {}", input_dir.display());

        // If a directory is specified, can it for FASTA files.
        let paths: Vec<PathBuf> = std::fs::read_dir(input_dir)?
            .filter_map(Result::ok)
            .filter(|entry| entry.path().extension().is_some_and(|ext| ext == "fa" || ext == "fasta" || ext == "fna"))
            .map(|entry| entry.path())
            .collect();

        if paths.is_empty() {
            return Err(anyhow::anyhow!("No FASTA files found in specified directory."));
        }
        
        let mut groups = HashMap::new();
        groups.insert("default".to_string(), paths);
        groups
    } else {
        return Err(anyhow::anyhow!("No input specified. Use --input_dir or --genome_group_table."));
    };

    // check that all genomic FASTA files exist
    if !args.skip_file_check {
        info!("Verifying all genomic FASTA files exist.");
        let num_genomes = groups.values().map(|v| v.len()).sum::<usize>();
        let progress_bar = progress_bar(num_genomes as u64);
        for (group, genome_paths) in &groups {
            for path in genome_paths {
                if !path.exists() {
                    return Err(anyhow::anyhow!("Genome file {} in group '{}' does not exist.", path.display(), group));
                }
                progress_bar.inc(1);
            }
        }
        progress_bar.finish();
    }

    // open output file for group k-hill and per genome entropy results
    std::fs::create_dir_all(&args.out_dir)?;
    let khill_out_file = File::create(args.out_dir.join("khill.tsv"))?;
    let mut khill_writer = BufWriter::new(khill_out_file);
    writeln!(khill_writer, "group_id\tnum_genomes\tk-hill")?;

    let genome_entropy_out_file = File::create(args.out_dir.join("genome_entropy.tsv"))?;
    let mut genome_entropy_writer = BufWriter::new(genome_entropy_out_file);
    writeln!(genome_entropy_writer, "genome_id\tbeta_entropy\tkl_divergence\tweight")?;

    // process each group of genomes
    let sketch_params = SketchParams::new(args.kmer_length, args.scale, true);
    info!("Processing {} genome groups:", groups.len());
    for (group, genome_paths) in &groups {
        let (k_hill, genome_stats) = khill(genome_paths, &sketch_params)?;
        info!(" - {}: {} genomes; {:.2} effective genomes", group, genome_paths.len(), k_hill);
        writeln!(khill_writer, "{}\t{}\t{}", group, genome_paths.len(), k_hill)?;

        for (genome_id, components) in genome_stats {
            writeln!(genome_entropy_writer, "{}\t{}\t{}\t{}", 
            genome_id, 
            components.weight * components.kl_divergence,
            components.kl_divergence, 
            components.weight)?;
        }
    }

    info!("Elapsed time (sec): {:.2}", start.elapsed().as_secs_f32());
    info!("Done.");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::write;
    use tempfile::NamedTempFile;
    
    #[test]
    fn test_parse_genome_groups_file() -> Result<()> {
        // Create a temporary test file
        let temp_file = NamedTempFile::new()?;
        let test_content = "# group_id\tpath\n\
                           group1\t/path/to/genome1.fna\n\
                           group1\t/path/to/genome2.fna\n\
                           group2\t/path/to/genome3.fna";
        write(temp_file.path(), test_content)?;
        
        let groups = parse_genome_groups_file(&temp_file.path().to_path_buf())?;
        
        assert_eq!(groups.len(), 2);
        assert_eq!(groups["group1"].len(), 2);
        assert_eq!(groups["group2"].len(), 1);
        assert_eq!(groups["group1"][0], PathBuf::from("/path/to/genome1.fna"));
        
        Ok(())
    }
}