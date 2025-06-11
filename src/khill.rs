//! This module implements the K-Hill method for calculating beta entropy (K-hill number) across a set of genomes.
//!
//! It provides functionality to:
//! - Sketch genome sequences into k-mer hashes in parallel.
//! - Aggregate k-mer counts across genomes.
//! - Compute the K-Hill number and per-genome KL-divergence and weights.
//!
//! The main entry point is the `khill` function, which returns the K-Hill number and detailed entropy components for each genome.

use std::fs::File;
use std::path::PathBuf;
use std::collections::{BTreeMap, HashMap};

use anyhow::{Context, Result};
use needletail::parse_fastx_reader;
use rayon::iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator};

use crate::hashing::ItemHash;
use crate::sketch_params::SketchParams;
use crate::frac_min_hash::Hashes;
use crate::io_utils::genome_id_from_filename;

/// Hill components.
#[derive(Clone, Debug, PartialEq)]
pub struct HillComponent {
    pub kl_divergence: f64,
    pub weight: f64,
}

/// Calculate beta entropy using the K-Hill method.
pub fn khill(genome_files: &Vec<PathBuf>, sketch_params: &SketchParams) -> Result<(f64, HashMap<String, HillComponent>)> {
    // calculate hashes for all genomes in parallel
    let genome_hashes: HashMap<String, Hashes> = genome_files
        .par_iter()
        .map(|genome_file| {
            let genome_id = genome_id_from_filename(genome_file);
            let hashes = sketch_file(genome_file, sketch_params).expect("Failed to create sketch");
            (genome_id, hashes)
        })
        .collect();

    // determine k-mers across all genomes in parallel using map-reduce
    let all_kmers = genome_hashes.values()
        .par_bridge()
        .map(|hashes| {
            // Create a local BTreeMap for each thread
            let mut local_map = BTreeMap::<ItemHash, u64>::new();
            for (hash, count) in hashes {
                *local_map.entry(*hash).or_insert(0) += *count as u64;
            }
            local_map
        })
        .reduce(
            // Initial empty map
            BTreeMap::<ItemHash, u64>::new,
            // Combine two maps
            |mut acc, map| {
                for (hash, count) in map {
                    *acc.entry(hash).or_insert(0) += count;
                }
                acc
            }
        );

    

    // calculate the K-hill number in parallel
    let total_num_hashes: u64 = all_kmers.values().sum();
    
    let genome_results: HashMap<String, HillComponent> = genome_hashes
        .par_iter()
        .map(|(genome_id, hashes)| {
            let num_genome_hashes: u64 = hashes.values().map(|&v| v as u64).sum();
            
            let kl_divergence = hashes.iter()
                .filter_map(|(hash, count)| {
                    all_kmers.get(hash).map(|total_count| {
                        let p_si = *count as f64 / num_genome_hashes as f64;
                        let p_i = *total_count as f64 / total_num_hashes as f64;
                        p_si * (p_si/p_i).ln()
                    })
                })
                .sum();

            let weight = num_genome_hashes as f64 / total_num_hashes as f64;
            (genome_id.clone(), HillComponent { kl_divergence, weight })
        })
        .collect();

    // calculate final K-hill value
    let khill = genome_results.values()
        .map(|comp| comp.weight * comp.kl_divergence)
        .sum::<f64>();

    Ok((khill.exp(), genome_results))
}

/// Create sketch from sequence file.
pub fn sketch_file(seq_file: &PathBuf, sketch_params: &SketchParams) -> Result<Hashes> {
    let mut sketcher = sketch_params.create_sketcher();
    let reader = File::open(seq_file)
        .context(format!("Failed to open {}", seq_file.display()))?;

    let mut fastx_reader = parse_fastx_reader(reader)?;
    while let Some(rec) = fastx_reader.next() {
        let record = rec?;
        sketcher.process_seq(&record);
    }

    Ok(sketcher.to_hashes())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::{tempdir, TempDir};

    // Helper to create a temporary FASTA file with given contents
    fn write_temp_fasta(contents: &str, filename: &str, dir: &TempDir) -> PathBuf {
        let file_path = dir.path().join(filename);
        let mut file = File::create(&file_path).unwrap();
        file.write_all(contents.as_bytes()).unwrap();
        file.sync_all().unwrap();
        file_path
    }

    #[test]
    fn test_khill_with_two_simple_genomes() {
        // Create a temporary directory that will stay alive for the whole test
        let temp_dir = tempdir().unwrap();

        // Create two simple FASTA files in the same directory
        let fasta1 = ">seq1\nACGTACGTACGT\n";
        let fasta2 = ">seq2\nACGTACGTACGA\n";
        let file1 = write_temp_fasta(fasta1, "genome1.fa", &temp_dir);
        let file2 = write_temp_fasta(fasta2, "genome2.fa", &temp_dir);

        // Use default sketch params
        let sketch_params = SketchParams::new(3, 1, true);

        let genome_files = vec![file1.clone(), file2.clone()];
        let result = khill(&genome_files, &sketch_params);

        assert!(result.is_ok());
        let (khill_value, genome_entropy) = result.unwrap();

        // Check for expected K-hill value
        assert!(khill_value == 1.0376237334557157);
        
        // There should be two entries in genome_entropy
        assert_eq!(genome_entropy.len(), 2);

        // Each genome should have a weight > 0 and KL divergence >= 0
        for comp in genome_entropy.values() {
            assert!(comp.weight > 0.0);
            assert!(comp.kl_divergence >= 0.0);
        }

        // temp_dir is automatically cleaned up when it goes out of scope
    }
}
