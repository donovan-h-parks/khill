
//! This module provides the `FracMinHash` struct and related types for computing
//! FracMinHash sketches from DNA sequences. FracMinHash is a probabilistic data structure
//! used for efficient similarity estimation between large sets, such as k-mer sets from
//! biological sequences. The implementation uses a scale factor to subsample hashes and
//! supports counting both unique and weighted k-mers. The module depends on the `needletail`
//! crate for sequence parsing and a custom hashing implementation for DNA k-mers.
//! 
//! See Hera et al., 2024: https://www.biorxiv.org/content/10.1101/2023.11.06.565843v3

use std::collections::BTreeMap;

use needletail::parser::SequenceRecord;

use crate::hashing::{dna_hashes, ItemHash};

pub type KmerCount = u16;
pub type Hashes = BTreeMap<ItemHash, KmerCount>;

#[derive(Clone, Debug)]
pub struct FracMinHash {
    hashes: Hashes,
    kmer_length: u8,
    max_hash: u64,
    kmer_total_count: u64,
    bp_count: u64,
}

impl FracMinHash {
    pub fn new(kmer_length: u8, scale: u64) -> Self {
        FracMinHash {
            hashes: BTreeMap::new(),
            kmer_length,
            max_hash: ItemHash::MAX / scale,
            kmer_total_count: 0,
            bp_count: 0,
        }
    }

    pub fn process_seq(&mut self, seq: &SequenceRecord) {
        self.bp_count += seq.num_bases() as u64;
        self.kmer_total_count += seq.num_bases() as u64 - self.kmer_length as u64 + 1;

        dna_hashes(
            &seq.seq(),
            &mut self.hashes,
            self.max_hash,
            self.kmer_length,
        );
    }

    pub fn unique_hash_count(&self) -> u64 {
        self.hashes.len() as u64
    }

    pub fn weighted_hash_count(&self) -> u64 {
        self.hashes.values().map(|v| *v as u64).sum()
    }

    pub fn kmer_total_count(&self) -> u64 {
        self.kmer_total_count
    }

    pub fn bp_count(&self) -> u64 {
        self.bp_count
    }

    pub fn to_hashes(self) -> Hashes {
        self.hashes
    }
}
