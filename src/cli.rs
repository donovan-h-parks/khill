
//! Command-line interface definition for the khill application.
//!
//! This file defines the `Cli` struct using the `clap` crate to parse and validate command-line arguments.
//! It includes options for specifying input directories or genome group tables, output directory, k-mer length,
//! sketch scaling factor, and number of threads. Custom value parsers are provided for k-mer length and thread count.
//! The CLI output is styled using the `anstyle` crate for improved readability.

use std::path::PathBuf;

use clap::Parser;

const DEFAULT_K: u8 = 19;
const DEFAULT_SCALE: u64 = 100;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(styles=get_styles())]
#[command(disable_help_subcommand = true)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    /// Directory of genomes to process
    #[arg(short = 'i', long, help_heading = "Inputs", group= "input", value_parser = clap::value_parser!(PathBuf))]
    pub input_dir: Option<PathBuf>,

    /// TSV file indicating groups of genomes to process (group_id, path to FASTA file)
    #[arg(short = 'g', long, help_heading = "Inputs", group = "input", value_parser = clap::value_parser!(PathBuf))]
    pub genome_group_table: Option<PathBuf>,

    /// Output directory
    #[arg(short = 'o', long, help_heading = "Output", value_parser = clap::value_parser!(PathBuf))]
    pub out_dir: PathBuf,

    /// Length of k-mers to use
    #[arg(short, long, help_heading = "Sketching parameters", default_value_t = DEFAULT_K, value_parser = validate_kmer_length)]
    pub kmer_length: u8,

    /// Sketch scaling factor (e.g. 100 will examine ~1% of k-mers)
    #[arg(short = 's', long, help_heading = "Sketching parameters", default_value_t = DEFAULT_SCALE)]
    pub scale: u64,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 1, value_parser = validate_threads)]
    pub threads: usize,

    /// Skip verification that genomic FASTA files exist
    #[arg(long, default_value_t = false)]
    pub skip_file_check: bool,
}

fn validate_kmer_length(k: &str) -> Result<u8, String> {
    let k: u8 = k
        .parse()
        .map_err(|_| format!("`{k}` isn't a valid k-mer length"))?;

    if !(1..=32).contains(&k) {
        return Err("k-mer length must be in the range [1, 32]".to_string());
    }

    Ok(k)
}

fn validate_threads(threads: &str) -> Result<usize, String> {
    let threads: usize = threads
        .parse()
        .map_err(|_| format!("`{threads}` isn't a valid value"))?;

    if !(1..=1024).contains(&threads) {
        return Err("Threads  must be in the range [1, 1024]".to_string());
    }

    Ok(threads)
}

fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
        .literal(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .invalid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .error(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .valid(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .placeholder(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
}

#[test]
fn test_verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
