//! # config.rs
//!
//! This module defines the configuration structures and logic for parsing input parameters
//! and building splice analysis configurations for HIV splicing assays.
//! It includes:
//! - `InputConfig`: Parses command-line arguments.
//! - `SpliceConfig`: Constructs a full splicing configuration from a config file and fasta data.
//! - `SpliceStep`: Represents a splice step pattern between donor and acceptor sites.
//! - `SpliceAssayType`: Enumerates types of splice assays.
//!
//! The configuration is based on a strain-specific query list loaded from a TOML file,
//! and includes logic to build all relevant splicing steps for downstream analysis.
//!
//! This module also includes test cases to verify correct behavior of configuration logic.

use crate::open_fasta_file;
use clap::{Parser, ValueEnum};
use config::Config;
use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use std::process;
use tap::Pipe;

///
/// Configuration parsed from CLI input arguments for initiating the splicing pipeline.
///TODO: Add output path
#[derive(Debug, Clone, Parser)]
///TODO review name/about information below.  Can reimplement version information from ENV here if we start using versioning
#[command(
    name = "viRust-splicing",
    version = "1.0",
    about = "High-performance 🦀Rust🦀-based tool that processes RNA sequencing data with unique molecular identifier (UMI) to identify and quantify complex HIV splicing patterns."
)]
pub struct InputConfig {
    #[arg(long, required = true)]
    pub query: String,
    #[arg(long, required = true)]
    pub distance: u8,
    #[arg(long, required = true)]
    pub filename_r1: String,
    #[arg(long, required = true)]
    pub filename_r2: String,
    #[arg(long, required = true, value_enum)]
    pub assay_type: SpliceAssayType,
    // pub output_path: String,
}

///
/// Full configuration for HIV splicing analysis, derived from input strain, query data,
/// and the reference genome FASTA file.
#[derive(Debug)]
pub struct SpliceConfig {
    pub strain: String,
    pub splice_assay_type: SpliceAssayType,
    pub distance: u8,
    pub query_list: HashMap<String, String>,
    pub d1: String,
    pub d2: String,
    pub d2b: String,
    pub d2b_unspliced: String,
    pub d3: String,
    pub d1_to_all: Vec<SpliceStep>,
    pub d2_to_all: Vec<SpliceStep>,
    pub d2b_to_all: Vec<SpliceStep>,
    pub d3_to_all: Vec<SpliceStep>,
    pub full_length_sequence: String,
    pub d4_breakpoint_position: usize,
    pub a7_breakpoint_position: usize,
}

///
/// Path to the TOML configuration file specifying splice form definitions.
pub const SPICE_FORM_CONFIG: &str = "data/splice_form_config.toml";

///
/// Represents a splicing step with donor, acceptor, and the combined sequence pattern.
#[derive(Debug)]
pub struct SpliceStep {
    pub donor: String,
    pub acceptor: String,
    pub pattern: String,
}

///
/// Enum for the type of assay used in splicing analysis.
#[derive(Debug, PartialEq, Clone, ValueEnum)]
pub enum SpliceAssayType {
    RandomReverse,
    KMer,
    SizeSpecific,
}

impl InputConfig {
    pub fn new(
        query: String,
        distance: u8,
        assay_type: SpliceAssayType,
        filename_r1: String,
        filename_r2: String,
    ) -> Self {
        InputConfig {
            query,
            distance,
            assay_type,
            filename_r1,
            filename_r2,
        }
    }
    ///
    /// Parses command-line arguments into an `InputConfig` structure.
    /// TODO: Output config to a file
    ///
    /// # Errors
    /// Returns an error if the arguments are malformed or invalid.
    pub fn build() -> Result<InputConfig, Box<dyn Error>> {
        let input_config: InputConfig = InputConfig::parse();

        if Path::new(&input_config.filename_r1).exists() == false {
            return Err(format!("File {} does not exist", input_config.filename_r1).into());
        }

        if Path::new(&input_config.filename_r2).exists() == false {
            return Err(format!("File {} does not exist", input_config.filename_r2).into());
        }

        Ok(input_config)
    }
}

impl SpliceConfig {
    ///
    /// Builds a `SpliceConfig` from a given strain, distance, and assay type.
    ///
    /// This function loads the splice form config and reference FASTA file,
    /// and prepares splicing steps accordingly.
    ///
    /// # Errors
    /// Returns an error if the config or FASTA file cannot be read.
    pub fn build(
        strain: String,
        distance: u8,
        splice_assay_type: SpliceAssayType,
    ) -> Result<SpliceConfig, Box<dyn Error>> {
        let query_list = from_splice_form_file_to_hashmap(&strain)?;

        // first step check D1
        let d1 = query_list.get("D1").unwrap().clone();

        // this is the most common acceptor sites list
        // there is an intrinsic order in the acceptor sites, so we need to keep the order

        let common_acceptor_sites = vec!["A3", "A4d", "A4c", "A4a", "A4b", "A5a", "A5b", "A7"];

        // second step check D1-unspliced and all the A sites

        let mut second_step_acceptor_sites = vec!["D1-unspliced", "A1", "A2"];

        second_step_acceptor_sites.append(&mut common_acceptor_sites.clone());

        let second_step = SpliceStep::build("D1", second_step_acceptor_sites, &query_list);

        // this is the third step, when A1 is found, check D2

        let d2 = query_list.get("D2").unwrap().clone();

        // this is parellel step to step 3, when A2 is found, check D3, followed by D3-unspliced and A3 to A7

        // this is the fourth step, when D2 is found, check D2-unspliced and A2 to A7

        let mut fouth_step_acceptor_sites = vec!["D2-unspliced", "A2"];

        fouth_step_acceptor_sites.append(&mut common_acceptor_sites.clone());

        let fourth_step = SpliceStep::build("D2", fouth_step_acceptor_sites, &query_list);

        // this is parellel step to step 4, when D2 is absent, check D2b

        let d2b = query_list.get("D2b").unwrap().clone();

        // this is step 5, only when D2-unspliced is found, check D2b-unspliced

        let d2b_unspliced = query_list.get("D2b-unspliced").unwrap().clone();

        // this is parellel step to step 5, when D2b is found is absent, check D2b-unspliced and A2 to A7

        let mut fifth_step_acceptor_sites = vec!["D2b-unspliced", "A2"];

        fifth_step_acceptor_sites.append(&mut common_acceptor_sites.clone());

        let fifth_step = SpliceStep::build("D2b", fifth_step_acceptor_sites, &query_list);

        // this is step 6, check D3, same as the parellel step to step 6

        let d3 = query_list.get("D3").unwrap().clone();

        // this is step 7, check D3-unspliced and A3 to A7, it will be used in several steps, will be called d3_unspliced_a3_a7
        // same as the parellel step to step 7

        let mut d3_unspliced_a3_a7 = vec!["D3-unspliced"];

        d3_unspliced_a3_a7.append(&mut common_acceptor_sites.clone());

        let seventh_step = SpliceStep::build("D3", d3_unspliced_a3_a7, &query_list);

        let strain_file_name = format!("data/{}.fasta", strain.to_lowercase());

        let fasta_reader = open_fasta_file(&strain_file_name)?;

        if let Some(record) = fasta_reader.records().next() {
            let record = record?;
            let full_length_sequence = record.seq().to_vec().pipe(String::from_utf8).unwrap();
            let d4_breakpoint_position = query_list
                .get("D4_breakpoint_position")
                .unwrap()
                .parse::<usize>()
                .unwrap();
            let a7_breakpoint_position = query_list
                .get("A7_breakpoint_position")
                .unwrap()
                .parse::<usize>()
                .unwrap();

            Ok(SpliceConfig {
                strain,
                splice_assay_type,
                distance,
                query_list,
                d1,
                d2,
                d2b,
                d2b_unspliced,
                d3,
                d1_to_all: second_step,
                d2_to_all: fourth_step,
                d2b_to_all: fifth_step,
                d3_to_all: seventh_step,
                full_length_sequence,
                d4_breakpoint_position,
                a7_breakpoint_position,
            })
        } else {
            Err("No record found in fasta file".into())
        }
    }

    ///
    /// Convenience method to build `SpliceConfig` from an `InputConfig`.
    pub fn build_from_input(input_config: InputConfig) -> Result<SpliceConfig, Box<dyn Error>> {
        let strain = input_config.query;
        let distance = input_config.distance;
        let splice_assay_type = input_config.assay_type;
        let splice_config = SpliceConfig::build(strain, distance, splice_assay_type)?;
        Ok(splice_config)
    }
}

pub fn from_splice_form_file_to_hashmap(
    strain: &str,
) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let splice_form_config: HashMap<String, HashMap<String, String>> = Config::builder()
        .add_source(config::File::with_name(SPICE_FORM_CONFIG))
        .build()?
        .try_deserialize::<HashMap<String, HashMap<String, String>>>()?;

    let query_list = splice_form_config
        .get(strain)
        .unwrap_or_else(|| {
            eprintln!("Strain {} not found in config file", strain);
            process::exit(1);
        })
        .clone();

    Ok(query_list)
}

impl SpliceStep {
    ///
    /// Constructs a list of `SpliceStep` entries from a donor and multiple acceptors.
    ///
    /// # Arguments
    /// * `donor` - The donor splice site key.
    /// * `acceptor_list` - A list of acceptor splice site keys.
    /// * `query_list` - A map from splice site keys to sequence strings.
    ///
    /// # Returns
    /// A vector of `SpliceStep` structs, each containing a donor-acceptor pattern.
    pub fn build(
        donor: &str,
        acceptor_list: Vec<&str>,
        query_list: &HashMap<String, String>,
    ) -> Vec<SpliceStep> {
        let mut splice_steps: Vec<SpliceStep> = Vec::new();
        let donor_pattern = query_list.get(donor).unwrap().to_string();
        for acceptor in acceptor_list {
            let pattern = donor_pattern.clone() + query_list.get(acceptor).unwrap();
            splice_steps.push(SpliceStep {
                donor: donor.to_string(),
                acceptor: acceptor.to_string(),
                pattern,
            });
        }
        splice_steps
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_splice_config() {
        let strain = "NL43".to_lowercase().to_string();
        let distance = 10;
        let splice_config =
            SpliceConfig::build(strain, distance, SpliceAssayType::RandomReverse).unwrap();
        assert_eq!(splice_config.strain, "nl43");
        assert_eq!(splice_config.distance, 10);
        assert_eq!(splice_config.query_list.get("D1").unwrap(), "GGGGCGGCGACTG");
        assert_eq!(
            splice_config.query_list.get("D1-unspliced").unwrap(),
            "GTGAGTACGCCAAAAA"
        );
        assert_eq!(
            splice_config.d1_to_all[0].pattern,
            "GGGGCGGCGACTGGTGAGTACGCCAAAAA"
        );
        assert_eq!(splice_config.d1_to_all[0].donor, "D1");
        assert_eq!(splice_config.d1_to_all[0].acceptor, "D1-unspliced");
        assert_eq!(
            splice_config.d1_to_all[1].pattern,
            "GGGGCGGCGACTGGGACAGCAGA"
        );
        assert_eq!(splice_config.d1_to_all[1].donor, "D1");
        assert_eq!(splice_config.d1_to_all[1].acceptor, "A1");
        assert_eq!(
            splice_config.d1_to_all[2].pattern,
            "GGGGCGGCGACTGAATCTGCTAT"
        );
        assert_eq!(splice_config.d1_to_all[2].donor, "D1");
        assert_eq!(splice_config.d1_to_all[2].acceptor, "A2");
        assert_eq!(splice_config.d2_to_all[0].pattern, "CTCTGGAAAGGTGAAGGGGC");
        assert_eq!(splice_config.d2_to_all[0].donor, "D2");
        assert_eq!(splice_config.d2_to_all[0].acceptor, "D2-unspliced");
        assert_eq!(splice_config.d2_to_all[1].pattern, "CTCTGGAAAGAATCTGCTAT");
        assert_eq!(splice_config.d2_to_all[1].donor, "D2");
        assert_eq!(splice_config.d2_to_all[1].acceptor, "A2");
        assert_eq!(splice_config.d4_breakpoint_position, 5301);
        assert_eq!(splice_config.a7_breakpoint_position, 7625);
    }

    #[test]
    fn test_build_splice_step() {
        let list = vec!["A1", "A2", "A3"];
        let mut query_list = HashMap::new();
        query_list.insert("D1".to_string(), "ACGT".to_string());
        query_list.insert("A1".to_string(), "TGCA".to_string());
        query_list.insert("A2".to_string(), "ACGT".to_string());
        query_list.insert("A3".to_string(), "CCCC".to_string());
        let splice_steps = SpliceStep::build("D1", list, &query_list);
        assert_eq!(splice_steps[0].pattern, "ACGTTGCA".to_string());
        assert_eq!(splice_steps[0].donor, "D1".to_string());
        assert_eq!(splice_steps[0].acceptor, "A1".to_string());
        assert_eq!(splice_steps[1].pattern, "ACGTACGT".to_string());
        assert_eq!(splice_steps[1].donor, "D1".to_string());
        assert_eq!(splice_steps[1].acceptor, "A2".to_string());
        assert_eq!(splice_steps[2].pattern, "ACGTCCCC".to_string());
        assert_eq!(splice_steps[2].donor, "D1".to_string());
        assert_eq!(splice_steps[2].acceptor, "A3".to_string());
    }

    #[test]
    fn test_args() {
        let invalid_reponse_missing_args = InputConfig::try_parse_from([
            "virust-splicing",
            "--query",
            "nl43",
            "--assay-type",
            "reverse",
        ]);

        assert!(
            invalid_reponse_missing_args.is_err(),
            "Expected an error, but parsing succeeded"
        );

        let valid_args = InputConfig::try_parse_from([
            "virust-splicing",
            "--query",
            "nl43",
            "--distance",
            "10",
            "--filename-r1",
            "./sim_data/mockseq_r1.fasta",
            "--filename-r2",
            "./sim_data/mockseq_r2.fasta",
            "--assay-type",
            "random-reverse",
        ]);

        if let Err(e) = valid_args {
            eprintln!("Error parsing arguments: {}", e);
            panic!("Failed to parse valid arguments");
        }

        assert!(
            valid_args.is_ok(),
            "Expected success, but parsing failed with error"
        );
    }
}
