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

use crate::ref_sequence::get_ref_seq_from_map;
use clap::{Parser, ValueEnum};
use config::Config;
use std::collections::HashMap;
use std::error::Error;
use std::fmt::Display;
use std::path::Path;
use std::process;

/// Configuration parsed from CLI input arguments for initiating the splicing pipeline.
#[derive(Debug, Clone, Parser)]
///TODO review name/about information below.  Can reimplement version information from ENV here if we start using versioning
#[command(
    name = "viRust-splicing",
    version = env!("CARGO_PKG_VERSION"),
    about = "High-performance 🦀Rust🦀-based tool that processes RNA sequencing data with unique molecular identifier (UMI) to identify and quantify complex HIV splicing patterns."
)]
pub struct InputConfig {
    #[arg(short, long, required = true)]
    pub query: String, // strain name, e.g., NL43
    #[arg(short, long, required = true)]
    pub distance: u8,
    #[arg(short = '1', long = "file1", required = true)]
    pub filename_r1: String,
    #[arg(short = '2', long = "file2", required = true)]
    pub filename_r2: String,
    #[arg(short = 'a', long = "assay", required = true, value_enum)]
    pub assay_type: SpliceAssayType,
    #[arg(long = "output", short = 'o', required = false)]
    pub output_path: Option<String>,
}

///
/// Full configuration for HIV splicing analysis, derived from input strain, query data,
/// and the reference genome FASTA file.
/// TODO: Consider adding steps for d1prime to all acceptors.
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
    pub d1prime_to_all: Vec<SpliceStep>, // TODO: add this step to the config
    pub alternative_d1_to_all: Vec<SpliceStep>, // this is for sequences missing D1, but looking for alternative D1 and acceptors
    pub d2_to_all: Vec<SpliceStep>,
    pub d2b_to_all: Vec<SpliceStep>,
    pub d3_to_all: Vec<SpliceStep>,
    pub full_length_sequence: String,
    pub d4_breakpoint_position: usize,
    pub a7_breakpoint_position: usize,
}

///
/// Path to the TOML configuration file specifying splice form definitions.
pub const SPICE_FORM_CONFIG_STR: &str = include_str!("../data/splice_form_config.toml");

/// Represents a splicing step with donor, acceptor, and the combined sequence pattern.
#[derive(Debug)]
pub struct SpliceStep {
    /// name of the donor site (e.g., "D1")
    pub donor: String,
    /// name of the acceptor site (e.g., "A1")
    pub acceptor: String,
    /// combined sequence pattern, not a name but a actual sequence as a string, eg. "GTAAGTAG"
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

impl Display for SpliceAssayType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let assay_str = match self {
            SpliceAssayType::RandomReverse => "random-reverse",
            SpliceAssayType::KMer => "kmer",
            SpliceAssayType::SizeSpecific => "size-specific",
        };
        write!(f, "{}", assay_str)
    }
}

impl InputConfig {
    pub fn new(
        query: String,
        distance: u8,
        assay_type: SpliceAssayType,
        filename_r1: String,
        filename_r2: String,
        output_path: Option<String>,
    ) -> Self {
        InputConfig {
            query,
            distance,
            assay_type,
            filename_r1,
            filename_r2,
            output_path,
        }
    }
    ///
    /// Parses command-line arguments into an `InputConfig` structure.
    ///
    /// # Errors
    /// Returns an error if the arguments are malformed or invalid.
    pub fn build() -> Result<InputConfig, Box<dyn Error>> {
        let mut input_config: InputConfig = InputConfig::parse();

        if Path::new(&input_config.filename_r1).exists() == false {
            return Err(format!("File {} does not exist", input_config.filename_r1).into());
        }

        if Path::new(&input_config.filename_r2).exists() == false {
            return Err(format!("File {} does not exist", input_config.filename_r2).into());
        }

        // If output path is set, then ensure it's not a file
        //    and create the directory if it does not exist
        // Else set output_path to filename_r1 directory
        if let Some(ref path) = input_config.output_path {
            let path_obj = Path::new(path);

            if path_obj.is_file() {
                return Err(format!("Output path is a file, please use a valid directory").into());
            }

            if !path_obj.exists() {
                std::fs::create_dir_all(path_obj)?;
            }
        } else {
            let r1_path = input_config.filename_r1.clone();
            let default_output_path = Path::new(&r1_path).parent().unwrap();

            let output_path: String = default_output_path
                .to_str()
                .ok_or("Failed to convert output directory path to string")?
                .to_string();

            input_config.output_path = Some(output_path);
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
        let query_list = from_splice_form_file_to_hashmap(&strain.to_lowercase())?;

        // first step check D1
        let d1 = query_list.get("D1").unwrap().clone();

        // this is the most common acceptor sites list
        // there is an intrinsic order in the acceptor sites, so we need to keep the order

        let common_acceptor_sites = vec!["A3", "A4d", "A4c", "A4a", "A4b", "A5a", "A5b", "A7"];

        // second step check D1-unspliced and all the A sites
        // we decide to put D1-unspliced and gag-AUG at the two ends of the list, because normally 70% of the sequences are unspliced.
        // reverse the order of the list may impact the performance, but we will see.

        let mut second_step_acceptor_sites = vec!["D1-unspliced", "gag-AUG", "A1", "A2"];

        second_step_acceptor_sites.append(&mut common_acceptor_sites.clone());

        let second_step = SpliceStep::build("D1", second_step_acceptor_sites, &query_list);

        // D1primer steps

        let mut d1prime_acceptor_list = vec!["D1-unspliced", "A1", "A2"];
        d1prime_acceptor_list.append(&mut common_acceptor_sites.clone());

        let d1prime_step =
            SpliceStep::build_d1prime(d1prime_acceptor_list, &mut query_list.clone())?;

        // Alternative D1 and acceptors, for sequences missing D1

        let mut alternative_d1_acceptor_sites = vec!["A1", "A2"];

        alternative_d1_acceptor_sites.append(&mut common_acceptor_sites.clone());

        // this order of the list matters. The search for pattern will end when the first match is found.
        // So we put gag-AUG at the end. We search all possible acceptors first.
        // If not, we look for intact gag-AUG (unspliced)
        alternative_d1_acceptor_sites.push("gag-AUG");

        let alternative_d1_to_all =
            SpliceStep::build("alternative_d1", alternative_d1_acceptor_sites, &query_list);

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

        let full_length_sequence = get_ref_seq_from_map(&strain.to_ascii_lowercase())
            //.ok_or_else(|| format!("Strain {} not found in config file", strain))?
            .or(get_ref_seq_from_map("nl43")) // default to NL43 if strain not found
            .ok_or("Reference sequence not found")?
            .to_string();

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
            d1prime_to_all: d1prime_step, // TODO: add this step to the config
            alternative_d1_to_all,
            d2_to_all: fourth_step,
            d2b_to_all: fifth_step,
            d3_to_all: seventh_step,
            full_length_sequence,
            d4_breakpoint_position,
            a7_breakpoint_position,
        })
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
        .add_source(config::File::from_str(
            SPICE_FORM_CONFIG_STR,
            config::FileFormat::Toml,
        ))
        .build()?
        .try_deserialize::<HashMap<String, HashMap<String, String>>>()?;
    #[cfg(debug_assertions)]
    dbg!(&splice_form_config);
    let query_list = splice_form_config
        .get(strain)
        .unwrap_or_else(|| {
            eprintln!("Strain {} not found in config file", strain);
            process::exit(1);
        })
        .clone();
    #[cfg(debug_assertions)]
    dbg!(&query_list);
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
        let donor_pattern = query_list.get(donor).unwrap_or(&"".to_string()).to_string();
        for acceptor in acceptor_list {
            let pattern = donor_pattern.clone() + query_list.get(acceptor).unwrap();
            if acceptor == "gag-AUG" {
                // this is a special case, we only look for intact gag-AUG if no other acceptor is found
                // we DO NOT combine D1 + gag-AUG, because it may match unspliced RNA when D1 is not present.

                splice_steps.push(SpliceStep {
                    donor: donor.to_string(),
                    acceptor: acceptor.to_string(),
                    pattern: query_list.get(acceptor).unwrap().to_string(),
                });
                continue;
            }
            splice_steps.push(SpliceStep {
                donor: donor.to_string(),
                acceptor: acceptor.to_string(),
                pattern,
            });
        }
        splice_steps
    }

    pub fn build_d1prime(
        acceptor_list: Vec<&str>,
        query_list: &mut HashMap<String, String>,
    ) -> Result<Vec<SpliceStep>, Box<dyn Error>> {
        let mut d1_sequence = query_list
            .get("D1")
            .cloned()
            .ok_or::<Box<dyn Error>>("Cannot find D1 sequence.".into())?;

        let d1_unspliced = query_list
            .get("D1-unspliced")
            .cloned()
            .ok_or::<Box<dyn Error>>("Cannot find D1-unspliced sequence.".into())?;

        d1_sequence = d1_sequence + &d1_unspliced[..4];

        query_list.insert("D1".to_string(), d1_sequence);
        query_list.insert("D1-unspliced".to_string(), d1_unspliced[4..].to_string());

        Ok(SpliceStep::build("D1", acceptor_list, query_list))
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
            "GTGAGTACGCC"
        );
        dbg!(&splice_config.d1_to_all);
        assert_eq!(
            splice_config.d1_to_all[0].pattern,
            "GGGGCGGCGACTGGTGAGTACGCC"
        );
        assert_eq!(splice_config.d1_to_all[0].donor, "D1");
        assert_eq!(splice_config.d1_to_all[0].acceptor, "D1-unspliced");
        // assert_eq!(
        //     splice_config.d1_to_all[2].pattern,
        //     "GGGGCGGCGACTGGGACAGCAGA"
        // );
        // assert_eq!(splice_config.d1_to_all[2].donor, "D1");
        // assert_eq!(splice_config.d1_to_all[2].acceptor, "A1");
        // assert_eq!(
        //     splice_config.d1_to_all[3].pattern,
        //     "GGGGCGGCGACTGAATCTGCTAT"
        // );

        // assert_eq!(splice_config.d1_to_all[3].donor, "D1");
        // assert_eq!(splice_config.d1_to_all[3].acceptor, "A2");

        // assert_eq!(splice_config.d1_to_all[1].acceptor, "gag-AUG");

        // assert_eq!(splice_config.d1_to_all[1].pattern, "GAGAGATGGGTGC");
        assert_eq!(splice_config.d2_to_all[0].pattern, "CTCTGGAAAGGTGAAGGGGC");
        assert_eq!(splice_config.d2_to_all[0].donor, "D2");
        assert_eq!(splice_config.d2_to_all[0].acceptor, "D2-unspliced");
        assert_eq!(
            splice_config.d2_to_all[1].pattern,
            "CTCTGGAAAGAATCTGCTATAAGAAA"
        );
        assert_eq!(splice_config.d2_to_all[1].donor, "D2");
        assert_eq!(splice_config.d2_to_all[1].acceptor, "A2");
        assert_eq!(splice_config.d4_breakpoint_position, 5301);
        assert_eq!(splice_config.a7_breakpoint_position, 7625);

        dbg!(&splice_config.alternative_d1_to_all);
        dbg!(&splice_config.d1prime_to_all);
        assert_eq!(
            splice_config.alternative_d1_to_all[0].pattern,
            "GGACAGCAGAGATCCA"
        );
        assert_eq!(
            splice_config.alternative_d1_to_all[0].donor,
            "alternative_d1"
        );
        assert_eq!(splice_config.alternative_d1_to_all[0].acceptor, "A1");
        assert_eq!(
            splice_config.alternative_d1_to_all[splice_config.alternative_d1_to_all.len() - 1]
                .acceptor,
            "gag-AUG"
        );
        assert_eq!(
            splice_config.alternative_d1_to_all[splice_config.alternative_d1_to_all.len() - 1]
                .pattern,
            "GAGAGATGGGTGC"
        );
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

        let valid_long_args = InputConfig::try_parse_from([
            "virust-splicing",
            "--query",
            "nl43",
            "--distance",
            "10",
            "--file1",
            "./sim_data/mockseq_r1.fasta",
            "--file2",
            "./sim_data/mockseq_r2.fasta",
            "--assay",
            "random-reverse",
        ]);

        assert!(
            valid_long_args.is_ok(),
            "Expected success, but parsing failed with error"
        );

        let valid_short_args = InputConfig::try_parse_from([
            "virust-splicing",
            "-q",
            "nl43",
            "-d",
            "10",
            "-1",
            "./sim_data/mockseq_r1.fasta",
            "-2",
            "./sim_data/mockseq_r2.fasta",
            "-a",
            "random-reverse",
        ]);

        assert!(
            valid_short_args.is_ok(),
            "Expected success, but parsing failed with error"
        );
    }
}
