use bio::io::fasta;
use rayon::prelude::*;
use runner::*;
use splice_events::find_umi_family_from_events;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;

pub mod config;
pub mod io;
pub mod joined_umi_sequence;
pub mod ref_sequence;
pub mod runner;
pub mod splice_events;
pub mod umi;

use crate::config::{InputConfig, SpliceConfig};
use crate::io::*;

/// Main function to run the viRust-Splicing analysis workflow.
/// It takes an InputConfig object as input and performs the following steps:
/// 1. Validates the input FASTA files.
/// 2. Opens the FASTA files and reads the sequences.
/// 3. Prepares the output file path and name based on the configuration.
/// 4. Builds the SpliceConfig from the InputConfig.
/// 5. Processes each pair of sequences in parallel to identify splice events,
///   filtering out sequences with homopolymers.
/// 6. Finds UMI families from the identified splice events.
/// 7. Writes the results to an output TSV file.
/// 8. Checks for R installation and required packages, and summarizes the data using R.
/// 9. Checks for Quarto and Python3 installation, and generates an HTML report using Quarto.
/// The function returns a Result indicating success or failure of the entire workflow.
pub fn run(config: InputConfig) -> Result<(), Box<dyn Error>> {
    let forward_n_size = 4; // TODO consider moving to master config
    let umi_size = 14; // TODO consider moving to master config

    let r1_file_path = &config.filename_r1;
    let r2_file_path = &config.filename_r2;

    let sequence_files = validate_input_files(r1_file_path, r2_file_path)?;
    let records = open_sequence_files(&sequence_files)?;

    let output_path_str = &config.output_path.clone().unwrap();
    let output_path = Path::new(output_path_str);
    let file_base_name = "output".to_owned()
        + "_strain_"
        + &config.query
        + "_distance_"
        + &config.distance.to_string()
        + "_assay_"
        + &config.assay_type.to_string();
    let file_name = file_base_name.clone() + ".tsv";
    let output_tsv_file = output_path.join(file_name);

    let splice_config = SpliceConfig::build_from_input(config.clone())?;

    #[cfg(debug_assertions)]
    dbg!(&splice_config);

    let splice_events: Vec<_> = records
        .par_iter()
        .filter(|(r1, r2)| !(homopolymer_check(r1.seq()) || homopolymer_check(r2.seq())))
        .map(|(r1_record, r2_record)| {
            let mut joined_umi_sequence = joined_umi_sequence::JoinedUmiSequnce::from_fasta_record(
                &r1_record,
                &r2_record,
                forward_n_size,
                umi_size,
            );

            joined_umi_sequence.join();

            joined_umi_sequence
                .check_splice_event(&splice_config)
                .unwrap() // not sure how to pass the error
        })
        .collect();

    #[cfg(debug_assertions)]
    dbg!(println!("{:#?}", splice_events));

    let splice_events_with_umi_family = find_umi_family_from_events(splice_events);

    let mut file = File::create(&output_tsv_file)?;
    println!(
        "\tWriting output to: {}",
        output_tsv_file
            .to_str()
            .unwrap_or("Error converting path to string...")
    );
    writeln!(
        file,
        "sequence_id\tumi\tumi_family\tsplice_category\tsize_class\tfinal_category\talternative_d1_used\tunknown_sequence_after_d1"
    )?;
    for res in splice_events_with_umi_family.iter() {
        writeln!(file, "{}", res.to_string())?;
    }

    if check_r_installed().is_err() {
        println!("R is not installed or not found in PATH. Skipping data summarization with R.");
        return Ok(());
    }

    if check_r_packages().is_err() {
        println!("Required R packages are not installed. Skipping data summarization with R.");
        return Ok(());
    }

    let output_summary_file = output_path.join(file_base_name.clone() + "_summary.csv");

    if let Err(e) = r_summarize_data(
        output_tsv_file.to_str().unwrap(),
        output_summary_file.to_str().unwrap(),
    ) {
        println!("Data summarization failed: {}", e);
        return Ok(());
    } else {
        println!(
            "✅ Data summarization completed. \n\tSummary file created at: {}",
            output_summary_file
                .to_str()
                .unwrap_or("Error converting path to string...")
        );
    }

    let check_quarto = check_quarto_installed();
    let check_python3 = check_python3_installed();

    if check_quarto.is_err() || check_python3.is_err() {
        println!(
            "Quarto cannot be executed. Skipping HTML report generation.\nQuarto check: {:?}\nPython3 check: {:?}",
            check_quarto, check_python3
        );
        return Ok(());
    }

    if let Err(e) = run_quarto_report(
        &config.assay_type.to_string(),
        &config.query,
        &config.distance,
        r1_file_path,
        r2_file_path,
        &(file_base_name.clone() + "_summary.csv"),
        output_path,
    ) {
        println!("Quarto report generation failed: {}", e);
    } else {
        println!(
            "\tHTML report generated at: {}",
            output_path
                .join("report.html")
                .to_str()
                .unwrap_or("Error converting path to string...")
        );
    }

    Ok(())
}

/// Opens a FASTA file and returns a FASTA reader.
/// # Arguments
/// * `file_path` - A string slice that holds the path to the FASTA file
/// # Returns
/// * `Result<fasta::Reader<BufReader<BufReader<File>>>, Box<dyn Error>>` - A result containing the FASTA reader or an error
pub fn open_fasta_file(
    file_path: &str,
) -> Result<fasta::Reader<BufReader<BufReader<File>>>, Box<dyn Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    Ok(fasta::Reader::new(reader))
}

fn homopolymer_check(seq: &[u8]) -> bool {
    let max_homopolymer_length = 10; // TODO consider moving to master config
    let mut current_char = 0;
    let mut current_count = 0;

    for &c in seq {
        if c == current_char {
            current_count += 1;
            if current_count > max_homopolymer_length {
                return true; // Found a homopolymer longer than the threshold
            }
        } else {
            current_char = c;
            current_count = 1;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use joined_umi_sequence::pattern_search;

    #[test]
    fn test_pattern_search() {
        let pattern = b"TTTTCCTAGGATATGGCTCCATAACTTAGGACAA";
        let nl43_file = "data/nl43.fasta";
        let fasta_reader = open_fasta_file(nl43_file).unwrap();
        let record = fasta_reader.records().next().unwrap().unwrap();
        let sequence = record.seq().to_vec();
        let aln = pattern_search(&sequence, pattern, 2).unwrap();
        assert_eq!((aln.ystart, aln.yend, aln.score), (4913, 4947, 0));
    }

    #[test]
    fn test_homopolymer_check() {
        let seq_with_homopolymer = b"AAACCCCTTTTTTTTTTTTGGG";
        let seq_without_homopolymer = b"AAACCCCTTTTGGG";
        assert!(homopolymer_check(seq_with_homopolymer));
        assert!(!homopolymer_check(seq_without_homopolymer));
    }
}
