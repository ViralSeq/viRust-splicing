use bio::io::fasta::{self, Record};
use rayon::prelude::*;
use splice_events::find_umi_family_from_events;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;

pub mod config;
pub mod joined_umi_sequence;
pub mod ref_sequence;
pub mod splice_events;
pub mod umi;

// pub mod splice_events;

use crate::config::{InputConfig, SpliceConfig};

// TODO: modify error handling
pub fn run(config: InputConfig) -> Result<(), Box<dyn Error>> {
    let forward_n_size = 4; // TODO consider moving to master config
    let umi_size = 14; // TODO consider moving to master config

    let r1_file_path = &config.filename_r1;
    let r2_file_path = &config.filename_r2;

    let fasta_reader_r1 = open_fasta_file(r1_file_path)?;
    let fasta_reader_r2 = open_fasta_file(r2_file_path)?;

    let output_path_str = &config.output_path.clone().unwrap();
    let output_path = Path::new(output_path_str);
    let file_name = "output".to_owned()
        + "_strain_"
        + &config.query
        + "_distance_"
        + &config.distance.to_string()
        + "_assay_"
        + &config.assay_type.to_string()
        + ".tsv";
    let output_tsv_file = output_path.join(file_name);

    // collect records
    // need to add homopolymer check here.
    // records with homopolymer (10 or more) in the either of the paired R1 or R2 will be skipped.
    let records: Vec<(Record, Record)> = fasta_reader_r1
        .records()
        .zip(fasta_reader_r2.records())
        .map(|(r1_res, r2_res)| {
            let r1 = r1_res?;
            let r2 = r2_res?;
            Ok::<_, Box<dyn Error>>((r1, r2))
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;

    let splice_config = SpliceConfig::build_from_input(config)?;

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

    // TODO: efficiency test needed
    let splice_events_with_umi_family = find_umi_family_from_events(splice_events);

    let mut file = File::create(&output_tsv_file)?;
    println!(
        "Writing output to: {}",
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

    Ok(())
}

// TODO: Currently only fasta reader. Consider adding fastq reader.
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
