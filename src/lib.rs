use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use bio::io::fasta;
use std::io::Write;
use rayon::prelude::*;
use splice_events::find_umi_family_from_events;

pub mod config;
pub mod joined_umi_sequence;
pub mod splice_events;
pub mod umi;

// pub mod splice_events;

use crate::config::{InputConfig, SpliceConfig};

// TODO: modify error handling
pub fn run(config: InputConfig) -> Result<(), Box<dyn Error>> {

    let forward_n_size = 4; // consider moving to master config
    let umi_size = 14; // consider moving to master config

    let r1_file_path = &config.filename_r1;
    let r2_file_path = &config.filename_r2;

    let fasta_reader_r1 = open_fasta_file(r1_file_path)?;
    let fasta_reader_r2 = open_fasta_file(r2_file_path)?;

    // collect records
    let records: Result<Vec<_>, Box<dyn Error>> = fasta_reader_r1
        .records()
        .zip(fasta_reader_r2.records())
        .map(|(r1, r2)| Ok((r1?, r2?)))
        .collect();

    let records = records?;

    let splice_config = SpliceConfig::build_from_input(config)?;

    let splice_events: Vec<_> = records.par_iter().map(|(r1_record, r2_record)| {
        let mut joined_umi_sequence = joined_umi_sequence::JoinedUmiSequnce::from_fasta_record(
            &r1_record,
            &r2_record,
            forward_n_size,
            umi_size);

        joined_umi_sequence.join();

        joined_umi_sequence.check_splice_event(&splice_config).unwrap() // not sure how to pass the error

    }).collect();

    // TODO: efficiency test needed
    let splice_events_with_umi_family = find_umi_family_from_events(splice_events);

    let mut file = File::create("output.tsv")?;
    writeln!(file,"sequence_id\tumi\tumi_family\tsplice_category\tsize_class\tfinal_category")?;
    for res in splice_events_with_umi_family.iter() {
        writeln!(file, "{}", res.to_string())?;
    }


    Ok(())
}

pub fn open_fasta_file(file_path: &str) -> Result<fasta::Reader<BufReader<BufReader<File>>>, Box<dyn Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    Ok(fasta::Reader::new(reader))
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
}
