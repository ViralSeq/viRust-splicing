use bio::io::fasta::Record;
use bio::io::fasta::Writer;
use csv;
use rand::Rng;
use rand::prelude::IndexedRandom;
use regex::Regex;
use std::error::Error;
use std::{collections::HashMap, env::args};
use virust_splicing::joined_umi_sequence::reverse_complement;

use virust_splicing::config::from_splice_form_file_to_hashmap;
use virust_splicing::open_fasta_file;
use virust_splicing::umi::simulate_umi_distribution;

fn main() {
    let number_to_generate = args().nth(1).unwrap().parse::<usize>().unwrap();

    // read csv file from /data/splice_forms_mock_ratio.csv
    let mut rdr = csv::Reader::from_path("data/splice_forms_mock_ratio.csv").unwrap();
    let mut ratio_hash: HashMap<String, f32> = HashMap::new();
    for result in rdr.records() {
        let record = result.unwrap();
        let key = record.get(0).unwrap().to_string();
        let value = record.get(1).unwrap().parse::<f32>().unwrap();
        ratio_hash.insert(key, value);
    }

    let weighted_seq_hash = ratio_hash
        .iter()
        .map(|(k, v)| (k.clone(), (*v * number_to_generate as f32).round() as i32))
        .collect::<HashMap<String, i32>>();

    let mut writer_r1 = Writer::to_file("sim_data/mockseq_r1.fasta").unwrap();
    let mut writer_r2 = Writer::to_file("sim_data/mockseq_r2.fasta").unwrap();

    // query_list is a hashmap of query name and query sequence for all D and A sequences
    let query_list = from_splice_form_file_to_hashmap("nl43").unwrap();
    generate_mock_sequence(&weighted_seq_hash, &query_list)
        .iter()
        .for_each(|record| {
            writer_r1.write_record(&record.0).unwrap();
            writer_r2.write_record(&record.1).unwrap();
        });
}

fn generate_mock_sequence(
    weighted_seq_hash: &HashMap<String, i32>,
    query_list: &HashMap<String, String>,
) -> Vec<(Record, Record)> {
    let mut mock_sequences: Vec<(Record, Record)> = Vec::new();
    let mut seq_id = 1;
    let mut rng = rand::rng();
    let bases = vec!['A', 'C', 'G', 'T'];

    let mut umi_hash: HashMap<String, Vec<String>> = HashMap::new();

    for (splice_event, seq_number_for_this_event) in weighted_seq_hash.iter() {
        let mut target_umi_numbers = if *seq_number_for_this_event < 10 {
            1
        } else {
            *seq_number_for_this_event / 10
        };

        target_umi_numbers = target_umi_numbers.min(100);

        let umi_dist = simulate_umi_distribution(
            *seq_number_for_this_event as u32,
            14,
            target_umi_numbers as u32,
            10,
            0.005,
        )
        .unwrap();

        let mut umi_list = Vec::new();
        for (umi, count) in umi_dist.iter() {
            for _i in 0..*count {
                umi_list.push(umi.clone());
            }
        }

        umi_hash.insert(splice_event.clone(), umi_list);
    }

    for (splice_event, v) in weighted_seq_hash.iter() {
        let umi_list = umi_hash.get(splice_event).unwrap();
        let mut r1_sequence = (0..4)
            .map(|_| {
                let idx = rng.random_range(0..bases.len());
                bases[idx]
            })
            .collect::<String>();

        r1_sequence += "A".repeat(20).as_str();

        let (r1_string, r2_sequence, size_tag) =
            from_splice_event_to_mock_r1_r2(splice_event, query_list).unwrap();

        r1_sequence += &r1_string;

        r1_sequence = fill_out_sequence(r1_sequence, "A", 300);

        for i in 0..*v {
            let seq_name_r1 = format!("mock_seq_{}_{}_{}_r1", splice_event, size_tag, seq_id);
            let seq_name_r2 = format!("mock_seq_{}_{}_{}_r2", splice_event, size_tag, seq_id);

            let r2_sequence_with_umi = umi_list[i as usize].clone() + &r2_sequence;
            mock_sequences.push((
                Record::with_attrs(&seq_name_r1, None, r1_sequence.as_bytes()),
                Record::with_attrs(&seq_name_r2, None, r2_sequence_with_umi.as_bytes()),
            ));
            seq_id += 1;
        }
    }

    mock_sequences
}

fn from_splice_event_to_mock_r1_r2(
    splice_event: &str,
    query_list: &HashMap<String, String>,
) -> Result<(String, String, String), Box<dyn Error>> {
    let mut mock_seq_r1 = String::new();

    let mut events = splice_event.split("_");

    let first_event = events.next().ok_or("No first part")?;

    if first_event == "noD1" {
        return Ok(("A".repeat(276), "G".repeat(286), "Unknown".to_string()));
    } else if first_event != "D1" {
        return Err("First event Not a D1 event".into());
    } else {
        mock_seq_r1 += query_list.get(first_event).unwrap();
    }

    let re = Regex::new(r"^(?:A.*|.*un.*)$").unwrap();

    for event in events {
        let pre_seq = if re.is_match(event) {
            "".to_string()
        } else {
            "A".repeat(20)
        };
        mock_seq_r1 += &pre_seq;
        let default_seq = "A".repeat(10);
        let event_seq = query_list.get(event).unwrap_or(&default_seq);
        mock_seq_r1 += event_seq;
    }
    let mock_seq_r2 = generate_r2_mock_seq();

    Ok((
        mock_seq_r1,
        reverse_complement(&mock_seq_r2.0),
        mock_seq_r2.1,
    ))
}

fn generate_r2_mock_seq() -> (String, String) {
    let nl43_file = "data/nl43.fasta";
    let fasta_reader = open_fasta_file(nl43_file).unwrap();
    let record = fasta_reader.records().next().unwrap().unwrap();
    let nl43_ref_seq = record.seq().to_vec();

    let seq_pool: Vec<(Vec<u8>, String)> = vec![
        (nl43_ref_seq[5000..5286].to_vec(), "Both".to_string()), // before D4
        (nl43_ref_seq[7000..7286].to_vec(), "4kb".to_string()),  // between D4 and A7
        (nl43_ref_seq[8000..8286].to_vec(), "1.8kb".to_string()), // after A7
        ("G".repeat(286).as_bytes().to_vec(), "Unknown".to_string()), // unknown
    ];

    let selected = seq_pool.choose(&mut rand::rng()).unwrap().clone();

    (String::from_utf8(selected.0).unwrap(), selected.1)
}

fn fill_out_sequence(seq: String, fill_base: &str, target_length: usize) -> String {
    let mut filled_seq = seq.clone();
    let current_length = filled_seq.len();
    if current_length < target_length {
        let fill_length = target_length - current_length;
        filled_seq += &fill_base.repeat(fill_length);
    }
    filled_seq
}
