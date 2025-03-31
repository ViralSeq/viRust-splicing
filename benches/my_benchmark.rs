#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]

use criterion::{criterion_group, criterion_main, Criterion};
use virust_splicing:: open_fasta_file;
use virust_splicing::joined_umi_sequence::pattern_search;
use tap::Pipe;
use bio::io::fasta;
use bio::io::fasta::Record;
use std::error::Error;
use rayon::prelude::*;
use virust_splicing::joined_umi_sequence::JoinedUmiSequnce;
use virust_splicing::config::InputConfig;
use virust_splicing::splice_events::find_umi_family_from_events;
use virust_splicing::config::SpliceConfig;
use virust_splicing::config::SpliceAssayType;
use std::sync::LazyLock;


// typical pattern matching senaria in our dataset
const PATTERN :&[u8] = b"TAGATGTAAAAGGCACCAAGGAAGCCTTAG";
const SEQUENCE :&[u8] = b"CTATTGTGTGCATCAAAGGATAGATGTAAAAGACACCAAGGAAGCCTTAGATAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAGGCACAGCAAGCAGCAGCTGACACAGGAAACAACAGCCAGGTCAGCCAAAATTACCCTATAGTGCAGAACCTCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTAATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAA";



static CONFIG: LazyLock<InputConfig> = LazyLock::new(|| InputConfig {
    query: "nl43".to_string(),
    distance: 1,
    filename_r1: "sim_data/mockseq_r1.fasta".to_string(),
    filename_r2: "sim_data/mockseq_r2.fasta".to_string(),
    assay_type: SpliceAssayType::RandomReverse, // Kmer and SizeSpecific are about 5 times faster than random reverse.
});

fn bench_fasta_reader() {
    let nl43_file = "data/nl43.fasta";
    let fasta_reader = open_fasta_file(nl43_file).unwrap();
    let record = fasta_reader.records().next().unwrap().unwrap();
    record.seq().to_vec().pipe(String::from_utf8).unwrap();
}

fn bench_end_join() {
    let r1 = fasta::Record::with_attrs("seq1_r1", None, b"NNNNTTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGA");
    let r2 = fasta::Record::with_attrs("seq1_r2", None, b"NNNNNNNNNNNNNNAAATCTACTAATTTTCTCCATTTAGTACTGTCTTTTTTCTTTATGGCAAATACTGGAGTATTGTATGGATTTTCAGGCCCAATTTTTGAAATTTTCCCTTC");
    let mut joined_umi_sequence = JoinedUmiSequnce::from_fasta_record(&r1, &r2, 4, 14);
    joined_umi_sequence.join();
}

fn run_rayon() -> Result<usize, Box<dyn Error>> {
    let r1_path = "./temp/R1.fasta";
    let r2_path = "./temp/R2.fasta";

    let fasta_reader_r1 = open_fasta_file(r1_path)?;
    let fasta_reader_r2 = open_fasta_file(r2_path)?;

    let records: Result<Vec<_>, Box<dyn Error>> = fasta_reader_r1
        .records()
        .zip(fasta_reader_r2.records())
        .map(|(r1, r2)| Ok((r1?, r2?)))
        .collect();

    let records = records?;

    let count: Vec<_> = records.par_iter().map(|(r1_record, r2_record)| {
        let mut joined_umi_sequence = JoinedUmiSequnce::from_fasta_record(&r1_record, &r2_record, 4, 11);
        joined_umi_sequence.join();
        joined_umi_sequence.find_sequence_for_search().len()

    }).collect();

    Ok(count.iter().sum())

}

fn run() -> Result<usize, Box<dyn Error>> {
    let r1_path = "./temp/R1.fasta";
    let r2_path = "./temp/R2.fasta";

    let fasta_reader_r1 = open_fasta_file(r1_path)?;
    let fasta_reader_r2 = open_fasta_file(r2_path)?;

    let records: Result<Vec<_>, Box<dyn Error>> = fasta_reader_r1
        .records()
        .zip(fasta_reader_r2.records())
        .map(|(r1, r2)| Ok((r1?, r2?)))
        .collect();

    let records = records?;

    let count: Vec<_> = records.iter().map(|(r1_record, r2_record)| {
        let mut joined_umi_sequence = JoinedUmiSequnce::from_fasta_record(&r1_record, &r2_record, 4, 11);
        joined_umi_sequence.join();
        joined_umi_sequence.find_sequence_for_search().len()
    }).collect();

    Ok(count.iter().sum())
}

fn run_chunk() -> Result<usize, Box<dyn Error>> {
    let r1_path = "./temp/R1.fasta";
    let r2_path = "./temp/R2.fasta";

    let fasta_reader_r1 = open_fasta_file(r1_path)?;
    let fasta_reader_r2 = open_fasta_file(r2_path)?;

    // Get mutable iterators over the FASTA records.
    let mut r1_iter = fasta_reader_r1.records();
    let mut r2_iter = fasta_reader_r2.records();

    // Define a chunk size. Adjust based on memory and performance trade-offs.
    let chunk_size = 1000;
    let mut count = Vec::new();

    loop{
        let mut chunk = Vec::with_capacity(chunk_size);

        for _ in 0..chunk_size {
            match (r1_iter.next(), r2_iter.next()) {
                (Some(Ok(r1_record)), Some(Ok(r2_record))) => {
                    chunk.push((r1_record, r2_record));
                }
                _ => break,
            }
        }

         // Break the outer loop if no records were collected.
         if chunk.is_empty() {
            break;
        }

        // Process the current chunk concurrently.
        let mut c: Vec<_> = chunk.par_iter().map(|(r1_record, r2_record)| {
            let mut joined_umi_sequence = JoinedUmiSequnce::from_fasta_record(&r1_record, &r2_record, 4, 11);
            joined_umi_sequence.join();
            joined_umi_sequence.find_sequence_for_search().len()
        }).collect();

        count.append(&mut c);
    }
    Ok(count.iter().sum())
}



fn integrated_bm(config: InputConfig, records: &Vec<(Record, Record)>) -> Result<(), Box<dyn Error>> {

    let forward_n_size = 4;
    let umi_size = 14;

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

    // dbg!(&splice_config);

    let splice_events: Vec<_> = records.par_iter().map(|(r1_record, r2_record)| {
        let mut joined_umi_sequence = JoinedUmiSequnce::from_fasta_record(
            &r1_record,
            &r2_record,
            forward_n_size,
            umi_size);

        joined_umi_sequence.join();

        joined_umi_sequence.check_splice_event(&splice_config).unwrap() // not sure how to pass the error

    }).collect();


    // find_umi_family_from_events(splice_events);

    Ok(())
}



fn criterion_benchmark(c: &mut Criterion) {

    let r1_file_path = &CONFIG.filename_r1;
    let r2_file_path = &CONFIG.filename_r2;

    let fasta_reader_r1 = open_fasta_file(r1_file_path).unwrap();
    let fasta_reader_r2 = open_fasta_file(r2_file_path).unwrap();

    // collect records
    let records: Result<Vec<_>, Box<dyn Error>> = fasta_reader_r1
        .records()
        .zip(fasta_reader_r2.records())
        .map(|(r1, r2)| Ok((r1?, r2?)))
        .collect();

    let records = records.unwrap();


    // c.bench_function("Pattern Search", |b| b.iter(|| pattern_search(SEQUENCE, PATTERN, 2).unwrap()));
    // c.bench_function("Read file", |b| b.iter(|| bench_fasta_reader()));
    // c.bench_function("End Join", |b| b.iter(|| bench_end_join()));
    c.bench_function("Integrated Benchmark", |b| b.iter(|| integrated_bm((*CONFIG).clone(), &records)));

    // let mut group = c.benchmark_group("expensive reading");
    // group.sample_size(10);
    // group.bench_function("Rayon", |b| b.iter(|| run_rayon().unwrap()));
    // group.bench_function("Sequential", |b| b.iter(|| run().unwrap()));
    // group.bench_function("Chunked", |b| b.iter(|| run_chunk().unwrap()));
}


criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
