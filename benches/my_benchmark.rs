use criterion::{criterion_group, criterion_main, Criterion};
use virust_splicing:: open_fasta_file;
use virust_splicing::joined_umi_sequence::pattern_search;
use tap::Pipe;
use bio::io::fasta;
use virust_splicing::joined_umi_sequence::JoinedUmiSequnce;


// typical pattern matching senaria in our dataset
const PATTERN :&[u8] = b"TAGATGTAAAAGGCACCAAGGAAGCCTTAG";
const SEQUENCE :&[u8] = b"CTATTGTGTGCATCAAAGGATAGATGTAAAAGACACCAAGGAAGCCTTAGATAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAGGCACAGCAAGCAGCAGCTGACACAGGAAACAACAGCCAGGTCAGCCAAAATTACCCTATAGTGCAGAACCTCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTAATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAA";



fn bench_fasta_reader() {
    let nl43_file = "tests/nl43.fasta";
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

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("Pattern Search", |b| b.iter(|| pattern_search(SEQUENCE, PATTERN, 2).unwrap()));
    c.bench_function("Read file", |b| b.iter(|| bench_fasta_reader()));
    c.bench_function("End Join", |b| b.iter(|| bench_end_join()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);