//! # joined_umi_sequence.rs
//!
//! This module provides utilities for parsing, joining, and analyzing paired-end FASTA reads with Unique Molecular Identifiers (UMIs)
//! in the context of HIV splicing analysis. It includes functionality to:
//! - Extract metadata and sequences from forward and reverse reads.
//! - Merge paired-end reads with configurable overlap detection.
//! - Represent joined reads and convert them into FASTA records.
//! - Detect splicing events through recursive search guided by a splicing configuration.
//!
//! The core structure in this module is `JoinedUmiSequence`, which stores parsed read data and facilitates downstream
//! processing and splice event detection.

use crate::config::{SpliceConfig, SpliceStep};
use crate::splice_events::{ProcessedSpliceRec, SpliceChain, SpliceEvents};
use bio::alignment::Alignment;
use bio::alphabets::dna;
///MARK: Library imports
use bio::io::fasta;
use bio::pattern_matching::myers;
use serde::{Deserialize, Serialize};
use std::error::Error;
use tap::Pipe;

///MARK: JoinedUmiSequence
/// Represents a joined UMI sequence derived from paired-end FASTA records.
///
/// This struct contains metadata and sequences extracted from the forward and reverse reads,
/// as well as the joined sequence if the reads are successfully merged.
#[derive(Debug, Deserialize, Serialize)]
pub struct JoinedUmiSequnce {
    /// The sequence ID derived from the forward read.
    pub sequence_id: String,
    /// The first few bases of the forward read, used as a unique identifier.
    pub forward_ns: String,
    /// The UMI (Unique Molecular Identifier) extracted from the reverse read.
    pub umi: String,
    /// The remaining sequence of the forward read after the forward Ns.
    pub forward_sequence: String,
    /// The remaining sequence of the reverse read after the UMI.
    pub reverse_sequence: String,
    /// The joined sequence obtained by merging the forward and reverse reads.
    pub joined_sequence: Option<String>,
}

impl JoinedUmiSequnce {
    /// Creates a `JoinedUmiSequnce` from paired-end FASTA records.
    ///
    /// # Arguments
    /// * `record_r1` - The forward read record.
    /// * `record_r2` - The reverse read record.
    /// * `forward_n_size` - The number of bases to extract as forward Ns.
    /// * `umi_size` - The number of bases to extract as the UMI.
    ///
    /// # Returns
    /// A new `JoinedUmiSequnce` instance.
    pub fn from_fasta_record(
        record_r1: &fasta::Record,
        record_r2: &fasta::Record,
        forward_n_size: usize,
        umi_size: usize,
    ) -> JoinedUmiSequnce {
        let sequence_id = record_r1.id().to_string();

        let forward_seq = record_r1
            .seq()
            .to_vec()
            .pipe(String::from_utf8)
            .unwrap()
            .to_uppercase();
        let reverse_seq = record_r2
            .seq()
            .to_vec()
            .pipe(String::from_utf8)
            .unwrap()
            .to_uppercase();

        let forward_ns = forward_seq[..forward_n_size].to_string();
        let forward_sequence = forward_seq[forward_n_size..].to_string();
        let umi = reverse_seq[..umi_size].to_string();
        let reverse_sequence = reverse_seq[umi_size..].to_string();
        let joined_sequence = None;

        JoinedUmiSequnce {
            sequence_id,
            forward_ns,
            umi,
            forward_sequence,
            reverse_sequence,
            joined_sequence,
        }
    }

    /// Attempts to join the forward and reverse sequences into a single contiguous sequence.
    ///
    /// This method uses a minimum overlap and error rate to determine if the reads can be merged.
    pub fn join(&mut self) {
        let min_overlap = 10; // TODO Minimum overlap length required for merging. Consider moving to master config.
        let error_rate = 0.02; // TODO Maximum error rate allowed in the overlap region. Consider moving to master config.

        if let Some(joined) = join_reads(
            &self.forward_sequence,
            &self.reverse_sequence,
            min_overlap,
            error_rate,
        ) {
            self.joined_sequence = Some(joined);
        }
    }

    /// Converts the joined sequence into a FASTA record.
    ///
    /// # Returns
    /// An `Option<fasta::Record>` containing the joined sequence, or `None` if the sequence is not joined.
    pub fn joined_seq_to_fasta_record(&self) -> Option<fasta::Record> {
        let id = format!("{}-joined", self.sequence_id);
        let desc = Some(format!("ForwardNs: {}, UMI: {}", self.forward_ns, self.umi));
        if let None = self.joined_sequence {
            return None;
        } else {
            let seq = self.joined_sequence.as_ref().unwrap();
            Some(fasta::Record::with_attrs(
                &id,
                desc.as_deref(),
                seq.as_bytes(),
            ))
        }
    }

    /// Finds the sequence to use for searching splice events.
    ///
    /// If the joined sequence exists, it is returned. Otherwise, the forward and reverse sequences are concatenated.
    ///
    /// # Returns
    /// A `String` containing the sequence for search.
    pub fn find_sequence_for_search(&self) -> String {
        match &self.joined_sequence {
            Some(joined_sequence) => joined_sequence.to_string(),
            None => self.forward_sequence.to_string() + &reverse_complement(&self.reverse_sequence),
        }
    }

    /// MARK: Splice event check
    /// Checks for splice events in the sequence using the provided splice configuration.
    ///
    /// # Arguments
    /// * `splice_config` - The configuration for splice event detection.
    ///
    /// # Returns
    /// A `Result` containing the detected `SpliceEvents` or an error.
    pub fn check_splice_event(
        &self,
        splice_config: &SpliceConfig,
    ) -> Result<SpliceEvents, Box<dyn Error>> {
        let seq = self.find_sequence_for_search();
        let seq = seq.as_bytes();
        let mut chain = SpliceChain::new();
        // Start processing at stage 1 (which handles D1).
        let initial_splice_rec = ProcessedSpliceRec::init(seq);
        let final_processed_splice_rec =
            process_splice_rec(initial_splice_rec, &mut chain, splice_config, 1);
        let final_seq = final_processed_splice_rec.seq;
        let alternative_d1 = final_processed_splice_rec.alternative_d1;
        let unknown_sequence_after_d1 = final_processed_splice_rec.unknown_sequence_after_d1;
        let mut event = SpliceEvents::from_joined_umi_with_event(self, chain);
        event.add_post_splice_sequence(String::from_utf8(final_seq.to_vec())?);
        event.add_alternative_d1(alternative_d1);
        event.add_unknown_sequence_after_d1(unknown_sequence_after_d1);
        event.find_size_class(splice_config)?;
        // println!("Splice events: {:?}", event.size_class);
        event.predict_final_category();
        Ok(event)
    }
}

/// MARK: Utility functions
/// Processes splicing steps recursively based on the stage.
///
/// The stages are defined as:
/// - **1:** Search for D1.
/// - **2:** Find acceptor downstream of D1.
///   - If acceptor is `"A1"`, continue to stage 3.
///   - If acceptor is `"A2"`, jump to stage 5.
/// - **3:** Search for D2 (A1 branch).
/// - **4:** Find acceptor downstream of D2.
///   - If acceptor is `"D2-unspliced"`, continue to stage 6.
///   - If acceptor is `"A2"`, jump to stage 5.
/// - **5:** Process the common D3 branch (for both A2 outcomes).
/// - **6:** Search for D2b (in the D2-unspliced branch).
///   - If an acceptor from D2b is found and it is `"A2"`, jump to stage 5.
/// - **7:** Look for an acceptor after D3.
///
/// # Returns
/// The remaining sequence after processing.
/// TODO！Fix all the issues and complete the TODOs in this function.
fn process_splice_rec<'a>(
    splice_rec: ProcessedSpliceRec<'a>,
    chain: &mut SpliceChain,
    config: &SpliceConfig,
    stage: u8,
) -> ProcessedSpliceRec<'a> {
    let distance = config.distance;
    match stage {
        // Stage 1: Process D1.
        1 => {
            if let Some(new_seq) =
                pattern_search_trim_seq(splice_rec.seq, config.d1.as_bytes(), distance)
            {
                chain.add_splice_event("D1".to_string());
                process_splice_rec(ProcessedSpliceRec::init(new_seq), chain, config, 2)
            } else {
                chain.add_splice_event("noD1".to_string());

                //TODO: We need to continue searching for acceptors even if D1 is not found.
                //If we find an acceptor, we should output the sequence before the acceptor for futher inspection.
                //And we need to map these D1 alternative donors on the sequence map.
                process_splice_rec(splice_rec, chain, config, 8) // Jump to stage 8 to search for acceptors without D1. And find out what is the alternative donor.
            }
        }
        // Stage 2: Process the acceptor immediately after D1.
        2 => {
            if let Some((acc, new_seq)) =
                pattern_search_trim_seq_batch(splice_rec.seq, &config.d1_to_all, distance)
            {
                chain.add_splice_event(acc.clone());

                let new_splice_rec = ProcessedSpliceRec::init_with_option_fields(
                    new_seq,
                    splice_rec.alternative_d1,
                    splice_rec.unknown_sequence_after_d1,
                );
                match acc.as_str() {
                    "A1" => process_splice_rec(new_splice_rec, chain, config, 3),
                    "A2" => process_splice_rec(new_splice_rec, chain, config, 5),
                    _ => new_splice_rec,
                }
            } else {
                chain.add_splice_event("unknown".to_string());

                //Need to know what sequence is left here for futher inspection.
                //Export the unknown sequence to the field of unknown_sequence_after_d1 in SpliceEvents struct.
                //may only need to keep up to 15 bases

                let query_length = splice_rec.seq.len().min(15); // Limit to the first 15 bases, if the query sequence is shorter, take all of it.
                if query_length > 0 {
                    let unknown_sequence_after_d1 = Some(&splice_rec.seq[..query_length]);

                    ProcessedSpliceRec::init_with_option_fields(
                        splice_rec.seq,
                        splice_rec.alternative_d1,
                        unknown_sequence_after_d1,
                    )
                } else {
                    splice_rec
                }
            }
        }
        // Stage 3: Process D2 (for the A1 branch).
        3 => {
            if let Some(new_seq) =
                pattern_search_trim_seq(splice_rec.seq, config.d2.as_bytes(), distance)
            {
                chain.add_splice_event("D2".to_string());
                let new_splice_rec = ProcessedSpliceRec::init_with_option_fields(
                    new_seq,
                    splice_rec.alternative_d1,
                    splice_rec.unknown_sequence_after_d1,
                );
                process_splice_rec(new_splice_rec, chain, config, 4)
            } else {
                chain.add_splice_event("noD2".to_string());
                splice_rec
            }
        }
        // Stage 4: Process acceptor after D2.
        4 => {
            if let Some((acc, new_seq)) =
                pattern_search_trim_seq_batch(splice_rec.seq, &config.d2_to_all, distance)
            {
                chain.add_splice_event(acc.clone());
                let new_splice_rec = ProcessedSpliceRec::init_with_option_fields(
                    new_seq,
                    splice_rec.alternative_d1,
                    splice_rec.unknown_sequence_after_d1,
                );
                match acc.as_str() {
                    "D2-unspliced" => process_splice_rec(new_splice_rec, chain, config, 6),
                    "A2" => process_splice_rec(new_splice_rec, chain, config, 5),
                    _ => new_splice_rec,
                }
            } else {
                chain.add_splice_event("unknown".to_string());
                splice_rec
            }
        }
        // Stage 5: Process D3 branch (common for any A2 outcome).
        5 => {
            if let Some(new_seq) =
                pattern_search_trim_seq(splice_rec.seq, config.d3.as_bytes(), distance)
            {
                chain.add_splice_event("D3".to_string());
                let new_splice_rec = ProcessedSpliceRec::init_with_option_fields(
                    new_seq,
                    splice_rec.alternative_d1,
                    splice_rec.unknown_sequence_after_d1,
                );
                process_splice_rec(new_splice_rec, chain, config, 7)
            } else {
                chain.add_splice_event("noD3".to_string());
                splice_rec
            }
        }
        // Stage 6: Process D2b for the "D2-unspliced" branch.
        6 => {
            if let Some(new_seq) =
                pattern_search_trim_seq(splice_rec.seq, config.d2b.as_bytes(), distance)
            {
                chain.add_splice_event("D2b".to_string());
                if let Some((acc, new_seq2)) =
                    pattern_search_trim_seq_batch(new_seq, &config.d2b_to_all, distance)
                {
                    chain.add_splice_event(acc.clone());
                    let new_splice_rec = ProcessedSpliceRec::init_with_option_fields(
                        new_seq2,
                        splice_rec.alternative_d1,
                        splice_rec.unknown_sequence_after_d1,
                    );

                    if acc == "A2" {
                        process_splice_rec(new_splice_rec, chain, config, 5)
                    } else {
                        new_splice_rec
                    }
                } else {
                    chain.add_splice_event("unknown".to_string());
                    ProcessedSpliceRec::init_with_option_fields(
                        new_seq,
                        splice_rec.alternative_d1,
                        splice_rec.unknown_sequence_after_d1,
                    )
                }
            } else {
                chain.add_splice_event("noD2b".to_string());
                splice_rec
            }
        }
        // Stage 7: Process the acceptor downstream of D3.
        7 => {
            if let Some((acc, new_seq)) =
                pattern_search_trim_seq_batch(splice_rec.seq, &config.d3_to_all, distance)
            {
                chain.add_splice_event(acc.clone());
                let new_splice_rec = ProcessedSpliceRec::init_with_option_fields(
                    new_seq,
                    splice_rec.alternative_d1,
                    splice_rec.unknown_sequence_after_d1,
                );
                new_splice_rec
            } else {
                chain.add_splice_event("unknown".to_string());
                splice_rec
            }
        }
        // Stage 8: Search for acceptors without D1 (for sequences missing D1).
        8 => {
            // "Need to know what sequence is left here for further inspection.
            // First thing would be redo the search function. The current search function will return trimmed the sequence.
            // But we want to export the leading sequence before the acceptor for further inspection.
            // It will export to the field of pre_splice_sequence in SpliceEvents struct."
            // in this step we we force the mistatch distance to be 1, because we want to find exact matches only.

            if let Some(matched_pattern) = pattern_search_trim_seq_with_donor_batch(
                splice_rec.seq,
                &config.alternative_d1_to_all,
                1,
            ) {
                // If we find an acceptor, we should output the sequence before the acceptor for further inspection.
                let acc = matched_pattern.acceptor_name.unwrap();
                let new_seq = matched_pattern.trimmed_sequence;
                let alternative_d1 = matched_pattern.donor_sequence;
                chain.add_splice_event(acc.clone());
                let new_splice_rec = ProcessedSpliceRec::init_with_option_fields(
                    new_seq,
                    Some(alternative_d1),
                    splice_rec.unknown_sequence_after_d1,
                );
                match acc.as_str() {
                    "A1" => process_splice_rec(new_splice_rec, chain, config, 3),
                    "A2" => process_splice_rec(new_splice_rec, chain, config, 5),
                    "gag-AUG" => {
                        // If we find intact gag-AUG (unspliced), we consider it as unspliced and stop further searching.
                        chain.add_splice_event("unspliced".to_string());
                        new_splice_rec
                    }
                    _ => new_splice_rec,
                }
            } else {
                chain.add_splice_event("no-acceptor".to_string());
                splice_rec
            }
        }
        _ => splice_rec,
    }
}

/// Computes the reverse complement of a DNA sequence.
///
/// # Arguments
/// * `seq` - The DNA sequence to reverse complement.
///
/// # Returns
/// A `String` containing the reverse complement of the input sequence.
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| c as u8)
        .map(|c| dna::complement(c))
        .collect::<Vec<_>>()
        .pipe(String::from_utf8)
        .unwrap()
}

/// Checks if the overlap between two sequences is acceptable given the error rate.
///
/// # Arguments
/// * `s1` - The first sequence.
/// * `s2` - The second sequence.
/// * `error_rate` - The maximum allowed error rate.
///
/// # Returns
/// `true` if the overlap is acceptable, `false` otherwise.
fn is_overlap_acceptable(s1: &str, s2: &str, error_rate: f64) -> bool {
    let overlap_length = s1.len();
    let mismatches = s1.chars().zip(s2.chars()).filter(|(a, b)| a != b).count();
    // The overlap is acceptable if the proportion of mismatches does not exceed error_rate.
    (mismatches as f64) <= (overlap_length as f64 * error_rate)
}

/// Finds the maximum overlap between two sequences that meets the error rate threshold.
///
/// # Arguments
/// * `r1` - The first sequence.
/// * `r2` - The second sequence.
/// * `min_overlap` - The minimum required overlap length.
/// * `error_rate` - The maximum allowed error rate.
///
/// # Returns
/// An `Option<usize>` containing the length of the acceptable overlap, or `None` if no overlap is found.
fn find_overlap(r1: &str, r2: &str, min_overlap: usize, error_rate: f64) -> Option<usize> {
    let max_overlap = r1.len().min(r2.len());
    // Check overlaps from largest possible to min_overlap
    for overlap in (min_overlap..=max_overlap).rev() {
        let r1_overlap = &r1[r1.len() - overlap..];
        let r2_overlap = &r2[..overlap];
        if is_overlap_acceptable(r1_overlap, r2_overlap, error_rate) {
            return Some(overlap);
        }
    }
    None
}

/// Joins two paired-end reads by reverse-complementing the second read and merging them.
///
/// # Arguments
/// * `r1` - The forward read.
/// * `r2` - The reverse read.
/// * `min_overlap` - The minimum required overlap length.
/// * `error_rate` - The maximum allowed error rate.
///
/// # Returns
/// An `Option<String>` containing the joined sequence, or `None` if the reads cannot be merged.
fn join_reads(r1: &str, r2: &str, min_overlap: usize, error_rate: f64) -> Option<String> {
    let r2_rc = reverse_complement(r2);
    if let Some(overlap) = find_overlap(r1, &r2_rc, min_overlap, error_rate) {
        // Merge by taking full R1 and appending the non-overlapping part of the reverse-complemented R2.
        let joined = format!("{}{}", r1, &r2_rc[overlap..]);
        Some(joined)
    } else {
        None
    }
}

/// Searches for a pattern in a sequence using the Myers algorithm.
///
/// # Arguments
/// * `sequence` - The sequence to search.
/// * `pattern` - The pattern to search for.
/// * `distance` - The maximum allowed edit distance.
///
/// # Returns
/// An `Option<Alignment>` containing the best alignment, or `None` if no match is found.
pub fn pattern_search(sequence: &[u8], pattern: &[u8], distance: u8) -> Option<Alignment> {
    let mut aln = Alignment::default();
    let mut myers = myers::Myers::<u64>::new(pattern);
    let mut lazy_matches = myers.find_all_lazy(sequence, distance);

    match lazy_matches.by_ref().min_by_key(|&(_, dist)| dist) {
        Some((best_end, _)) => {
            lazy_matches.alignment_at(best_end, &mut aln);

            Some(aln)
        }
        None => None,
    }
}

/// Searches for a pattern in a sequence and trims the sequence up to the match, aka, removing everything before the pattern, but retaining the pattern and everything after it.
///
/// # Arguments
/// * `sequence` - The sequence to search.
/// * `pattern` - The pattern to search for.
/// * `distance` - The maximum allowed edit distance.
///
/// # Returns
/// An `Option<&[u8]>` containing the trimmed sequence, or `None` if no match is found.
pub fn pattern_search_trim_seq<'a>(
    sequence: &'a [u8],
    pattern: &[u8],
    distance: u8,
) -> Option<&'a [u8]> {
    let mut aln = Alignment::default();
    let mut myers = myers::Myers::<u64>::new(pattern);
    let mut lazy_matches = myers.find_all_lazy(sequence, distance);

    match lazy_matches.by_ref().min_by_key(|&(_, dist)| dist) {
        Some((best_end, _)) => {
            lazy_matches.alignment_at(best_end, &mut aln);

            Some(&sequence[aln.ystart..])
        }
        None => None,
    }
}

#[derive(Debug, Clone)]
pub struct MatchedPattern<'a> {
    pub acceptor_name: Option<String>, // Name of the matched acceptor, if applicable
    pub donor_sequence: &'a [u8],
    pub trimmed_sequence: &'a [u8],
}

pub fn pattern_search_trim_with_donor<'a>(
    sequence: &'a [u8],
    pattern: &[u8],
    distance: u8,
) -> Option<MatchedPattern<'a>> {
    let mut aln = Alignment::default();
    let mut myers = myers::Myers::<u64>::new(pattern);
    let mut lazy_matches = myers.find_all_lazy(sequence, distance);

    match lazy_matches.by_ref().min_by_key(|&(_, dist)| dist) {
        Some((best_end, _)) => {
            lazy_matches.alignment_at(best_end, &mut aln);

            let mut donor_sequence = &sequence[..aln.ystart];
            let trimmed_sequence = &sequence[aln.ystart..];

            // only keep up to 15 bases of donor sequence
            if donor_sequence.len() > 15 {
                donor_sequence = &donor_sequence[donor_sequence.len() - 15..];
            }

            Some(MatchedPattern {
                acceptor_name: None, // Placeholder, can be set if needed
                donor_sequence,
                trimmed_sequence,
            })
        }
        None => None,
    }
}

/// Searches for multiple patterns in a sequence and trims the sequence up to the first match.
///
/// # Arguments
/// * `sequence` - The sequence to search.
/// * `list` - A list of splice steps containing patterns and acceptors.
/// * `distance` - The maximum allowed edit distance.
///
/// # Returns
/// An `Option<(String, &[u8])>` containing the matched acceptor and the trimmed sequence, or `None` if no match is found.
pub fn pattern_search_trim_seq_batch<'a>(
    sequence: &'a [u8],
    list: &Vec<SpliceStep>,
    distance: u8,
) -> Option<(String, &'a [u8])> {
    for step in list {
        let pattern = step.pattern.as_bytes();

        if let Some(trimmed_sequence) = pattern_search_trim_seq(sequence, pattern, distance) {
            // println!("matched acceptor: {:?}", step.acceptor);
            let matched_acceptor = step.acceptor.clone();
            return Some((matched_acceptor, trimmed_sequence));
        }
    }
    None
}

pub fn pattern_search_trim_seq_with_donor_batch<'a>(
    sequence: &'a [u8],
    list: &Vec<SpliceStep>,
    distance: u8,
) -> Option<MatchedPattern<'a>> {
    for step in list {
        // be careful here, the order of the list matters.
        let pattern = step.pattern.as_bytes();

        if let Some(matched_pattern) = pattern_search_trim_with_donor(sequence, pattern, distance) {
            // println!("matched acceptor: {:?}", step.acceptor);
            let mut matched_pattern = matched_pattern;
            matched_pattern.acceptor_name = Some(step.acceptor.clone());
            return Some(matched_pattern);
        }
    }
    None
}

#[cfg(test)]

mod tests {

    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("aaccggtt"), "aaccggtt");
        assert_eq!(reverse_complement("AGGTACCA"), "TGGTACCT");
        assert_eq!(reverse_complement("ARWN"), "NWYT");
    }

    #[test]
    fn test_find_overlap() {
        let seq1 = &format!(
            "{}{}{}{}",
            "A".repeat(6),
            "C".repeat(6),
            "T".repeat(6),
            "G".repeat(6)
        );
        let seq2 = "GGGGGGTTTTTCCCCCCTTTTT";
        let seq3 = "CCCCCCTTTTTTGGGGGGCCCAAAGG";
        let seq4 = "GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let seq5 = "AAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT";
        let min_overlap = 10;
        let error_rate = 0.01;
        assert_eq!(find_overlap(seq1, seq2, min_overlap, error_rate), None);
        assert_eq!(find_overlap(seq1, seq3, min_overlap, error_rate), Some(18));
        assert_eq!(find_overlap(seq4, seq5, min_overlap, error_rate), Some(100));
        assert_eq!(find_overlap(seq4, seq5, min_overlap, 0.0), None);
    }

    #[test]

    fn test_joined_umi_sequence() {
        let r1 = fasta::Record::with_attrs("seq1_r1", None, b"NNNNTTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGA");
        let r2 = fasta::Record::with_attrs("seq1_r2", None, b"NNNNNNNNNNNNNNAAATCTACTAATTTTCTCCATTTAGTACTGTCTTTTTTCTTTATGGCAAATACTGGAGTATTGTATGGATTTTCAGGCCCAATTTTTGAAATTTTCCCTTC");
        let mut joined_umi_sequence = JoinedUmiSequnce::from_fasta_record(&r1, &r2, 4, 14);
        joined_umi_sequence.join();
        assert_eq!(joined_umi_sequence.joined_sequence, Some("TTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT".to_string()));
        assert_eq!(joined_umi_sequence.joined_seq_to_fasta_record(), Some(fasta::Record::with_attrs("seq1_r1-joined", Some("ForwardNs: NNNN, UMI: NNNNNNNNNNNNNN"), b"TTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT")));
        assert_eq!(joined_umi_sequence.find_sequence_for_search(), "TTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT".to_string());
    }

    #[test]
    fn test_pattern_search() {
        let sequence = b"TTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT";
        let pattern = b"AGAGATGGAAAAGGAAGGGAAA";
        let distance = 2;
        let aln = pattern_search(sequence, pattern, distance);
        dbg!(&aln);
        assert!(aln.is_some());
    }

    #[test]
    fn test_pattern_search_trim_seq() {
        let sequence = b"TTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT";
        let pattern = b"AGAGCTGGAAAGGAAGGGAAA";
        let distance = 2;
        let trimmed_seq = pattern_search_trim_seq(sequence, pattern, distance);
        dbg!(&trimmed_seq);
        assert_eq!(trimmed_seq, Some(b"AGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT".as_ref()));

        let pattern2 = b"AGAGCTGGAAACGAAGGGAAA"; // 3 mismatches, should return None
        let trimmed_seq2 = pattern_search_trim_seq(sequence, pattern2, distance);
        dbg!(&trimmed_seq2);
        assert_eq!(trimmed_seq2, None);
    }

    #[test]
    fn test_pattern_search_trim_with_donor() {
        let sequence = b"TTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT";
        let pattern = b"AGAGATGGAAAAGGAAGGGAAA";
        let distance = 2;
        let matched_pattern = pattern_search_trim_with_donor(sequence, pattern, distance);
        dbg!(&matched_pattern);
        assert!(matched_pattern.is_some());
        let matched_pattern = matched_pattern.unwrap();
        assert_eq!(matched_pattern.donor_sequence, b"AGTAGAAATTTGTAC");
        assert_eq!(matched_pattern.trimmed_sequence, b"AGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT".as_ref());
    }

    #[test]
    fn test_pattern_search_trim_with_donor2() {
        let seq = b"TTATGAACAGGGACTTGAAAGCGAAAGTAAAGCCAGAGGAGATCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGAATTGGGTGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGGTCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTCCATTTCCAAGTTTGTTTCA";
        let pattern = b"AATCTGCTAT";
        let distance = 0;
        let matched_pattern = pattern_search_trim_with_donor(seq, pattern, distance);
        dbg!(&matched_pattern);
        assert!(matched_pattern.is_none());
    }
}
