use bio::io::fasta;
use tap::Pipe;
use bio::alphabets::dna;
use serde::{Deserialize, Serialize};
use bio::pattern_matching::myers;
use bio::alignment::Alignment;
use crate::config::{SpliceConfig, SpliceStep};
use crate::splice_events::{SpliceEvents, SpliceChain};

#[derive(Debug, Deserialize, Serialize)]
pub struct JoinedUmiSequnce {
    pub sequence_id: String,
    pub forward_ns: String,
    pub umi: String,
    pub forward_sequence: String,
    pub reverse_sequence: String,
    pub joined_sequence: Option<String>,
}

impl JoinedUmiSequnce {
    pub fn from_fasta_record (record_r1: &fasta::Record, record_r2: &fasta::Record, forward_n_size: usize, umi_size: usize) -> JoinedUmiSequnce {
        let sequence_id = record_r1.id()[..(record_r1.id().len() - 3)].to_string();

        let forward_seq = record_r1.seq().to_vec().pipe(String::from_utf8).unwrap().to_uppercase();
        let reverse_seq = record_r2.seq().to_vec().pipe(String::from_utf8).unwrap().to_uppercase();

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

    pub fn join(&mut self) {
        let min_overlap = 10; // Minimum overlap length required for merging. Consider moving to master config.
        let error_rate = 0.02; // Maximum error rate allowed in the overlap region. Consider moving to master config.
        if let Some(joined) = join_reads(&self.forward_sequence, &self.reverse_sequence, min_overlap, error_rate) {
            self.joined_sequence = Some(joined);
        }
    }

    pub fn joined_seq_to_fasta_record(&self) -> Option<fasta::Record> {
        let id = format!("{}-joined", self.sequence_id);
        let desc = Some(format!("ForwardNs: {}, UMI: {}", self.forward_ns, self.umi));
        if let None = self.joined_sequence {
            return None;
        } else {
            let seq = self.joined_sequence.as_ref().unwrap();
            Some(fasta::Record::with_attrs(&id, desc.as_deref(), seq.as_bytes()))
        }
    }

    pub fn find_sequence_for_search(&self) -> String {
        match &self.joined_sequence {
            Some(joined_sequence) => joined_sequence.to_string(),
            None => self.forward_sequence.to_string() + &self.reverse_sequence,
        }
    }

    pub fn check_splice_event(&self, splice_config: &SpliceConfig) -> SpliceEvents {
        let seq = self.find_sequence_for_search();
        let mut sequence_to_search = seq.as_bytes();
        let distance = splice_config.distance;
        let mut splice_chain = SpliceChain::new();

        let d1 = splice_config.d1.clone();
        let d2 = splice_config.d2.clone();
        let d2b = splice_config.d2b.clone();
        let d3 = splice_config.d3.clone();

        let d1_to_all = &splice_config.d1_to_all;
        let d2_to_all = &splice_config.d2_to_all;
        let d2b_to_all = &splice_config.d2b_to_all;
        let d3_to_all = &splice_config.d3_to_all;


        // step 1, check D1

        match pattern_search_trim_seq(sequence_to_search, &d1.as_bytes(), distance) {
            Some(trimmed) => {
                splice_chain.add_splice_event("D1".to_string());
                sequence_to_search = trimmed;
                
                // step2, check after D1

                match pattern_search_trim_seq_batch(sequence_to_search, d1_to_all, distance) {
                    Some((matched_acceptor, trimmed)) => {
                        splice_chain.add_splice_event(matched_acceptor.clone());
                        sequence_to_search = trimmed;

                        match matched_acceptor.as_str() {
                            "A1" => 
                            {
                                // step 3, check if D2 is present
                                
                                match pattern_search_trim_seq(sequence_to_search, &d2.as_bytes(), distance) {
                                    Some(trimmed) =>
                                    {
                                        splice_chain.add_splice_event("D2".to_string());
                                        sequence_to_search = trimmed;

                                        // step 4, when D2 is present check after D2

                                        match pattern_search_trim_seq_batch(sequence_to_search, d2_to_all, distance) {

                                            Some((matched_acceptor, trimmed)) => {
                                                splice_chain.add_splice_event(matched_acceptor.clone());
                                                sequence_to_search = trimmed;

                                                match matched_acceptor.as_str() {
                                                    "D2-unspliced" => 
                                                    {
                                                        // step 5, when D2 is not spliced, check if D2b is present
                                                        
                                                        match pattern_search_trim_seq(sequence_to_search, &d2b.as_bytes(), distance) {
                                                            Some(trimmed) => {
                                                                splice_chain.add_splice_event("D2b".to_string());
                                                                sequence_to_search = trimmed;

                                                                // step 6, check downstream of D2b

                                                                match pattern_search_trim_seq_batch(sequence_to_search, d2b_to_all, distance) {
                                                                    Some((matched_acceptor, trimmed)) => {
                                                                        splice_chain.add_splice_event(matched_acceptor.clone());
                                                                        sequence_to_search = trimmed;

                                                                        match matched_acceptor.as_str() {
                                                                            "A2" => 
                                                                            {
                                                                                // step 7, when A2 is present, check D3
                                                                                match pattern_search_trim_seq(sequence_to_search, &d3.as_bytes(), distance) {
                                                                                    Some(trimmed) => {
                                                                                        splice_chain.add_splice_event("D3".to_string());
                                                                                        sequence_to_search = trimmed;
                                                                                        
                                                                                        // step 8 check downstream of D3, end of all regardless of the match
                                                                            
                                                                                        match pattern_search_trim_seq_batch(sequence_to_search, d3_to_all, distance) {
                                                                                            Some((matched_acceptor, trimmed)) => {
                                                                                                splice_chain.add_splice_event(matched_acceptor.clone());
                                                                                                sequence_to_search = trimmed;
                                                                                            },
                                                                                            None => {
                                                                                                // end of the splicing chain at step 8
                                                                                                splice_chain.add_splice_event("unknown".to_string());
                                                                                            }
                                                                                        }                       
                                                                                    },
                                                                                    None =>
                                                                                    {
                                                                                        //D3 absent, end of the splicing chain at step 7
                                                                                        splice_chain.add_splice_event("noD3".to_string());
                                                                        
                                                                                    }
                                                                                }
                                                                            },
                                                                            // other downstream acceptors (A3-A7), end of the splicing chain at step 7
                                                                            _ => {}
                                                                        }

                                                                    },
                                                                    None => {
                                                                        // end of the splicing chain at step 6
                                                                        splice_chain.add_splice_event("unknown".to_string());
                                                                    }
                                                                }

                                                            },
                                                            None => {
                                                                // end of the splicing chain at step 5 when D2b is absent
                                                                splice_chain.add_splice_event("noD2b".to_string());
                                                            }
                                                        }

                                                    },
                                                   "A2" =>
                                                    {
                                                        // step 5, when D2 is spliced, check D3
                                                        match pattern_search_trim_seq(sequence_to_search, &d3.as_bytes(), distance) {
                                                            Some(trimmed) => {
                                                                splice_chain.add_splice_event("D3".to_string());
                                                                sequence_to_search = trimmed;
                                                                
                                                                // step 6 check downstream of D3, end of all regardless of the match
                        
                                                                match pattern_search_trim_seq_batch(sequence_to_search, d3_to_all, distance) {
                                                                    Some((matched_acceptor, trimmed)) => {
                                                                        splice_chain.add_splice_event(matched_acceptor.clone());
                                                                        sequence_to_search = trimmed;
                                                                    },
                                                                    None => {
                                                                        // end of the splicing chain at step 4
                                                                        splice_chain.add_splice_event("unknown".to_string());
                                                                    }
                                                                }                       
                                                            },
                                                            None =>
                                                            {
                                                                //D3 absent, end of the splicing chain at step 6
                                                                splice_chain.add_splice_event("noD3".to_string());
                        
                                                            }
                                                        }

                                                    },
                                                  // other downstream acceptors (A3-A7), end of the splicing chain at step 4
                                                    _ => {

                                                    }
                                                };
                                            },
                                            None => {
                                                // end of the splicing chain at step 4
                                                splice_chain.add_splice_event("unknown".to_string());
                                            }
                                        };

                                    },
                                    None => 
                                    {
                                        // D2 absent, end of the splicing chain at step 3
                                        splice_chain.add_splice_event("noD2".to_string());
                                        
                                    }
                                }
                            },
                            "A2" => 
                            {
                                // step 3, parrell, check D3
                                match pattern_search_trim_seq(sequence_to_search, &d3.as_bytes(), distance) {
                                    Some(trimmed) => {
                                        splice_chain.add_splice_event("D3".to_string());
                                        sequence_to_search = trimmed;
                                        
                                        // step 4 check downstream of D3, end of all regardless of the match

                                        match pattern_search_trim_seq_batch(sequence_to_search, d3_to_all, distance) {
                                            Some((matched_acceptor, trimmed)) => {
                                                splice_chain.add_splice_event(matched_acceptor.clone());
                                                sequence_to_search = trimmed;
                                            },
                                            None => {
                                                // end of the splicing chain at step 4
                                                splice_chain.add_splice_event("unknown".to_string());
                                            }
                                        }
                                    },
                                    None =>
                                    {
                                        //D3 absent, end of the splicing chain at step 3
                                        splice_chain.add_splice_event("noD3".to_string());

                                    }
                                }
                            },
                            // matching with D1-unspliced, A3-A7, end of splicing chain at step 2
                            _ => {}
                        }                        
                    },
                    None => {

                        // unknown acceptor after D1, end of splicing chain at step 2
                        splice_chain.add_splice_event("unknown".to_string());
                    }
                }
            },
            None => {

                // no D1, end of splicing chain at step 1
                splice_chain.add_splice_event("noD1".to_string());
            }
        }

        let mut event = SpliceEvents::from_joined_umi_with_event(&self, splice_chain);
        event.add_post_splice_sequence(String::from_utf8(sequence_to_search.to_vec()).unwrap());

        event
    }

}

/// Compute the reverse complement of a DNA sequence.
/// Use the `dna` alphabet from the `bio` crate for maximum performance.
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| c as u8)
        .map(|c| dna::complement(c))
        .collect::<Vec<_>>()
        .pipe(String::from_utf8)
        .unwrap()
}

/// Check if the overlap between two sequence slices is acceptable given the error rate.
/// The overlap is acceptable if the fraction of mismatched bases is less than or equal to `error_rate`.
fn is_overlap_acceptable(s1: &str, s2: &str, error_rate: f64) -> bool {
    let overlap_length = s1.len();
    let mismatches = s1
        .chars()
        .zip(s2.chars())
        .filter(|(a, b)| a != b)
        .count();
    // The overlap is acceptable if the proportion of mismatches does not exceed error_rate.
    (mismatches as f64) <= (overlap_length as f64 * error_rate)
}

/// Find the maximum overlap between the suffix of `r1` and the prefix of `r2` that meets the error rate threshold.
/// Returns the length of the acceptable overlap if at least `min_overlap` bases match (within error rate), or None otherwise.
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

/// Join two paired-end reads by reverse-complementing R2, finding an overlapping region
/// (allowing for a sequencing error rate), and merging them into one contiguous sequence.
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

pub fn pattern_search(sequence: &[u8], pattern: &[u8], distance: u8) -> Option<Alignment> {

    let mut aln = Alignment::default();
    let mut myers = myers::Myers::<u64>::new(pattern);
    let mut lazy_matches = myers.find_all_lazy(sequence, distance);

    match lazy_matches.by_ref().min_by_key(|&(_, dist)| dist) {
        Some((best_end, _)) => {
            lazy_matches.alignment_at(best_end, &mut aln);

            Some(aln)
        }
        None => {
            None
        }
    }
}

pub fn pattern_search_trim_seq<'a>(sequence: &'a [u8], pattern: &[u8], distance: u8) -> Option<&'a [u8]> {

    let mut aln = Alignment::default();
    let mut myers = myers::Myers::<u64>::new(pattern);
    let mut lazy_matches = myers.find_all_lazy(sequence, distance);

    match lazy_matches.by_ref().min_by_key(|&(_, dist)| dist) {
        Some((best_end, _)) => {
            lazy_matches.alignment_at(best_end, &mut aln);

            Some(&sequence[aln.ystart..])
        }
        None => {
            None
        }
    }
}

pub fn pattern_search_trim_seq_batch<'a>(sequence: &'a [u8], list: &Vec<SpliceStep>, distance: u8) -> Option<(String, &'a [u8])> {

    for step in list {
        let pattern = step.pattern.as_bytes();

        if let Some(trimmed_sequence) = pattern_search_trim_seq(sequence, pattern, distance) {
            let matched_acceptor = step.acceptor.clone();
            Some((matched_acceptor, trimmed_sequence));
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
    fn test_find_overlap(){
        let seq1 = &format!("{}{}{}{}", "A".repeat(6), "C".repeat(6), "T".repeat(6), "G".repeat(6));
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
        println!("{:#?}", joined_umi_sequence);
        assert_eq!(joined_umi_sequence.joined_sequence, Some("TTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT".to_string()));
        assert_eq!(joined_umi_sequence.joined_seq_to_fasta_record(), Some(fasta::Record::with_attrs("seq1-joined", Some("ForwardNs: NNNN, UMI: NNNNNNNNNNNNNN"), b"TTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT")));
        assert_eq!(joined_umi_sequence.find_sequence_for_search(), "TTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTT".to_string());
    }





}






