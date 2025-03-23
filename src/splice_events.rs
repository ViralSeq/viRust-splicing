use serde::{Deserialize, Serialize};
use crate::config::SpliceConfig;
use std::error::Error;
use crate::joined_umi_sequence::JoinedUmiSequnce;
use crate::joined_umi_sequence::pattern_search;
use crate::config::SpliceAssayType;

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct SpliceEvents {
    pub sequence_id: String,
    pub search_sequence: String,
    pub post_splice_sequence: Option<String>, // sequence post last splicing site that can be identified. used to check size class.
    pub umi: String,
    pub umi_family: Option<String>,
    pub splice_category: SpliceChain,
    pub size_class: Option<SizeClass>,

}

#[derive(Debug, Deserialize, Serialize, PartialEq, Clone)]
pub struct SpliceChain {
    pub splice_event: Vec<String>,
}

#[derive(Debug, Deserialize, Serialize, PartialEq, Clone)]
pub enum SizeClass {
    OnePointEightKb,
    FourKb,
    BothClass,
    Unknown,
}

impl SpliceEvents {
    pub fn from_joined_umi_with_event(joined_umi: &JoinedUmiSequnce, splice_category: SpliceChain) -> SpliceEvents {
        let sequence_id = joined_umi.sequence_id.clone();
        let umi = joined_umi.umi.clone();
        let umi_family = None;
        let search_sequence = joined_umi.find_sequence_for_search().clone();
        let post_splice_sequence = None;
        let size_class = None;

        SpliceEvents {
            sequence_id,
            search_sequence,
            post_splice_sequence,
            umi,
            umi_family,
            splice_category,
            size_class,
        }
    }

    pub fn add_post_splice_sequence(&mut self, sequence: String) {
        self.post_splice_sequence = Some(sequence);
    }

    pub fn find_size_class(&mut self, config: &SpliceConfig) -> Result<(), Box<dyn Error>> {
        for event in self.splice_category.splice_event.iter() {
            if event == "A7" { // if we find an A7 event, we know it is 1.8kb, nef transcript
                self.size_class = Some(SizeClass::OnePointEightKb);
                return Ok(());
            } else if event == "unknown" || event == "noD1" || event == "D1-unspliced" {
                self.size_class = Some(SizeClass::Unknown); 
                return Ok(());
            }
        };

        let assay_type = &config.splice_assay_type;
        let post_splice_sequence = match &self.post_splice_sequence {
            Some(seq) => seq,
            None => {
                self.size_class = Some(SizeClass::Unknown);
                return Ok(());
            },
        };

        match assay_type {
            SpliceAssayType::SizeSpecific | SpliceAssayType::KMer => {
                if post_splice_sequence.len() < 30 {
                    self.size_class = Some(SizeClass::Unknown);
                    return Ok(());
                }

                let search_pattern = &post_splice_sequence[post_splice_sequence.len() - 30..post_splice_sequence.len() - 15];
                if let Some(aln) = pattern_search(config.full_length_sequence.as_bytes(), search_pattern.as_bytes(), config.distance) {

                    self.size_class = match assay_type {
                        SpliceAssayType::SizeSpecific => {
                            if aln.ystart < config.d4_breakpoint_position {
                                Some(SizeClass::FourKb)
                            } else if aln.ystart > config.a7_breakpoint_position {
                                Some(SizeClass::OnePointEightKb)
                            } else {
                                Some(SizeClass::Unknown)
                            }
                        }
                        SpliceAssayType::KMer => {
                            if aln.ystart > config.a7_breakpoint_position {
                                Some(SizeClass::OnePointEightKb)
                            } else {
                                Some(SizeClass::FourKb)
                            }
                        }
                        _ => unreachable!(),
                    };
                } else {
                    self.size_class = Some(SizeClass::Unknown);
                }
            }
            SpliceAssayType::RandomReverse => {
                let slice_of_patterns = slice_from_end(post_splice_sequence, 15, 15, 10);
                if slice_of_patterns.is_empty() {
                    self.size_class = Some(SizeClass::Unknown);
                    return Ok(());
                }

                let mut found_1_8k = false;
                let mut found_both = false;

                for pattern in slice_of_patterns.iter() {
                    if let Some(aln) = pattern_search(config.full_length_sequence.as_bytes(), pattern.as_bytes(), config.distance) {
                        if aln.ystart < config.a7_breakpoint_position && aln.ystart > config.d4_breakpoint_position {
                            self.size_class = Some(SizeClass::FourKb); // if we find a 4kb pattern (between D4 and A7), we can stop
                            return Ok(());
                        } else if aln.ystart > config.a7_breakpoint_position {
                            found_1_8k = true;
                        } else if aln.ystart < config.d4_breakpoint_position {
                            found_both = true;
                        }
                    }
                }

                if found_1_8k  {
                    self.size_class = Some(SizeClass::OnePointEightKb); // if we find a 1.8kb pattern (after A7), we can stop, it doesn't matter if we found a both pattern (before D4)
                } else if found_both {
                    self.size_class = Some(SizeClass::BothClass); // if the only thing we can find is the the pattern before D4, it can be both size classes. 
                } else {
                    self.size_class = Some(SizeClass::Unknown); // if no alignment is found, we don't know the size class
                }
            
            }
        }

        Ok(())
    }


}

impl SpliceChain {
    pub fn new() -> SpliceChain {
        SpliceChain {
            splice_event: Vec::new(),
        }
    }
    pub fn add_splice_event(&mut self, event: String) {
        self.splice_event.push(event);
    }
}


// TODO: Work on this tomorrow. 
pub fn find_umi_family_from_events(mut events: Vec<SpliceEvents>) -> Vec<SpliceEvents> {
    let mut outcome_events: Vec<SpliceEvents> = Vec::new();

    outcome_events
}

fn slice_from_end(input: &str, chunk_size: u8, offset: u8, min_length: u8) -> Vec<String> {

    let mut chunks = Vec::new();

    let mut end = input.len() - offset as usize;
    while end > (min_length as usize - 1) {
        let start = if end >= chunk_size as usize { end - chunk_size as usize } else { 0 };
        chunks.push(input[start..end].to_string());
        end = start;
    }
    chunks
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::open_fasta_file;
    use tap::Pipe;

    #[test]
    fn test_slice_from_end() {
        let input = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        let chunk_size = 10;
        let offset = 15;
        let expected = vec![
            "BCDEFGHIJK".to_string(),
            // "A".to_string(),
        ];
        let result = slice_from_end(input, chunk_size, offset, 5);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_size_class() {

        let nl43_file = "data/nl43.fasta";
        let fasta_reader = open_fasta_file(nl43_file).unwrap();
        let record = fasta_reader.records().next().unwrap().unwrap();
        let nl43_ref_seq = record.seq().to_vec();
        let search_seq_before_d4 = nl43_ref_seq[5100..5160].to_vec().pipe(String::from_utf8).unwrap();
        let search_seq_after_a7 = nl43_ref_seq[8000..8060].to_vec().pipe(String::from_utf8).unwrap();
        let search_seq_between_d4_a7 = nl43_ref_seq[7000..7060].to_vec().pipe(String::from_utf8).unwrap();
        let search_seq_cross_a7 = nl43_ref_seq[7580..7700].to_vec().pipe(String::from_utf8).unwrap();

        let mut splice_event1 = SpliceEvents {
            sequence_id: "test".to_string(),
            search_sequence: "TTTTCCTAGGATATGGCTCCATAACTTAGGACAA".to_string(),
            post_splice_sequence: None,
            umi: "TTTTCCTAGGATATGGCTCCATAACTTAGGACAA".to_string(),
            umi_family: None,
            splice_category: SpliceChain {
                splice_event: vec!["test".to_string()],
            },
            size_class: None,
        };

        let config = SpliceConfig::build("nl43".to_string(), 2, SpliceAssayType::SizeSpecific).unwrap();
        let config2 = SpliceConfig::build("nl43".to_string(), 2, SpliceAssayType::KMer).unwrap();   
        let config3 = SpliceConfig::build("nl43".to_string(), 2, SpliceAssayType::RandomReverse).unwrap();

        splice_event1.find_size_class(&config).unwrap();
        assert_eq!(splice_event1.size_class, Some(SizeClass::Unknown));
        splice_event1.find_size_class(&config2).unwrap();
        assert_eq!(splice_event1.size_class, Some(SizeClass::Unknown));
        splice_event1.find_size_class(&config3).unwrap();
        assert_eq!(splice_event1.size_class, Some(SizeClass::Unknown));


        let mut splice_event2 = splice_event1.clone();
        splice_event2.add_post_splice_sequence(search_seq_before_d4.clone());
        assert_eq!(splice_event2.post_splice_sequence, Some("TCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTG".to_string()));
        splice_event2.find_size_class(&config).unwrap(); // should be 4kb if it is size specific
        assert_eq!(splice_event2.size_class, Some(SizeClass::FourKb));

        splice_event2.find_size_class(&config2).unwrap(); // should be 4kb if it is kmer
        assert_eq!(splice_event2.size_class, Some(SizeClass::FourKb));

        splice_event2.find_size_class(&config3).unwrap(); // could be both class if it is random reverse
        assert_eq!(splice_event2.size_class, Some(SizeClass::BothClass));

        let mut splice_event3 = splice_event1.clone();
        splice_event3.add_post_splice_sequence(search_seq_after_a7.clone());
        splice_event3.find_size_class(&config).unwrap(); // should be 1.8kb if it is size specific
        assert_eq!(splice_event3.size_class, Some(SizeClass::OnePointEightKb));
        splice_event3.find_size_class(&config2).unwrap(); // should be 1.8kb if it is kmer
        assert_eq!(splice_event3.size_class, Some(SizeClass::OnePointEightKb));
        splice_event3.find_size_class(&config3).unwrap(); // should be 1.8kb if it is random reverse
        assert_eq!(splice_event3.size_class, Some(SizeClass::OnePointEightKb));

        let mut splice_event4 = splice_event1.clone();
        splice_event4.add_post_splice_sequence(search_seq_between_d4_a7.clone());
        splice_event4.find_size_class(&config).unwrap(); // should be unknown if it is size specific
        assert_eq!(splice_event4.size_class, Some(SizeClass::Unknown));
        splice_event4.find_size_class(&config2).unwrap(); // should be 4kb if it is kmer
        assert_eq!(splice_event4.size_class, Some(SizeClass::FourKb));
        splice_event4.find_size_class(&config3).unwrap(); // should be 4kb class if it is random reverse    
        assert_eq!(splice_event4.size_class, Some(SizeClass::FourKb));

        let mut splice_event5 = splice_event1.clone();
        splice_event5.add_post_splice_sequence(search_seq_cross_a7.clone());
        splice_event5.find_size_class(&config).unwrap(); // should be 1.8kb if it is size specific
        assert_eq!(splice_event5.size_class, Some(SizeClass::OnePointEightKb));
        splice_event5.find_size_class(&config2).unwrap(); // should be 4kb if it is kmer
        assert_eq!(splice_event5.size_class, Some(SizeClass::OnePointEightKb));
        splice_event5.find_size_class(&config3).unwrap(); // should be 4kb class if it is random reverse
        assert_eq!(splice_event5.size_class, Some(SizeClass::FourKb));

        let mut splice_event6 = splice_event1.clone();
        splice_event6.splice_category.splice_event = vec!["A7".to_string()];
        splice_event6.find_size_class(&config).unwrap(); // should be 1.8kb if it is size specific
        assert_eq!(splice_event6.size_class, Some(SizeClass::OnePointEightKb));
        splice_event6.find_size_class(&config2).unwrap(); // should be 1.8kb if it is kmer
        assert_eq!(splice_event6.size_class, Some(SizeClass::OnePointEightKb));
        splice_event6.find_size_class(&config3).unwrap(); // should be 1.8kb class if it is random reverse
        assert_eq!(splice_event6.size_class, Some(SizeClass::OnePointEightKb));


    }
}