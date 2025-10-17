/// The `SpliceEvents` struct represents a splicing event with associated metadata.
///
/// # Fields
/// - `sequence_id` (`String`): The identifier for the sequence.
/// - `search_sequence` (`String`): The sequence used for searching.
/// - `post_splice_sequence` (`Option<String>`): The sequence after the last splicing site, used to check size class.
/// - `umi` (`String`): The unique molecular identifier.
/// - `umi_family` (`Option<String>`): The UMI family, if identified.
/// - `splice_category` (`SpliceChain`): The category of splicing events.
/// - `size_class` (`Option<SizeClass>`): The size class of the splicing event.
/// - `final_category` (`Option<String>`): The final predicted category of the splicing event.
///
/// # Methods
/// - `from_joined_umi_with_event`: Creates a `SpliceEvents` instance from a `JoinedUmiSequnce` and a `SpliceChain`.
/// - `add_post_splice_sequence`: Adds a post-splice sequence to the event.
/// - `find_size_class`: Determines the size class of the splicing event based on the configuration.
/// - `predict_final_category`: Predicts the final category of the splicing event based on its size class and splice events.
///
/// The struct also implements the `Display` trait for easy string formatting.
///
/// # Associated Types
/// - `SpliceChain`: Represents a chain of splicing events.
/// - `SizeClass`: Enum representing the size class of the splicing event.
///
/// # Utility Functions
/// - `find_umi_family_from_events`: Groups splicing events by their final category and assigns UMI families.
/// - `group_by_final_category`: Groups splicing events into a `HashMap` by their final category.
/// - `slice_from_end`: Extracts slices from the end of a string based on chunk size, offset, and minimum length.
///
/// # Testing
/// The module includes unit tests for:
/// - `slice_from_end`: Validates slicing logic.
/// - `find_size_class`: Tests size class determination under various configurations.
/// - `Display` implementations for `SizeClass` and `SpliceEvents`.
/// - `find_umi_family_from_events`: Ensures correct grouping and UMI family assignment.
///
/// # Notes
/// - The `find_umi_family_from_events` function uses multiple `HashMap` and `Vec` structures, which may not be optimal.
/// - Consider optimizing the implementation for better performance.
///
/// MARK: Imports
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fmt::{self, Display, Formatter};

use crate::config::{SpliceAssayType, SpliceConfig};
use crate::joined_umi_sequence::{JoinedUmiSequnce, pattern_search};
use crate::umi::find_umi_family;

/// MARK: SpliceEvents
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct SpliceEvents {
    pub sequence_id: String,
    pub search_sequence: String,
    pub post_splice_sequence: Option<String>, // sequence post last splicing site that can be identified. used to check size class.
    pub umi: String,
    pub umi_family: Option<String>,
    pub splice_category: SpliceChain,
    pub size_class: Option<SizeClass>,
    pub final_category: Option<String>,
    // in case that we do not find D1 but we find the downstream acceptor (A1, etc),
    // we should have a field to store the sequence information before the first splicing site.
    // This is the alternative d1 donor site.
    pub alternative_d1: Option<String>,
    // sometimes, we do not know what happened after D1, they are not the known events (:d1_to_all),
    // we should populate the sequence right after D1 (but known) for future analysis.
    pub unknown_sequence_after_d1: Option<String>,
}

/// MARK: SpliceChain
/// Represents a chain of splicing events.
/// for example: D1_A1_D2_A2_D3_A3 (the struct itself is a vector of strings, but can be joined by _ to make a string)
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
    Unspliced,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct ProcessedSpliceRec<'a> {
    pub seq: &'a [u8],
    pub alternative_d1: Option<&'a [u8]>, // the alternative donor site if D1 is not found. Only store the sequence when we find an acceptor without D1.
    pub unknown_sequence_after_d1: Option<&'a [u8]>,
}

impl ProcessedSpliceRec<'_> {
    pub fn init<'a>(seq: &'a [u8]) -> ProcessedSpliceRec<'a> {
        ProcessedSpliceRec {
            seq,
            alternative_d1: None,
            unknown_sequence_after_d1: None,
        }
    }
    pub fn init_with_option_fields<'a>(
        seq: &'a [u8],
        alternative_d1: Option<&'a [u8]>,
        unknown_sequence_after_d1: Option<&'a [u8]>,
    ) -> ProcessedSpliceRec<'a> {
        ProcessedSpliceRec {
            seq,
            alternative_d1,
            unknown_sequence_after_d1,
        }
    }
}

/// MARK: Impl SpliceEvents
impl SpliceEvents {
    pub fn from_joined_umi_with_event(
        joined_umi: &JoinedUmiSequnce,
        splice_category: SpliceChain,
    ) -> SpliceEvents {
        let sequence_id = joined_umi.sequence_id.clone();
        let umi = joined_umi.umi.clone();
        let umi_family = None;
        let search_sequence = joined_umi.find_sequence_for_search().clone();
        let post_splice_sequence = None;
        let size_class = None;
        let final_category = None;
        let alternative_d1 = None;
        let unknown_sequence_after_d1 = None;

        SpliceEvents {
            sequence_id,
            search_sequence,
            post_splice_sequence,
            umi,
            umi_family,
            splice_category,
            size_class,
            final_category,
            alternative_d1,
            unknown_sequence_after_d1,
        }
    }

    pub fn add_post_splice_sequence(&mut self, sequence: String) {
        self.post_splice_sequence = Some(sequence);
    }

    pub fn add_alternative_d1(&mut self, sequence: Option<&[u8]>) {
        self.alternative_d1 = sequence.map(|s| String::from_utf8_lossy(s).into_owned());
    }

    pub fn add_unknown_sequence_after_d1(&mut self, sequence: Option<&[u8]>) {
        self.unknown_sequence_after_d1 = sequence.map(|s| String::from_utf8_lossy(s).into_owned());
    }

    /// MARK: find_size_class
    pub fn find_size_class(&mut self, config: &SpliceConfig) -> Result<(), Box<dyn Error>> {
        for event in self.splice_category.splice_event.iter() {
            if event == "A7" {
                // if we find an A7 event, we know it is 1.8kb, nef transcript
                self.size_class = Some(SizeClass::OnePointEightKb);
                return Ok(());
            } else if event == "unknown" || event == "no-acceptor" {
                self.size_class = Some(SizeClass::Unknown);
                return Ok(());
            } else if event == "D1-unspliced" || event == "gag-AUG" {
                self.size_class = Some(SizeClass::Unspliced);
                return Ok(());
            }
        }

        let assay_type = &config.splice_assay_type;
        let post_splice_sequence = match &self.post_splice_sequence {
            Some(seq) => seq,
            None => {
                self.size_class = Some(SizeClass::Unknown);
                return Ok(());
            }
        };

        match assay_type {
            SpliceAssayType::SizeSpecific | SpliceAssayType::KMer => {
                if post_splice_sequence.len() < 30 {
                    self.size_class = Some(SizeClass::Unknown);
                    return Ok(());
                }

                let search_pattern = &post_splice_sequence
                    [post_splice_sequence.len() - 30..post_splice_sequence.len() - 15];
                if let Some(aln) = pattern_search(
                    config.full_length_sequence.as_bytes(),
                    search_pattern.as_bytes(),
                    config.distance,
                ) {
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
                let slice_of_patterns = slice_from_end(post_splice_sequence, 20, 15, 10);
                // dbg!(&slice_of_patterns);
                if slice_of_patterns.is_empty() {
                    self.size_class = Some(SizeClass::Unknown);
                    return Ok(());
                }

                let mut found_1_8k = false;
                let mut found_both = false;

                for pattern in slice_of_patterns[..slice_of_patterns.len().min(10)].iter() {
                    // at most 10 patterns will be used to reduce run time.
                    if let Some(aln) = pattern_search(
                        config.full_length_sequence.as_bytes(),
                        pattern.as_bytes(),
                        config.distance,
                    ) {
                        // println!("{:?}", aln.ystart);
                        // println!("{:?}", pattern);
                        if aln.ystart < config.a7_breakpoint_position
                            && aln.ystart > config.d4_breakpoint_position
                        {
                            self.size_class = Some(SizeClass::FourKb); // if we find a 4kb pattern (between D4 and A7), we can stop
                            return Ok(());
                        } else if aln.ystart > config.a7_breakpoint_position {
                            // println!("found 1.8kb");
                            found_1_8k = true;
                        } else if aln.ystart < config.d4_breakpoint_position {
                            // println!("found both");
                            found_both = true;
                        }
                    }
                }

                if found_1_8k {
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

    pub fn predict_final_category(&mut self) {
        let mut key_event = self.splice_category.splice_event.join("_"); // TODO: add a refining function before joining
        let key_class = match self.size_class.as_ref() {
            Some(class) => class.to_string(),
            None => "Unknown".to_string(),
        };
        key_event = key_event + "_" + &key_class;
        self.final_category = Some(key_event);
    }
}

/// MARK: Impl SizeClass
/// Implements the `Display` trait for the `SpliceEvents` struct, allowing it to be
/// formatted as a string. The output is a tab-separated string containing the following fields:
///
/// - `sequence_id`: The identifier for the sequence.
/// - `umi`: The unique molecular identifier.
/// - `umi_family`: The UMI family, or "None" if not present.
/// - `splice_category.splice_event`: A joined string of splice events separated by underscores.
/// - `size_class`: The size class, or "Unknown" if not specified.
/// - `final_category`: The final category, or "Unknown" if not specified.
///
/// This implementation ensures that `SpliceEvents` can be easily converted to a human-readable
/// string representation for logging or debugging purposes.
impl Display for SpliceEvents {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.sequence_id,
            self.umi,
            self.umi_family.as_deref().unwrap_or("None"),
            self.splice_category.splice_event.join("_"), // TODO: add a refining function before joining
            self.size_class
                .as_ref()
                .unwrap_or(&SizeClass::Unknown)
                .to_string(),
            self.final_category.as_deref().unwrap_or("Unknown"),
            self.alternative_d1.as_deref().unwrap_or("None"),
            self.unknown_sequence_after_d1.as_deref().unwrap_or("None"),
        )
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

    pub fn remove_splice_event(&mut self, event: String) {
        self.splice_event.retain(|e| e != &event);
    }
}

impl Display for SizeClass {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            SizeClass::OnePointEightKb => write!(f, "1.8kb"),
            SizeClass::FourKb => write!(f, "4kb"),
            SizeClass::BothClass => write!(f, "Both"),
            SizeClass::Unknown => write!(f, "Unknown"),
            SizeClass::Unspliced => write!(f, "Unspliced"),
        }
    }
}

/// MARK: FIND UMI FAMILY
// TODO: This function uses multiple HashMaps and Vecs to store the data. It is not very efficient. Consider optimizing it.
// add a way to find the common alternative_d1 and unknown_sequence_after_d1 for each umi family.
pub fn find_umi_family_from_events(events: Vec<SpliceEvents>) -> Vec<SpliceEvents> {
    let events_by_final_category = group_by_final_category(events);
    let mut outcome_events = Vec::new();

    for events in events_by_final_category.values() {
        let umis: Vec<&str> = events.iter().map(|event| event.umi.as_str()).collect();

        let families = find_umi_family(umis);

        outcome_events.extend(events.iter().map(|event| {
            let mut new_event = event.clone();
            new_event.umi_family = families
                .iter()
                .find(|&&x| x == event.umi)
                .map(|_| event.umi.clone());
            new_event
        }));
    }

    outcome_events
}

fn group_by_final_category(events: Vec<SpliceEvents>) -> HashMap<String, Vec<SpliceEvents>> {
    let mut category_map: HashMap<String, Vec<SpliceEvents>> = HashMap::new();

    for event in events.iter() {
        let category = event.final_category.as_ref().unwrap();
        if category_map.contains_key(category) {
            category_map.get_mut(category).unwrap().push(event.clone());
        } else {
            category_map.insert(category.clone(), vec![event.clone()]);
        }
    }

    category_map
}

fn slice_from_end(input: &str, chunk_size: u8, offset: u8, min_length: u8) -> Vec<String> {
    let mut chunks = Vec::new();

    let mut end = input.len() - offset as usize;
    while end > (min_length as usize - 1) {
        let start = if end >= chunk_size as usize {
            end - chunk_size as usize
        } else {
            0
        };
        if let Some(substr) = input.get(start..end) {
            chunks.push(substr.to_string());
            end = start;
        }
    }
    chunks
}

/// MARK: Tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::open_fasta_file;
    use itertools::Itertools;
    use tap::Pipe;
    // use std::fs::File;
    // use std::io::Write;

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

        let input2 = "ABCDEFGHIJKLMNO";

        let result2 = slice_from_end(input2, chunk_size, offset, 5);
        assert!(result2.is_empty());
    }

    #[test]
    fn test_size_class() {
        let nl43_file = "data/nl43.fasta";
        let fasta_reader = open_fasta_file(nl43_file).unwrap();
        let record = fasta_reader.records().next().unwrap().unwrap();
        let nl43_ref_seq = record.seq().to_vec();
        let search_seq_before_d4 = nl43_ref_seq[5100..5160]
            .to_vec()
            .pipe(String::from_utf8)
            .unwrap();
        let search_seq_after_a7 = nl43_ref_seq[8000..8060]
            .to_vec()
            .pipe(String::from_utf8)
            .unwrap();
        let search_seq_between_d4_a7 = nl43_ref_seq[7000..7060]
            .to_vec()
            .pipe(String::from_utf8)
            .unwrap();
        let search_seq_cross_a7 = nl43_ref_seq[7580..7700]
            .to_vec()
            .pipe(String::from_utf8)
            .unwrap();

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
            final_category: None,
            alternative_d1: None,
            unknown_sequence_after_d1: None,
        };

        let config =
            SpliceConfig::build("nl43".to_string(), 2, SpliceAssayType::SizeSpecific).unwrap();
        let config2 = SpliceConfig::build("nl43".to_string(), 2, SpliceAssayType::KMer).unwrap();
        let config3 =
            SpliceConfig::build("nl43".to_string(), 2, SpliceAssayType::RandomReverse).unwrap();

        splice_event1.find_size_class(&config).unwrap();
        assert_eq!(splice_event1.size_class, Some(SizeClass::Unknown));
        splice_event1.find_size_class(&config2).unwrap();
        assert_eq!(splice_event1.size_class, Some(SizeClass::Unknown));
        splice_event1.find_size_class(&config3).unwrap();
        assert_eq!(splice_event1.size_class, Some(SizeClass::Unknown));

        let mut splice_event2 = splice_event1.clone();
        splice_event2.add_post_splice_sequence(search_seq_before_d4.clone());
        assert_eq!(
            splice_event2.post_splice_sequence,
            Some("TCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTG".to_string())
        );
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

    #[test]
    fn test_display_size_class() {
        assert_eq!(SizeClass::OnePointEightKb.to_string(), "1.8kb");
        assert_eq!(SizeClass::FourKb.to_string(), "4kb");
        assert_eq!(SizeClass::BothClass.to_string(), "Both");
        assert_eq!(SizeClass::Unknown.to_string(), "Unknown");
    }

    #[test]
    fn test_display_splice_events() {
        let mut splice_event = SpliceEvents {
            sequence_id: "test".to_string(),
            search_sequence: "TTTTCCTAGGATATGGCTCCATAACTTAGGACAA".to_string(),
            post_splice_sequence: None,
            umi: "TTTTCCTAGGATATGGCTCCATAACTTAGGACAA".to_string(),
            umi_family: None,
            splice_category: SpliceChain {
                splice_event: vec![
                    "D1".to_string(),
                    "A1".to_string(),
                    "D2".to_string(),
                    "A2".to_string(),
                    "D3".to_string(),
                    "A3".to_string(),
                ],
            },
            size_class: Some(SizeClass::FourKb),
            final_category: None,
            alternative_d1: None,
            unknown_sequence_after_d1: None,
        };
        splice_event.predict_final_category();

        assert_eq!(
            splice_event.to_string(),
            "test\tTTTTCCTAGGATATGGCTCCATAACTTAGGACAA\tNone\tD1_A1_D2_A2_D3_A3\t4kb\tD1_A1_D2_A2_D3_A3_4kb\tNone\tNone"
        );
    }

    #[test]

    fn test_find_umi_family_from_events() {
        let mut events: Vec<SpliceEvents> = Vec::new();
        let mut events_template_1 = SpliceEvents {
            sequence_id: "test".to_string(),
            search_sequence: "TTTTCCTAGGATATGGCTCCATAACTTAGGACAA".to_string(),
            post_splice_sequence: None,
            umi: "AAAAA".to_string(),
            umi_family: None,
            splice_category: SpliceChain {
                splice_event: vec![
                    "D1".to_string(),
                    "A1".to_string(),
                    "D2".to_string(),
                    "A2".to_string(),
                    "D3".to_string(),
                    "A3".to_string(),
                ],
            },
            size_class: Some(SizeClass::FourKb),
            final_category: Some("D1_A1_D2_A2_D3_A3_4kb".to_string()),
            alternative_d1: None,
            unknown_sequence_after_d1: None,
        };
        let mut events_template_2 = SpliceEvents {
            sequence_id: "test".to_string(),
            search_sequence: "TTTTCCTAGGATATGGCTCCATAACTTAGGACAA".to_string(),
            post_splice_sequence: None,
            umi: "AAAAA".to_string(),
            umi_family: Some("AAAAA".to_string()),
            splice_category: SpliceChain {
                splice_event: vec![
                    "D1".to_string(),
                    "A2".to_string(),
                    "D3".to_string(),
                    "A7".to_string(),
                ],
            },
            size_class: Some(SizeClass::OnePointEightKb),
            final_category: Some("D1_A2_D3_A7_1.8kb".to_string()),
            alternative_d1: None,
            unknown_sequence_after_d1: None,
        };

        for _ in 0..1000 {
            events.push(events_template_1.clone());
            events.push(events_template_2.clone());
        }

        events_template_1.umi = "TTTTT".to_string();
        events_template_2.umi = "TTTTT".to_string();
        for _ in 0..50 {
            events.push(events_template_1.clone());
            events.push(events_template_2.clone());
        }

        events_template_1.umi = "CCCCC".to_string();
        events_template_2.umi = "CCCCC".to_string();
        for _ in 0..5 {
            events.push(events_template_1.clone());
            events.push(events_template_2.clone());
        }

        events_template_1.umi = "GGGGG".to_string();
        events_template_2.umi = "GGGGG".to_string();
        for _ in 0..500 {
            events.push(events_template_1.clone());
            events.push(events_template_2.clone());
        }

        let result = find_umi_family_from_events(events.clone());

        let family: Vec<String> = result
            .iter()
            .map(|event| {
                event
                    .umi_family
                    .as_ref()
                    .unwrap_or(&"None".to_string())
                    .to_string()
            })
            .collect();
        let uniq_family: Vec<String> = family.into_iter().unique().collect();

        assert_eq!(uniq_family.len(), 4);
        assert_eq!(uniq_family, vec!["AAAAA", "TTTTT", "None", "GGGGG",]);

        let group = group_by_final_category(events);

        let mut group_keys = group.keys().collect::<Vec<_>>();
        group_keys.sort();
        assert_eq!(
            group_keys,
            vec!["D1_A1_D2_A2_D3_A3_4kb", "D1_A2_D3_A7_1.8kb"]
        );

        // TODO: was trying to write the lines to a file but run into warnings. test file write later.
        // let mut file = File::create("output.tsv").expect("Unable to create file");
        // writeln!(file,"sequence_id\tumi\tumi_family\tsplice_category\tsize_class\tfinal_category").expect("Unable to write to file");
        // for res in result.iter() {
        //     writeln!(file, "{}", res.to_string()).expect("Unable to write to file");
        // }

        assert_eq!(result.len(), 3110);
    }

    // #[test]
    // fn test_remove_events() {
    //     let mut events = vec![
    //         "D1".to_string(),
    //         "unknown".to_string(),
    //         "A1".to_string(),
    //         "D2".to_string(),
    //         "A2".to_string(),
    //     ];
    //     let refined = refine_events(&mut events);
    //     assert_eq!(
    //         refined,
    //         &vec![
    //             "D1".to_string(),
    //             "A1".to_string(),
    //             "D2".to_string(),
    //             "A2".to_string()
    //         ]
    //     );

    //     let mut events2 = vec![
    //         "D1".to_string(),
    //         "A1".to_string(),
    //         "unknown".to_string(),
    //         "D2".to_string(),
    //         "A2".to_string(),
    //     ];
    //     let refined2 = refine_events(&mut events2);
    //     assert_eq!(
    //         refined2,
    //         &vec![
    //             "D1".to_string(),
    //             "A1".to_string(),
    //             "unknown".to_string(),
    //             "D2".to_string(),
    //             "A2".to_string()
    //         ]
    //     );

    //     let mut events3 = vec!["unknown".to_string(), "D1".to_string()];
    //     let refined3 = refine_events(&mut events3);
    //     assert_eq!(refined3, &vec!["unknown".to_string(), "D1".to_string(),]);

    //     let mut events4 = vec!["D1".to_string()];
    //     let refined4 = refine_events(&mut events4);
    //     assert_eq!(refined4, &vec!["D1".to_string(),]);
    // }
}
