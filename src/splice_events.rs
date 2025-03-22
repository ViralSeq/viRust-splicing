use serde::{Deserialize, Serialize};
use crate::config::SpliceConfig;
use crate::joined_umi_sequence::JoinedUmiSequnce;
use crate::joined_umi_sequence::pattern_search;
use crate::config::SpliceAssayType;

#[derive(Debug, Deserialize, Serialize)]
pub struct SpliceEvents {
    pub sequence_id: String,
    pub search_sequence: String,
    pub post_splice_sequence: Option<String>, // sequence post last splicing site that can be identified. used to check size class.
    pub umi: String,
    pub umi_family: Option<String>,
    pub splice_category: SpliceChain,
    pub size_class: Option<SizeClass>,

}

#[derive(Debug, Deserialize, Serialize)]
pub struct SpliceChain {
    pub splice_event: Vec<String>,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub enum SizeClass {
    OnePointEightKb,
    FourKb ,
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

    pub fn find_size_class(&mut self, config: &SpliceConfig) -> Result<(), Box<dyn std::error::Error>> {
        let assay_type = &config.splice_assay_type;
        if let Some(post_splice_sequence) = &self.post_splice_sequence {

            if assay_type == &SpliceAssayType::SizeSpecific || assay_type == &SpliceAssayType::KMer {
                if post_splice_sequence.len() < 30 {
                    self.size_class = Some(SizeClass::Unknown);
                } else {
                    // last 15 bases are primer sequences. take 15 bases before that for locating.
                    let search_pattern = post_splice_sequence[post_splice_sequence.len() - 30..post_splice_sequence.len() - 15].to_string();
                    let aln = pattern_search(config.full_length_sequence.as_bytes(), search_pattern.as_bytes(), config.distance);
                    if let Some(aln) = aln {

                        if assay_type == &SpliceAssayType::SizeSpecific {
                            if aln.ystart < config.d4_breakpoint_position  {
                                self.size_class = Some(SizeClass::FourKb);
                            } else if aln.ystart > config.a7_breakpoint_position {
                                self.size_class = Some(SizeClass::OnePointEightKb);
                            } else {
                                self.size_class = Some(SizeClass::Unknown);}
                        } else { // KMer
                            if aln.ystart > config.a7_breakpoint_position {
                                self.size_class = Some(SizeClass::OnePointEightKb);
                            } else {
                                self.size_class = Some(SizeClass::FourKb);
                            }
                        }

                    } else {
                        self.size_class = Some(SizeClass::Unknown);
                    }
                }
            } else if assay_type == &SpliceAssayType::RandomReverse {
                self.size_class = Some(SizeClass::Unknown); // not implemented
            } else {
                return Err("Unknown assay type".into());
            }

        } else {
            return Err("No post splice sequence found".into());
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
