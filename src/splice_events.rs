use serde::{Deserialize, Serialize};
use crate::joined_umi_sequence::JoinedUmiSequnce;

#[derive(Debug, Deserialize, Serialize)]
pub struct SpliceEvents {
    pub sequence_id: String,
    pub search_sequence: String,
    pub post_splice_sequence: String, // sequence post last splicing site that can be identified. used to check size class.
    pub umi: String,
    pub umi_family: Option<String>,
    pub splice_category: SpliceChain,

}

#[derive(Debug, Deserialize, Serialize)]
pub struct SpliceChain {
    pub splice_event: Vec<String>,
}

impl SpliceEvents {
    pub fn from_joined_umi_with_event(joined_umi: &JoinedUmiSequnce, splice_category: SpliceChain) -> SpliceEvents {
        let sequence_id = joined_umi.sequence_id.clone();
        let umi = joined_umi.umi.clone();
        let umi_family = None;
        let search_sequence = joined_umi.find_sequence_for_search().clone();
        let post_splice_sequence = String::new();

        SpliceEvents {
            sequence_id,
            search_sequence,
            post_splice_sequence,
            umi,
            umi_family,
            splice_category,
        }
    }

    pub fn add_post_splice_sequence(&mut self, sequence: String) {
        self.post_splice_sequence = sequence;
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
