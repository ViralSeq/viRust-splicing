impl JoinedUmiSequnce {
    pub fn check_splice_event(&self, splice_config: &SpliceConfig) -> SpliceEvents {
        let seq = self.find_sequence_for_search();
        let mut sequence_to_search = seq.as_bytes();
        let distance = splice_config.distance;
        let mut splice_chain = SpliceChain::new();

        // Recursive function to process downstream events dynamically
        fn process_downstream<'a>(
            sequence: &'a [u8],
            splice_chain: &mut SpliceChain,
            splice_config: &SpliceConfig,
            distance: u8,
            current_to_all: &Vec<SpliceStep>,
        ) -> &'a [u8] {
            let mut current_sequence = sequence;

            while let Some((matched_acceptor, trimmed)) =
                pattern_search_trim_seq_batch(current_sequence, current_to_all, distance)
            {
                splice_chain.add_splice_event(matched_acceptor.clone());
                current_sequence = trimmed;

                match matched_acceptor.as_str() {
                    "A1" => {
                        // Check D2
                        if let Some(trimmed) = pattern_search_trim_seq(
                            current_sequence,
                            splice_config.d2.as_bytes(),
                            distance,
                        ) {
                            splice_chain.add_splice_event("D2".to_string());
                            current_sequence = trimmed;

                            // Process downstream of D2
                            current_sequence = process_downstream(
                                current_sequence,
                                splice_chain,
                                splice_config,
                                distance,
                                &splice_config.d2_to_all,
                            );
                        } else {
                            splice_chain.add_splice_event("noD2".to_string());
                        }
                    }
                    "D2-unspliced" => {
                        // Check D2b
                        if let Some(trimmed) = pattern_search_trim_seq(
                            current_sequence,
                            splice_config.d2b.as_bytes(),
                            distance,
                        ) {
                            splice_chain.add_splice_event("D2b".to_string());
                            current_sequence = trimmed;

                            // Process downstream of D2b
                            current_sequence = process_downstream(
                                current_sequence,
                                splice_chain,
                                splice_config,
                                distance,
                                &splice_config.d2b_to_all,
                            );
                        } else {
                            splice_chain.add_splice_event("noD2b".to_string());
                        }
                    }
                    "A2" => {
                        // Check D3
                        if let Some(trimmed) = pattern_search_trim_seq(
                            current_sequence,
                            splice_config.d3.as_bytes(),
                            distance,
                        ) {
                            splice_chain.add_splice_event("D3".to_string());
                            current_sequence = trimmed;

                            // Process downstream of D3
                            current_sequence = process_downstream(
                                current_sequence,
                                splice_chain,
                                splice_config,
                                distance,
                                &splice_config.d3_to_all,
                            );
                        } else {
                            splice_chain.add_splice_event("noD3".to_string());
                        }
                    }
                    _ => {
                        // Handle other downstream acceptors (A3-A7)
                        break;
                    }
                }
            }

            current_sequence
        }

        // Step 1: Check D1
        if let Some(trimmed) =
            pattern_search_trim_seq(sequence_to_search, splice_config.d1.as_bytes(), distance)
        {
            splice_chain.add_splice_event("D1".to_string());
            sequence_to_search = trimmed;

            // Step 2: Process downstream of D1
            sequence_to_search = process_downstream(
                sequence_to_search,
                &mut splice_chain,
                splice_config,
                distance,
                &splice_config.d1_to_all,
            );
        } else {
            splice_chain.add_splice_event("noD1".to_string());
        }

        // Finalize the splice event
        let mut event = SpliceEvents::from_joined_umi_with_event(&self, splice_chain);
        event.add_post_splice_sequence(String::from_utf8(sequence_to_search.to_vec()).unwrap());

        event
    }
}