//! IUPAC-aware sequence filters used before splice-event classification.

/// Returns whether `pattern` contains at least one valid IUPAC nucleotide code and no invalid
/// characters.
pub fn is_iupac_sequence(pattern: &str) -> bool {
    !pattern.is_empty() && pattern.bytes().all(|base| iupac_mask(base).is_some())
}

/// Matches an IUPAC pattern against the beginning of a sequence with a mismatch allowance.
///
/// Two positions match when their IUPAC nucleotide sets overlap. `N` matches every nucleotide,
/// which makes it suitable for UMI positions in a read header. The pattern length determines the
/// prefix length; sequences shorter than the pattern do not match.
pub fn matches_iupac_prefix(sequence: &str, pattern: &str, max_mismatches: u8) -> bool {
    if !is_iupac_sequence(pattern) || sequence.len() < pattern.len() {
        return false;
    }

    pattern
        .bytes()
        .zip(sequence.bytes())
        .filter(|(pattern_base, sequence_base)| {
            let pattern_mask = iupac_mask(*pattern_base).expect("pattern was validated above");
            let sequence_mask = iupac_mask(*sequence_base).unwrap_or(0);
            pattern_mask & sequence_mask == 0
        })
        .take(max_mismatches as usize + 1)
        .count()
        <= max_mismatches as usize
}

fn iupac_mask(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0b0001),
        b'C' => Some(0b0010),
        b'G' => Some(0b0100),
        b'T' | b'U' => Some(0b1000),
        b'R' => Some(0b0101),
        b'Y' => Some(0b1010),
        b'S' => Some(0b0110),
        b'W' => Some(0b1001),
        b'K' => Some(0b1100),
        b'M' => Some(0b0011),
        b'B' => Some(0b1110),
        b'D' => Some(0b1101),
        b'H' => Some(0b1011),
        b'V' => Some(0b0111),
        b'N' => Some(0b1111),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_matches_iupac_prefix() {
        assert!(matches_iupac_prefix("ACGTACGT", "ACGN", 0));
        assert!(matches_iupac_prefix("AGCTACGT", "ARCN", 0));
        assert!(matches_iupac_prefix("ACGTACGT", "ACAT", 1));
        assert!(!matches_iupac_prefix("ACGTACGT", "ACAT", 0));
        assert!(!matches_iupac_prefix("ACG", "ACGT", 0));
    }

    #[test]
    fn test_iupac_sequence_validation() {
        assert!(is_iupac_sequence("ARYN"));
        assert!(is_iupac_sequence("acgu"));
        assert!(!is_iupac_sequence("ACGX"));
        assert!(!is_iupac_sequence(""));
    }
}
