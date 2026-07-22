use itertools::Itertools;
use plotters::prelude::*;
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;
use rand_distr::Normal;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::error::Error;

use crate::config::{DEFAULT_UMI_AUTO_NEIGHBOR_RISK, UMI_AUTO_COVERAGE_RATIO, UmiFamilyModel};

/// Per-category diagnostics and the UMI-family model selected for that category.
#[derive(Debug, Clone)]
pub struct UmiFamilyModelDecision {
    pub effective_model: UmiFamilyModel,
    pub unique_umi_count: usize,
    pub unique_umi_fraction: f64,
    pub neighbor_risk: Option<f64>,
}

/// Selects the effective UMI-family model for one final category.
///
/// Auto mode uses max-model when unique UMIs account for at most 10% of category reads, or when
/// the expected number of random equal-length UMI neighbors reaches the configured threshold.
/// Otherwise it uses the standard cluster-model at `max_mismatches`.
pub fn select_umi_family_model(
    umis: &[&str],
    requested_model: &UmiFamilyModel,
    max_mismatches: u8,
    auto_neighbor_risk: f64,
) -> UmiFamilyModelDecision {
    let mut unique_umis_by_length: HashMap<usize, HashSet<&str>> = HashMap::new();
    let mut all_umis_are_acgt = true;

    for &umi in umis {
        all_umis_are_acgt &= umi
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'));
        unique_umis_by_length
            .entry(umi.len())
            .or_default()
            .insert(umi);
    }

    let unique_umi_count = unique_umis_by_length.values().map(HashSet::len).sum();
    let unique_umi_fraction = if umis.is_empty() {
        0.0
    } else {
        unique_umi_count as f64 / umis.len() as f64
    };

    if *requested_model != UmiFamilyModel::AutoModel {
        return UmiFamilyModelDecision {
            effective_model: requested_model.clone(),
            unique_umi_count,
            unique_umi_fraction,
            neighbor_risk: None,
        };
    }

    // The random-neighbor model assumes each UMI base is one of A, C, G, or T.
    let neighbor_risk = all_umis_are_acgt.then(|| {
        unique_umis_by_length
            .iter()
            .map(|(&length, umis)| random_neighbor_risk(umis.len(), length, max_mismatches))
            .fold(0.0_f64, f64::max)
    });
    let effective_model = if unique_umi_fraction <= UMI_AUTO_COVERAGE_RATIO
        || neighbor_risk.is_none_or(|risk| risk >= auto_neighbor_risk)
    {
        UmiFamilyModel::MaxModel
    } else {
        UmiFamilyModel::ClusterModel
    };

    UmiFamilyModelDecision {
        effective_model,
        unique_umi_count,
        unique_umi_fraction,
        neighbor_risk,
    }
}

fn random_neighbor_risk(unique_umi_count: usize, umi_length: usize, max_mismatches: u8) -> f64 {
    if unique_umi_count <= 1 || max_mismatches == 0 {
        return 0.0;
    }
    if max_mismatches as usize >= umi_length {
        return (unique_umi_count - 1) as f64;
    }

    let mut combinations = 1.0_f64;
    let mut substitutions = 1.0_f64;
    let mut neighbor_count = 0.0_f64;
    for mismatches in 1..=max_mismatches as usize {
        combinations *= (umi_length - mismatches + 1) as f64 / mismatches as f64;
        substitutions *= 3.0;
        neighbor_count += combinations * substitutions;
    }

    let sequence_space = 4.0_f64.powi(umi_length as i32);
    (unique_umi_count - 1) as f64 * neighbor_count / (sequence_space - 1.0)
}

/// Finds UMI families based on their frequency in the input list.
///
/// # Arguments
///
/// * `umis` - A vector of UMI strings.
///
/// # Returns
///
/// A vector of UMI strings that belong to families based on a frequency cutoff.
pub fn find_umi_family(mut umis: Vec<&str>) -> Vec<&str> {
    let mut families: Vec<&str> = Vec::new();

    umis.par_sort_unstable();

    let umi_freq = umis.into_iter().dedup_with_count();

    let mut freq_count: Vec<usize> = Vec::new();
    umi_freq.clone().into_iter().for_each(|(count, _)| {
        freq_count.push(count);
    });

    let max_freq = freq_count.iter().max().unwrap();

    let umi_cut_off = umi_cut_off(*max_freq);

    umi_freq.into_iter().for_each(|(count, umi)| {
        if count > umi_cut_off as usize {
            families.push(umi);
        }
    });

    families
}

/// Assigns each non-singleton UMI family to a representative.
///
/// `MaxModel` preserves the existing frequency-cutoff behavior. `ClusterModel` and
/// `ClusterModelPlusone` build abundance-ordered clusters from equal-length UMIs within the
/// configured mismatch threshold or that threshold plus one, respectively. Each UMI must be
/// directly within the threshold of its family parent. Each family must meet
/// `umi_min_fraction` of reads in its final category, and ties use lexical order. `AutoModel`
/// should normally be resolved with `select_umi_family_model` before this function is called.
pub fn find_umi_family_assignments(
    umis: Vec<&str>,
    umi_family_model: &UmiFamilyModel,
    max_mismatches: u8,
    umi_min_fraction: f64,
) -> HashMap<String, String> {
    if umis.is_empty() {
        return HashMap::new();
    }

    match umi_family_model {
        UmiFamilyModel::MaxModel => find_umi_family(umis)
            .into_iter()
            .map(|umi| (umi.to_owned(), umi.to_owned()))
            .collect(),
        UmiFamilyModel::ClusterModel => {
            cluster_umi_families(umis, max_mismatches, umi_min_fraction)
        }
        UmiFamilyModel::ClusterModelPlusone => {
            cluster_umi_families(umis, max_mismatches.saturating_add(1), umi_min_fraction)
        }
        UmiFamilyModel::AutoModel => {
            let decision = select_umi_family_model(
                &umis,
                umi_family_model,
                max_mismatches,
                DEFAULT_UMI_AUTO_NEIGHBOR_RISK,
            );
            find_umi_family_assignments(
                umis,
                &decision.effective_model,
                max_mismatches,
                umi_min_fraction,
            )
        }
    }
}

fn cluster_umi_families(
    umis: Vec<&str>,
    max_mismatches: u8,
    umi_min_fraction: f64,
) -> HashMap<String, String> {
    let mut counts: Vec<(String, usize)> = umis
        .into_iter()
        .fold(HashMap::new(), |mut counts, umi| {
            *counts.entry(umi.to_owned()).or_insert(0) += 1;
            counts
        })
        .into_iter()
        .collect();
    counts.sort_unstable_by(|left, right| left.0.cmp(&right.0));

    let mut indices_by_length: HashMap<usize, Vec<usize>> = HashMap::new();
    for (index, (umi, _)) in counts.iter().enumerate() {
        indices_by_length.entry(umi.len()).or_default().push(index);
    }

    let mut candidates_by_chunk: HashMap<(usize, usize, Vec<u8>), Vec<usize>> = HashMap::new();
    for (&length, indices) in &indices_by_length {
        if max_mismatches as usize >= length {
            continue;
        }

        // With at most d substitutions, two sequences must share at least one of d + 1 chunks.
        let chunk_count = max_mismatches as usize + 1;
        for &index in indices {
            let umi = counts[index].0.as_bytes();
            for chunk_index in 0..chunk_count {
                let (start, end) = chunk_bounds(length, chunk_index, chunk_count);
                candidates_by_chunk
                    .entry((length, chunk_index, umi[start..end].to_vec()))
                    .or_default()
                    .push(index);
            }
        }
    }

    let mut parent_order: Vec<usize> = (0..counts.len()).collect();
    parent_order.sort_unstable_by(|left, right| {
        counts[*right]
            .1
            .cmp(&counts[*left].1)
            .then_with(|| counts[*left].0.cmp(&counts[*right].0))
    });

    let mut assigned_parents = vec![None; counts.len()];
    let mut family_sizes = vec![0; counts.len()];
    let total_reads = counts.iter().map(|(_, count)| count).sum::<usize>();
    let minimum_family_size = (umi_min_fraction * total_reads as f64).ceil() as usize;

    for parent in parent_order {
        if assigned_parents[parent].is_some() {
            continue;
        }

        assigned_parents[parent] = Some(parent);
        let mut family_size = counts[parent].1;
        let parent_umi = counts[parent].0.as_bytes();
        let length = parent_umi.len();
        let candidate_indices: Vec<usize> = if max_mismatches as usize >= length {
            indices_by_length[&length].clone()
        } else {
            let chunk_count = max_mismatches as usize + 1;
            let mut candidates = HashSet::new();
            for chunk_index in 0..chunk_count {
                let (start, end) = chunk_bounds(length, chunk_index, chunk_count);
                if let Some(indices) =
                    candidates_by_chunk.get(&(length, chunk_index, parent_umi[start..end].to_vec()))
                {
                    candidates.extend(indices.iter().copied());
                }
            }
            candidates.into_iter().collect()
        };

        for candidate in candidate_indices {
            if assigned_parents[candidate].is_none()
                && hamming_distance_within(
                    parent_umi,
                    counts[candidate].0.as_bytes(),
                    max_mismatches,
                )
            {
                assigned_parents[candidate] = Some(parent);
                family_size += counts[candidate].1;
            }
        }
        family_sizes[parent] = family_size;
    }

    counts
        .iter()
        .enumerate()
        .filter_map(|(index, (umi, _))| {
            let parent = assigned_parents[index].expect("every UMI is assigned to a parent");
            if family_sizes[parent] < minimum_family_size {
                return None;
            }
            Some((umi.clone(), counts[parent].0.clone()))
        })
        .collect()
}

fn chunk_bounds(length: usize, chunk_index: usize, chunk_count: usize) -> (usize, usize) {
    let base_length = length / chunk_count;
    let remainder = length % chunk_count;
    let start = chunk_index * base_length + chunk_index.min(remainder);
    let end = start + base_length + usize::from(chunk_index < remainder);
    (start, end)
}

fn hamming_distance_within(left: &[u8], right: &[u8], max_mismatches: u8) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .filter(|(left_base, right_base)| left_base != right_base)
            .take(max_mismatches as usize + 1)
            .count()
            <= max_mismatches as usize
}

/// Generates a single UMI of the specified length using a random seed.
///
/// # Arguments
///
/// * `umi_length` - The length of the UMI to generate.
/// * `seed` - A seed value for reproducible random generation.
///
/// # Returns
///
/// A string representing the generated UMI.
pub fn generate_one_umi(umi_length: u32, seed: u64) -> String {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let bases = ['A', 'T', 'C', 'G'];
    let mut umi = String::new();

    // Generate a single UMI
    for _ in 0..umi_length {
        let base = bases
            .choose(&mut rng)
            .expect("Failed to choose a base")
            .to_owned(); // Choose a random base from the array
        umi.push(base);
    }

    umi
}

/// Generates a list of UMIs of the specified length.
///
/// # Arguments
///
/// * `umi_length` - The length of each UMI.
/// * `num_umis` - The number of UMIs to generate.
///
/// # Returns
///
/// A vector of strings representing the generated UMIs.
pub fn generate_umis(umi_length: u32, num_umis: u32) -> Vec<String> {
    let mut umis: Vec<String> = Vec::new();
    // Generate UMIs
    for i in 0..num_umis {
        let umi = generate_one_umi(umi_length, (i + 1) as u64); // Use i+1 as seed for reproducibility
        umis.push(umi);
    }
    umis
}

/// Simulates the distribution of UMIs across sequences with optional mutations.
///
/// # Arguments
///
/// * `num_sequences` - The number of sequences to simulate.
/// * `umi_length` - The length of each UMI.
/// * `umi_number` - The number of unique UMIs to generate.
/// * `seed` - A seed value for reproducible random generation.
/// * `mutation_rate` - The probability of mutating each base in a UMI.
///
/// # Returns
///
/// A `HashMap` where keys are UMIs and values are their frequencies, or an error if validation fails.
pub fn simulate_umi_distribution(
    num_sequences: u32,
    umi_length: u32,
    umi_number: u32,
    seed: u64,
    mutation_rate: f64,
) -> Result<HashMap<String, u32>, Box<dyn Error>> {
    // Validate the number of UMIs by length
    if !validate_umis_by_length(umi_number, umi_length) {
        return Err(
            "The number of UMIs exceeds the maximum possible unique UMIs for the given length."
                .into(),
        );
    }
    let mut umi_distribution: HashMap<String, u32> = HashMap::new();

    // Generate UMIs

    let umis = generate_umis(umi_length, umi_number);
    // Simulate the distribution of UMIs across the sequences

    // need to introduce more skewness in the distribution to simulate real-world data
    // For example, we can randomly select UMIs based on a weighted distribution
    let normal = Normal::new(umi_number as f64 / 2.0, umi_number as f64 / 6.0).unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(seed);

    for _ in 0..num_sequences {
        // Randomly select a UMI from the generated UMIs
        // Generate a random number from the normal distribution
        let mut random_number = normal.sample(&mut rng).round() as isize;
        random_number = random_number
            .max(0) // Ensure the index is non-negative
            .min(umi_number as isize - 1); // Ensure the index is within bounds

        let selected_umi = &umis[random_number as usize];

        // Introduce some mutation to the selected UMI

        let mutated_umi = mutate_umi(selected_umi, mutation_rate);
        // Update the distribution count for the mutated UMI
        *umi_distribution.entry(mutated_umi).or_insert(0) += 1;
    }
    Ok(umi_distribution)
}

/// Introduces mutations to a UMI based on a specified mutation rate.
///
/// # Arguments
///
/// * `umi` - The original UMI string.
/// * `mutation_rate` - The probability of mutating each base in the UMI.
///
/// # Returns
///
/// A string representing the mutated UMI.
pub fn mutate_umi(umi: &str, mutation_rate: f64) -> String {
    let mut rng = ChaCha8Rng::from_rng(&mut rand::rng());
    let bases = ['A', 'T', 'C', 'G'];
    let mut mutated_umi = String::new();

    for base in umi.chars() {
        let mut new_base_pool = bases.to_vec();
        new_base_pool.retain(|c| c != &base); // Exclude the original base
        if rng.random::<f64>() < mutation_rate {
            // Mutate the base
            let new_base = new_base_pool
                .choose(&mut rng)
                .expect("Failed to choose a base")
                .to_owned(); // Choose a random base from the vec of bases excluding the original base
            mutated_umi.push(new_base);
        } else {
            mutated_umi.push(base);
        }
    }

    mutated_umi
}

/// Validates whether the requested number of UMIs can be generated for a given length.
///
/// # Arguments
///
/// * `num_umis` - The number of UMIs to validate.
/// * `umi_length` - The length of each UMI.
///
/// # Returns
///
/// `true` if the number of UMIs is valid, `false` otherwise.
pub fn validate_umis_by_length(num_umis: u32, umi_length: u32) -> bool {
    let maxi_umis = 4.0_f64.powi(umi_length as i32); // Maximum number of unique UMIs possible
    // Check if the number of UMIs exceeds the maximum possible unique UMIs
    if num_umis as f64 > maxi_umis {
        false
    } else {
        true
    }
}

/// Plots the distribution of UMIs as a histogram and saves it to a file.
///
/// # Arguments
///
/// * `umi_distribution` - A `HashMap` where keys are UMIs and values are their frequencies.
/// * `path` - The file path to save the histogram image.
///
/// # Returns
///
/// An empty result if successful, or an error if plotting fails.
pub fn plot_umi_distribution(
    umi_distribution: &HashMap<String, u32>,
    path: &str,
) -> Result<(), Box<dyn Error>> {
    let mut dist = umi_distribution.values().cloned().collect::<Vec<u32>>();

    dist.sort();
    dist.reverse();

    // Plotting code can be added here using a plotting library
    let mut expanded_umi_list: Vec<String> = Vec::new();
    for (umi, count) in umi_distribution {
        for _ in 0..*count {
            expanded_umi_list.push(umi.clone());
        }
    }

    let root = BitMapBackend::new(path, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_bin = dist.iter().max().unwrap_or(&0) + 1;
    let bins = (0..=max_bin).step_by(1).collect::<Vec<_>>();

    let mut histogram = vec![0; bins.len() - 1];
    for &value in &dist {
        for i in 0..bins.len() - 1 {
            if value >= bins[i] && value < bins[i + 1] {
                histogram[i] += 1;
                break;
            }
        }
    }

    let mut chart = ChartBuilder::on(&root)
        .caption("UMI Distribution Histogram", ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(
            dist.iter().map(|&x| x as u32).min().unwrap()
                ..dist.iter().map(|&x| x as u32).max().unwrap(),
            0..*histogram.iter().max().unwrap(),
        )?;

    chart
        .configure_mesh()
        .x_desc("UMI Bin Size")
        .y_desc("Number of Unique UMIs")
        .draw()?;

    chart.draw_series(dist.iter().map(|&value| {
        Rectangle::new(
            [
                (value as u32, 0),
                (value as u32 + 1, histogram[value as usize]),
            ],
            BLUE.filled(),
        )
    }))?;

    Ok(())
}

/// Calculates the cutoff frequency for UMI families based on the maximum frequency.
///
/// # Arguments
///
/// * `m` - The maximum frequency of UMIs.
///
/// # Returns
///
/// An integer representing the cutoff frequency.
fn umi_cut_off(m: usize) -> i32 {
    let n: f64;
    if m <= 10 {
        n = 2.0;
    } else {
        n = -9.59e-27 * (m as f64).powi(6) + 3.27e-21 * (m as f64).powi(5)
            - 3.05e-16 * (m as f64).powi(4)
            + 1.2e-11 * (m as f64).powi(3)
            - 2.19e-7 * (m as f64).powi(2)
            + 0.004044 * (m as f64)
            + 2.273;
    }

    let mut n_rounded = n.round() as i32;
    if n_rounded < 3 {
        n_rounded = 2;
    }

    n_rounded
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::iter::repeat;

    #[test]
    fn test_umi_cut_off() {
        assert_eq!(umi_cut_off(9), 2);
        assert_eq!(umi_cut_off(100), 3);
        assert_eq!(umi_cut_off(400), 4);
        assert_eq!(umi_cut_off(1000), 6);
        assert_eq!(umi_cut_off(5000), 18);
    }

    #[test]
    fn test_find_umi_family() {
        let mut umis: Vec<&str> = Vec::new();
        umis.extend(repeat("AAAAA").take(1000));
        umis.extend(repeat("TTTTT").take(50));
        umis.extend(repeat("CCCCC").take(5));
        umis.extend(repeat("GGGGG").take(500));
        let mut families = find_umi_family(umis.clone());
        families.sort();
        assert_eq!(families.len(), 3);
        assert_eq!(families, vec!["AAAAA", "GGGGG", "TTTTT"]);

        let mut filtered_umis: Vec<Option<String>> = Vec::new();

        filtered_umis.extend(
            umis.iter()
                .map(|x| families.iter().find(|&y| x == y).map(|y| y.to_string())),
        );

        let n = filtered_umis.iter().filter(|umi| umi.is_some()).count();

        assert_eq!(n, 1550);
    }

    #[test]
    fn test_find_umi_family_2() {
        let mut umis: Vec<&str> = Vec::new();
        umis.extend(repeat("AAAAA").take(2));
        umis.extend(repeat("TTTTT").take(2));
        umis.extend(repeat("CCCCC").take(2));
        umis.extend(repeat("GGGGG").take(2));
        let mut families = find_umi_family(umis);
        families.sort();
        assert_eq!(families.len(), 0);
    }

    #[test]
    fn test_cluster_umi_family_assignments() {
        let umis = vec![
            "AAAAA", "AAAAA", "AAAAA", "AAAAT", "AAAAT", "AAATT", "GGGGG", "AAAA",
        ];

        let assignments = find_umi_family_assignments(umis, &UmiFamilyModel::ClusterModel, 1, 0.2);

        // AAAAT is directly within one mismatch of AAAAA, but AAATT is not.
        assert_eq!(assignments["AAAAA"], "AAAAA");
        assert_eq!(assignments["AAAAT"], "AAAAA");
        assert!(!assignments.contains_key("AAATT"));
        assert!(!assignments.contains_key("GGGGG"));
        assert!(!assignments.contains_key("AAAA"));
    }

    #[test]
    fn test_cluster_umi_family_assignments_do_not_chain_distant_umis() {
        let umis = vec![
            "AAAAA", "AAAAA", "AAAAA", "AAAAT", "AAAAT", "AAATT", "AAATT",
        ];

        let assignments =
            find_umi_family_assignments(umis, &UmiFamilyModel::ClusterModel, 1, 0.005);

        assert_eq!(assignments["AAAAA"], "AAAAA");
        assert_eq!(assignments["AAAAT"], "AAAAA");
        assert_eq!(assignments["AAATT"], "AAATT");
        for (umi, family) in assignments {
            assert!(hamming_distance_within(
                umi.as_bytes(),
                family.as_bytes(),
                1
            ));
        }
    }

    #[test]
    fn test_cluster_umi_family_assignments_breaks_frequency_ties_lexically() {
        let assignments = find_umi_family_assignments(
            vec!["AAAAT", "AAAAA"],
            &UmiFamilyModel::ClusterModel,
            1,
            0.005,
        );

        assert_eq!(assignments["AAAAA"], "AAAAA");
        assert_eq!(assignments["AAAAT"], "AAAAA");
    }

    #[test]
    fn test_cluster_umi_family_assignments_respects_configured_mismatches() {
        let assignments = find_umi_family_assignments(
            vec!["AAAAA", "AAAAA", "AAATT"],
            &UmiFamilyModel::ClusterModel,
            1,
            0.5,
        );

        assert_eq!(assignments["AAAAA"], "AAAAA");
        assert!(!assignments.contains_key("AAATT"));
    }

    #[test]
    fn test_cluster_umi_family_assignments_allow_one_extra_mismatch() {
        let assignments = find_umi_family_assignments(
            vec!["AAAAA", "AAAAA", "AAATT"],
            &UmiFamilyModel::ClusterModelPlusone,
            1,
            0.5,
        );

        // AAATT differs from AAAAA by two substitutions: --distance (1) plus one for clustering.
        assert_eq!(assignments["AAAAA"], "AAAAA");
        assert_eq!(assignments["AAATT"], "AAAAA");
    }

    #[test]
    fn test_cluster_umi_family_assignments_use_relative_cutoff() {
        let umis = vec!["AAAAA", "AAAAA", "CCCCC", "GGGGG", "TTTTT", "ATATA"];

        let assignments =
            find_umi_family_assignments(umis.clone(), &UmiFamilyModel::ClusterModel, 0, 0.005);
        assert_eq!(assignments.len(), 5);

        let assignments = find_umi_family_assignments(umis, &UmiFamilyModel::ClusterModel, 0, 0.2);
        assert_eq!(assignments.len(), 1);
        assert_eq!(assignments["AAAAA"], "AAAAA");
    }

    #[test]
    fn test_auto_model_uses_max_model_for_sufficient_coverage() {
        let umis = vec!["AAAAAAAAAAAAAA"; 10];
        let decision = select_umi_family_model(
            &umis,
            &UmiFamilyModel::AutoModel,
            2,
            DEFAULT_UMI_AUTO_NEIGHBOR_RISK,
        );

        assert_eq!(decision.unique_umi_fraction, UMI_AUTO_COVERAGE_RATIO);
        assert_eq!(decision.neighbor_risk, Some(0.0));
        assert_eq!(decision.effective_model, UmiFamilyModel::MaxModel);
    }

    #[test]
    fn test_auto_model_uses_cluster_model_for_sparse_low_risk_category() {
        let umis = vec![
            "AAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAA",
            "CAAAAAAAAAAAAA",
            "GAAAAAAAAAAAA",
            "TAAAAAAAAAAAA",
            "ACAAAAAAAAAAAA",
        ];
        let decision = select_umi_family_model(
            &umis,
            &UmiFamilyModel::AutoModel,
            2,
            DEFAULT_UMI_AUTO_NEIGHBOR_RISK,
        );

        assert_eq!(decision.unique_umi_count, 5);
        assert_eq!(decision.unique_umi_fraction, 5.0 / 6.0);
        assert!(decision.neighbor_risk.expect("risk should be available") < 0.05);
        assert_eq!(decision.effective_model, UmiFamilyModel::ClusterModel);
    }

    #[test]
    fn test_auto_model_uses_max_model_when_random_neighbor_risk_is_high() {
        let umis = vec![
            "AAAAAA", "AAAAAC", "AAAAAG", "AAAAAT", "AAAACA", "AAAACC", "AAAACG", "AAAACT",
            "AAAAGA", "AAAAGC", "AAAAGG", "AAAAGT", "AAAATA", "AAAATC", "AAAATG", "AAAATT",
            "AAACAA", "AAACAC", "AAACAG", "AAACAT",
        ];
        let decision = select_umi_family_model(
            &umis,
            &UmiFamilyModel::AutoModel,
            2,
            DEFAULT_UMI_AUTO_NEIGHBOR_RISK,
        );

        assert!(random_neighbor_risk(20, 6, 2) > 0.05);
        assert!(decision.neighbor_risk.expect("risk should be available") > 0.05);
        assert_eq!(decision.effective_model, UmiFamilyModel::MaxModel);

        let permissive_decision =
            select_umi_family_model(&umis, &UmiFamilyModel::AutoModel, 2, 0.8);
        assert_eq!(
            permissive_decision.effective_model,
            UmiFamilyModel::ClusterModel
        );
    }

    #[test]
    fn test_auto_model_uses_max_model_when_umi_is_not_acgt() {
        let umis = vec!["AAAAAN", "CCCCCC"];
        let decision = select_umi_family_model(
            &umis,
            &UmiFamilyModel::AutoModel,
            1,
            DEFAULT_UMI_AUTO_NEIGHBOR_RISK,
        );

        assert_eq!(decision.neighbor_risk, None);
        assert_eq!(decision.effective_model, UmiFamilyModel::MaxModel);
    }

    #[test]
    fn test_generate_umi() {
        let umi = generate_one_umi(10, 1);

        assert_eq!(umi, "CTGAACTAGT".to_string()); // fixed seed for reproducibility (seed = 1)
    }

    #[test]
    fn test_generate_umis() {
        let umi_length = 10;
        let num_umis = 5;
        let umis = generate_umis(umi_length, num_umis);

        // Check if the number of UMIs generated is correct
        assert_eq!(umis.len(), num_umis as usize);

        // Check if each UMI has the correct length
        for umi in &umis {
            assert_eq!(umi.len(), umi_length as usize);
        }

        // Since the random generation is based on a fixed seed, we can check for known output
        assert_eq!(umis[0], "CTGAACTAGT"); // This is based on the fixed seed used in generate_one_umi
        assert_eq!(umis[1], "AGCCAGAGAA");
        assert_eq!(umis[2], "ACCAGTGAGC");
        assert_eq!(umis[3], "GCCATGGGTA");
        assert_eq!(umis[4], "AATCTGTGGT");
    }

    #[test]
    fn test_validate_umis_by_length() {
        let num_umis = 1000;
        let umi_length = 5;

        // This test should pass because 1000 UMIs is less than the maximum possible unique UMIs for length 5
        let is_valid = validate_umis_by_length(num_umis, umi_length);
        assert!(
            is_valid,
            "The number of UMIs should be valid for the given length"
        );
        // Check for a case where the number of UMIs exceeds the maximum possible unique UMIs
        let num_umis_invalid = 200000; // This is more than the maximum possible unique UMIs for length 5
        let is_valid_invalid = validate_umis_by_length(num_umis_invalid, umi_length);
        assert!(
            !is_valid_invalid,
            "The number of UMIs should not be valid for the given length"
        );
    }

    #[test]
    fn test_simulate_umi_distribution() {
        let num_sequences = 10000;
        let umi_length = 10;
        let umi_number = 100; // Number of unique UMIs to generate

        // Simulate the UMI distribution
        let result = simulate_umi_distribution(num_sequences, umi_length, umi_number, 1, 0.05);

        // Check if the result is Ok
        assert!(result.is_ok(), "The simulation should succeed");

        // Get the distribution map
        let umi_distribution = result.unwrap();

        let mut dist = umi_distribution.values().cloned().collect::<Vec<u32>>();

        dist.sort();
        dist.reverse();
        // Print the UMI distribution
        println!("UMI Distribution: {:?}", &dist);

        plot_umi_distribution(&umi_distribution, "temp/umi_distribution_histogram.png").unwrap();

        assert_eq!(umi_distribution.values().sum::<u32>(), num_sequences);
        // Check if the counts are within the expected range
    }

    #[test]
    #[should_panic]
    fn test_simulate_umi_distribution_panic() {
        let num_sequences = 10000;
        let umi_length = 5;
        let umi_number = 200000; // This exceeds the maximum possible unique UMIs for length 5

        // Simulate the UMI distribution
        let result = simulate_umi_distribution(num_sequences, umi_length, umi_number, 1, 0.05);

        // Check if the result is Ok
        assert!(result.is_ok(), "The simulation should succeed");

        // Get the distribution map
        let _umi_distribution = result.unwrap();
    }
}
