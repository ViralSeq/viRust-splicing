use rayon::prelude::*;
use itertools::Itertools;
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;
use std::collections::HashMap;
use std::error::Error;
use rand_distr::Normal;
use plotters::prelude::*;

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
        if count > umi_cut_off as usize{
            families.push(umi);
        }
    });

    families
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
pub fn generate_one_umi(
    umi_length: u32,
    seed:u64
) -> String {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let bases = ['A', 'T', 'C', 'G'];
    let mut umi = String::new();

    // Generate a single UMI
    for _ in 0..umi_length {
        let base = bases.choose(&mut rng)
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
pub fn generate_umis(
    umi_length: u32,
    num_umis: u32
) -> Vec<String> {
    let mut umis: Vec<String> = Vec::new();
    // Generate UMIs
    for i in 0..num_umis {
        let umi = generate_one_umi(umi_length, (i+1) as u64); // Use i+1 as seed for reproducibility
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
    seed:u64,
    mutation_rate: f64
) -> Result<HashMap<String, u32>, Box<dyn Error>> {
    // Validate the number of UMIs by length
    if !validate_umis_by_length(umi_number, umi_length) {
        return Err("The number of UMIs exceeds the maximum possible unique UMIs for the given length.".into());
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
        random_number = random_number.max(0) // Ensure the index is non-negative
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
pub fn mutate_umi(
    umi: &str,
    mutation_rate: f64,
) -> String {
    let mut rng = ChaCha8Rng::from_os_rng();
    let bases = ['A', 'T', 'C', 'G'];
    let mut mutated_umi = String::new();

    for base in umi.chars() {
        let mut new_base_pool = bases.to_vec();
        new_base_pool.retain(|c| c != &base); // Exclude the original base
        if rng.random::<f64>() < mutation_rate {
            // Mutate the base
            let new_base = new_base_pool.choose(&mut rng)
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
pub fn validate_umis_by_length(
    num_umis: u32,
    umi_length: u32,
) -> bool {

    let maxi_umis = 4.0_f64.powi(umi_length as i32); // Maximum number of unique UMIs possible
    // Check if the number of UMIs exceeds the maximum possible unique UMIs
    if num_umis as f64  > maxi_umis {
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
    path: &str
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
       
        let root = BitMapBackend::new(path, (640, 480))
            .into_drawing_area();
        root.fill(&WHITE)?;

        let max_bin = dist.iter().max().unwrap_or(&0) +1;
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
            dist.iter().map(|&x| x as u32).min().unwrap()..dist.iter().map(|&x| x as u32).max().unwrap(),
            0..*histogram.iter().max().unwrap(),
            )?;

        chart
            .configure_mesh()
            .x_desc("UMI Bin Size")
            .y_desc("Number of Unique UMIs")
            .draw()?;
 

        chart
            .draw_series(
            dist.iter()
            .map(|&value| {
                Rectangle::new([(value as u32, 0), (value as u32 + 1, histogram[value as usize])], BLUE.filled())
            }),
            )?;
   

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
    let n :f64;
    if m <= 10 {
        n = 2.0;
    } else {
        n = -9.59e-27 * (m as f64).powi(6)
            + 3.27e-21 * (m as f64).powi(5)
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
            umis
            .iter()
            .map(
                |x| {families.iter().find(|&y|x == y).map(|y|y.to_string())
                }
            )
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
    fn test_generate_umi() {
        let umi = generate_one_umi(10, 1);

        assert_eq!(umi, "CTGAACTAGT".to_string()); // fixed seed for reproducibility (seed = 1)
    }

    #[test]
    fn test_generate_umis(){
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
        assert!(is_valid, "The number of UMIs should be valid for the given length");
        // Check for a case where the number of UMIs exceeds the maximum possible unique UMIs
        let num_umis_invalid = 200000; // This is more than the maximum possible unique UMIs for length 5
        let is_valid_invalid = validate_umis_by_length(num_umis_invalid, umi_length);
        assert!(!is_valid_invalid, "The number of UMIs should not be valid for the given length");
    }

    #[test]
    fn test_simulate_umi_distribution() {
        let num_sequences = 10000;
        let umi_length = 10;
        let umi_number = 100; // Number of unique UMIs to generate

        // Simulate the UMI distribution
        let result = simulate_umi_distribution(
            num_sequences, umi_length, umi_number, 1, 0.05);

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
        let result = simulate_umi_distribution(
            num_sequences, umi_length, umi_number, 1, 0.05);

        // Check if the result is Ok
        assert!(result.is_ok(), "The simulation should succeed");

        // Get the distribution map
        let _umi_distribution = result.unwrap();
    }
}
