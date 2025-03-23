use rayon::prelude::*;
use itertools::Itertools;

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
}