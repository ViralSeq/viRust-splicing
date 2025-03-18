use config::Config;
use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::process;

#[derive(Debug)]
pub struct InputConfig {
    pub query: String,
    pub distance: u8,
    pub filename_r1: String,
    pub filename_r2: String,
}

#[derive(Debug)]
pub struct SpliceConfig {
    pub strain: String,
    pub distance: u8,
    pub query_list: HashMap<String, String>,
    pub d1: String,
    pub d2: String,
    pub d2b: String,
    pub d2b_unspliced: String,
    pub d3: String,
    pub d1_to_all: Vec<SpliceStep>,
    pub d2_to_all: Vec<SpliceStep>,
    pub d2b_to_all: Vec<SpliceStep>,
    pub d3_to_all: Vec<SpliceStep>,
}

#[derive(Debug)]
pub struct SpliceStep {
    pub donor: String,
    pub acceptor: String,
    pub pattern: String,
}


impl InputConfig {
    pub fn build() -> Result<InputConfig, Box<dyn Error>> {
        let args: Vec<String> = env::args().collect();
        if args.len() != 5 {
            return Err("Wrong number of Arguments\n Usage: cargo run <query_tag> <distance> <r1_file_path> <r1_file_path>".into());
        }

        let query = args[1].clone();

        let distance = 0;

        match args[2].parse::<u8>() {
            Ok(distance) => distance,
            Err(_) => {
                return Err("Distance must be a number".into());
            }
        };

        let filename_r1 = args[3].clone();
        let filename_r2 = args[4].clone();

        Ok(InputConfig {
            query,
            distance,
            filename_r1,
            filename_r2,
        })
    }
}

impl SpliceConfig {
    pub fn build(strain: String, distance: u8) -> Result<SpliceConfig, Box<dyn Error>> {
        let splice_form_config: HashMap<String, HashMap<String, String>> = Config::builder()
            .add_source(config::File::with_name("data/splice_form_config.toml"))
            .build()?
            .try_deserialize::<HashMap<String, HashMap<String, String>>>()?;

        let query_list = splice_form_config
            .get(&strain)
            .unwrap_or_else(|| {
                eprintln!("Strain {} not found in config file", strain);
                process::exit(1);
            })
            .clone();

        // first step check D1
        let d1 = query_list.get("D1").unwrap().clone();

        // this is the most common acceptor sites list

        let mut common_acceptor_sites = vec!["A3", "A4a", "A4b", "A4c", "A4d", "A5a", "A5b", "A7"];

        // second step check D1-unspliced and all the A sites

        let mut second_step_acceptor_sites = vec!["D1-unspliced", "A1", "A2"];

        second_step_acceptor_sites.append(&mut common_acceptor_sites);

        let second_step = SpliceStep::build("D1", second_step_acceptor_sites, &query_list);

        // this is the third step, when A1 is found, check D2

        let d2 = query_list.get("D2").unwrap().clone();

        // this is parellel step to step 3, when A2 is found, check D3, followed by D3-unspliced and A3 to A7

        // this is the fourth step, when D2 is found, check D2-unspliced and A2 to A7

        let mut fouth_step_acceptor_sites = vec!["D2-unspliced", "A2"];

        fouth_step_acceptor_sites.append(&mut common_acceptor_sites);

        let fourth_step = SpliceStep::build("D2", fouth_step_acceptor_sites, &query_list);

        // this is parellel step to step 4, when D2 is absent, check D2b

        let d2b = query_list.get("D2b").unwrap().clone();

        // this is step 5, only when D2-unspliced is found, check D2b-unspliced

        let d2b_unspliced = query_list.get("D2b-unspliced").unwrap().clone();

        // this is parellel step to step 5, when D2b is found is absent, check D2b-unspliced and A2 to A7

        let mut fifth_step_acceptor_sites = vec!["D2b-unspliced", "A2"];

        fifth_step_acceptor_sites.append(&mut common_acceptor_sites);

        let fifth_step = SpliceStep::build("D2b", fifth_step_acceptor_sites, &query_list);

        // this is step 6, check D3, same as the parellel step to step 6

        let d3 = query_list.get("D3").unwrap().clone();

        // this is step 7, check D3-unspliced and A3 to A7, it will be used in several steps, will be called d3_unspliced_a3_a7
        // same as the parellel step to step 7

        let mut d3_unspliced_a3_a7 = vec!["D3-unspliced"];

        d3_unspliced_a3_a7.append(&mut common_acceptor_sites);

        let seventh_step = SpliceStep::build("D3", d3_unspliced_a3_a7, &query_list);

        Ok(SpliceConfig {
            strain,
            distance,
            query_list,
            d1,
            d2,
            d2b,
            d2b_unspliced,
            d3,
            d1_to_all: second_step,
            d2_to_all: fourth_step,
            d2b_to_all: fifth_step,
            d3_to_all: seventh_step,
        })
    }
}



impl SpliceStep {
    pub fn build(
        donor: &str,
        acceptor_list: Vec<&str>,
        query_list: &HashMap<String, String>,
    ) -> Vec<SpliceStep> {
        let mut splice_steps: Vec<SpliceStep> = Vec::new();
        let donor_pattern = query_list.get(donor).unwrap().to_string();
        for acceptor in acceptor_list {
            let pattern = donor_pattern.clone() + query_list.get(acceptor).unwrap();
            splice_steps.push(SpliceStep {
                donor: donor.to_string(),
                acceptor: acceptor.to_string(),
                pattern,
            });
        }
        splice_steps
    }
}

#[cfg(test)]

mod tests {
    use super::*;
    #[test]
    fn test_splice_config() {
        let strain = "NL43".to_lowercase().to_string();
        let distance = 10;
        let splice_config = SpliceConfig::build(strain, distance).unwrap();
        assert_eq!(splice_config.strain, "nl43");
        assert_eq!(splice_config.distance, 10);
        assert_eq!(splice_config.query_list.get("D1").unwrap(), "GGGGCGGCGACTG");
        assert_eq!(splice_config.query_list.get("D1-unspliced").unwrap(), "GTGAGTACGCCAAAAA");
        assert_eq!(splice_config.d1_to_all[0].pattern, "GGGGCGGCGACTGGTGAGTACGCCAAAAA");
        assert_eq!(splice_config.d1_to_all[0].donor, "D1");
        assert_eq!(splice_config.d1_to_all[0].acceptor, "D1-unspliced");
        assert_eq!(splice_config.d1_to_all[1].pattern, "GGGGCGGCGACTGGGACAGCAGAGA");
        assert_eq!(splice_config.d1_to_all[1].donor, "D1");
        assert_eq!(splice_config.d1_to_all[1].acceptor, "A1");
        assert_eq!(splice_config.d1_to_all[2].pattern, "GGGGCGGCGACTGAATCTGCTATAA");
        assert_eq!(splice_config.d1_to_all[2].donor, "D1");
        assert_eq!(splice_config.d1_to_all[2].acceptor, "A2");   
        assert_eq!(splice_config.d2_to_all[0].pattern, "AG");
        assert_eq!(splice_config.d2_to_all[0].donor, "D2");
        assert_eq!(splice_config.d2_to_all[0].acceptor, "D2-unspliced");
        assert_eq!(splice_config.d2_to_all[1].pattern, "AAATCTGCTATAA");
        assert_eq!(splice_config.d2_to_all[1].donor, "D2");
        assert_eq!(splice_config.d2_to_all[1].acceptor, "A2");
    }

    #[test]

    fn test_build_splice_step() {
        let list = vec!["A1", "A2", "A3"];
        let mut query_list = HashMap::new();
        query_list.insert("D1".to_string(), "ACGT".to_string());
        query_list.insert("A1".to_string(), "TGCA".to_string());
        query_list.insert("A2".to_string(), "ACGT".to_string());
        query_list.insert("A3".to_string(), "CCCC".to_string());
        let splice_steps = SpliceStep::build("D1", list, &query_list);
        assert_eq!(splice_steps[0].pattern, "ACGTTGCA".to_string());
        assert_eq!(splice_steps[0].donor, "D1".to_string());
        assert_eq!(splice_steps[0].acceptor, "A1".to_string());
        assert_eq!(splice_steps[1].pattern, "ACGTACGT".to_string());
        assert_eq!(splice_steps[1].donor, "D1".to_string());
        assert_eq!(splice_steps[1].acceptor, "A2".to_string());
        assert_eq!(splice_steps[2].pattern, "ACGTCCCC".to_string());
        assert_eq!(splice_steps[2].donor, "D1".to_string());
        assert_eq!(splice_steps[2].acceptor, "A3".to_string());
    }


}
