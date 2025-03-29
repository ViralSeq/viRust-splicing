use std::{collections::HashMap, env::args};
use csv;



fn main() {

    let number_to_generate = args().nth(1).unwrap().parse::<usize>().unwrap();
    
    // read csv file from /data/splice_forms_mock_ratio.csv
    let mut rdr = csv::Reader::from_path("data/splice_forms_mock_ratio.csv").unwrap();
    let mut ratio_hash: HashMap<String, f32> = HashMap::new();
    for result in rdr.records() {
        let record = result.unwrap();
        let key = record.get(0).unwrap().to_string();
        let value = record.get(1).unwrap().parse::<f32>().unwrap();
        ratio_hash.insert(key, value);
    }
    println!("{:#?}", ratio_hash);

    println!("{}", ratio_hash.get("noD1").unwrap());

    println!("Number to generate: {}", number_to_generate);
}
