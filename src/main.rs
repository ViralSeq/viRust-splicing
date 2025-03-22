use std::process;
use colored::Colorize;

use virust_splicing::config::InputConfig;

fn main() {
    println!("{}", "ViRust-Splicing for processing HIV splicing data.\n Placeholder for more information.".cyan().bold());
    let config = InputConfig::build()
        .unwrap_or_else(|err| {
            eprintln!("Problem parsing arguments: {}", err.to_string().red().bold());
            process::exit(1);
    });


    if let Err(e) = virust_splicing::run(config) {
        eprintln!("Application error: {}", e.to_string().red().bold());
        process::exit(1);
    }

}
