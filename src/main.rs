use colored::Colorize;
use std::process;

use virust_splicing::config::InputConfig;

// TODO: Intergrated test needed.
fn main() {
    //TODO: Add a help message to the config struct.
    println!(
        "{}",
        "ViRust-Splicing for processing HIV splicing data.\n Placeholder for more information."
            .cyan()
            .bold()
    );
    let config = InputConfig::build().unwrap_or_else(|err| {
        eprintln!(
            "Problem parsing arguments: {}",
            err.to_string().red().bold()
        );
        process::exit(1);
    });

    if let Err(e) = virust_splicing::run(config) {
        eprintln!("Application error: {}", e.to_string().red().bold());
        process::exit(1);
    }
}
