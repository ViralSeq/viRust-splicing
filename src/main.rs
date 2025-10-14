use colored::Colorize;
use std::process;

use virust_splicing::config::InputConfig;

fn main() {
    println!(
        "{} {}",
        "ViRust-Splicing for processing HIV splicing data.\n By Shuntai Zhou and Michael Clark @ Zhou lab, 2025\n Version:"
            .cyan()
            .bold(),
        env!("CARGO_PKG_VERSION").cyan().bold()
    );
    let config = InputConfig::build().unwrap_or_else(|err| {
        eprintln!(
            "Problem parsing arguments: {}",
            err.to_string().red().bold()
        );
        process::exit(1);
    });
    println!("✅ Configurations validated, starting analysis...");
    #[cfg(debug_assertions)]
    dbg!(&config);
    if let Err(e) = virust_splicing::run(config) {
        eprintln!("Application error: {}", e.to_string().red().bold());
        process::exit(1);
    }
}
