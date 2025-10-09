use std::error::Error;
use std::process::Command;
use std::time::Instant;

pub fn check_r_installed() -> Result<(), Box<dyn Error>> {
    // In a real implementation, you would check if R is available in the system PATH
    let output = Command::new("R").arg("--version").output();

    match output {
        Ok(out) => {
            if out.status.success() {
                Ok(())
            } else {
                Err("R command failed to execute properly.".into())
            }
        }
        Err(_) => Err("R is not installed or not found in PATH.".into()),
    }
}

pub fn check_quarto_installed() -> Result<(), Box<dyn Error>> {
    // In a real implementation, you would check if Quarto is available in the system PATH
    let output = Command::new("quarto").arg("--version").output();

    match output {
        Ok(out) => {
            if out.status.success() {
                Ok(())
            } else {
                Err("Quarto command failed to execute properly.".into())
            }
        }
        Err(_) => Err("Quarto is not installed or not found in PATH.".into()),
    }
}

pub fn check_r_packages() -> Result<(), Box<dyn Error>> {
    let check_env_script = include_str!("../scripts/R/check_env.R");

    let output = Command::new("Rscript")
        .arg("-e")
        .arg(check_env_script)
        .output()?;

    if output.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        Err(format!("R packages check failed: {}", stderr).into())
    }
}

pub fn r_summarize_data(tsv_path: &str, output_path: &str) -> Result<(), Box<dyn Error>> {
    let summarize_script = include_str!("../scripts/R/summarize_data.R");

    let output = Command::new("Rscript")
        .arg("-e")
        .arg(summarize_script)
        .arg(tsv_path)
        .arg(output_path)
        .output()?;

    if output.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        Err(format!("R summarize data failed: {}", stderr).into())
    }
}

pub fn run_quarto_report(
    input: &str,
    spliceform_file: &str,
    output_html: &str,
) -> std::io::Result<()> {
    let qmd = include_str!("../scripts/quarto/report.qmd");
    let version = env!("CARGO_PKG_VERSION"); // compile-time version
    let start = Instant::now();

    // Run Quarto
    let status = Command::new("quarto")
        .arg("render")
        .arg(qmd)
        .arg("--execute")
        .arg("--to")
        .arg("html")
        .arg("--output")
        .arg(output_html)
        .arg("--")
        .arg(input)
        .arg(version)
        .arg(spliceform_file)
        .status()?;

    let elapsed = start.elapsed();
    println!("✅ Report generated in {:.2?}", elapsed);

    if !status.success() {
        eprintln!("⚠️ Quarto render failed");
    }

    Ok(())
}

pub fn write_reference_csv(spliceform_file: &str) -> std::io::Result<()> {
    let reference_csv = include_str!("../data/spliceforms_short.csv");
    std::fs::write(spliceform_file, reference_csv)
}
