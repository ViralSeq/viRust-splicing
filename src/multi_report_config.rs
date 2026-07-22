use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use anyhow::Result;
use clap::Parser;
use csv::Writer;
use getset::{Getters, Setters};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Parser)]
#[command(
    name = "viRust-splicing multi-report",
    about = "Pull and analyze multiple splicing reports at once.\n
    IMPORTANT: Splicing reports from different assays should not be analyzed together.\n
    IMPORTANT: Pooling reports from more than 10 samples may lead to visualization issues in the HTML report.\n
    Please consider analyzing in smaller batches.", 
    version = env!("CARGO_PKG_VERSION")
)]
pub struct MultiReportConfig {
    /// Input directory containing multiple splicing output subdirectories from viRust-splicing.
    /// Each subdirectory should contain the output files for a single sample.
    /// The program will mainly check the *_summary.csv file in each subdirectory.
    /// The name of the subdirectories will be used to label the samples in the combined report.
    /// Example structure:
    /// input_dir/
    /// ├── sample1/
    /// │   ├── sample1_summary.csv
    /// │   └── other_output_files...
    /// ├── sample2/
    /// │   ├── sample2_summary.csv
    /// │   └── other_output_files...
    /// └── sampleN/
    ///     ├── sampleN_summary.csv
    ///     └── other_output_files...
    #[arg(short, long, required = true)]
    pub input_dir: String,
    /// Output directory to save the combined report and analysis results.
    /// The directory will be created if it does not exist.
    #[arg(short, long, required = true)]
    pub output_dir: String,
}

#[derive(Debug, Serialize, Deserialize, Clone, Getters, Setters)]
pub struct ValidatedMultiReportConfig {
    #[getset(get = "pub")]
    input_dir: String,
    #[getset(get = "pub")]
    output_dir: String,
    #[getset(get = "pub")]
    summary_file_map: Vec<SummaryFileFeatures>, // mapping from subdirectory name to summary CSV file path
}

#[derive(Debug, Serialize, Deserialize, Clone, Getters, Setters)]
pub struct SummaryFileFeatures {
    #[getset(get = "pub")]
    sample_name: String,
    #[getset(get = "pub")]
    reference_strain: String,
    #[getset(get = "pub")]
    distance: String,
    #[getset(get = "pub")]
    assay: String,
    #[getset(get = "pub")]
    file_path: PathBuf,
}

impl MultiReportConfig {
    /// Validate the input directory and output directory.
    /// Check if the input directory exists and is a directory.
    /// Check if the input directory contains at least one subdirectory with a summary CSV file.
    /// Create the output directory if it does not exist.
    /// Return a HashMap mapping subdirectory names to their summary CSV file paths.
    /// Return an error if any validation fails.
    pub fn validate(&self) -> Result<ValidatedMultiReportConfig> {
        let input_path = Path::new(&self.input_dir);
        if !input_path.exists() || !input_path.is_dir() {
            return Err(anyhow::anyhow!(
                "Input directory '{}' does not exist or is not a directory.",
                self.input_dir
            ));
        }
        let mut summary_file_map: Vec<SummaryFileFeatures> = Vec::new();
        let mut assays = Vec::new();
        for entry in std::fs::read_dir(input_path)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                // use regex to find the summary csv file
                // the summary file should have the format as: output_strain_<strain>_distance_<distance>_assay_<assay>_summary.csv
                // where <strain>, <distance>, and <assay> can be any string without spaces

                let re = regex::Regex::new(
                    r"^output_strain_(.*)_distance_(.*)_assay_(.*)_summary\.csv$",
                )
                .unwrap();

                for sub_entry in std::fs::read_dir(path)? {
                    let sub_entry = sub_entry?;
                    let sub_path = sub_entry.path();
                    if sub_path.is_file() {
                        if let Some(file_name) = sub_path.file_name() {
                            if let Some(file_name_str) = file_name.to_str() {
                                if re.is_match(file_name_str) {
                                    // capture the features from the file name
                                    if let Some(caps) = re.captures(file_name_str) {
                                        if let (
                                            Some(reference_strain),
                                            Some(distance),
                                            Some(assay),
                                        ) = (caps.get(1), caps.get(2), caps.get(3))
                                        {
                                            summary_file_map.push(SummaryFileFeatures {
                                                sample_name: entry
                                                    .file_name()
                                                    .to_string_lossy()
                                                    .to_string(),
                                                reference_strain: reference_strain
                                                    .as_str()
                                                    .to_string(),
                                                distance: distance.as_str().to_string(),
                                                assay: assay.as_str().to_string(),
                                                file_path: sub_path.clone(),
                                            });
                                            assays.push(assay.as_str().to_string());
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        let has_valid_subdir = !summary_file_map.is_empty();
        if !has_valid_subdir {
            return Err(anyhow::anyhow!(
                "Input directory '{}' does not contain any subdirectories with a summary CSV file.",
                self.input_dir
            ));
        }
        let output_path = Path::new(&self.output_dir);

        if !output_path.exists() {
            std::fs::create_dir_all(output_path)?;
        } else if !output_path.is_dir() {
            return Err(anyhow::anyhow!(
                "Output path '{}' exists but it is a file.",
                self.output_dir
            ));
        }
        assays.sort();
        assays.dedup();
        if assays.len() > 1 {
            println!(
                "Warning: Multiple assay types found in the input directory: {:?}. \
                It is recommended to analyze reports from the same assay type together.",
                assays
            );
        }
        Ok(ValidatedMultiReportConfig {
            input_dir: self.input_dir.clone(),
            output_dir: self.output_dir.clone(),
            summary_file_map,
        })
    }
}

const MULTI_QMD: &str = include_str!("../scripts/quarto/multi_report.qmd");

impl ValidatedMultiReportConfig {
    pub fn write_feature_summary_csv(&self) -> Result<()> {
        let feature_summary_path = Path::new(&self.output_dir).join("run_features.csv");
        let mut wtr = Writer::from_path(feature_summary_path)?;
        // Implementation to write the feature summary CSV goes here
        for record in self.summary_file_map() {
            wtr.serialize(record)?;
        }
        wtr.flush()?;

        // copy all the *_summary.csv files into the output directory with a new name
        for record in self.summary_file_map() {
            let new_file_name = format!("{}_summary.csv", record.sample_name());
            let new_file_path = Path::new(&self.output_dir).join(new_file_name);
            std::fs::copy(record.file_path(), new_file_path)?;
        }
        Ok(())
    }

    pub fn run_multi_report(&self) -> Result<()> {
        let start = Instant::now();
        // Write the multi_report.qmd file to the output directory
        let qmd_path = Path::new(&self.output_dir).join("multi_report.qmd");
        std::fs::write(&qmd_path, MULTI_QMD)?;

        // write spliceforms_short.csv to the output directory
        let splice_isoform_csv = Path::new(&self.output_dir).join("spliceforms_short.csv");
        crate::runner::write_reference_csv(&splice_isoform_csv)?;

        // Run Quarto to generate the multi-sample report, no args needed.
        let output = Command::new("quarto")
            .arg("render")
            .arg(qmd_path)
            .output()?;

        let elapsed = start.elapsed();
        println!("✅ Report generated in {:.2?}", elapsed);

        if !output.status.success() {
            eprintln!("⚠️ Quarto render failed.");
            eprintln!("--- stdout ---");
            eprintln!("{}", String::from_utf8_lossy(&output.stdout));
            eprintln!("--- stderr ---");
            eprintln!("{}", String::from_utf8_lossy(&output.stderr));
            Err(anyhow::anyhow!("Quarto render command failed"))
        } else {
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_input_config() {
        let config = MultiReportConfig {
            input_dir: "tests/multi_report/test_input".to_string(),
            output_dir: "tests/multi_report/test_output".to_string(),
        };
        let result = config.validate();
        dbg!(&result);
        assert!(result.is_ok());
        let summary_map = result.unwrap();
        dbg!(summary_map.clone());
        assert_eq!(summary_map.summary_file_map.len(), 3);
        assert!(
            summary_map
                .summary_file_map
                .iter()
                .any(|s| s.sample_name == "d3" && s.assay == "random-reverse")
        );
    }

    #[test]
    fn test_write_feature_summary_csv() {
        let config = MultiReportConfig {
            input_dir: "tests/multi_report/test_input".to_string(),
            output_dir: "tests/multi_report/test_output".to_string(),
        };
        let validated_config = config.validate().unwrap();
        let result = validated_config.write_feature_summary_csv();
        assert!(result.is_ok());
        let feature_summary_path = Path::new(&config.output_dir).join("run_features.csv");
        // read and print out csv content for debugging
        let mut rdr = csv::Reader::from_path(&feature_summary_path).unwrap();
        rdr.deserialize::<SummaryFileFeatures>().for_each(|result| {
            let record = result.unwrap();
            println!("{:?}", record);
        });
        assert!(feature_summary_path.exists());
        // Clean up (optional)
        // std::fs::remove_file(feature_summary_path).unwrap();
        // std::fs::remove_dir_all(&config.output_dir).unwrap();
    }

    #[test]
    fn test_run_multi_report() {
        let config = MultiReportConfig {
            input_dir: "tests/multi_report/test_input".to_string(),
            output_dir: "tests/multi_report/test_output".to_string(),
        };
        let validated_config = config.validate().unwrap();
        validated_config.write_feature_summary_csv().unwrap(); // ensure the summary CSVs are in place
        let result = validated_config.run_multi_report();
        assert!(result.is_ok());
        let report_html_path = Path::new(&config.output_dir).join("multi_report.html");
        assert!(report_html_path.exists());
    }
}
