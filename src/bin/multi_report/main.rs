//TODO: Implement this binary for pull and analyze multiple splicing reports at once.

use anyhow::Result;
use clap::Parser;

use virust_splicing::multi_report_config::MultiReportConfig;

fn main() -> Result<()> {
    let config = MultiReportConfig::parse();

    let validated_config = config.validate()?;

    validated_config.write_feature_summary_csv()?;

    validated_config.run_multi_report()?;

    Ok(())
}
