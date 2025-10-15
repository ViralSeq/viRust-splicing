# HIV Splicing analysis

![License](https://img.shields.io/badge/license-MIT-green)
![Rust](https://img.shields.io/badge/language-Rust-orange)
![GitHub last commit](https://img.shields.io/github/last-commit/ViralSeq/viRust-splicing)
![GitHub repo size](https://img.shields.io/github/repo-size/ViralSeq/viRust-splicing)
![GitHub Repo stars](https://img.shields.io/github/stars/ViralSeq/viRust-splicing?style=social)

**viRust-splicing** is a high-performance command-line program written in Rust for analyzing HIV-1 splicing patterns from UMI/Primer-ID–based next-generation sequencing (NGS) data.

It parses paired-end FASTQ/FASTA files (also support .gz files), maps donor–acceptor sites, quantifies splice isoforms, and generates structured outputs suitable for downstream statistical analysis and visualization.Rust App to analyze NGS data to detect HIV splice forms.

## Key Features

- ⚡ Ultra-fast processing using Rust parallelism (Rayon)

- 🧬 UMI/Primer-ID aware read deduplication and merging

- 🧭 Donor/acceptor mapping to HXB2 reference coordinates

- 📊 Isoform quantification and junction usage statistics

- 🪄 TSV/CSV/JSON output for downstream R/Python/Quarto visualization

- 🧰 Seamless integration with Plotly, ggplot2, and D3.js

- 🌐 Self-contained HTML report with downloadable data.

## Installation

1. From GitHub (lastest dev)

```bash
git clone https://github.com/ViralSeq/viRust-splicing.git
cd viRusts-splicing
cargo build --release
./target/release/virust-splicing
```

### Dependencies:

- Rust ≥ 1.7.5

- R(≥ 4.0.0) Python3 and Quarto are required for report generation.

## Usage

```
Usage: virust-splicing [OPTIONS] --query <QUERY> --distance <DISTANCE> --file1 <FILENAME_R1> --file2 <FILENAME_R2> --assay <ASSAY_TYPE>

Options:
  -q, --query <QUERY>
  -d, --distance <DISTANCE>
  -1, --file1 <FILENAME_R1>
  -2, --file2 <FILENAME_R2>
  -a, --assay <ASSAY_TYPE>    [possible values: random-reverse, k-mer, size-specific]
  -o, --output <OUTPUT_PATH>
  -h, --help                  Print help
  -V, --version               Print version

```

## Example usage

### Run from repo

```bash
cargo run --release -- --query nl43 --distance 1 --assay random-reverse --file1 ./sample_data/mock_data_r1.fasta.gz --file2 ./sample_data/mock_data_r2.fasta.gz -o ./
```

### Run using binary

```bash
./virust-splicing -q nl43 -d 1 -a random-reverse -1 ./sample_data/mock_data_r1.fasta.gz -2 ./sample_data/mock_data_r2.fasta.gz -o ./
```

<small>Note that inputs use '-' instead of '\_' , use assay "random-reverse" instead of "random_reverse"</small>

### Example Output

```
results/
 ├─ output_strain_nl43_distance_1_assay_random-reverse.tsv
 ├─ output_strain_nl43_distance_1_assay_random-reverse_summary.csv
 ├─ report.html
 └─ report.qmd
```

## Benchmark & Performance

    •	Processes 10M+ reads in minutes using multi-threading
    •	Memory efficient (<2 GB typical usage)
    •	Scales linearly with read count and cores

## Biological Context

HIV-1 relies on extensive alternative splicing to generate over 40 mRNAs. This tool identifies splice junction usage from NGS data, enabling:

- Quantification of donor–acceptor site usage (D1→A1–A7, D2→A3–A7, etc.)

- Detection of cryptic junctions

- Longitudinal splicing dynamics analysis

## Contributing

We welcome pull requests and bug reports!

    1.	Fork the repo
    2.	Create a feature branch (git checkout -b feature/foo)
    3.	Commit and push
    4.	Open a PR

## License

This project is licensed under the MIT License.

## Achnowledgements

- NIH/NIAID for funding support.

- Center for Structural Biology of HIV RNA (CRNA) (U54 AI170660-01)

- UNC Center for AIDS Research (CFAR) (P30 AI050410)

- Dr. Ann Emery

## Splice Event Idenfitication Workflow

```mermaid

%%{
  init: {
    'theme': 'base',
    'themeVariables': {
      'primaryColor': '#1966C9',
      'primaryTextColor': '#FEFFE2',
      'primaryBorderColor': '#3A8CE5',
      'lineColor': '#3A8CE5',
      'secondaryColor': '#F87C07',
      'secondaryBorderColor': '#FFC78F',
      'tertiaryColor': '#FF8BA0'
    }
  }
}%%


flowchart TD
    A(Start) --> B(Check for D1)
    B -- Found --> C(Add D1)
    B -- Not Found --> End1(No D1 and End)

    C --> D(Check acceptor after D1)
    D -- A1 --> E(Add A1 event)
    D -- A2 --> F(Add A2 event)
    D -- Other --> End2(Found A3,A4a,A4b,A4c,A4d,A5a,A5b,A7; End)
    D -- Not Found --> End3(Add unknown and End)

    E --> G(Check for D2)
    G -- Found --> H(Add D2)
    G -- Not Found --> End4(No D2 and End)

    H --> I(Check acceptor after D2)
    I -- D2-unspliced --> J(Add D2-unspliced event)
    I -- A2 --> K(Add A2 event)
    I -- Other --> End5(Add unknown and End)

    J --> L(Check for D2b)
    L -- Found --> M(Add D2b)
    L -- Not Found --> End6(No D2b and End)

    M --> N(Check acceptor after D2b)
    N -- A2 --> O(Add A2 event - Process D3)
    N -- Other --> End7(Add unknown and End)

    %% Convergence to the D3 branch:
    K --> P(Check for D3)
    F --> P
    O --> P

    P -- Found --> Q(Add D3)
    P -- Not Found --> End8(No D3 and End)

    Q --> R(Check acceptor after D3 )
    R -- Found --> End9(Add acceptor event: A3,A4a, A4b, A4c, A4d, A5, A7; and End)
    R -- Not Found --> End10(Add unknown and End)


```
