# Changelog

All notable changes to this project will be documented in this file.

## [0.3.0] - 2026-07-21

- Updated the project for newer Rust/rand compatibility and resolved stale rust-analyzer diagnostics.
- Added internal-only homopolymer filtering, keeping terminal homopolymers allowed.
- Added optional IUPAC-aware --headerMatch filtering on the beginning of R1.
- Expanded UMI family handling:
  - max-model
  - cluster-model
  - cluster-model-plusone
- auto-model
- Reworked clustering to use direct, abundance-ordered parent assignments, preventing transitive UMI chains.
- Added --umiMinFraction, UMI family size/fraction TSV fields, and Quarto parameter reporting.
- Added auto selection per final_category: max model for well-covered categories or collision-prone UMI pools; cluster model only for sparse, low-risk categories.
- Added --umiAutoNeighborRisk, default 0.05, plus TSV diagnostics for actual model, unique UMI count/fraction, and lambda.
- Expanded the README with the method background, formula, assumptions, and the Primer ID validation paper.-

## [0.2.1] - 2025-10-17

- A helper bin `multi_report` added to perform comparisons across samples

- Fixed a bug in counting D1 vs. D1prime as the first splicing donor

- Fixed a bug to skip size class accessment for sequences using alternative D1

## [0.2.0] - 2025-10-15

- Initial stable build
