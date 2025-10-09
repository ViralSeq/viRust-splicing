use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};

use bio::io::fasta;
use bio::io::fastq;
use flate2::read::MultiGzDecoder;

#[derive(Debug, PartialEq)]
pub enum DataType {
    Fasta,
    Fastq,
    FastaGz,
    FastqGz,
}

pub fn determine_data_type(file_path: &Path) -> Result<DataType, Box<dyn Error>> {
    let file_name = file_path
        .file_name()
        .ok_or("Invalid file path: no file name found")?
        .to_string_lossy()
        .to_string();

    if file_name.ends_with(".fasta") || file_name.ends_with(".fa") {
        Ok(DataType::Fasta)
    } else if file_name.ends_with(".fastq") || file_name.ends_with(".fq") {
        Ok(DataType::Fastq)
    } else if file_name.ends_with(".fasta.gz") || file_name.ends_with(".fa.gz") {
        Ok(DataType::FastaGz)
    } else if file_name.ends_with(".fastq.gz") || file_name.ends_with(".fq.gz") {
        Ok(DataType::FastqGz)
    } else {
        Err(format!("Unsupported file format of input file {}", file_name).into())
    }
}

#[derive(Debug)]
pub struct SequenceFiles {
    pub r1_file: PathBuf,
    pub r2_file: PathBuf,
    pub data_type: DataType,
}

// pub fn read_sequence_files()
pub fn validate_input_files(r1_path: &str, r2_path: &str) -> Result<SequenceFiles, Box<dyn Error>> {
    let r1_file = PathBuf::from(r1_path);
    let r2_file = PathBuf::from(r2_path);

    let data_type_r1 = determine_data_type(&r1_file)?;
    let data_type_r2 = determine_data_type(&r2_file)?;

    if data_type_r1 != data_type_r2 {
        return Err(format!(
            "Input files are of different formats, r1 is {:?}, r2 is {:?}. Program expects both files to be of the same format.",
            data_type_r1, data_type_r2
        )
        .into());
    }

    Ok(SequenceFiles {
        r1_file,
        r2_file,
        data_type: data_type_r1,
    })
}

pub fn open_sequence_files(
    files: &SequenceFiles,
) -> std::io::Result<Vec<(fasta::Record, fasta::Record)>> {
    let r1_file = File::open(&files.r1_file)?;
    let r2_file = File::open(&files.r2_file)?;

    let (r1_stream, r2_stream): (Box<dyn std::io::Read>, Box<dyn std::io::Read>) =
        match files.data_type {
            DataType::Fastq | DataType::Fasta => (
                Box::new(BufReader::new(r1_file)),
                Box::new(BufReader::new(r2_file)),
            ),
            DataType::FastqGz | DataType::FastaGz => (
                Box::new(MultiGzDecoder::new(BufReader::new(r1_file))),
                Box::new(MultiGzDecoder::new(BufReader::new(r2_file))),
            ),
        };
    if files.data_type == DataType::Fastq || files.data_type == DataType::FastqGz {
        let r1_reader = fastq::Reader::new(r1_stream);
        let r2_reader = fastq::Reader::new(r2_stream);

        let pairs: Vec<(fasta::Record, fasta::Record)> = r1_reader
            .records()
            .zip(r2_reader.records())
            .map(|(r1_res, r2_res)| {
                let r1 = r1_res.map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
                let r2 = r2_res.map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
                Ok::<_, io::Error>((fastq_to_fasta(&r1), fastq_to_fasta(&r2)))
            })
            .collect::<Result<Vec<_>, io::Error>>()?;
        Ok(pairs)
    } else {
        let r1_reader = fasta::Reader::new(r1_stream);
        let r2_reader = fasta::Reader::new(r2_stream);

        let pairs: Vec<(fasta::Record, fasta::Record)> = r1_reader
            .records()
            .zip(r2_reader.records())
            .map(|(r1_res, r2_res)| {
                let r1 = r1_res?;
                let r2 = r2_res?;
                Ok::<_, io::Error>((r1, r2))
            })
            .collect::<Result<Vec<_>, io::Error>>()?;
        Ok(pairs)
    }
}

fn fastq_to_fasta(record: &fastq::Record) -> fasta::Record {
    fasta::Record::with_attrs(record.id(), record.desc(), record.seq())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_determine_data_type() {
        let fasta_path = Path::new("sample.fasta");
        let fastq_path = Path::new("sample.fastq");
        let fasta_gz_path = Path::new("sample.fasta.gz");
        let fastq_gz_path = Path::new("sample.fastq.gz");
        let invalid_path = Path::new("sample.txt");

        assert_eq!(determine_data_type(fasta_path).unwrap(), DataType::Fasta);
        assert_eq!(determine_data_type(fastq_path).unwrap(), DataType::Fastq);
        assert_eq!(
            determine_data_type(fasta_gz_path).unwrap(),
            DataType::FastaGz
        );
        assert_eq!(
            determine_data_type(fastq_gz_path).unwrap(),
            DataType::FastqGz
        );
        assert!(determine_data_type(invalid_path).is_err());
    }

    #[test]
    fn test_validate_input_files() {
        let r1_path = "sample_R1.fastq";
        let r2_path = "sample_R2.fastq";
        let result = validate_input_files(r1_path, r2_path).unwrap();
        assert_eq!(result.r1_file, PathBuf::from(r1_path));
        assert_eq!(result.r2_file, PathBuf::from(r2_path));
        assert_eq!(result.data_type, DataType::Fastq);

        let r1_path_mismatch = "sample_R1.fastq";
        let r2_path_mismatch = "sample_R2.fasta";
        let result_mismatch = validate_input_files(r1_path_mismatch, r2_path_mismatch);
        assert!(result_mismatch.is_err());
    }
}
