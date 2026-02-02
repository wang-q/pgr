use std::fs::File;
use std::io::Read;
use std::path::Path;

pub fn is_fq<P: AsRef<Path>>(path: P) -> bool {
    let path = path.as_ref();

    // Create a buffer to store the first two bytes
    let mut buffer = [0; 2];
    {
        let mut file = match File::open(path) {
            Err(why) => panic!("could not open {}: {}", path.display(), why),
            Ok(file) => file,
        };
        file.read_exact(&mut buffer).unwrap();
    }

    // Check if the file is in Gzip format
    let is_fq;
    if buffer[0] == 0x1f && buffer[1] == 0x8b {
        let mut decoder = flate2::read::GzDecoder::new(File::open(path).unwrap());
        let mut buffer = [0; 2]; // Recreate the buffer
        decoder.read_exact(&mut buffer).unwrap();

        // Determine the format of the decompressed file
        match buffer[0] as char {
            '>' => is_fq = false,
            '@' => is_fq = true,
            _ => unreachable!("Unknown file format"),
        }
    } else {
        // The file is in plain text format, determine the format
        match buffer[0] as char {
            '>' => is_fq = false,
            '@' => is_fq = true,
            _ => unreachable!("Unknown file format"),
        }
    }

    is_fq
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_is_fq_plain_text() {
        let dir = tempdir().unwrap();

        // Create a plain text FASTQ file
        let fq_file_path = dir.path().join("test.fq");
        {
            let mut file = File::create(&fq_file_path).unwrap();
            writeln!(file, "@SEQ_ID").unwrap(); // FASTQ format
        }
        assert!(is_fq(&fq_file_path));

        // Create a plain text FASTA file
        let fasta_file_path = dir.path().join("test.fasta");
        {
            let mut file = File::create(&fasta_file_path).unwrap();
            writeln!(file, ">SEQ_ID").unwrap(); // FASTA format
        }
        assert!(!is_fq(&fasta_file_path));
    }

    #[test]
    fn test_is_fq_gzip() {
        let dir = tempdir().unwrap();

        // Create a Gzip FASTQ file
        let fq_file_path = dir.path().join("test.fq.gz");
        {
            let file = File::create(&fq_file_path).unwrap();
            let mut encoder = GzEncoder::new(file, flate2::Compression::default());
            writeln!(encoder, "@SEQ_ID").unwrap(); // FASTQ format
            encoder.finish().unwrap();
        }
        assert!(is_fq(&fq_file_path));

        // Create a Gzip FASTA file
        let fasta_file_path = dir.path().join("test.fasta.gz");
        {
            let file = File::create(&fasta_file_path).unwrap();
            let mut encoder = GzEncoder::new(file, flate2::Compression::default());
            writeln!(encoder, ">SEQ_ID").unwrap(); // FASTA format
            encoder.finish().unwrap();
        }
        assert!(!is_fq(&fasta_file_path));
    }
}
