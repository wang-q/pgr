use clap::*;
use noodles_bgzf as bgzf;
use std::io::{Read, Seek};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("gz")
        .about("Compressing a file using the BGZF format")
        .after_help(
            r###"
This command compresses a file using BGZF (Blocked Gzip Format).

Features:
* Parallel compression with multiple threads
* Creates index file (.gzi) for random access
* Supports stdin as input
* Preserves original file

Output files:
* <infile>.gz: Compressed file
* <infile>.gz.gzi: Index file

Notes:
* Cannot compress already gzipped files
* Default thread count is 1
* Index creation is automatic

Examples:
1. Compress a file with default settings, and the outfile is input.fa.gz:
   pgr fa gz input.fa

2. Multi-threaded compression:
   pgr fa gz input.fa -p 4

3. Set compression level (0-9, default -1):
   pgr fa gz input.fa -l 9

4. Create index for existing .gz file (reindex):
   pgr fa gz input.fa.gz -r

5. From stdin with custom output:
   cat input.fa | pgr fa gz stdin -o output.fa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input FA file to compress"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .value_parser(value_parser!(std::num::NonZeroUsize))
                .num_args(1)
                .default_value("1")
                .help("Number of threads for parallel compression"),
        )
        .arg(
            Arg::new("compress-level")
                .long("compress-level")
                .short('l')
                .value_parser(value_parser!(i32))
                .num_args(1)
                .default_value("-1")
                .help("Compression level (0-9, or -1 for default)"),
        )
        .arg(
            Arg::new("reindex")
                .long("reindex")
                .short('r')
                .action(ArgAction::SetTrue)
                .help("Create BGZF index (.gzi) for an existing .gz file"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .help("Output filename (default: <infile>.gz)"),
        )
}

// command implementation
fn is_bgzf(path: &str) -> bool {
    if let Ok(mut file) = std::fs::File::open(path) {
        let mut header = [0u8; 4];
        if file.read_exact(&mut header).is_ok() {
            // Check GZIP magic (0x1f 0x8b) and FEXTRA flag (4) in FLG (byte 3)
            return header[0] == 0x1f && header[1] == 0x8b && (header[3] & 4) != 0;
        }
    }
    false
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();

    if args.get_flag("reindex") {
        if !std::path::Path::new(infile).exists() {
            return Err(anyhow::anyhow!("Input file not found: {}", infile));
        }
        if !is_bgzf(infile) {
            return Err(anyhow::anyhow!("Input file is not a valid BGZF file: {}", infile));
        }
        build_gzi_index(infile)?;
        return Ok(());
    }

    let opt_parallel = *args.get_one::<std::num::NonZeroUsize>("parallel").unwrap();
    let compress_level = *args.get_one::<i32>("compress-level").unwrap();

    let outfile = if args.contains_id("outfile") {
        args.get_one::<String>("outfile").unwrap().to_string()
    } else {
        format!("{}.gz", infile)
    };

    //----------------------------
    // Input
    //----------------------------
    let mut reader: Box<dyn std::io::BufRead> = if infile == "stdin" {
        Box::new(std::io::BufReader::new(std::io::stdin()))
    } else {
        let path = std::path::Path::new(infile);
        let file = match std::fs::File::open(path) {
            Err(why) => panic!("could not open {}: {}", path.display(), why),
            Ok(file) => file,
        };

        Box::new(std::io::BufReader::new(file))
    };

    let inner_writer = Box::new(std::io::BufWriter::new(
        std::fs::File::create(&outfile).unwrap(),
    ));

    let mut builder = bgzf::io::multithreaded_writer::Builder::default()
        .set_worker_count(opt_parallel);

    if compress_level >= 0 && compress_level <= 9 {
        use noodles_bgzf::io::writer::CompressionLevel;
        builder = builder.set_compression_level(CompressionLevel::new(compress_level as u8).unwrap());
    }

    let mut writer = builder.build_from_writer(inner_writer);

    //----------------------------
    // Output
    //----------------------------
    std::io::copy(&mut reader, &mut writer)?;
    writer.finish()?;

    // Generate GZI index
    build_gzi_index(&outfile)?;

    Ok(())
}

/// Build a .gzi index for a BGZF file
///
/// The GZI format is defined by `bgzip` and used for random access.
/// It consists of:
/// 1. A header (u64): Number of entries
/// 2. Entries (pairs of u64): (compressed_offset, uncompressed_offset)
///
/// Note:
/// * The format is Little-Endian.
/// * The first BGZF block (offset 0, 0) is implicitly skipped and NOT included in the index.
/// * Empty blocks (like EOF markers with ISIZE=0) are also skipped.
fn build_gzi_index(path: &str) -> anyhow::Result<()> {
    let mut file = std::fs::File::open(path)?;
    let mut index_data = Vec::new();
    let mut uncompressed_offset = 0;
    let mut compressed_offset = 0;

    loop {
        // Seek to start of current block
        file.seek(std::io::SeekFrom::Start(compressed_offset))?;
        
        // 1. Read fixed header (12 bytes)
        // [0-1]   ID1, ID2: 0x1f, 0x8b (GZIP magic)
        // [2]     CM: 8 (Deflate)
        // [3]     FLG: FEXTRA (4) must be set
        // [4-7]   MTIME
        // [8]     XFL
        // [9]     OS
        // [10-11] XLEN: Length of extra fields
        let mut header_fixed = [0u8; 12];
        match file.read_exact(&mut header_fixed) {
            Ok(_) => {},
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break, // EOF reached cleanly
            Err(e) => return Err(e.into()),
        }

        // Verify GZIP magic bytes
        // Verified against htslib/bgzf.c lines 64-78 (header definition)
        // and lines 2382-2415 (bgzf_index_dump_hfile logic).
        if header_fixed[0] != 0x1f || header_fixed[1] != 0x8b {
            break; // Not a GZIP block
        }
        
        // Verify BGZF flag (FEXTRA = 4)
        let flg = header_fixed[3];
        if (flg & 4) == 0 {
             break; // Standard GZIP, not BGZF
        }

        // Get Extra Field Length (XLEN)
        let xlen = u16::from_le_bytes([header_fixed[10], header_fixed[11]]) as u64;
        if xlen == 0 {
             break; // Should not happen in valid BGZF
        }

        // 2. Read Extra Fields
        let mut extra = vec![0u8; xlen as usize];
        file.read_exact(&mut extra)?;

        // 3. Find 'BC' subfield to get Block Size (BSIZE)
        // Subfield format: SI1(1), SI2(1), SLEN(2), DATA(SLEN)
        let mut bsize = 0u16;
        let mut cursor = 0;
        let mut found_bc = false;
        
        while cursor + 4 <= extra.len() {
             let si1 = extra[cursor];
             let si2 = extra[cursor+1];
             let slen = u16::from_le_bytes([extra[cursor+2], extra[cursor+3]]);
             
             // BGZF block size is stored in 'BC' subfield with SLEN=2
             if si1 == b'B' && si2 == b'C' && slen == 2 {
                 if cursor + 6 <= extra.len() {
                     bsize = u16::from_le_bytes([extra[cursor+4], extra[cursor+5]]);
                     found_bc = true;
                 }
                 break;
             }
             // Move to next subfield
             cursor += 4 + slen as usize;
        }

        if !found_bc {
            return Err(anyhow::anyhow!("Missing BC subfield in BGZF block at offset {}", compressed_offset));
        }

        // BSIZE is total block size - 1
        let block_size = bsize as u64 + 1;

        // 4. Read ISIZE (Input SIZE / Uncompressed Size)
        // Located at the last 4 bytes of the BGZF block
        file.seek(std::io::SeekFrom::Start(compressed_offset + block_size - 4))?;
        let mut isize_buf = [0u8; 4];
        file.read_exact(&mut isize_buf)?;
        let isize = u32::from_le_bytes(isize_buf) as u64;

        // 5. Record index entry
        // Rules compatible with `bgzip`:
        // - Skip the first block (offset 0)
        // - Skip empty blocks (ISIZE = 0), e.g., EOF marker
        if compressed_offset > 0 && isize > 0 {
            index_data.push((compressed_offset, uncompressed_offset));
        }

        // Advance offsets
        compressed_offset += block_size;
        uncompressed_offset += isize;
    }
    
    // Write the GZI index
    let index = bgzf::gzi::Index::from(index_data);
    let index_path = format!("{}.gzi", path);
    let mut writer = std::fs::File::create(index_path)?;
    bgzf::gzi::io::Writer::new(&mut writer).write_index(&index)?;

    Ok(())
}
