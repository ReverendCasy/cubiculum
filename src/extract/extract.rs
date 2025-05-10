//! # cubiculum::extract
//! 
//! 
//!
//! Author: Yury V.Malovichko
//!
//! Year: 2025

#[allow(dead_code)]

// use anyhow::{Error, Result};
use std::cmp;
use std::fmt::Display;
use std::ops;

use crate::structs::structs::{BedEntry, Coordinates};

#[derive(Debug)]
pub enum CubiculumError {
    ParseError(String),
    MissingTraitError(String),
    FormattingError(String),
}

impl Display for CubiculumError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            CubiculumError::ParseError(x) => {write!(f, "ParseError: {}", x)},
            CubiculumError::MissingTraitError(x) => {write!(f, "MissingTraitError: {}", x)},
            CubiculumError::FormattingError(x) => {write!(f, "FormattingError: {}", x)},
        }
    }
}

impl std::error::Error for CubiculumError {}

#[derive(PartialEq)]
pub enum BedFractionMode {
    All,
    Cds,
    Utr,
    Utr5,
    Utr3
}


/// Basic BED file line parser
/// 
/// # Arguments
/// 
pub fn parse_bed(
    line: String, format: usize, skip_blank: bool
) -> Option<BedEntry> {
    // BED file cannot contain less than three fields, and BED12+ are not currently accepted
    if format < 3 || format > 12 {
        panic!("Illegal BED file format specification! Accepted formats are BED3 through BED12");
    }
    if format == 10 || format == 11 {
        panic!(
            "BED10 and BED11 formats contain incomplete data on the sequence block structure. 
If you want to parse an incomplete BED entry, consider BED9 format instead"
        );
    } 
    let data: Vec<&str>  = line
        .trim()
        .split("\t")
        .collect::<Vec<&str>>();
    if data.len() == 0 {
        // panic if skip_blank was not set
        return None;
    }

    if (data.len() as usize) < format {
        // panic here
    }

    let chrom: String = data[0].to_string();
    let thin_start: u64 = data[1]
        .parse::<u64>()
        .expect("ThickStart is not a valid positive integer");
    let thin_end: u64 = data[2]
        .parse::<u64>()
        .expect("ThickEnd is not a valid positive integer");
    assert!(thin_start <= thin_end);

    if format == 3 {
        return Some(BedEntry::bed3(chrom, thin_start, thin_end));
    }

    let name: String = data[3].to_string();
    if format == 4 {
        return Some(BedEntry::bed4(chrom, thin_start, thin_end, name));
    }

    let score: String = data[4].to_string();
    if format == 5 {
        return Some(BedEntry::bed5(chrom, thin_start, thin_end, name, score));
    }

    let strand: bool = data[5] == "+";
    if format == 6 {
        return Some(BedEntry::bed6(chrom, thin_start, thin_end, name, score, strand));
    }

    let thick_start: u64 = data[6]
        .parse::<u64>()
        .expect("thinStart is not a valid positive integer");
    if thick_start < thin_start {
        panic!("thickStart value ({}) cannot be smaller than thinStart ({})", thick_start, thin_start)
    }
    let thick_end: u64 = data[7]
        .parse::<u64>()
        .expect("thinEnd is not a valid positive integer");
    if thick_end > thin_end {
        panic!("thickEnd value ({}) cannot be larger than thinEnd ({})", thick_end, thin_end)
    }
    if thick_start > thick_end {
        panic!("thickStart value ({}) cannot be larger than thickEnd ({})", thick_start, thick_end)
    }

    if format == 8 {
        return Some(BedEntry::bed8(chrom, thin_start, thin_end, name, score, strand, thick_start, thick_end))
    }

    let rgb: String = data[8].to_string();
    if format == 9 {
        return Some(
            BedEntry::bed9(chrom, thin_start, thin_end, name, score, strand, thick_start, thick_end, rgb)
        )
    }

    let ex_num: u16 = data[9]
        .parse::<u16>()
        .expect("Exon number is not a valid positive integer");
    let exon_sizes: Vec<u64> = data[10]
        .split(',')
        .filter(|x|
            !x.is_empty()
        )
        .map(|x|
            x.parse::<u64>().expect("Invalid exon size value")
        )
        .collect::<Vec<u64>>();
    let exon_starts: Vec<u64> = data[11]
        .split(',') 
        .filter(|x|
            !x.is_empty()
        )
        .map(|x|
            x.parse::<u64>().expect("Invalid exon start position")
        )
        .collect::<Vec<u64>>();
    return Some(
        BedEntry::bed12(
            chrom, thin_start, thin_end, name, score, strand, thick_start, thick_end, rgb, 
            ex_num, exon_sizes, exon_starts 
        )
    )
}

// pub fn extract_fraction(input: &BedEntry, mode: BedFractionMode, intron: bool) -> BedEntry {
        // let mut output
// }

/// Format a BedEntry object into a tab-separated BED file line
/// 
/// # Arguments
/// `bed_entry`: a BedEntry object to convert
/// `format`: number of columns in the output line, three through twelve.
/// (WARNING: BED12+ files are currently not accepted)
/// 
/// # Returns
/// A String representation of the input BedEntry
/// 
pub fn to_line(bed_entry: BedEntry, format: u8) -> Result<String, CubiculumError> {
    let entry_format = match bed_entry.format() {
        0 => {return Err(CubiculumError::MissingTraitError("Undefined BED format for the entry".to_string()))}
        x  => {x},
    };
    if entry_format < format {
        return Err(
            CubiculumError::FormattingError(
                format!("Cannot format BED{} entry into a BED{} line", format, entry_format)
            )
        );
    }
    if format < 3 || format == 7 || (format > 9 && format < 12) || format > 12 {
        return Err(
            CubiculumError::FormattingError(
                format!("Provided format BED{} is not supported. Accepted formats are : BED3,4,5,6,8,9,12", format)
            )
        );
    }
    let out = String::new();
    let chrom = match bed_entry.chrom() {
        Some(x) => {x},
        None => {return Err(CubiculumError::MissingTraitError("Undefined chromosome field".to_string()))}
    };
    let thin_start = match bed_entry.thin_start() {
        Some(x) => {x},
        None => {return Err(CubiculumError::MissingTraitError("Undefined thinStart field".to_string()))}
    };
    let thin_end = match bed_entry.thin_end() {
        Some(x) => {x},
        None => {return Err(CubiculumError::MissingTraitError("Undefined thinEnd field".to_string()))}
    };
    if format == 3 {
        return Ok(format!("{}\t{}\t{}", chrom, thin_start, thin_end));
    }
    let name = match bed_entry.name() {
        Some(x) => {x},
        None => {return Err(CubiculumError::MissingTraitError("Undefined name field".to_string()))}
    };
    if format == 4 {
        return Ok(format!("{}\t{}\t{}\t{}", chrom, thin_start, thin_end, name));
    }
    let score = match bed_entry.score() {
        Some(x) => {x},
        None => {return Err(CubiculumError::MissingTraitError("Undefined score field".to_string()))}
    };
    if format == 5 {
        return Ok(format!("{}\t{}\t{}\t{}\t{}", chrom, thin_start, thin_end, name, score));
    }
    let strand = match bed_entry.strand() {
        Some(x) => {
            match x {
                true => {'+'},
                false => {'-'}
            }
        },
        None => {return Err(CubiculumError::MissingTraitError("Undefined strand field".to_string()))}
    };
    if format == 6 {
        return Ok(format!("{}\t{}\t{}\t{}\t{}\t{}", chrom, thin_start, thin_end, name, score, strand));
    }
    let thick_start = match bed_entry.thick_start() {
        Some(x) => {x},
        None => {return Err(CubiculumError::MissingTraitError("Undefined thickStart field".to_string()))}
    };
    let thick_end = match bed_entry.thick_end() {
        Some(x) => {x},
        None => {return Err(CubiculumError::MissingTraitError("Undefined thickEnd field".to_string()))}
    };
    if format == 8 {
        return Ok(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                chrom, thin_start, thin_end, name, score, strand, thick_start, thick_end
            )
        );
    }
    let rgb = match bed_entry.rgb() {
        Some(x) => {x},
        None => {return Err(CubiculumError::MissingTraitError("Undefined Rgb field".to_string()))}
    };
    if format == 9 {
        return Ok(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                chrom, thin_start, thin_end, name, score, strand, thick_start, thick_end, rgb
            )
        );
    }
    let exon_num = match bed_entry.exon_num() {
        Some(x) => {x},
        None => {return Err(CubiculumError::MissingTraitError("Undefined exonNumber field".to_string()))}
    };
    let exon_sizes = match bed_entry.exon_sizes() {
        Some(x) => {
            x
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
                .join(",") + ","
        },
        None => {return Err(CubiculumError::MissingTraitError("Undefined exonSizes field".to_string()))}
    };
    let exon_starts = match bed_entry.exon_starts() {
        Some(x) => {
            x
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
                .join(",") + ","
        },
        None => {return Err(CubiculumError::MissingTraitError("Undefined exonStarts field".to_string()))}
    };
    return Ok(
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
            chrom, thin_start, thin_end, name, score, strand, 
            thick_start, thick_end, rgb,
            exon_num, exon_sizes, exon_starts
        )
    );
}

/// A highly optimized version for bed12ToFraction command line utility
/// 
/// # Arguments
/// 
/// * `line`: a BED12 format line string to parse
/// * `mode`: fraction of annotated blocks to report [accepted values: "all", "cds", "utr", "3utr", "5utr"]
/// * `intron`: boolean value specifying whether introns should be reported instead of exons
/// * `bed6`: boolean value specifying whether the resulting fraction should be split into separate BED6 records
/// 
/// 

pub fn bed_to_fraction(
    line: String, mode: &str, intron: bool, bed6: bool
) -> Option<String> {
    let mode: BedFractionMode = match mode {
        "all" => { BedFractionMode::All },
        "cds" => { BedFractionMode::Cds },
        "utr" => { BedFractionMode::Utr },
        "5utr" => { BedFractionMode::Utr5 },
        "3utr" => { BedFractionMode::Utr3 },
        _ => {
            panic!("Invalid 'mode' has been provided: {}. Valid modes are: all, cds, utr, 3utr, 5utr", mode)
        }
    };

    let data: Vec<&str>  = line
        .trim()
        .split("\t")
        .collect::<Vec<&str>>();
    if data.len() == 0 {
        return None;
    }
    if data.len() != 12 {
        panic!("Error: File contains improperly formatted lines. Make sure all lines in the file are in BED12 format");
    }
    let chrom: &str = data[0];
    let mut thin_start: u64 = data[1]
        .parse::<u64>()
        .expect("ThickStart is not a valid positive integer");
    let mut thin_end: u64 = data[2]
        .parse::<u64>()
        .expect("ThickEnd is not a valid positive integer");
    assert!(thin_start <= thin_end);
    let name: &str = data[3];
    let score: &str = data[4];
    let strand_line: &str = data[5];
    let strand: bool = strand_line == "+";
    let mut thick_start: u64 = data[6]
        .parse::<u64>()
        .expect("thinStart is not a valid positive integer");
    if thick_start < thin_start {
        panic!("thickStart value ({}) cannot be smaller than thinStart ({})", thick_start, thin_start)
    }
    let mut thick_end: u64 = data[7]
        .parse::<u64>()
        .expect("thinEnd is not a valid positive integer");
    if thick_end > thin_end {
        panic!("thickEnd value ({}) cannot be larger than thinEnd ({})", thick_end, thin_end)
    }
    if thick_start > thick_end {
        panic!("thickStart value ({}) cannot be larger than thickEnd ({})", thick_start, thick_end)
    }
    let rgb: &str = data[8];
    let ex_num: u64 = data[9]
        .parse::<u64>()
        .expect("Exon number is not a valid positive integer");
    let exon_sizes: Vec<u64> = data[10]
        .split(',')
        .filter(|x|
            !x.is_empty()
        )
        .map(|x|
            x.parse::<u64>().expect("Invalid exon size value")
        )
        .collect::<Vec<u64>>();
    let exon_starts: Vec<u64> = data[11]
        .split(',') 
        .filter(|x|
            !x.is_empty()
        )
        .map(|x|
            x.parse::<u64>().expect("Invalid exon start position")
        )
        .collect::<Vec<u64>>();

    // create shortcuts to control behaviour in UTR-targeted modes
    // the definition of 5' and 3' depends on the strand
    // the transcript is located on
    let report_up: bool = strand && mode == BedFractionMode::Utr5 || !strand && mode == BedFractionMode::Utr3;
    let report_down: bool = strand && mode == BedFractionMode::Utr3 || !strand && mode == BedFractionMode::Utr5;
    let noncoding: bool = (thick_end - thick_start) == 0;
    let report_coding: bool = !noncoding & (mode == BedFractionMode::Cds || mode == BedFractionMode::All);

    // infer the new sequence's start position
    let mut seq_start: u64 = match mode {
        BedFractionMode::Cds => thick_start, // will not change down the road
        // "3utr" => thick_end, // can be further set to the first 3'-UTR exon start
        _ => thin_start 
        // for "intron", will be set to the end of the first coding exon;
        // for utr, can be set to the start of the 3'-UTR
        // set in stone for 5utr 
    };

    // create storage objects for updated block coordinates
    let mut upd_block_starts: Vec<u64> = Vec::new();
    let mut upd_block_sizes: Vec<u64> = Vec::new();

    let range: ops::Range<u64>= if intron {0..ex_num-1} else {0..ex_num};
    if range.is_empty() {return None};
    for i in range {
        let i: usize = i as usize;

        let block_start: u64 = exon_starts[i] + thin_start;
        let block_end: u64 = block_start + exon_sizes[i];

        // first, check if the current block lies entirely within UTR
        // save block coordinates if respective mode is set, continue otherwise
        let upstream_to_cds: bool = block_end <= thick_start;
        let downstream_to_cds: bool = block_start >= thick_end;
        let upstream_and_report: bool = upstream_to_cds & (
            mode == BedFractionMode::Utr || report_up || noncoding
        );
        let downstream_and_report: bool = downstream_to_cds & (
            mode == BedFractionMode::Utr || report_down || noncoding
        );
        // current block is either completely upstream or completely downstream to CDS 
        if upstream_to_cds || downstream_to_cds {
            if upstream_and_report || downstream_and_report || mode == BedFractionMode::All {
                if intron {
                    // update the starting position for the intron track
                    if upd_block_starts.len() == 0 {seq_start = block_end};
                    upd_block_starts.push(block_end - seq_start);
                    upd_block_sizes.push(exon_starts[i+1] + thin_start - block_end);
                } else {
                    // update the starting position if the needed sequence portion started only with 3'-UTR
                    if downstream_to_cds && upd_block_starts.len() == 0 {
                        seq_start = block_start - thin_start;
                    }
                    upd_block_starts.push(block_start - seq_start);
                    upd_block_sizes.push(block_end - block_start);
                }
            };
            continue;
        }

        // if we reached this point, the block belongs to
        // the coding sequence, at least partially

        // for 3'-UTR/5'-UTR on the negative strand, we can safely skip blocks up until the CDS end
        if report_down && block_end <= thick_end {continue}; 

        // for introns, boundaries are block end and next block's start
        if intron & report_coding {
            if upd_block_starts.len() == 0 {seq_start = block_end};
            if block_end >= thick_end {break};
            upd_block_starts.push(block_end - seq_start);
            upd_block_sizes.push(exon_starts[i+1] + thin_start - block_end);
            continue;
        };

        if mode == BedFractionMode::All {
            upd_block_starts.push(block_start - seq_start);
            upd_block_sizes.push(block_end - block_start);
            continue
        }

        // for blocks overlapping with the coding sequence, assess their boundaries
        let upd_block_start: u64 = cmp::max(block_start, thick_start);
        let upd_block_end: u64 = cmp::min(block_end, thick_end);

        // for 'cds' mode, save the updated coordinates
        if mode == BedFractionMode::Cds {
            upd_block_starts.push(upd_block_start - seq_start);
            upd_block_sizes.push(upd_block_end - upd_block_start);
            continue
        }

        // for UTR-related modes, clip the coding part and save the UTR bases
        if (upd_block_start > block_start) & (mode == BedFractionMode::Utr || report_up) & !intron{
                upd_block_starts.push(block_start - seq_start);
                upd_block_sizes.push(upd_block_start - block_start);
                // for 5'-UTR/3'-UTR on the negative strand, the loop can be safely exited
                if report_up {break};
                continue
        }
        if (upd_block_end < block_end) & (mode == BedFractionMode::Utr || report_down) {
            if intron {
                if upd_block_starts.len() == 0 {seq_start = block_end};
                upd_block_starts.push(block_end - seq_start);
                upd_block_sizes.push(exon_starts[i+1] + thin_start - block_end);
            } else {
                if upd_block_starts.len() == 0 {seq_start = upd_block_end};
                upd_block_starts.push(upd_block_end - seq_start);
                upd_block_sizes.push(block_end - upd_block_end);    
            }
        }
        // // if the mode was set to 'utr' but no upstream UTR was found so far,
        // // updating the starting point
        // if (mode == "utr" || mode == "3utr") && upd_block_sizes.len() == 0 {
        //     seq_start = thick_end
        // }
    };

    assert!(upd_block_starts.len() == upd_block_sizes.len());
    let upd_block_count: usize = upd_block_sizes.len();
    if upd_block_count == 0 {return None};

    // set the start and end points
    match (mode, intron) {
        (BedFractionMode::All, false) => {
            // pass
        }
        (BedFractionMode::Cds, false) => {
            thin_start = thick_start;
            thin_end = thick_end;
        },
        (BedFractionMode::Utr, false) => {
            thick_end = thin_end;
            thick_start = thick_end;
        },
        _ => {
            thin_start = seq_start;
            thin_end = seq_start + upd_block_starts[..].last().unwrap() + upd_block_sizes[..].last().unwrap();
            thick_end = thin_end;
            thick_start = thick_end;
        }
    };

    // if bed6 output is expected, modify the lines
    if bed6 {
        let mut bed6_line: String = String::new();
        for i in 0..upd_block_count {
            let block_start: u64 = seq_start + upd_block_starts[i];
            let block_end: u64 = block_start + upd_block_sizes[i];
            let block_num: u64 = if strand {i as u64 + 1} else {(upd_block_count - i) as u64};
            let out_line: String = format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                chrom, block_start, block_end, name, block_num, strand_line
            );
            bed6_line.push_str(&out_line);
            if i < upd_block_count - 1 {
                bed6_line.push('\n');
            }
        }
        return Some(bed6_line);
    }
    let size_line: String = upd_block_sizes
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
        .join(",") + ",";
    let start_line: String = upd_block_starts
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
        .join(",") + ",";

    let result: String = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        chrom, thin_start, thin_end, name, score, strand_line, 
        thick_start, thick_end, rgb, upd_block_count, size_line, start_line
    );
    return Some(result)

    // let result: String = String::new();
    // return Some(result);
}

// //////////////
// UNIT TESTS
// //////////////

#[cfg(test)]
mod test {
    use super::*;

    // PANIC TESTS
    #[test]
    #[should_panic]
    fn invalid_mode_test() {
        // should panic due to an unknown mode name
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        bed_to_fraction(input, "hello", false, false);
    }

    #[test]
    #[should_panic]
    fn truncated_line_test(){
        // should panic due to the input line containing the number of columns different from twelve
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,");
        bed_to_fraction(input, "cds", false, false);
    }

    #[test]
    #[should_panic]
    fn invalid_value() {
        // should panic due to one of the numeric fields occupied by a non-numeric value
        let input: String = String::from("chr9	AAA	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        bed_to_fraction(input, "cds", false, false);
    }

    #[test]
    #[should_panic]
    fn negative_length(){
        // should panic due to start value exceeding end value
        let input: String = String::from("chr9	101385006	101360416	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        bed_to_fraction(input, "cds", false, false);
    }

    #[test]
    #[should_panic]
    fn out_of_boundary_cds() {
        // should panic due to thickEnd exceeding thinEnd
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101385008	0	4	2599,203,525,152,	0,7703,10522,24438,");
        bed_to_fraction(input, "cds", false, false);
    }

    // PERFORMANCE TESTS
    #[test]
    fn basic_test() {
        // tests whether the bed_to_fraction can return the same line as provided as input
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        let expected: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        assert_eq!(expected, bed_to_fraction(input, "all", false, false).unwrap());
    }

    #[test]
    fn basic_cds_test() {
        // tests whether "cds" mode returns the same line as input if a record contains coding sequence only
        let input: String = String::from("chr1	149156055	149163998	XM_047439510.1#LOC124904581	0	+	149156055	149163998	0	4	36,75,602,112,	0,2796,6686,7831,");
        let expected: String = String::from("chr1	149156055	149163998	XM_047439510.1#LOC124904581	0	+	149156055	149163998	0	4	36,75,602,112,	0,2796,6686,7831,");
        assert_eq!(expected, bed_to_fraction(input, "cds", false, false).unwrap());
    }

    #[test]
    fn cds_exon_test() {
        // tests the cds mode for a sequence with both merged and intron-separated UTR exons present
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        let expected: String = String::from("chr9	101362427	101371404	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	3	588,203,466,	0,5692,8511,");
        assert_eq!(expected, bed_to_fraction(input, "cds", false, false).unwrap());
    }

    #[test]
    fn single_exon_cds_with_utrs() {
        // tests cds mode on a single exon sequence with arbitrary UTRs
        let input: String = String::from("chr9	129489948	129513686	XM_047424327.1#LINC00963	0	+	129490480	129491083	0	4	1180,177,350,268,	0,3470,13374,23470,");
        let expected: String = String::from("chr9	129490480	129491083	XM_047424327.1#LINC00963	0	+	129490480	129491083	0	1	603,	0,");
        assert_eq!(expected, bed_to_fraction(input, "cds", false, false).unwrap());
    }

    #[test]
    fn intron_test() {
        // tests the intron mode
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        let expected: String = String::from("chr9	101363015	101370938	ENST00000259407.7#BAAT	0	-	101370938	101370938	0	2	5104,2616,	0,5307,");
        assert_eq!(expected, bed_to_fraction(input, "cds", true, false).unwrap());
    }

    #[test]
    fn uuu() {
        // tests the intron mode
        let input: String = String::from("chr19	47403123	47422233	NM_001346148.2#MEIS3	0	-	47406476	47422191	0	13	430,67,84,59,77,149,112,150,51,51,160,173,45,	0,3336,3764,3955,4228,5975,6312,11593,11927,13528,13680,14054,19065,");
        println!("{}", bed_to_fraction(input, "cds", true, false).unwrap());
    }

    #[test]
    fn zero_intron_test(){
        // tests intron mode for single-exon transcripts; must return None
        let input: String = String::from("chr9	129490480	129491083	XM_047424327.1#LINC00963	0	+	129490480	129491083	0	1	603,	0,");
        // assert_eq!(None, bed_to_fraction(input, "intron", false));
        assert!(bed_to_fraction(input, "cds", true, false).is_none());
    }

    #[test]
    fn zero_intron_with_utrs_test() {
        // the same as the test above but in the presence of UTR blocks
        let input: String = String::from("chr9	129489948	129513686	XM_047424327.1#LINC00963	0	+	129490480	129491083	0	4	1180,177,350,268,	0,3470,13374,23470,");
        assert!(bed_to_fraction(input, "cds", true, false).is_none());
    }

    #[test]
    fn full_utr_test() {
        // tests full UTR mode
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        let expected: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101385006	101385006	0	3	2011,59,152,	0,10988,24438,");
        assert_eq!(expected, bed_to_fraction(input, "utr", false, false).unwrap());
    }

    #[test]
    fn no_utr_test(){
        // tests utr mode for pure CDS record; must return None
        let input: String = String::from("chr9	101362427	101371404	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	3	588,203,466,	0,5692,8511,");
        assert!(bed_to_fraction(input, "utr", false, false).is_none());
    }

    #[test]
    fn fiveprime_utr_forward_test() {
        // tets 5utr mode for a plus-strand transcript
        let input: String = String::from("chr19	45692665	45703987	NM_001163377.2#QPCTL	0	+	45692703	45703049	0	6	245,144,153,100,117,1084,	0,747,5881,6135,9132,10238,");
        let expected: String = String::from("chr19	45692665	45692703	NM_001163377.2#QPCTL	0	+	45692703	45692703	0	1	38,	0,");
        assert_eq!(expected, bed_to_fraction(input, "5utr", false, false).unwrap());
    }

    #[test]
    fn fiveprime_utr_reverse_test(){
        // tets 5utr mode for a minus-strand transcript
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        let expected: String = String::from("chr9	101371404	101385006	ENST00000259407.7#BAAT	0	-	101385006	101385006	0	2	59,152,	0,13450,");
        assert_eq!(expected, bed_to_fraction(input, "5utr", false, false).unwrap());
    }

    #[test]
    fn threeprime_utr_forward_test(){
        // tets 3utr mode for a plus-strand transcript
        let input: String = String::from("chr19	45692665	45703987	NM_001163377.2#QPCTL	0	+	45692703	45703049	0	6	245,144,153,100,117,1084,	0,747,5881,6135,9132,10238,");
        let expected: String = String::from("chr19	45703049	45703987	NM_001163377.2#QPCTL	0	+	45703987	45703987	0	1	938,	0,");
        assert_eq!(expected, bed_to_fraction(input, "3utr", false, false).unwrap());
    }

    #[test]
    fn threeprime_utr_reverse_test(){
        // tets 3utr mode for a minus-strand transcript
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        let expected: String = String::from("chr9	101360416	101362427	ENST00000259407.7#BAAT	0	-	101362427	101362427	0	1	2011,	0,");
        assert_eq!(expected, bed_to_fraction(input, "3utr", false, false).unwrap());
    }

    #[test]
    fn full_utr_noncoding_test() {
        // tests utr mode performance on pseudogenes; must return the same line as input
        let input: String = String::from("chr1	3205900	3216344	ENSMUST00000162897	0	-	3216344	3216344	0	2	1417,2736,	0,7708,");
        assert_eq!(input, bed_to_fraction(input.clone(), "utr", false, false).unwrap());
    }

    #[test]
    fn fiveprime_utr_noncoding_test(){
        // tests side-bound utr mode for pseudogenes; returns the same string as input
        let input: String = String::from("chr1	3205900	3216344	ENSMUST00000162897	0	-	3216344	3216344	0	2	1417,2736,	0,7708,");
        assert_eq!(input, bed_to_fraction(input.clone(), "5utr", false, false).unwrap());
    }

    #[test]
    fn utr_intron_test() {
        // tests intron mode for UTRs only
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        let expected: String = String::from("chr9	101371463	101384854	ENST00000259407.7#BAAT	0	-	101384854	101384854	0	1	13391,	0,");
        assert_eq!(expected, bed_to_fraction(input, "utr", true, false).unwrap());
    }

    #[test]
    fn fiveprime_utr_intron_test() {
        // tests intron mode for 5'-UTRs only; returns the same as the test above
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        let expected: String = String::from("chr9	101371463	101384854	ENST00000259407.7#BAAT	0	-	101384854	101384854	0	1	13391,	0,");
        assert_eq!(expected, bed_to_fraction(input, "5utr", true, false).unwrap());
    }

    #[test]
    fn fiveprime_utr_multiintron_test(){
        // tests intron mode for 5'-UTR in case of multiple 5'-UTR introns
        let input: String = String::from("chr19	46746056	46758575	ENST00000318584.10#FKRP	0	+	46755450	46756938	0	4	34,62,151,3164,	0,1970,2458,9355,");
        let expected: String = String::from("chr19	46746090	46755411	ENST00000318584.10#FKRP	0	+	46755411	46755411	0	3	1936,426,6746,	0,1998,2575,");
        assert_eq!(expected, bed_to_fraction(input, "5utr", true, false).unwrap());
    }

    #[test]
    fn threeprime_utr_intron_test() {
        // tests intron mode for 3'-UTRs only
        let input: String = String::from("chr19	47403123	47422233	NM_001346148.2#MEIS3	0	-	47406476	47422191	0	13	430,67,84,59,77,149,112,150,51,51,160,173,45,	0,3336,3764,3955,4228,5975,6312,11593,11927,13528,13680,14054,19065,");
        let expected: String = String::from("chr19	47403553	47406459	NM_001346148.2#MEIS3	0	-	47406459	47406459	0	1	2906,	0,");
        assert_eq!(expected, bed_to_fraction(input, "3utr", true, false).unwrap());
    }

    #[test]
    fn threeprime_utr_positive_strand_test() {
        // same as the test above, but this time the strand is positive
        let input: String = String::from("chr19	47778702	47784682	ENST00000601048.6#SELENOW	0	+	47778785	47781370	0	6	112,25,54,75,99,393,	0,2022,2161,2405,2587,5587,");
        let expected: String = String::from("chr19	47781388	47784289	ENST00000601048.6#SELENOW	0	+	47784289	47784289	0	1	2901,	0,");
        assert_eq!(expected, bed_to_fraction(input, "3utr", true, false).unwrap());
    }

    #[test]
    fn no_utr_intron_test()
    {
        // tests intron for UTRs for a transcript with no UTR introns; None is expected
        let input: String = String::from("chr19	45692665	45703987	NM_001163377.2#QPCTL	0	+	45692703	45703049	0	6	245,144,153,100,117,1084,	0,747,5881,6135,9132,10238,");
        assert!(bed_to_fraction(input, "utr", true, false).is_none());
    }

    #[test]
    fn all_bed6_test() {
        let input: String = String::from("chr19	47778702	47784682	ENST00000601048.6#SELENOW	0	+	47778785	47781370	0	6	112,25,54,75,99,393,	0,2022,2161,2405,2587,5587,");
        let expected: String = String::from(
            "chr19	47778702	47778814	ENST00000601048.6#SELENOW	1	+
chr19	47780724	47780749	ENST00000601048.6#SELENOW	2	+
chr19	47780863	47780917	ENST00000601048.6#SELENOW	3	+
chr19	47781107	47781182	ENST00000601048.6#SELENOW	4	+
chr19	47781289	47781388	ENST00000601048.6#SELENOW	5	+
chr19	47784289	47784682	ENST00000601048.6#SELENOW	6	+"
        );
        assert_eq!(expected, bed_to_fraction(input, "all", false, true).unwrap());
    }

    #[test]
    fn cds_bed6_test() {
        let input: String = String::from("chr19	47778702	47784682	ENST00000601048.6#SELENOW	0	+	47778785	47781370	0	6	112,25,54,75,99,393,	0,2022,2161,2405,2587,5587,");
        let expected: String = String::from(
            "chr19	47778785	47778814	ENST00000601048.6#SELENOW	1	+
chr19	47780724	47780749	ENST00000601048.6#SELENOW	2	+
chr19	47780863	47780917	ENST00000601048.6#SELENOW	3	+
chr19	47781107	47781182	ENST00000601048.6#SELENOW	4	+
chr19	47781289	47781370	ENST00000601048.6#SELENOW	5	+"
        );
        assert_eq!(expected, bed_to_fraction(input, "cds", false, true).unwrap());
    }

    #[test]
    fn cds_intron_bed6_test() {
        let input: String = String::from("chr9	101360416	101385006	ENST00000259407.7#BAAT	0	-	101362427	101371404	0	4	2599,203,525,152,	0,7703,10522,24438,");
        let expected: String = String::from(
            "chr9	101363015	101368119	ENST00000259407.7#BAAT	2	-
chr9	101368322	101370938	ENST00000259407.7#BAAT	1	-"
        );
        assert_eq!(expected, bed_to_fraction(input, "cds", true, true).unwrap());
    }

    #[test]
    fn all_intron_bed6_test() {
        // tests bed6 mode for all introns in the sequence
        let input: String = String::from("chr19	47403123	47422233	NM_001346148.2#MEIS3	0	-	47406476	47422191	0	13	430,67,84,59,77,149,112,150,51,51,160,173,45,	0,3336,3764,3955,4228,5975,6312,11593,11927,13528,13680,14054,19065,");
        let expected: String = String::from(
            "chr19	47403553	47406459	NM_001346148.2#MEIS3	12	-
chr19	47406526	47406887	NM_001346148.2#MEIS3	11	-
chr19	47406971	47407078	NM_001346148.2#MEIS3	10	-
chr19	47407137	47407351	NM_001346148.2#MEIS3	9	-
chr19	47407428	47409098	NM_001346148.2#MEIS3	8	-
chr19	47409247	47409435	NM_001346148.2#MEIS3	7	-
chr19	47409547	47414716	NM_001346148.2#MEIS3	6	-
chr19	47414866	47415050	NM_001346148.2#MEIS3	5	-
chr19	47415101	47416651	NM_001346148.2#MEIS3	4	-
chr19	47416702	47416803	NM_001346148.2#MEIS3	3	-
chr19	47416963	47417177	NM_001346148.2#MEIS3	2	-
chr19	47417350	47422188	NM_001346148.2#MEIS3	1	-"
        );
        assert_eq!(expected, bed_to_fraction(input, "all", true, true).unwrap());
    }

}
