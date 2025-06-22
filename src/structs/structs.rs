use std::cmp::{min, max};

use crate::extract::extract::{parse_bed, to_line};
use crate::merge::merge::{intersection, merge_multiple};

/// Contains data on storage structures for annotation manipulations in Cubiculum and associated packages

#[derive(Clone, Debug)]
pub struct Interval {
    chrom: Option<String>,
    start: Option<u64>,
    end: Option<u64>,
    name: Option<String>
}

impl Interval {
    pub fn new() -> Interval {
        Interval { chrom: None, start: None, end: None, name: None }
    }

    pub fn from(chrom: Option<String>, start: Option<u64>, end: Option<u64>, name: Option<String>) -> Interval {
        Interval {chrom: chrom, start: start, end: end, name: name}
    }

    pub fn update_name(&mut self, name: String) {
        self.name = Some(name);
    }

    pub fn update_chrom(&mut self, chrom: String) {
        self.chrom = Some(chrom);
    }

    pub fn update_start(&mut self, start: u64) {
        self.start = Some(start);
    }

    pub fn update_end(&mut self, end: u64) {
        self.end = Some(end);
    }
}
#[derive(Clone, Debug)]
pub struct BedEntry{
    format: Option<u8>,
    chrom: Option<String>,
    thin_start: Option<u64>,
    thin_end: Option<u64>,
    name: Option<String>,
    score: Option<String>,
    strand: Option<bool>,
    thick_start: Option<u64>,
    thick_end: Option<u64>,
    rgb: Option<String>,
    exon_num: Option<u16>,
    exon_sizes: Option<Vec<u64>>,
    exon_starts: Option<Vec<u64>>
}

impl BedEntry{
    pub fn empty() -> BedEntry{
        BedEntry{
            format: None,
            chrom: None, 
            thin_start: None, 
            thin_end: None, 
            name: None, 
            score: None, 
            strand: None, 
            thick_start: None, 
            thick_end: None, 
            rgb: None, 
            exon_num: None, 
            exon_sizes: None, 
            exon_starts: None
        }
    }

    pub fn from_interval<T>(inter: T) -> Option<BedEntry> 
    where T:
        Coordinates
    {
        let mut output: BedEntry = BedEntry::empty();
        let mut format: u8 = 0;
        let chrom = match inter.chrom() {
            Some(x) => {x.clone()},
            None => {return None}
        };
        let thin_start = match inter.start() {
            Some(x) => {*x},
            None => {return None}
        };
        let thin_end = match inter.end() {
            Some(x) => {*x},
            None => {return None}
        };
        Some(BedEntry::bed3(chrom, thin_start, thin_end))
    }

    pub fn bed3(chrom: String, start: u64, end: u64) -> BedEntry {
        BedEntry{
            format: Some(3), 
            chrom: Some(chrom), 
            thin_start: Some(start), 
            thin_end: Some(end), 
            name: None, 
            score: None, 
            strand: None, 
            thick_start: None, 
            thick_end: None, 
            rgb: None, 
            exon_num: None, 
            exon_sizes: None, 
            exon_starts: None
        }
    }

    pub fn bed4(chrom: String, start: u64, end: u64, name: String) -> BedEntry {
        BedEntry{
            format: Some(4), 
            chrom: Some(chrom), 
            thin_start: Some(start), 
            thin_end: Some(end), 
            name: Some(name), 
            score: None, 
            strand: None, 
            thick_start: None, 
            thick_end: None, 
            rgb: None, 
            exon_num: None, 
            exon_sizes: None, 
            exon_starts: None
        }
    }

    pub fn bed5(chrom: String, start: u64, end: u64, name: String, score: String) -> BedEntry {
        BedEntry{
            format: Some(5), 
            chrom: Some(chrom), 
            thin_start: Some(start), 
            thin_end: Some(end), 
            name: Some(name), 
            score: Some(score), 
            strand: None, 
            thick_start: None, 
            thick_end: None, 
            rgb: None, 
            exon_num: None, 
            exon_sizes: None, 
            exon_starts: None
        }
    }

    pub fn bed6(chrom: String, start: u64, end: u64, name: String, score: String, strand: bool) -> BedEntry {
        BedEntry{
            format: Some(6), 
            chrom: Some(chrom), 
            thin_start: Some(start), 
            thin_end: Some(end), 
            name: Some(name), 
            score: Some(score), 
            strand: Some(strand), 
            thick_start: None, 
            thick_end: None, 
            rgb: None, 
            exon_num: None, 
            exon_sizes: None, 
            exon_starts: None
        }
    }

    pub fn bed8(
        chrom: String, start: u64, end: u64, name: String, score: String, strand: bool, 
        thick_start: u64, thick_end: u64 
    ) -> BedEntry {
        BedEntry{
            format: Some(8), 
            chrom: Some(chrom), 
            thin_start: Some(start), 
            thin_end: Some(end), 
            name: Some(name), 
            score: Some(score), 
            strand: Some(strand), 
            thick_start: Some(thick_start), 
            thick_end: Some(thick_end), 
            rgb: None, 
            exon_num: None, 
            exon_sizes: None, 
            exon_starts: None
        }
    }

    pub fn bed9(
        chrom: String, start: u64, end: u64, name: String, score: String, strand: bool, 
        thick_start: u64, thick_end: u64, rgb: String
    ) -> BedEntry {
        BedEntry{
            format: Some(9), 
            chrom: Some(chrom), 
            thin_start: Some(start), 
            thin_end: Some(end), 
            name: Some(name), 
            score: Some(score), 
            strand: Some(strand), 
            thick_start: Some(thick_start), 
            thick_end: Some(thick_end), 
            rgb: Some(rgb), 
            exon_num: None, 
            exon_sizes: None, 
            exon_starts: None
        }
    }

    pub fn bed12(
        chrom: String, start: u64, end: u64, name: String, score: String, strand: bool, 
        thick_start: u64, thick_end: u64, rgb: String, 
        exon_num: u16, exon_sizes: Vec<u64>, exon_starts: Vec<u64>
    ) -> BedEntry {
        BedEntry{
            format: Some(12), 
            chrom: Some(chrom), 
            thin_start: Some(start), 
            thin_end: Some(end), 
            name: Some(name), 
            score: Some(score), 
            strand: Some(strand), 
            thick_start: Some(thick_start), 
            thick_end: Some(thick_end), 
            rgb: Some(rgb), 
            exon_num: Some(exon_num), 
            exon_sizes: Some(exon_sizes), 
            exon_starts: Some(exon_starts)
        }
    }

    pub fn format(&self) -> u8 {
        if let None = self.format {return 0};
        self.format.unwrap()
    }

    pub fn thin_start(&self) -> Option<u64> {
        self.thin_start
    }

    pub fn thin_end(&self) -> Option<u64> {
        self.thin_end
    }

    pub fn name(&self) -> Option<&String> {
        self.name.as_ref()
    }

    pub fn score(&self) -> Option<&String> {
        self.score.as_ref()
    }
    pub fn strand(&self) -> Option<bool> {
        self.strand
    }

    pub fn thick_start(&self) -> Option<u64> {
        self.thick_start
    }

    pub fn thick_end(&self) -> Option<u64> {
        self.thick_end
    }

    pub fn rgb(&self) -> Option<&String> {
        self.rgb.as_ref()
    }

    pub fn exon_num(&self) -> Option<u16> {
        self.exon_num
    }

    pub fn exon_sizes(&self) -> Option<&Vec<u64>> {
        self.exon_sizes.as_ref()
    }

    pub fn exon_starts(&self) -> Option<&Vec<u64>> {
        self.exon_starts.as_ref()
    }

    pub fn update_thin_start(&mut self, thin_start: u64) {
        self.thin_start = Some(thin_start)
    }

    /// Returns the length sum for all the blocks
    /// 
    pub fn block_length(&self) -> u64 {
        if self.format() < 12 {
            return match (self.thin_start(), self.thin_end()) {
                (Some(x), Some(y)) => {y - x},
                _ => {0}
            };
        }
        let mut length_sum: u64 = 0;
        for x in self.exon_sizes().unwrap() {
            length_sum += *x
        }
        return length_sum;
    }

    pub fn to_interval(&mut self) -> Interval {
        Interval::from(
            self.chrom.clone(),
            self.thin_start, 
            self.thin_end, 
            self.name.clone()
        )
    }

    pub fn to_blocks(&self) -> Option<Vec<BedEntry>> {
        if self.format.unwrap() != 12 {
            return None;
        }
        let chrom: &str = match &self.chrom {
            Some(x) => {x},
            None => {return None}
        };
        let thin_start: u64 = match self.thin_start {
            Some(x) => {x},
            None => {return None}
        };
        let name: &str = match &self.name {
            Some(x) => {x},
            None => {return None}
        };
        let ex_num = match self.exon_num {
            Some(x) => {x as usize},
            None => {return None}
        };
        let score: &str = match &self.score {
            Some(x) => {x},
            None => {"0"}
        };
        let strand: bool = match self.strand {
            Some(x) => {x},
            None => {return None}
        };
        let mut blocks: Vec<BedEntry> = Vec::with_capacity(ex_num);
        for i in 0..ex_num {
            let start: u64 = thin_start + self.exon_starts.as_ref().unwrap()[i];
            let end: u64 = start + self.exon_sizes.as_ref().unwrap()[i];
            blocks.push(
                BedEntry::bed6(
                    chrom.to_string(),
                    start, end, 
                    name.to_string(), 
                    score.to_string(), 
                    strand
                )
            );
        }
        Some(blocks)
    }

    pub fn clip_by(&mut self, start: Option<u64>, end: Option<u64>, inplace: bool) -> Option<BedEntry> {
        let chrom: &str = match &self.chrom {
            Some(x) => {x},
            None => {return None}
        };
        let thin_start: u64 = match self.thin_start {
            Some(x) => {x},
            None => {return None}
        };
        let name: &str = match &self.name {
            Some(x) => {x},
            None => {return None}
        };
        let mut new_thin_start: u64 = match start {
            Some(x) => {max(self.thin_start.unwrap(), x)},
            None => self.thin_start.unwrap()
        };
        let mut new_thin_end: u64 = match end {
            Some(x) => {min(self.thin_end.unwrap(), x)},
            None => {self.thin_end.unwrap()}
        };
        let mut new_thick_start = match self.thick_start {
            Some(x) => {max(x, new_thin_start)},
            None => {new_thin_start}
        };
        let mut new_thick_end = match self.thick_end {
            Some(x) => {min(x, new_thin_end)},
            None => {new_thin_end}
        };
        let (new_ex_num, new_ex_sizes, new_ex_starts) = match (self.exon_num, &self.exon_sizes, &self.exon_starts) {
            (Some(x), Some(y), Some(z)) => {
                let mut ex_counter: u16 = 0;
                let mut _sizes: Vec<u64> = Vec::new();
                let mut _starts: Vec<u64> = Vec::new();
                for i in 0..y.len() {
                    let ex_start: u64= z[i] + self.thin_start.unwrap();
                    if ex_start > new_thin_end {break};
                    let ex_end: u64 = y[i] + ex_start;
                    if ex_end < new_thin_start {continue};
                    let new_ex_start = max(new_thin_start, ex_start);
                    if new_ex_start > new_thin_start {
                        new_thin_start = new_ex_start;
                        new_thick_start = max(new_thin_start, new_thick_start);
                    }
                    let new_ex_size = match min(ex_end, new_thin_end).checked_sub(new_ex_start){
                        Some(x) => {if x == 0 {continue} else {x}},
                        None => {continue} 
                    };
                    let new_ex_end = new_ex_start + new_ex_size;
                    if new_ex_end < new_thin_end {
                        new_thin_end = new_ex_end;
                        new_thick_end = min(new_thin_end, new_thick_end);
                    }
                    _sizes.push(new_ex_size);
                    _starts.push(new_ex_start - new_thin_start);
                    ex_counter += 1;
                }
                (Some(ex_counter), Some(_sizes), Some(_starts))
            },
            _ => {(None, None, None)}
        };
        if inplace {
            self.thin_start = Some(new_thin_start);
            self.thin_end = Some(new_thin_end);
            self.thick_start = Some(new_thick_start);
            self.thick_end = Some(new_thick_end);
            if let Some(x) = Some(new_ex_num) {
                self.exon_num = new_ex_num
            };
            if let Some(x) = new_ex_sizes {
                self.exon_sizes = Some(x)
            };
            if let Some(x) = new_ex_starts {
                self.exon_starts = Some(x)
            };
            return None;
        };
        let mut clipped_bed = BedEntry::empty();
        // TODO: rewrite with if-lets
        clipped_bed.format = self.format;
        clipped_bed.chrom = match &self.chrom {
            Some(x) => {Some(x.clone())},
            None => None
        };
        clipped_bed.thin_start = Some(new_thin_start);
        clipped_bed.thin_end = Some(new_thin_end);
        clipped_bed.name = match &self.name {
            Some(x) => {Some(x.clone())},
            None => {None}
        };
        clipped_bed.score = match &self.score {
            Some(x) => {Some(x.clone())},
            None => {None}
        };
        clipped_bed.strand = match &self.strand {
            Some(x) => Some(*x),
            None => {None}
        };
        clipped_bed.thick_start = Some(new_thick_start);
        clipped_bed.thick_end = Some(new_thick_end);
        clipped_bed.exon_num = new_ex_num;
        clipped_bed.exon_sizes = new_ex_sizes;
        clipped_bed.exon_starts = new_ex_starts;
        Some(clipped_bed)

    }
    
    pub fn to_cds(&mut self, inplace: bool)  -> Option<BedEntry> {
        if self.format.unwrap() < 8 {return None};
        self.clip_by(self.thick_start, self.thick_end, inplace)
    }

    pub fn graft<T>(
        &mut self, graft: T, inplace: bool,
        chrom_compatible: bool,
        allow_overlaps: bool, 
        coding: bool,
        append_upstream: bool, 
        append_downstream: bool,
    ) -> Option<BedEntry> 
    where
        T: Coordinates + Clone
    {
        if append_upstream && append_downstream {
            panic!("Cannot append from both up- and downstream sides");
        }
        if self.format() != 12 {
            panic!("Cannot graft to a non-BED12 object");
        }
        if chrom_compatible {
            match (self.chrom(), graft.chrom()) {
                (Some(x), Some(y)) => {
                    if x != y {
                        panic!("BED12 and graft are located on different chromosomes ({} and {})", x, y)
                    }
                },
                _ => {panic!("Undefined chromosome for either BED12 or graft when `chrom_compatible` was set")}
            }
        }

        let mut thin_start = match self.thin_start {
            Some(x) => {x},
            None => {panic!("Undefined thinStart value for BED12")}
        };
        let mut thick_start = match self.thick_start {
            Some(x) => {x},
            None => {panic!("CRITICAL: Undefined thickStart value for BED12")}
        };
        let mut thin_end = match self.thin_end {
            Some(x) => {x},
            None => {panic!("CRITICAL: Undefined thinEnd value for BED12")}
        };
        let mut thick_end = match self.thick_end {
            Some(x) => {x},
            None => {panic!("CRITICAL: Undefined thickEnd value for BED12")}
        };
        
        let mut exon_num = match self.exon_num {
            Some(x) => {x},
            None => {panic!("CRITICAL: Exon number is not defined for the BED12 object")}
        };

        let mut exon_sizes = match &mut self.exon_sizes {
            Some(x) => {x.clone()},
            None => {panic!("CRITICAL: Exon sizes are not defined for the BED12 object")}
        };
        let mut exon_starts = match &mut self.exon_starts {
            Some(x) => {x.clone()},
            None => {panic!("CRITICAL: Exon starts are not defined for the BED12 object")}
        };

        let graft_start = match graft.start() {
            Some(x) => {*x},
            None => {panic!("CRITICAL: Undefined start coordinate for a grafted interval")}
        };
        let graft_end = match graft.end() {
            Some(x) => {*x},
            None => {panic!("CRITICAL: Undefined end coordinate for a grafted interval")}
        };
        let mut graft_len = graft.length().unwrap();

        // keep track on whether the final block should be merged
        let mut to_merge = false;

        // for appending upstream, only the start coordinate actually matters
        if append_upstream {
            if coding && thin_start != thick_start {
                panic!("CRITICAL: Attempting to graft a coding block to a sequence with non-coding upstream fraction")
            }
            if !coding && graft_start > thick_start {
                println!("WARNING: Graft start coordinate lies within the coding sequence");
                return None;
            };
            // update the start coordinate(s)
            let updated_start: bool = graft_start < thin_start;
            // thin_start = graft_start;
            if coding {thick_start = graft_start};
            // append the graft to the first block
            // first, increase the size of the first block
            // exon_sizes[0] += graft.length().unwrap();
            let mut grafted = false;
            for i in 0..exon_sizes.len() {
                let exon_start = thin_start + exon_starts[i];
                let exon_end =  exon_start + exon_sizes[i];
                if exon_start <= graft_start && graft_start <= exon_end {
                    if allow_overlaps {to_merge = true} else {return None}
                }
                // println!("graft_len={}, graft_start={}, graft_end={}, exon_start={}, exon_end={}, thick_start={}, thick_end={}", graft_len, graft_start, graft_end, exon_start, exon_end, thick_start, thick_end);
                if exon_end > thick_start && !grafted {
                    // check if exon has a non-coding fraction
                    if exon_start < thick_start {
                        if graft_start < exon_start {
                            exon_sizes[i] += exon_start - min(graft_start, exon_start);
                            if i != 0 {exon_starts[i] = min(graft_start, exon_start) - thin_start;}
                            graft_len = exon_start - min(graft_start, exon_start);
                        }
                    } else {
                        // just tilt the exon start and update its size
                        exon_starts[i] = if i == 0 {0} else {
                            min(graft_start, thick_start) - thin_start
                        };
                        exon_sizes[i] += exon_start - graft_start;
                        graft_len = exon_start - min(graft_start, exon_start)
                    }
                    grafted = true;
                    // potentially, nothing will happen further; break the loop
                    if !updated_start {break}
                } else {
                    if updated_start {
                        exon_starts[i] += graft_len//thin_start - graft_start;
                    }
                }
                // exon_starts[i] += graft_len;
            }
            thin_start = min(thin_start, graft_start)
        } else if append_downstream {
        // the reverse is true for downstream appending
            if coding && thin_end != thick_end {
                panic!("CRITICAL: Attempting to graft a coding block to a sequence with non-coding downstream fraction")
            }
            if !coding && graft_end < thick_end {
                println!("WARNING: Graft end coordinate lies within the coding sequence");
                return None;
            };
            // update the start coordinate(s)
            if coding {thick_end = graft_end};
            // append the graft to the last block
            // for that, find the last coding block first
            for mut i in 0..exon_sizes.len() {
                i = exon_sizes.len() - i - 1;
                let exon_start = thin_start + exon_starts[i];
                let exon_end =  exon_start + exon_sizes[i];
                if exon_start <= graft_end && graft_end <= exon_end {
                    if allow_overlaps {to_merge = true} else {return None}
                }
                if exon_start < thick_end {
                    // first (last) coding exon caught
                    if exon_end > thick_end {
                        if graft_end > thin_end {
                            exon_sizes[i] += graft_end - max(thick_end, graft_start);
                        }
                    } else {
                        exon_sizes[i] += graft_end - thick_end;
                    }
                    // further exons will not be affected; feel free to break
                    break
                }
            }
            // exon start positions will not change in this case though
            thin_end = max(graft_end, thin_end);
        } else { // grafting to the exact position requires invoking the merge_multiple() function 
            to_merge = true
        }
        if to_merge{
            // if graft_start > thin_start {
            //     println!("Graft start coordinate lies within the coding sequence");
            //     return None;
            // };
            // if graft_end < thin_end {
            //     println!("Graft end coordinate lies within the coding sequence");
            //     return None;
            // };
            let mut blocks = self.to_blocks().unwrap();
            // println!("blocks={:#?}", blocks);
            blocks.push(BedEntry::from_interval(graft).unwrap());
            let unmerged_block_num = blocks.len();
            blocks.sort_by(
                |a, b| if a.start().unwrap() == b.start().unwrap() {
                    a.end().unwrap().cmp(&b.end().unwrap())
                } else {
                    a.start().unwrap().cmp(&b.start().unwrap())
                }
            );
            let merged_blocks = merge_multiple(&mut blocks);
            if merged_blocks.len() < unmerged_block_num as usize && !allow_overlaps {
                // println!("Grafted interval overlaps some of the existing blocks. Consider setting allow overlap to allow merging blocks");
                return None;
            }
            // println!("merged_blocks={:#?},\nmerged_blocks.len()={}", merged_blocks, merged_blocks.len());
            // println!("blocks.len()={}, merged_blocks.len()={}", blocks.len(), merged_blocks.len());
            exon_sizes.clear();
            exon_starts.clear();
            thin_start = min(thin_start, graft_start);
            thin_end = max(thin_end, graft_end);
            if coding {
                thick_start = min(thick_start, graft_start);
                thick_end = max(thick_end, graft_end);
            }
            for i in 0..merged_blocks.len() {
                let inter = &merged_blocks[i];
                let start = inter.start().unwrap();
                let end = inter.end().unwrap();
                exon_sizes.push(end - start);
                exon_starts.push(start - thin_start);
            }
            exon_num = merged_blocks.len() as u16;
        }

        if inplace{
            self.thin_start = Some(thin_start);
            self.thin_end = Some(thin_end);
            self.thick_start = Some(thick_start);
            self.thick_start = Some(thick_start);
            self.exon_num = Some(exon_num as u16);
            self.exon_sizes = Some(exon_sizes);
            self.exon_starts = Some(exon_starts);
            return None;
        }
        let mut grafted_bed = BedEntry::empty();
        grafted_bed.format = Some(12);
        grafted_bed.chrom = self.chrom.clone();
        grafted_bed.thin_start = Some(thin_start);
        grafted_bed.thin_end = Some(thin_end);
        grafted_bed.name = self.name.clone();
        grafted_bed.score = self.score.clone();
        grafted_bed.strand = self.strand;
        grafted_bed.thick_start = Some(thick_start);
        grafted_bed.thick_end = Some(thick_end);
        grafted_bed.rgb = self.rgb.clone();
        grafted_bed.exon_num = Some(exon_num);
        grafted_bed.exon_sizes = Some(exon_sizes);
        grafted_bed.exon_starts = Some(exon_starts);
        Some(grafted_bed)
    }
}

#[cfg(test)]
mod test_graft {
    use super::*;

    #[test]
    fn test_graft_upstream(){
        // let input = BedEntry::bed12(
        //     String::from("chr1"),
        //     53298978,
        //     53308962,
        //     String::from("XM_047446425.1#ORMDL1#78"),
        //     String::from("0"),
        //     true,
        //     53298978,
        //     53308962,
        //     String::from("0,0,100"),
        //     3,
        //     vec![174,152,136,],
        //     vec![0,6476,9848,]
        // );
        let mut input = parse_bed(
            String::from("chr1	53298978	53308962	XM_047446425.1#ORMDL1#78	0	+	53298978	53308962	0,0,100	3	174,152,136,	0,6476,9848,"),
            12,
            false
        ).unwrap();
        let graft1 = parse_bed(
            String::from("chr1	53297131	53298145	XM_047446425.1#ORMDL1#78	1	+"),
            6,
            false
        ).unwrap();
        println!("Adding graft1");
        let grafted = input.graft(
            graft1, 
            true, 
            true, 
            false, 
            false, 
            false, 
            false
        );
        let graft2 = parse_bed(
            String::from("chr1	53298971	53298978	XM_047446425.1#ORMDL1#78	2	+"),
            6,
            false
        ).unwrap();
        println!("Adding graft2");
        let grafted = input.graft(
            graft2, 
            true, 
            true, 
            false, 
            false, 
            true, 
            false
        );
        let graft3 = parse_bed(
            String::from("chr1	53308962	53310298	XM_047446425.1#ORMDL1#78	3	+"),
            6,
            false
        ).unwrap();
        println!("Adding graft3");
        let grafted = input.graft(
            graft3, 
            true, 
            true, 
            false, 
            false, 
            false, 
            true
        );
        println!("{:#?}", input);
        println!("{}", to_line(&input, 12).unwrap());
    }

    #[test]
    fn test_graft_both() {
        let input_line = String::from(
            "chr4	49489819	49503120	ENST00000259407.7#BAAT#20	0	-	49489819	49503120	0,0,100	3	594,203,466,	0,9816,12835,"
        );
        let mut input = parse_bed(input_line, 12, false).unwrap();
        let graft1 = parse_bed(
            String::from(
                "chr4	49472245	49489819	ENST00000259407.7#BAAT|0	0	-"
            ),
            6,
            false
        ).unwrap();
        let grafted = input.graft(
            graft1, 
            true, 
            true, 
            false, 
            false, 
            true, 
            false
        );
        println!("{}", to_line(&input, 12).unwrap());
        
        let graft2 = parse_bed(
            String::from("chr4	49503120	49503179	ENST00000259407.7#BAAT|1	0	-"),
            6,
            false
        ).unwrap();
        let grafted = input.graft(
            graft2, 
            true, 
            true, 
            false, 
            false, 
            false, 
            true
        );
        println!("{}", to_line(&input, 12).unwrap());

        let graft3 = parse_bed(
            String::from("chr4	49506738	49510808	ENST00000259407.7#BAAT|2	0	-"),
            6,
            false
        ).unwrap();
        let grafted = input.graft(
            graft3, 
            true, 
            true, 
            false, 
            false, 
            false, 
            false
        );
        println!("{}", to_line(&input, 12).unwrap());
    }


    #[test]
    fn test_graft_downstream() {
        let mut input = parse_bed(
            String::from("chr4	136609684	136613103	ENST00000566855.4#TEX46#5	0	+	136609684	136613103	0,0,200	3	2,160,210,	0,843,3209,"),
            12,
            false
        ).unwrap();
        let graft = parse_bed(
            String::from("chr4	136613095	136613132	ENST00000566855.4#TEX46|1	0	+"),
            6,
            false
        ).unwrap();
        let result = input.graft(
            graft,
            false,
            true,
            false,
            false,
            false,
            true
        ).unwrap();
        println!(
            "{}", to_line(&result, 12).unwrap()
        );
    }

    #[test]
    fn graft_both_with_overlaps() {
        let mut input = parse_bed(
            String::from("chr10	81321231	81325954	ENST00000248420.9#CACTIN#261	0	+	81321231	81325954	0,0,100	9	215,568,146,163,115,193,120,287,491,	0,1206,1889,2387,2678,3071,3637,3858,4232,"),
            12,
            false
        ).unwrap();
        let graft_up = parse_bed(
            String::from("chr10\t81321176\t81321231\tENST00000248420.9#CACTIN\t0\t+"),
            6,
            false
        ).unwrap();
        let grafted_up = input.graft(
            graft_up, 
            false, 
            true, 
            false, 
            false, 
            true, 
            false
        ).unwrap();
        println!(
            "{}", to_line(&grafted_up, 12).unwrap()
        );

        let graft_down1 = parse_bed(
            String::from("chr10\t81325913\t81326232\tENST00000248420.9#CACTIN|2\t0\t+"),
            6,
            false
        ).unwrap();
        let mut grafted_down1 = input.graft(
            graft_down1, 
            false, 
            true, 
            true, 
            false, 
            false, 
            false
        ).unwrap();
        println!(
            "{}", to_line(&grafted_down1, 12).unwrap()
        );

        let graft_down2 = parse_bed(
            String::from("chr10\t81325954\t81325973\tENST00000248420.9#CACTIN|1\t0\t+"), 
            6, 
            false
        ).unwrap();
        let grafted_down2 = grafted_down1.graft(
            graft_down2, 
            false, 
            true, 
            true, 
            false, 
            false, 
            true
        ).unwrap();
        println!(
            "{}", to_line(&grafted_down2, 12).unwrap()
        );
    }

    #[test]
    pub fn graft_problematic1() {
        let mut input = parse_bed(
            String::from("chr5	33379227	33414891	A	0	+	33379227	33414891	0,0,100	13	98,331,121,396,113,129,106,172,123,184,112,175,94,	0,8428,9488,10414,12516,12828,23266,29581,30199,31631,32333,34671,35570,"),
            12,
            false
        ).unwrap();
        let graft1 = parse_bed(
            String::from("chr5\t33378734\t33379277\t0\t0\t+"),
            6,
            false
        ).unwrap();
        let _ = input.graft(
            graft1,
            true,
            true, 
            true, 
            false, 
            false, 
            false
        );
        println!(
            "{}", to_line(&input, 12).unwrap()
        );
        let graft2 = parse_bed(
            String::from("chr5\t33379225\t33379227\t1\t0\t+"),
            6,
            false
        );
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum UtrSide {
    FivePrime,
    ThreePrime
}

#[derive(Clone, Debug)]
pub struct UtrBlock {
    chrom: Option<String>,
    start: Option<u64>,
    end: Option<u64>,
    name: Option<String>,
    strand: Option<bool>,
    side: Option<UtrSide>,
    adjacent: Option<bool>
}

impl UtrBlock {
    pub fn new() -> UtrBlock {
        UtrBlock {chrom: None, start: None, end: None, name: None, strand: None, side: None, adjacent: None}
    }

    pub fn from_bed(source: &BedEntry) -> UtrBlock {
        let mut result = UtrBlock::new();
        if let Some(x) = source.chrom() {
            result.chrom = Some(x.clone())
        }
        if let Some(x) = source.start() {
            result.start = Some(*x)
        }
        if let Some(x) = source.end() {
            result.end = Some(*x)
        }
        if let Some(x) = source.name() {
            result.name = Some(x.clone())
        }
        if let Some(x) = source.strand() {
            result.strand = Some(x)
        }
        result
    }

    pub fn set_side(&mut self, side: UtrSide) {
        self.side = Some(side)
    }

    pub fn set_adjacency(&mut self, is_adjacent: bool) {
        self.adjacent = Some(is_adjacent)
    }
}

pub trait Coordinates{
    fn chrom(&self) -> Option<&String>;

    fn start(&self) -> Option<&u64>;

    fn end(&self) -> Option<&u64>;

    fn reset_start(&mut self);

    fn reset_end(&mut self);

    fn length(&self) -> Option<u64>;
}

impl Coordinates for Interval {
    fn chrom(&self) -> Option<&String> {
        self.chrom.as_ref()
    }

    fn start(&self) -> Option<&u64> {
        self.start.as_ref()
    }

    fn end(&self) -> Option<&u64> {
        self.end.as_ref()
    }

    fn reset_start(&mut self) {
        self.start = None;
    }

    fn reset_end(&mut self) {
        self.end = None;
    }

    fn length(&self) -> Option<u64> {
        match (self.start, self.end) {
            (Some(a), Some(b)) => {b.checked_sub(a)},
            _ => None
        }
    }
}

impl<'a> Coordinates for  &'a Interval {
// impl<'a, T> Coordinates for T 
// where 
//     &'a T: Coordinates
// {
    fn chrom(&self) -> Option<&String> {
        self.chrom.as_ref()
    }

    fn start(&self) -> Option<&u64> {
        self.start.as_ref()
    }

    fn end(&self) -> Option<&u64> {
        self.end.as_ref()
    }

    fn reset_start(&mut self) {
        // self.start = None;
    }

    fn reset_end(&mut self) {
        // self.end = None;
    }

    fn length(&self) -> Option<u64> {
        match (self.start, self.end) {
            (Some(a), Some(b)) => {b.checked_sub(a)},
            _ => None
        }
    }
}

impl<'a> Coordinates for &'a BedEntry {
    fn chrom(&self) -> Option<&String> {
        self.chrom.as_ref()
    }

    fn start(&self) -> Option<&u64> {
        self.thin_start.as_ref()
    }

    fn end(&self) -> Option<&u64> {
        self.thin_end.as_ref()
    }

    fn reset_start(&mut self) {
        // self.start = None;
    }

    fn reset_end(&mut self) {
        // self.end = None;
    }

    fn length(&self) -> Option<u64> {
        match (self.thin_start, self.thin_end) {
            (Some(a), Some(b)) => {b.checked_sub(a)},
            _ => None
        }
    }
}

impl Coordinates for BedEntry {
    fn chrom(&self) -> Option<&String> {
        self.chrom.as_ref()
    }

    fn start(&self) -> Option<&u64> {
        self.thin_start.as_ref()
    }

    fn end(&self) -> Option<&u64> {
        self.thin_end.as_ref()
    }

    fn reset_start(&mut self) {
        // self.start = None;
    }

    fn reset_end(&mut self) {
        // self.end = None;
    }

    fn length(&self) -> Option<u64> {
        match (self.thin_start, self.thin_end) {
            (Some(a), Some(b)) => {b.checked_sub(a)},
            _ => None
        }
    }
}

impl<'a> Coordinates for  &'a UtrBlock {

    fn chrom(&self) -> Option<&String> {
        self.chrom.as_ref()
    }

    fn start(&self) -> Option<&u64> {
        self.start.as_ref()
    }

    fn end(&self) -> Option<&u64> {
        self.end.as_ref()
    }

    fn reset_start(&mut self) {
        //
    }

    fn reset_end(&mut self) {
        // self.end = None;
    }

    fn length(&self) -> Option<u64> {
        match (self.start, self.end) {
            (Some(a), Some(b)) => {b.checked_sub(a)},
            _ => None
        }
    }
}

impl Coordinates for UtrBlock {
    fn chrom(&self) -> Option<&String> {
        self.chrom.as_ref()
    }

    fn start(&self) -> Option<&u64> {
        self.start.as_ref()
    }

    fn end(&self) -> Option<&u64> {
        self.end.as_ref()
    }

    fn reset_start(&mut self) {
        self.start = None;
    }

    fn reset_end(&mut self) {
        self.end = None;
    }

    fn length(&self) -> Option<u64> {
        match (self.start, self.end) {
            (Some(a), Some(b)) => {b.checked_sub(a)},
            _ => None
        }
    }
}

pub trait Stranded {
    fn strand(&self) -> bool;

    fn update_strand(&mut self, strand: bool);
}

impl Stranded for UtrBlock {
    fn strand(&self) -> bool {
        self.strand.unwrap()
    }

    fn update_strand(&mut self, strand: bool) {
        self.strand = Some(strand)
    }
}

pub trait Named {
    fn name(&self) -> Option<&str>;

    fn update_name(&mut self, new_name: &str);
}

impl Named for Interval {
    fn name(&self) -> Option<&str> {
        match self.name.as_ref() {
            Some(x) => Some(x),
            None => None
        }
    }

    fn update_name(&mut self, new_name: &str ) {
        self.name = Some(new_name.to_string());
    }
}

impl<'a> Named for &'a Interval {
    fn name(&self) -> Option<&str> {
        match self.name.as_ref() {
            Some(x) => Some(x),
            None => None
        }
    }

    fn update_name(&mut self, new_name: &str ) {
        // self.name = Some(new_name.to_string());
    }
}

impl Named for BedEntry{
    fn name(&self) -> Option<&str> {
        match self.name.as_ref() {
            Some(x) => Some(x),
            None => None
        }
    }

    fn update_name(&mut self, new_name: &str ) {
        self.name = Some(new_name.to_string());
    }
}

impl<'a> Named for &'a BedEntry{
    fn name(&self) -> Option<&str> {
        match self.name.as_ref() {
            Some(x) => Some(x),
            None => None
        }
    }

    fn update_name(&mut self, new_name: &str ) {
        // self.name = Some(new_name.to_string());
    }
}

impl Named for UtrBlock{
    fn name(&self) -> Option<&str> {
        match self.name.as_ref() {
            Some(x) => Some(x),
            None => None
        }
    }

    fn update_name(&mut self, new_name: &str ) {
        self.name = Some(new_name.to_string());
    }
}

impl<'a> Named for &'a UtrBlock{
    fn name(&self) -> Option<&str> {
        match self.name.as_ref() {
            Some(x) => Some(x),
            None => None
        }
    }

    fn update_name(&mut self, new_name: &str ) {
        // self.name = Some(new_name.to_string());
    }
}