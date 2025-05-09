use std::cmp::{min, max};
use std::ops::Sub;

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
        let new_thin_start: u64 = match start {
            Some(x) => {max(self.thin_start.unwrap(), x)},
            None => self.thin_start.unwrap()
        };
        let new_thin_end: u64 = match end {
            Some(x) => {min(self.thin_end.unwrap(), x)},
            None => {self.thin_end.unwrap()}
        };
        let new_thick_start = match self.thick_start {
            Some(x) => {max(x, new_thin_start)},
            None => {new_thin_end}
        };
        let new_thick_end = match self.thick_end {
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
                    if ex_start > new_thin_end {continue};
                    let ex_end: u64 = y[i] + ex_start;
                    if ex_end < new_thin_start {continue};
                    let new_ex_start = max(new_thin_start, ex_start) - new_thin_start;
                    let new_ex_size = match min(ex_end, new_thin_end).checked_sub(new_ex_start){
                        Some(x) => {if x == 0 {continue} else {x}},
                        None => {continue} 
                    };

                    _sizes.push(new_ex_size);
                    _starts.push(new_ex_start);
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


}

pub trait Coordinates{
    fn chrom(&self) -> Option<&String>;

    fn start(&self) -> Option<&u64>;

    fn end(&self) -> Option<&u64>;

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

    fn length(&self) -> Option<u64> {
        match (self.start, self.end) {
            (Some(a), Some(b)) => {Some(a + b)},
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

    fn length(&self) -> Option<u64> {
        match (self.thin_start, self.thin_end) {
            (Some(a), Some(b)) => {Some(a + b)},
            _ => None
        }
    }
}

pub trait Named {
    fn name(&self) -> Option<&str>;
}

impl Named for Interval {
    fn name(&self) -> Option<&str> {
        match self.name.as_ref() {
            Some(x) => Some(x),
            None => None
        }
    }
}

impl Named for BedEntry{
    fn name(&self) -> Option<&str> {
        match self.name.as_ref() {
            Some(x) => Some(x),
            None => None
        }
    }
}