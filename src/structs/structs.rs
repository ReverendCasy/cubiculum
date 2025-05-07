/// Contains data on storage structures for annotation manipulations in Cubiculum and associated packages

#[derive(Debug)]
pub struct Interval {
    chrom: Option<String>,
    start: Option<u32>,
    end: Option<u32>,
    name: Option<String>
}

impl Interval {
    pub fn new() -> Interval {
        Interval { chrom: None, start: None, end: None, name: None }
    }

    pub fn update_name(&mut self, name: String) {
        self.name = Some(name);
    }

    pub fn update_chrom(&mut self, chrom: String) {
        self.chrom = Some(chrom);
    }

    pub fn update_start(&mut self, start: u32) {
        self.start = Some(start);
    }

    pub fn update_end(&mut self, end: u32) {
        self.end = Some(end);
    }
}
#[derive(Debug)]
pub struct BedEntry{
    format: Option<u8>,
    chrom: Option<String>,
    thin_start: Option<u32>,
    thin_end: Option<u32>,
    name: Option<String>,
    score: Option<String>,
    strand: Option<bool>,
    thick_start: Option<u32>,
    thick_end: Option<u32>,
    rgb: Option<String>,
    exon_num: Option<u32>,
    exon_sizes: Option<Vec<u32>>,
    exon_starts: Option<Vec<u32>>
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

    pub fn bed3(chrom: String, start: u32, end: u32) -> BedEntry {
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

    pub fn bed4(chrom: String, start: u32, end: u32, name: String) -> BedEntry {
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

    pub fn bed5(chrom: String, start: u32, end: u32, name: String, score: String) -> BedEntry {
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

    pub fn bed6(chrom: String, start: u32, end: u32, name: String, score: String, strand: bool) -> BedEntry {
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
        chrom: String, start: u32, end: u32, name: String, score: String, strand: bool, 
        thick_start: u32, thick_end: u32 
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
        chrom: String, start: u32, end: u32, name: String, score: String, strand: bool, 
        thick_start: u32, thick_end: u32, rgb: String
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
        chrom: String, start: u32, end: u32, name: String, score: String, strand: bool, 
        thick_start: u32, thick_end: u32, rgb: String, 
        exon_num: u32, exon_sizes: Vec<u32>, exon_starts: Vec<u32>
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

    pub fn to_blocks(&self) -> Option<Vec<BedEntry>> {
        if self.format.unwrap() != 12 {
            return None;
        }
        let chrom: &str = match &self.chrom {
            Some(x) => {x},
            None => {return None}
        };
        let thin_start: u32 = match self.thin_start {
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
            let start: u32 = thin_start + self.exon_starts.as_ref().unwrap()[i];
            let end: u32 = start + self.exon_sizes.as_ref().unwrap()[i];
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


}

pub trait Coordinates{
    fn start(&self) -> Option<&u32>;

    fn end(&self) -> Option<&u32>;

    fn length(&self) -> Option<u32>;
}

impl Coordinates for Interval {
    fn start(&self) -> Option<&u32> {
        self.start.as_ref()
    }

    fn end(&self) -> Option<&u32> {
        self.end.as_ref()
    }

    fn length(&self) -> Option<u32> {
        match (self.start, self.end) {
            (Some(a), Some(b)) => {Some(a + b)},
            _ => None
        }
    }
}

impl Coordinates for BedEntry {
    fn start(&self) -> Option<&u32> {
        self.thin_start.as_ref()
    }

    fn end(&self) -> Option<&u32> {
        self.thin_end.as_ref()
    }

    fn length(&self) -> Option<u32> {
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