/// Contains data on storage structures for annotation manipulations in Cubiculum and associated packages

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
    exon_num: Option<u32>,
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
            format: Some(3), 
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
}

pub trait Coordinates{
    fn start(&self) -> Option<&u64>;

    fn end(&self) -> Option<&u64>;

    fn length(&self) -> Option<u64>;
}

impl Coordinates for Interval {
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