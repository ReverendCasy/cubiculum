use num_traits::CheckedSub;
use std::cmp::{Ord, PartialOrd, min, max};
use std::ops::Sub;

use crate::structs::structs::{Coordinates,  Interval};

/// Assess intersection between the two numeric intervals
/// 
/// # Arguments
/// `start1` = start coordinate of the first interval;
/// `end1` - end coordinate of the first interval;
/// `start2` - start coordinate of the first interval;
/// `end2` - start coordinate of the first interval;
/// 
/// # Returns
/// Numeric value of the intersection size if the intervals intersect, zero otherwise
/// 
/// # Usage
/// 
/// ```
/// use cubiculum;
/// let s1: u8 = 10;
/// let e1: u8 = 20;
/// let s2: u8 = 15;
/// let e2: u8 = 23;
/// let inter_size: u8 = intersection(s1, e1, s2, e2);
/// assert_eq!(inter_size, 5);
/// ```
pub fn intersection<T>(
    start1: T, end1: T, start2: T, end2: T
) -> Option<T> 
where T: Ord + PartialOrd + Sub<Output = T> + CheckedSub<Output = T>//<T: cmp::PartialOrd + ops::Sub<Output = T>>
{
    
    let min_end: T = min(end1, end2); //if end1 > end2 { end2 } else { end1 };
    let max_start: T = max(start1, start2); //if start1 > start2 { start1 } else { start2 };
    min_end.checked_sub(&max_start)
}


/// Merge two Coordinates objects into a single Interval object 
/// 
/// # Arguments
/// `inter1` - the first Coordinates object 
/// `inter2` - the second Cordinates object
/// 
/// # Returns
/// An Option containing the merged interval if the objects overlap, None otherwise
/// 
/// # Usage
/// ```
/// use cubiculum::merge::merge;
/// let inter1 = Interval::from(String::from("chr1"), 100, 200, String::from("inter1"));
/// let inter2 = Interval::from(String::from("chr1"), 170, 300, String::from("inter1"));
/// let merged = merge(inter1, inter2);
/// assert_eq!(merged, Interval::from(None, 100, 300, None));
/// ```
pub fn merge<T>(inter1: T, inter2: T) -> Option<Interval> 
where
    T: Coordinates
{
    let s1 = *inter1.start().expect("Cannot merge intervals with undefined coordinates");
    let e1 = *inter1.end().expect("Cannot merge intervals with undefined coordinates");
    let s2 = *inter2.start().expect("Cannot merge intervals with undefined coordinates");
    let e2 = *inter2.end().expect("Cannot merge intervals with undefined coordinates");
    match intersection(s1, e1, s2, e2) {
        None => {return None},
        Some(_) => {
            let mut merged: Interval = Interval::new();
            let merged_start = min(s1, s2);
            merged.update_start(merged_start);
            let merged_end = max(e1, e2);
            merged.update_end(merged_end);
            return Some(merged);
        }
    };
}


// merge all the overlapping intervals in the vector
pub fn merge_multiple<T>(intervals: &mut Vec<T>) -> Vec<Interval> 
where 
    T: Coordinates
{
    let mut out_vec: Vec<Interval> = Vec::new();
    if intervals.len() == 0 {return out_vec}
    let mut prev_start: u64 = 0;
    let mut prev_end: u64 = 0;
    for el in intervals {
        let curr_start = *el.start().unwrap();
        let curr_end = *el.end().unwrap();
        match intersection(prev_start, prev_end, curr_start, curr_end) {
            Some(_) => {
                // current item intersects the last interval in the output vector;
                // create a single intersecting item out of them 
                let _ = out_vec.pop();
                prev_start = min(prev_start, curr_start);
                prev_end = max(prev_end, curr_end);
                let mut merged: Interval = Interval::new();
                merged.update_chrom(el.chrom().unwrap().clone());
                merged.update_start(prev_start);
                merged.update_end(prev_end);
                out_vec.push(merged);
            },
            None => {
                // no intersection to the previous item; create a new interval, add it to the output vector
                prev_start = curr_start;
                prev_end = curr_end;
                // since the output value is the vector of Intervals, create an Interval decoy for this element
                let mut out_interval = Interval::new();
                out_interval.update_chrom(el.chrom().unwrap().clone());
                out_interval.update_start(prev_start);
                out_interval.update_end(prev_end);
            }
        };
    }
    out_vec
}

/// create an interval spanning over all the Coordinates objects in the vector
///
/// # Arguments
/// `intervals`: Vec collection containing the intervals
/// 
/// # Returns
pub fn total_span<T>(intervals: &mut Vec<T>) -> Interval 
where 
    T: Coordinates
{
    intervals.sort_by(
        |a, b| if a.start().unwrap() == b.start().unwrap() {
            a.end().unwrap().cmp(&b.end().unwrap())
        } else {
            a.start().unwrap().cmp(&b.start().unwrap())
        }
        );
    let chrom: String = intervals[0]
        .chrom()
        .expect("Intervals for total span inference must have a defined")
        .clone();
    let start: u64 = *intervals[0].start().unwrap();
    let end: u64 = *intervals[intervals.len() - 1].end().unwrap();
    let name: String = String::from(format!("{}:{}-{}", chrom, start, end));
    Interval::from(Some(chrom), Some(start), Some(end), Some(name))
}
