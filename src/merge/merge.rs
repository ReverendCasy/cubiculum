use num_traits::CheckedSub;
use std::cmp::{Ord, PartialOrd, min, max};

use crate::structs::{Coordinates, Named, Interval};

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
pub fn merge<T>(inter1: T, inter2: T) -> Option<Interval> 
where
    T: Coordinates
{
    let s1 = inter1.start().expect("Cannot merge intervals with undefined coordinates");
    let e1 = inter1.end().expect("Cannot merge intervals with undefined coordinates");
    let s2 = inter2.start().expect("Cannot merge intervals with undefined coordinates");
    let e2 = inter2.end().expect("Cannot merge intervals with undefined coordinates");
    if intersection(s1, e1, s2, e2) == 0 {return None};
    let merged: Interval = Interval::new();
    let merged_start = min(s1, s2);
    merged.update_start(merged_start);
    let merged_end = max(e1, e2);
    merged.update_end();
    merged
}