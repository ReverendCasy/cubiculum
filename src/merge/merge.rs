use anyhow::Error;
use fxhash::{FxHashMap, FxHashSet};
use num_traits::CheckedSub;
use std::cmp::{Ord, PartialOrd, min, max};
use std::ops::Sub;

use crate::structs::structs::{Coordinates,  Interval, Named};

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
where T: Ord + PartialOrd + Sub<Output = T> + CheckedSub<Output = T>
{
    
    let min_end: T = min(end1, end2); //if end1 > end2 { end2 } else { end1 };
    let max_start: T = max(start1, start2); //if start1 > start2 { start1 } else { start2 };
    min_end.checked_sub(&max_start)
}

/// ALPHA: Split interval into two by given gap coordinates
/// # Arguments
/// `block` - an interval object with defined coordinates (BedEntry, Interval, etc.)
/// `split_start` - start of the split interval 
/// `split_end` - end of the split interval
/// 
/// # Returns
/// A Vector of two value - left and right subintervals
/// 
/// # Example
/// ```
/// use cubiculum;
/// let a = Interval::from(None, Some(10), Some(20), None);
/// let split_start: u64 = 14;
/// let split_end: u64 = 16;
/// let subs = split_by_interval(a, split_start, split_end);
/// assert_eq!(subs.len(), 2);
/// assert_eq!(subs[0], Interval::from(None, Some(10), Some(14), None));
/// assert_eq!(subs[1], Interval::from(None, Some(16), Some(20), None));
/// ```
pub fn split_by_interval<T>(block: T, split_start: u64, split_end: u64) -> Vec<Interval>//Result<Vec<Interval>, Error>
where
    T: Coordinates {
    // 
    // let chrom = block.chrom().expect("No chromosome value recorded for the interval");
    let block_start = block.start().expect("No start value recorded for block start");
    let block_end = block.end().expect("No start value recorded for block end");
    let left = Interval::from(
        // Some(chrom.clone()),
        None, 
        Some(*block_start), 
        Some(split_start), 
        None,
    );
    let right = Interval::from(
        // Some(chrom.clone()),
        None,
        Some(split_end),
        Some(*block_end),
        None,
    );
    vec![left, right]
    // Ok(vec![left, right])s
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
                out_vec.push(out_interval);
            }
        };
    }
    out_vec.sort_by(
        |a, b| if a.start().unwrap() == b.start().unwrap() {
            a.end().unwrap().cmp(&b.end().unwrap())
        } else {
            a.start().unwrap().cmp(&b.start().unwrap())
        }
    );
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

/// split a vector of potentially overlapping intervals into discrete, non-overlapping ones,
/// and map the resulting intervals to names of original items overlapping the respective interval
/// 
/// 
pub fn discrete_interval_map<T>(intervals: &mut Vec<T>) -> (Vec<Interval>, FxHashMap<String, Vec<&str>>)
where 
    T: Coordinates + Named
{
    let mut interval_vec: Vec<Interval> = Vec::new();
    let mut out_map: FxHashMap<String, Vec<&str>> = FxHashMap::default();
    if intervals.len() == 0 {
        return (interval_vec, out_map);
    }
    intervals.sort_by(
        |a, b| if a.start().unwrap() == b.start().unwrap() {
            a.end().unwrap().cmp(&b.end().unwrap())
        } else {
            a.start().unwrap().cmp(&b.start().unwrap())
        }
    );
    
    let mut curr: usize = 0;
    let mut next: usize = 1;

    let mut start_points: Vec<u64> = Vec::new();
    let mut start2trs: FxHashMap<u64, Vec<&str>> = FxHashMap::default();
    let chrom: Option<String> = match intervals[0].chrom() {
        Some(x) => {Some(x.clone())},
        None => {None}
    };
    let mut curr_interval: u64 = 0;
    // let end2trs: FxHashMap<u64, Vec<&str>> = FxHashMap::default();
    while curr < intervals.len() {
        let first_start: u64 = match intervals[curr].start() {
            Some(x) => {*x},
            None => {
                panic!(
                    "Cannot discretize intervals with undefined coordinates; found an undefined start coordinate for interval {}", curr
                )
            }
        };
        let first_end: u64 = match intervals[curr].end() {
            Some(x) => {*x},
            None => {
                panic!(
                    "Cannot discretize intervals with undefined coordinates; found an undefined end coordinate for interval {}", curr
                )
            }
        };
        if let None = intervals[curr].name() {
            panic!("Cannot discretize unnamed intervals");
        }
        let mut curr_end = first_end;
        start_points.push(first_start);
        start_points.push(first_end);
        start2trs.entry(first_start).or_insert(Vec::new()).push(intervals[curr].name().unwrap());

        while next < intervals.len() {
            let next_start: u64 = match intervals[next].start() {
                Some(x) => {*x},
                None => {
                    panic!(
                        "Cannot discretize intervals with undefined coordinates; found an undefined start coordinate for interval {}", next
                    )
                }
            };
            let next_end: u64 = match intervals[next].end() {
                Some(x) => {*x},
                None => {
                    panic!(
                        "Cannot discretize intervals with undefined coordinates; found an undefined end coordinate for interval {}", next
                    )
                }
            };

            if let None = intervals[next].name() {
                panic!("Cannot discretize unnamed intervals");
            }

            if next_start > curr_end {
                break
            }

            curr_end = max(curr_end, next_end);
            if !start_points.contains(&next_start) {start_points.push(next_start)};
            if !start_points.contains(&next_end) {start_points.push(next_end)};
            for i in &start_points {
                if *i < next_start {continue};
                if *i >= next_end {continue};
                start2trs
                    .entry(*i)
                    .and_modify(|x| 
                        if x.contains(&intervals[next].name().unwrap()) {} else {x.push(intervals[next].name().unwrap())}
                    )
                    .or_insert(vec![intervals[next].name().unwrap()]);
            }

            // assess whether any of the previous intervals cover terminal coordinates for the current interval
            for i in curr..next+1 {
                // every interval that does not end before this point is attributed to this discrete interval
                let i_end: u64 = *intervals[i].end().unwrap();
                // println!("i_end={}, next_start={}, next_end={}", i_end, next_start, next_end);
                if i_end > next_start {
                    start2trs
                        .entry(next_start)
                        .and_modify(|x| 
                            if x.contains(&intervals[i].name().unwrap()) {} else {x.push(intervals[i].name().unwrap())}
                        )
                        .or_insert(vec![intervals[i].name().unwrap()]);
                };

                if i_end > next_end {
                    start2trs
                        .entry(next_end)
                        .and_modify(|x| 
                            if x.contains(&intervals[i].name().unwrap()) {} else {x.push(intervals[i].name().unwrap())}
                        )
                        .or_insert(vec![intervals[i].name().unwrap()]);
                };
            }
            next += 1;
        }
        // discretize the overlapping interval
        start_points.sort();
        for i in 1..start_points.len() {
            // define interval boundaries
            let inter_start: u64 = start_points[i-1];
            let inter_end: u64 = start_points[i];
            // define which transcripts correspond to this interval
            let tr_names: &Vec<&str>  = start2trs.get(&inter_start).unwrap_or_else(||
                {
                    println!("{:#?}", start2trs);
                    println!("{:#?}", start_points);
                    panic!("No transcripts overlapping this value: {}!", inter_start);
                }
            );
            // create an interval object and add the resulting values to the output collections
            let interval_name: String = curr_interval.to_string();
            out_map.insert(interval_name.clone(), tr_names.clone());
            let discrete_interval: Interval = Interval::from(
                chrom.clone(), Some(inter_start), Some(inter_end), Some(interval_name)
            );
            interval_vec.push(discrete_interval); 
            curr_interval += 1;
        }
        // clear the storage structures
        start_points.clear();
        start2trs.clear();
        // next iteration starts from the break point
        curr = next;
    }
    (interval_vec, out_map)
}

#[cfg(test)]
mod discretizer_test{
    use super::*;

    #[test]
    fn discretizer_identical(){
        let mut input: Vec<Interval> = vec![
            Interval::from(Some(String::from("chr1")), Some(100), Some(200), Some(String::from("one"))),
            Interval::from(Some(String::from("chr1")), Some(100), Some(200), Some(String::from("two")))
        ];
        let (vec, map) = discrete_interval_map(&mut input);
        println!("{:#?}", vec);
        println!("{:#?}", map);
    }

    #[test]
    fn discretizer_simple_overlap(){
        let mut input: Vec<Interval> = vec![
            Interval::from(Some(String::from("chr1")), Some(100), Some(200), Some(String::from("one"))),
            Interval::from(Some(String::from("chr1")), Some(150), Some(220), Some(String::from("two")))
        ];
        let (vec, map) = discrete_interval_map(&mut input);
        println!("{:#?}", vec);
        println!("{:#?}", map);
    }

    #[test]
    fn discretizer_nested_overlap(){
        let mut input: Vec<Interval> = vec![
            Interval::from(Some(String::from("chr1")), Some(100), Some(200), Some(String::from("one"))),
            Interval::from(Some(String::from("chr1")), Some(150), Some(180), Some(String::from("two")))
        ];
        let (vec, map) = discrete_interval_map(&mut input);
        println!("{:#?}", vec);
        println!("{:#?}", map);
    }

    #[test]
    fn discretizer_shared_start(){
        let mut input: Vec<Interval> = vec![
            Interval::from(Some(String::from("chr1")), Some(100), Some(200), Some(String::from("one"))),
            Interval::from(Some(String::from("chr1")), Some(100), Some(220), Some(String::from("two")))
        ];
        let (vec, map) = discrete_interval_map(&mut input);
        println!("{:#?}", vec);
        println!("{:#?}", map);
    }

    #[test]
    fn discretizer_three_intervals(){
        let mut input: Vec<Interval> = vec![
            Interval::from(Some(String::from("chr1")), Some(100), Some(200), Some(String::from("one"))),
            Interval::from(Some(String::from("chr1")), Some(100), Some(220), Some(String::from("two"))),
            Interval::from(Some(String::from("chr1")), Some(230), Some(250), Some(String::from("three")))
        ];
        let (vec, map) = discrete_interval_map(&mut input);
        println!("{:#?}", vec);
        println!("{:#?}", map);
    }

    #[test]
    fn real_life_test(){
        let mut input: Vec<Interval> = vec![
            Interval::from(Some(String::from("chr9")), Some(113042724), Some(113044268), Some(String::from("ENST00000374227.8#ZFP37_1"))),
            Interval::from(Some(String::from("chr9")), Some(113049361), Some(113049496), Some(String::from("ENST00000374227.8#ZFP37_2"))),
            Interval::from(Some(String::from("chr9")), Some(113049790), Some(113049872), Some(String::from("ENST00000374227.8#ZFP37_3"))),
            Interval::from(Some(String::from("chr9")), Some(113056556), Some(113056688), Some(String::from("ENST00000374227.8#ZFP37_4"))),
            Interval::from(Some(String::from("chr9")), Some(113042724), Some(113044268), Some(String::from("NM_001282515.2#ZFP37_1"))),
            Interval::from(Some(String::from("chr9")), Some(113049361), Some(113049496), Some(String::from("NM_001282515.2#ZFP37_2"))),
            Interval::from(Some(String::from("chr9")), Some(113049790), Some(113049917), Some(String::from("NM_001282515.2#ZFP37_3"))),
            Interval::from(Some(String::from("chr9")), Some(113056556), Some(113056688), Some(String::from("NM_001282515.2#ZFP37_4"))),
            Interval::from(Some(String::from("chr9")), Some(113042724), Some(113044268), Some(String::from("NM_001282518.2#ZFP37_1"))),
            Interval::from(Some(String::from("chr9")), Some(113049361), Some(113049496), Some(String::from("NM_001282518.2#ZFP37_2"))),
            Interval::from(Some(String::from("chr9")), Some(113049790), Some(113049875), Some(String::from("NM_001282518.2#ZFP37_3"))),
            Interval::from(Some(String::from("chr9")), Some(113056556), Some(113056688), Some(String::from("NM_001282518.2#ZFP37_4"))),
        ];
        let (vec, map) = discrete_interval_map(&mut input);
        println!("{:#?}", vec);
        println!("{:#?}", map);
    }
}


/// Given a Vector of Interval-like objects, 
/// returns a Vector of sorted unique intervals from the input collection
/// and interval mapping to names from the original collections
/// 
/// # Arguments
/// `intervals`: 
pub fn unique_interval_map<T>(intervals: &mut Vec<T>) -> (Vec<Interval>, FxHashMap<String, Vec<&str>>) 
where T: Coordinates + Named {
    let mut interval_vec: Vec<Interval> = Vec::new();
    let mut out_map: FxHashMap<String, Vec<&str>> = FxHashMap::default();
    let mut coord2name: FxHashMap<(&u64, &u64), String> = FxHashMap::default();

    // sort the input intervals
    intervals.sort_by(
        |a, b| if a.start().unwrap() == b.start().unwrap() {
            a.end().unwrap().cmp(&b.end().unwrap())
        } else {
            a.start().unwrap().cmp(&b.start().unwrap())
        }
    );

    // create a coordinate set to keep track of the already encountered coordinates

    let mut _counter: u16 = 1;
    for el in intervals {
        let name = el
            .name()
            .expect("Unnamed interval encountered in `unique_interval_map()`");
        let chrom = el
            .chrom()
            .expect(
                &format!("Interval {} has no defined chromosome", name)
        ); // do we really need this one?
        let start = el
            .start()
            .expect(
                &format!("Interval {} has no defined start coordinate", name)
        );
        let end = el
            .end()
            .expect(
                &format!("Interval {} has no defined end coordinate", name)
        );
        let coords = (start, end);
        if coord2name.contains_key(&coords) {
            // add the current interval's name to the map
            let new_inter_name = coord2name.get(&coords).unwrap();
            out_map
                .entry(new_inter_name.clone())
                .and_modify(|x| x.insert(x.len(), name));
        } else {
            // create a new interval
            let new_inter_name = format!("interval_{}", _counter);
            let new_inter = Interval::from(
                Some(chrom.clone()), 
                Some(*start), 
                Some(*end),
                Some(new_inter_name.clone())
            );
            interval_vec.insert(interval_vec.len(), new_inter);
            out_map.insert(new_inter_name.clone(), vec![name]);
            coord2name.insert((start, end), new_inter_name);
            _counter += 1;
        }

    }

    interval_vec.sort_by(
        |a, b| {
            if (a.start().unwrap() == b.start().unwrap()) {
                a.end().unwrap().cmp(b.end().unwrap())
            } else {
                a.start().unwrap().cmp(b.end().unwrap())
            }
        }
    );

    (interval_vec, out_map)
}

#[cfg(test)]
mod unique_interval_map_test {
    use fxhash::FxHashMap;
    use super::*;

    #[test]
    fn test1(){
        let mut input: Vec<Interval> = vec![
            Interval::from(Some(String::from("chr1")), Some(100), Some(200), Some(String::from("one"))),
            Interval::from(Some(String::from("chr1")), Some(100), Some(200), Some(String::from("two"))),
            Interval::from(Some(String::from("chr1")), Some(230), Some(250), Some(String::from("three")))
        ];
        let (vec, map) =unique_interval_map(&mut input);
        println!("{:#?}", vec);
        println!("{:#?}", map);
        let vec_to_compare = vec![
            Interval::from(Some(String::from("chr1")), Some(100), Some(200), Some(String::from("interval_1"))),
            Interval::from(Some(String::from("chr1")), Some(230), Some(250), Some(String::from("interval_2")))
        ];
        assert_eq!(
            vec, 
            vec_to_compare,
        );
        let mut map_to_compare: FxHashMap<String, Vec<&str>> = FxHashMap::default();
        map_to_compare.insert(String::from("interval_1"), vec!["one", "two"]);
        map_to_compare.insert(String::from("interval_2"), vec!["three"]);
        assert_eq!(
            map,
            map_to_compare,
        )
    }
}
