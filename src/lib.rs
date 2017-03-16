use std::ops::Range;
use std::cmp::Ordering;


extern crate cgmath;

type Point = cgmath::Point2<f64>;
type Vector = cgmath::Vector2<f64>;


/// Represent point id in delaunay data
#[derive(Debug, Clone, Copy, PartialEq)]
struct PointIndex(usize);
/// Represent triangle id in delaunay data
trait TrIndex{
    fn id(&self)->usize;
}

/// Represent any triangle index (further check required)
#[derive(Debug, Clone, Copy, PartialEq)]
struct TriangleIndex(usize);
impl TrIndex for TriangleIndex{
    #[inline]
    fn id(&self)->usize{
        self.0
    }
}

/// Represent ghost triangle index
#[derive(Debug, Clone, Copy, PartialEq)]
struct EdgeIndex(usize);
impl EdgeIndex{
    /// Create new EdgeIndex checking for correctness
    ///
    /// check appear only in debug build (no runtime cost in release)
    /// this check does not guarantee correctness, since referenced TriangleLike could be changed after
    /// index creation
    fn new(d: &Delaunay, id: TriangleIndex)->EdgeIndex{
        let edge = EdgeIndex(id.id());
        edge.test(d);
        edge
    }

    /// Check if EdgeIndex represent correct ghost triangle
    /// if not panic
    /// simply ommited in release
    #[inline]
    fn test(&self, d: &Delaunay){
        debug_assert!(d.tr(*self).points.2.is_none());
    }
}

impl TrIndex for EdgeIndex{
    #[inline]
    fn id(&self)->usize{
        self.0
    }
}

/// every EdgeIndex could be used as TriangleIndex
impl Into<TriangleIndex> for EdgeIndex{
    fn into(self)->TriangleIndex{
        TriangleIndex(self.id())
    }
}


/// Represent triangle in delaunay data
///
/// Contain 2 cases:
/// Regular triangle
///     points.3 is Some(PointIndex)
///     points kept in counterclockwise order (any points could be first)
///     neighbors[i] - triangle having same edge as opposite to points[i]
/// Ghost triangle
///     points.3 is None
///     represent hull edge
///     neighbors[0] - triangle on the opposite of edge
///     neighbors[1] - next counterclockwise edge on hull
///     neighbors[2] - previous counterclockwise edge on hull (next clockwise)
///     points kept in counterclockwise order on hull
#[derive(Debug, Clone, PartialEq)]
struct TriangleLike{
    neighbors: [TriangleIndex; 3],
    points: (PointIndex, PointIndex, Option<PointIndex>)
}
impl TriangleLike{
    /// for resize fashion
    fn default()->TriangleLike{
        TriangleLike{
            neighbors: [TriangleIndex(0), TriangleIndex(0), TriangleIndex(0)],
            points: (PointIndex(0), PointIndex(0), None)
        }
    }
}


/// Represent Delaunay triangulation
#[derive(Debug)]
pub struct Delaunay{
    points: Vec<Point>,
    triangles: Vec<TriangleLike>
}


impl Delaunay{
    /// Construct Delaunay triangulation from given set of points
    pub fn new(points: Vec<Point>)->Delaunay{
        let mut points = points;
        Delaunay::sort_points(&mut points);
        let len = points.len();
        let mut d = Delaunay{
            points: points, 
            triangles: Vec::with_capacity(len*2)
        };
        d.triangles.resize(len*2, TriangleLike::default());
        d.build(0..len);
        d
    }


    /// sort points in from left to right order
    ///
    /// if several points has similar x-coordinate, then sorting done by y-coordinate
    /// removes duplicate points
    fn sort_points(points: &mut Vec<Point>){
        points.sort_by(|a, b|{
            debug_assert!(
                a.x.is_finite()
                && a.y.is_finite()
                && b.x.is_finite()
                && b.y.is_finite()
                , "Delaunay do not support infinite and NaN point coordinates."
            );
            // since none coordinate is NaN unwrap is correct
            let x_cmp = a.x.partial_cmp(&b.x).unwrap();
            if x_cmp == Ordering::Equal{
                a.y.partial_cmp(&b.y).unwrap()
            }else{
                x_cmp
            }
        });
        points.dedup();
    }

    /// Construct triangulation on selected slice of points
    ///
    /// space for triangles must be preallocated
    /// triangles kept in range.start*2..range..(end*2-2)
    /// triangles (end*2-2)..end*2 reserved
    /// 
    /// TODO: leftmost and rightmost edges (edges containing leftmost point in .0 and rightmost in .1)
    /// stored as in reserved triangle end*2-2 as neighbors[0] and neighbors[1] correspondingly
    /// (merge)
    ///
    /// Divide and Conquer Gulbah-Stolfi algorithm
    fn build(&mut self, range: Range<usize>){
        match range.len(){
            2 => self.build_2points(range.start),
            3 => self.build_3points(range.start),
            _ => {
                let middle = (range.start + range.end) / 2;
                self.build(range.start..middle);
                self.build(middle..range.end);

                self.merge(range.start, middle, range.end);
            }
        }
    }

    /// Construct triangulation on 2 points
    ///
    /// consist of one edge convex hull
    /// through 'ghost' triangles
    fn build_2points(&mut self, start: usize){
        let i1 = TriangleIndex(start * 2); // index to t1
        let i2 = TriangleIndex(start * 2 + 1); // index to t2

        let t1 = TriangleLike{
            neighbors: [i2, i2, i2],
            points: (PointIndex(start), PointIndex(start+1), None)
        };
        let t2 = TriangleLike{
            neighbors: [i1, i1, i1],
            points: (PointIndex(start+1), PointIndex(start), None)
        };

        *self.tr_mut(i1) = t1;
        *self.tr_mut(i2) = t2;

        // store rightmost and left most edges
        let reversed_id = TriangleIndex(start*2 + 2); // same as end*2 - 2, since end == start+2
        self.tr_mut(reversed_id).neighbors[0] = i1; // t1 contain leftmost point as .0
        self.tr_mut(reversed_id).neighbors[1] = i1; // t1 contain rightmost point as .1
    }

    /// Construct triangulation on 3 points
    ///
    /// consist of triangle and
    /// hull through 'ghost' triangles
    ///     or
    /// of hull through 'ghost' triangles
    fn build_3points(&mut self, start: usize){
        let p_ids = [PointIndex(start), PointIndex(start+1), PointIndex(start+2)];
        let tr_ids = [TriangleIndex(start*2), TriangleIndex(start*2 + 1), TriangleIndex(start*2 + 2), TriangleIndex(start*2 + 3)];
        if self.is_on_line_index(&p_ids){
            // 4 ghost triangles by 2 edges
            *self.tr_mut(tr_ids[0]) = TriangleLike{
                points: (p_ids[0], p_ids[1], None),
                neighbors: [tr_ids[3], tr_ids[1], tr_ids[3]]
            };
            *self.tr_mut(tr_ids[1]) = TriangleLike{
                points: (p_ids[1], p_ids[2], None),
                neighbors: [tr_ids[2], tr_ids[2], tr_ids[0]]
            };
            *self.tr_mut(tr_ids[2]) = TriangleLike{
                points: (p_ids[2], p_ids[1], None),
                neighbors: [tr_ids[1], tr_ids[3], tr_ids[1]]
            };
            *self.tr_mut(tr_ids[3]) = TriangleLike{
                points: (p_ids[1], p_ids[0], None),
                neighbors: [tr_ids[0], tr_ids[0], tr_ids[2]]
            };

            // store rightmost and left most edges
            let reversed_id = TriangleIndex(start*2 + 4); // same as end*2 - 2, since end == start+3
            self.tr_mut(reversed_id).neighbors[0] = tr_ids[0]; // tr_ids[0] contain leftmost point as .0
            self.tr_mut(reversed_id).neighbors[1] = tr_ids[1]; // tr_ids[1] contain rightmost point as .1
        }else{
            // real triangle and 3 ghost triangles by its edges
            
            // make points appear in counterclockwise order
            let p_ids = if self.is_counterclockwise_index(&p_ids){
                p_ids
            }else{
                [p_ids[0], p_ids[2], p_ids[1]]
            };

            *self.tr_mut(tr_ids[0]) = TriangleLike{
                points: (p_ids[0], p_ids[1], Some(p_ids[2])),
                neighbors: [tr_ids[1], tr_ids[2], tr_ids[3]]
            };
            *self.tr_mut(tr_ids[1]) = TriangleLike{
                points: (p_ids[1], p_ids[2], None),
                neighbors: [tr_ids[0], tr_ids[2], tr_ids[3]]
            };
            *self.tr_mut(tr_ids[2]) = TriangleLike{
                points: (p_ids[2], p_ids[0], None),
                neighbors: [tr_ids[0], tr_ids[3], tr_ids[1]]
            };
            *self.tr_mut(tr_ids[3]) = TriangleLike{
                points: (p_ids[0], p_ids[1], None),
                neighbors: [tr_ids[0], tr_ids[1], tr_ids[2]]
            };

            // store rightmost and left most edges
            let reversed_id = TriangleIndex(start*2 + 4); // same as end*2 - 2, since end == start+3
            self.tr_mut(reversed_id).neighbors[0] = tr_ids[3]; // tr_ids[3] contain leftmost point as .0
            if self.p(p_ids[1]).x <= self.p(p_ids[2]).x{
                self.tr_mut(reversed_id).neighbors[1] = tr_ids[1]; // tr_ids[2] contain leftmost point as .1
            }else{
                self.tr_mut(reversed_id).neighbors[1] = tr_ids[3]; // tr_ids[1] contain leftmost point as .1
            }
        }
    }

    /// Check if points at given indexes lies on line
    fn is_on_line_index(&mut self, ids: &[PointIndex; 3])->bool{
        is_on_line(self.p(ids[0]), self.p(ids[1]), self.p(ids[2]))
    }

    /// Check if points at given indexes lies in counterclockwise order
    fn is_counterclockwise_index(&mut self, ids: &[PointIndex; 3])->bool{
        is_counterclockwise(self.p(ids[0]), self.p(ids[1]), self.p(ids[2]))
    }

    /// Get mutable reference on triangle by index
    fn tr_mut<T: TrIndex>(&mut self, id: T)->&mut TriangleLike{
        &mut self.triangles[id.id()]
    }

    /// Get immutable reference on triangle by index
    fn tr<T: TrIndex>(&self, id: T)->&TriangleLike{
        &self.triangles[id.id()]
    }

    /// Get immutable reference on Point by index
    fn p(&self, id: PointIndex)->&Point{
        &self.points[id.0]
    }

    /// merge to partial triangulation
    ///
    /// left triangulation consist from points in from..sep
    /// right from points in sep..to
    fn merge(&mut self, from: usize, sep: usize, to: usize){
        let left = self.find_rightmost_edge(from..sep);
        let right = self.find_leftmost_edge(sep..to);

        let (left, right) = self.find_lower_tangent(left, right);

        // use reserved space for lower tangent and upper tangent edges
        let lower_tangent_id = TriangleIndex(sep*2);
        let upper_tangent_id = TriangleIndex(sep*2 + 1); // for time of merge will be for moving edge there new triangles will appear

        let lower_tangent = TriangleLike{
            points: (self.tr(left).points.0, self.tr(right).points.1, None),
            neighbors: [upper_tangent_id, self.edge_clockwise(right).into(), self.edge_counterclockwise(left).into()]
        };

        let mut merge_edge = TriangleLike{
            points: (self.tr(right).points.1, self.tr(left).points.0, None),
            neighbors: [lower_tangent_id, left.into(), right.into()]
        };

        loop{
            while self.try_merge_right(&mut merge_edge){}

            let mut done = true;
            while self.try_merge_left(&mut merge_edge){
                done = false;
            }

            if done{
                break;
            }
        }

        *self.tr_mut(lower_tangent_id) = lower_tangent;
        *self.tr_mut(upper_tangent_id) = merge_edge; // since we added all we could merge_edge became upper_tangent
    }

    /// find index of ghost triangle, containing rightmost point at index 0
    fn find_rightmost_edge(&self, range: Range<usize>)->EdgeIndex{
        // see build() invariants to details
        EdgeIndex::new(self, self.triangles[range.end*2 - 2].neighbors[1])
    }

    /// find index of ghost triangle, containing leftmost point at index 1
    fn find_leftmost_edge(&self, range: Range<usize>)->EdgeIndex{
        // see build() invariants to details
        EdgeIndex::new(self, self.triangles[range.end*2 - 2].neighbors[0])
    }

    /// find lower tangent between two convex polygones
    ///
    /// left and right - 'ghost' triangles, containting points with clear line sight
    /// between them. left containt this point at .0, right at .1
    ///
    /// result is left and right 'ghost' triangles, containing left point in .0 and right in .1 correspondingly.
    fn find_lower_tangent(&self, left: EdgeIndex, right: EdgeIndex)->(EdgeIndex, EdgeIndex){
        // TODO
        // TODO DO DO
        (EdgeIndex(0), EdgeIndex(0))
    }

    /// find next counterclockwise ghost triangle on hull
    fn edge_counterclockwise(&self, id: EdgeIndex)->EdgeIndex{
        EdgeIndex::new(self, self.tr(id).neighbors[1])
    }

    /// find next clockwise ghost triangle on hull
    fn edge_clockwise(&self, id: EdgeIndex)->EdgeIndex{
        EdgeIndex::new(self, self.tr(id).neighbors[1])
    }

    /// Trying to merge on given edge with some point from right triangulation
    ///
    /// if merge appear return true, otherwise false
    /// if merge appear will move change merge_edge, to be upper edge of created triangle
    fn try_merge_right(&mut self, merge_edge: &mut TriangleLike)->bool{
        // TODO
        false
    }

    /// Trying to merge on given edge with some point from left triangulation
    ///
    /// if merge appear return true, otherwise false
    /// if merge appear will move change merge_edge, to be upper edge of created triangle
    fn try_merge_left(&mut self, merge_edge: &mut TriangleLike)->bool{
        // TODO
        false
    }
}

/// check if 3 point lies on one line
fn is_on_line(a: &Point, b: &Point, c: &Point)->bool{
    cross_product(&(b-a), &(c-a)) == 0.0
}

/// check if 3 points counterclockwise
fn is_counterclockwise(a: &Point, b: &Point, c: &Point)->bool{
    cross_product(&(b-a), &(c-a)) > 0.0
}

/// calculate z-component of cross product of vectors extended with z=0
fn cross_product(a: &Vector, b: &Vector) -> f64{
    a.x * b.y - a.y * b.x
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sort_points(){
        let mut points = vec![
            Point::new(1.0, 2.0),
            Point::new(0.0, 2.0),
            Point::new(3.0, 6.0),
            Point::new(0.0, 0.0),
            Point::new(-1.0, 2.0),
            Point::new(0.0, 2.0),
        ];
        let expected = vec![
            Point::new(-1.0, 2.0),
            Point::new(0.0, 0.0),
            Point::new(0.0, 2.0),
            Point::new(1.0, 2.0),
            Point::new(3.0, 6.0),
        ];

        Delaunay::sort_points(&mut points);
        assert_eq!(points, expected);
    }

    #[test]
    fn test_build_2point(){
        let points = vec![
            Point::new(0.0, 2.0),
            Point::new(1.0, 2.0),
            Point::new(3.0, 2.0),
            Point::new(4.0, 2.0),
        ];

        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(8, TriangleLike::default());
        d.build_2points(0);
        d.build_2points(2);

        let expected1 = [
            TriangleLike{
                neighbors: [TriangleIndex(1), TriangleIndex(1), TriangleIndex(1)], 
                points: (PointIndex(0), PointIndex(1), None)
            },
            TriangleLike{
                neighbors: [TriangleIndex(0), TriangleIndex(0), TriangleIndex(0)], 
                points: (PointIndex(1), PointIndex(0), None)
            },
        ];

        let expected2 = [
            TriangleLike{
                neighbors: [TriangleIndex(5), TriangleIndex(5), TriangleIndex(5)], 
                points: (PointIndex(2), PointIndex(3), None)
            },
            TriangleLike{
                neighbors: [TriangleIndex(4), TriangleIndex(4), TriangleIndex(4)], 
                points: (PointIndex(3), PointIndex(2), None)
            },
        ];

        assert_eq!(&d.triangles[0..2], &expected1);
        assert_eq!(&d.triangles[4..6], &expected2);
    }
    
    #[test]
    fn test_is_on_line(){
        let p1 = Point::new(1.0, 1.0);
        let p2 = Point::new(2.0, 3.0);
        let p3 = Point::new(3.0, 5.0);
        assert!(is_on_line(&p1, &p2, &p3));
        let p1 = Point::new(1.0, 1.0);
        let p2 = Point::new(2.0, 3.0);
        let p3 = Point::new(3.0, 3.0);
        assert!(!is_on_line(&p1, &p2, &p3));
    }

    #[test]
    fn test_is_counterclockwise(){
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(1.0, 0.0);
        let p3 = Point::new(0.0, 1.0);
        assert!(is_counterclockwise(&p1, &p2, &p3));
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(0.0, 1.0);
        let p3 = Point::new(1.0, 0.0);
        assert!(!is_counterclockwise(&p1, &p2, &p3));
    }

    #[test]
    fn test_build_3points(){
        // triangles test
        let points = vec![
            // start in clockwise
            Point::new(0.0, 0.0),
            Point::new(0.0, 1.0),
            Point::new(1.0, 0.0),

            // start in counterclockwise
            Point::new(1.0, 1.0),
            Point::new(2.0, 0.0),
            Point::new(3.0, 1.0)
        ];

        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(12, TriangleLike::default());
        d.build_3points(0);
        d.build_3points(3);

        let expected1 = &[
            TriangleLike{
                points: (PointIndex(0), PointIndex(2), Some(PointIndex(1))),
                neighbors: [TriangleIndex(1), TriangleIndex(2), TriangleIndex(3)]
            },
            TriangleLike{
                points: (PointIndex(2), PointIndex(1), None),
                neighbors: [TriangleIndex(0), TriangleIndex(2), TriangleIndex(3)]
            },
            TriangleLike{
                points: (PointIndex(1), PointIndex(0), None),
                neighbors: [TriangleIndex(0), TriangleIndex(3), TriangleIndex(1)]
            },
            TriangleLike{
                points: (PointIndex(0), PointIndex(2), None),
                neighbors: [TriangleIndex(0), TriangleIndex(1), TriangleIndex(2)]
            }
        ];

        let expected2 = &[
            TriangleLike{
                points: (PointIndex(3), PointIndex(4), Some(PointIndex(5))),
                neighbors: [TriangleIndex(7), TriangleIndex(8), TriangleIndex(9)]
            },
            TriangleLike{
                points: (PointIndex(4), PointIndex(5), None),
                neighbors: [TriangleIndex(6), TriangleIndex(8), TriangleIndex(9)]
            },
            TriangleLike{
                points: (PointIndex(5), PointIndex(3), None),
                neighbors: [TriangleIndex(6), TriangleIndex(9), TriangleIndex(7)]
            },
            TriangleLike{
                points: (PointIndex(3), PointIndex(4), None),
                neighbors: [TriangleIndex(6), TriangleIndex(7), TriangleIndex(8)]
            },
        ];

        assert_eq!(&d.triangles[0..4], expected1);
        assert_eq!(&d.triangles[6..10], expected2);

        // test 3 points on line
        let points = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 2.0),
            Point::new(3.0, 6.0),

            Point::new(4.0, 0.0),
            Point::new(5.0, -1.0),
            Point::new(6.0, -2.0),
        ];

        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(12, TriangleLike::default());
        d.build_3points(0);
        d.build_3points(3);

        let expected1 = &[
            TriangleLike{
                points: (PointIndex(0), PointIndex(1), None),
                neighbors: [TriangleIndex(3), TriangleIndex(1), TriangleIndex(3)]
            },
            TriangleLike{
                points: (PointIndex(1), PointIndex(2), None),
                neighbors: [TriangleIndex(2), TriangleIndex(2), TriangleIndex(0)]
            },
            TriangleLike{
                points: (PointIndex(2), PointIndex(1), None),
                neighbors: [TriangleIndex(1), TriangleIndex(3), TriangleIndex(1)]
            },
            TriangleLike{
                points: (PointIndex(1), PointIndex(0), None),
                neighbors: [TriangleIndex(0), TriangleIndex(0), TriangleIndex(2)]
            }
        ];

        let expected2 = &[
            TriangleLike{
                points: (PointIndex(3), PointIndex(4), None),
                neighbors: [TriangleIndex(9), TriangleIndex(7), TriangleIndex(9)]
            },
            TriangleLike{
                points: (PointIndex(4), PointIndex(5), None),
                neighbors: [TriangleIndex(8), TriangleIndex(8), TriangleIndex(6)]
            },
            TriangleLike{
                points: (PointIndex(5), PointIndex(4), None),
                neighbors: [TriangleIndex(7), TriangleIndex(9), TriangleIndex(7)]
            },
            TriangleLike{
                points: (PointIndex(4), PointIndex(3), None),
                neighbors: [TriangleIndex(6), TriangleIndex(6), TriangleIndex(8)]
            }
        ];
        assert_eq!(&d.triangles[0..4], expected1);
        assert_eq!(&d.triangles[6..10], expected2);
    }

    #[test]
    fn test_find_leftmost_rightmost_invariant(){
        // test for build2
        let points = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(1.0, 1.0),
            Point::new(1.0, 2.0)
        ];

        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(8, TriangleLike::default());
        d.build_2points(0);
        d.build_2points(2);

        assert_eq!(d.triangles[2].neighbors[0], TriangleIndex(0));
        assert_eq!(d.triangles[2].neighbors[1], TriangleIndex(0));

        assert_eq!(d.triangles[6].neighbors[0], TriangleIndex(4));
        assert_eq!(d.triangles[6].neighbors[1], TriangleIndex(4));

        // test for build3 triangle
        let points = vec![
            Point::new(0.0, 0.0),
            Point::new(0.0, 1.0),
            Point::new(1.0, 0.0),

            Point::new(1.0, 1.0),
            Point::new(2.0, 0.0),
            Point::new(2.0, 1.0)
        ];
        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(12, TriangleLike::default());
        d.build_3points(0);
        d.build_3points(3);

        assert_eq!(d.triangles[4].neighbors[0], TriangleIndex(3));
        assert_eq!(d.triangles[4].neighbors[1], TriangleIndex(3));

        assert_eq!(d.triangles[10].neighbors[0], TriangleIndex(9));
        assert_eq!(d.triangles[10].neighbors[1], TriangleIndex(7));

        // TODO: test following invariant for merge once it finished
    }
}
