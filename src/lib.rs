use std::ops::Range;
use std::cmp::Ordering;


extern crate cgmath;
use std::fmt;
use cgmath::prelude::*;

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
        self.tr_mut(reversed_id).neighbors[0] = i2; // t2 contain leftmost point as .1
        self.tr_mut(reversed_id).neighbors[1] = i2; // t2 contain rightmost point as .0
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
        if self.is_on_line_index(p_ids[0], p_ids[1], p_ids[2]){
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
            self.tr_mut(reversed_id).neighbors[0] = tr_ids[3]; // tr_ids[3] contain leftmost point as .1
            self.tr_mut(reversed_id).neighbors[1] = tr_ids[2]; // tr_ids[2] contain rightmost point as .0
        }else{
            // real triangle and 3 ghost triangles by its edges
            
            // make points appear in counterclockwise order
            let p_ids = if self.is_counterclockwise_index(p_ids[0], p_ids[1], p_ids[2]){
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
            self.tr_mut(reversed_id).neighbors[0] = tr_ids[2]; // tr_ids[2] contain leftmost point as .1
            if self.p(p_ids[1]).x <= self.p(p_ids[2]).x{
                self.tr_mut(reversed_id).neighbors[1] = tr_ids[2]; // tr_ids[2] contain rightmost point as .0
            }else{
                self.tr_mut(reversed_id).neighbors[1] = tr_ids[1]; // tr_ids[1] contain rightmost point as .0
            }
        }
    }

    /// Check if points at given indexes lies on line
    fn is_on_line_index(&self, a: PointIndex, b: PointIndex, c: PointIndex)->bool{
        is_on_line(self.p(a), self.p(b), self.p(c))
    }

    /// Check if points at given indexes lies in counterclockwise order
    fn is_counterclockwise_index(&self, a: PointIndex, b: PointIndex, c: PointIndex)->bool{
        is_counterclockwise(self.p(a), self.p(b), self.p(c))
    }

    /// Check if points at given indexes lies in counterclockwise order
    fn is_clockwise_index(&self, a: PointIndex, b: PointIndex, c: PointIndex)->bool{
        is_clockwise(self.p(a), self.p(b), self.p(c))
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

        println!("left: {:?}", self.tr(left));
        println!("right: {:?}", self.tr(right));

        // save new leftmost edge
        self.triangles[to*2 - 2].neighbors[0] = self.triangles[sep*2 - 2].neighbors[0];

        let (left, right) = self.find_lower_tangent(left, right);

        // use reserved space for lower tangent and upper tangent edges
        let lower_tangent_id = EdgeIndex(sep*2 - 2);
        let merge_edge_id = EdgeIndex(sep*2 - 1); // for time of merge will be for moving edge there new triangles will appear

        let lower_tangent = TriangleLike{
            points: (self.tr(left).points.1, self.tr(right).points.0, None),
            neighbors: [merge_edge_id.into(), right.into(), left.into()]
        };

        let merge_edge = TriangleLike{
            points: (self.tr(right).points.0, self.tr(left).points.1, None),
            neighbors: [lower_tangent_id.into(), self.edge_counterclockwise(left).into(), self.edge_clockwise(right).into()]
        };

        *self.tr_mut(lower_tangent_id) = lower_tangent;
        *self.tr_mut(merge_edge_id) = merge_edge; // once we merged all we could, merge_edge will became upper_tangent

        let merge_neig_left = self.edge_counterclockwise(left);
        let merge_neig_right = self.edge_clockwise(right);
        self.tr_mut(merge_neig_left).neighbors[2] = merge_edge_id.into();
        self.tr_mut(merge_neig_right).neighbors[1] = merge_edge_id.into();

        self.tr_mut(left).neighbors[1] = lower_tangent_id.into();
        self.tr_mut(right).neighbors[2] = lower_tangent_id.into();


        loop{
            let right_candidat = self.next_right_candidate(merge_edge_id);
            let left_candidat = self.next_left_candidate(merge_edge_id);

            match (right_candidat, left_candidat){
                (Some(r), Some(l)) => {
                    if self.circumcircle_contain( (self.tr(merge_edge_id).points.0, self.tr(merge_edge_id).points.1, self.tr(r).points.0), self.tr(l).points.1){
                        self.merge_candidate_l(merge_edge_id, l);
                    }else if self.circumcircle_contain( (self.tr(merge_edge_id).points.0, self.tr(merge_edge_id).points.1, self.tr(l).points.1), self.tr(r).points.0) {
                        self.merge_candidate_r(merge_edge_id, r);
                    }else{
                        panic!("Delaunay::merge: cound not add either of candidates!");
                    }
                },
                (Some(r), None) => self.merge_candidate_r(merge_edge_id, r),
                (None, Some(l)) => self.merge_candidate_l(merge_edge_id, l),
                (None, None) => break
            }
        }

        // repair rightmost and leftmost edges since they could have been transformed into
        // triangles
        if self.tr(self.triangles[to*2-2].neighbors[1]).points.2.is_some(){
            // rightmost replaced by upper tanget
            self.triangles[to*2-2].neighbors[1] = merge_edge_id.into();
        }
        if self.tr(self.triangles[to*2-2].neighbors[0]).points.2.is_some(){
            println!("lower tangent: {:?}", self.tr(lower_tangent_id));
            self.triangles[to*2-2].neighbors[0] = merge_edge_id.into();
        }
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
    /// result is left and right 'ghost' triangles, containing left point in .1 and right in .0 correspondingly.
    fn find_lower_tangent(&self, left: EdgeIndex, right: EdgeIndex)->(EdgeIndex, EdgeIndex){
        let mut left = self.edge_clockwise(left);
        let mut right = self.edge_counterclockwise(right);
        loop{
            while self.is_clockwise_index(self.tr(left).points.0, self.tr(left).points.1, self.tr(right).points.0){
                left = self.edge_clockwise(left);
            }
            let mut done = true;
            while self.is_clockwise_index(self.tr(left).points.1, self.tr(right).points.0, self.tr(right).points.1){
                right = self.edge_counterclockwise(right);
                done = false;
            }
            if done{
                break;
            }
        }
        (left, right)
    }

    /// find next counterclockwise ghost triangle on hull
    fn edge_counterclockwise(&self, id: EdgeIndex)->EdgeIndex{
        EdgeIndex::new(self, self.tr(id).neighbors[1])
    }

    /// find next clockwise ghost triangle on hull
    fn edge_clockwise(&self, id: EdgeIndex)->EdgeIndex{
        EdgeIndex::new(self, self.tr(id).neighbors[2])
    }

    /// Searching for next candidate from right to append some edge to it
    fn next_right_candidate(&mut self, merge_edge_id: EdgeIndex)->Option<EdgeIndex>{
        let l = self.edge_counterclockwise(merge_edge_id);
        let r = self.edge_clockwise(merge_edge_id);
        
        if !self.is_counterclockwise_index(self.tr(l).points.0, self.tr(r).points.1, self.tr(r).points.0){
            None
        }else{
            while self.circumcircle_contain_next_candidate(r, self.tr(l).points.0){
                self.remove_edge_r(r);
            }
            
            Some(r)
        }
    }

    fn next_left_candidate(&mut self, merge_edge_id: EdgeIndex)->Option<EdgeIndex>{
        let l = self.edge_counterclockwise(merge_edge_id);
        let r = self.edge_clockwise(merge_edge_id);
        
        if !self.is_counterclockwise_index(self.tr(l).points.1, self.tr(l).points.0, self.tr(r).points.1){
            None
        }else{
            while self.circumcircle_contain_next_candidate(l, self.tr(r).points.1){
                self.remove_edge_l(l);
            }
            
            Some(l)
        }    
    }

    fn merge_candidate_r(&mut self, merge_edge_id: EdgeIndex, candidate: EdgeIndex){
        let b = self.tr(merge_edge_id).neighbors[0];
        let b_id = self.find_neighbor_id(b, merge_edge_id);

        let prev = self.tr(candidate).neighbors[2];
        self.tr_mut(merge_edge_id).neighbors[2] = prev;
        self.tr_mut(prev).neighbors[1] = merge_edge_id.into();

        self.tr_mut(merge_edge_id).neighbors[0] = candidate.into();
        let prev_point = self.tr(candidate).points.0;
        self.tr_mut(merge_edge_id).points.0 = prev_point;

        let points = (self.tr(merge_edge_id).points.1, self.tr(candidate).points.1, Some(self.tr(candidate).points.0));
        self.tr_mut(candidate).points = points;

        self.tr_mut(candidate).neighbors[1] = merge_edge_id.into();
        self.tr_mut(candidate).neighbors[2] = b;

        self.tr_mut(b).neighbors[b_id] = candidate.into();
    }

    fn merge_candidate_l(&mut self, merge_edge_id: EdgeIndex, candidate: EdgeIndex){
        let b = self.tr(merge_edge_id).neighbors[0];
        let b_id = self.find_neighbor_id(b, merge_edge_id);

        let next = self.tr(candidate).neighbors[1];
        self.tr_mut(merge_edge_id).neighbors[1] = next;
        self.tr_mut(next).neighbors[2] = merge_edge_id.into();

        self.tr_mut(merge_edge_id).neighbors[0] = candidate.into();
        let next_point = self.tr(candidate).points.1;
        self.tr_mut(merge_edge_id).points.1 = next_point;

        let points = (self.tr(merge_edge_id).points.0, self.tr(candidate).points.1, Some(self.tr(candidate).points.0));
        self.tr_mut(candidate).points = points;

        self.tr_mut(candidate).neighbors[1] = b;
        self.tr_mut(candidate).neighbors[2] = merge_edge_id.into();

        self.tr_mut(b).neighbors[b_id] = candidate.into();
    }

    fn circumcircle_contain_next_candidate(&self, edge: EdgeIndex, third_point: PointIndex)->bool{
        let opposite = self.tr(edge).neighbors[0];
        if self.tr(opposite).points.2.is_none(){
            debug_assert!(self.tr(opposite).neighbors[0].id() == edge.id());
            false
        }else{
            let n = self.tr(edge).neighbors[0];
            let n = self.tr(n);
            let next_candidat_pid = if n.neighbors[0].id() == edge.id(){
                n.points.0
            }else if n.neighbors[1].id() == edge.id(){
                n.points.1
            }else{
                debug_assert!(n.neighbors[2].id() == edge.id());
                n.points.2.unwrap() // safe, since was checked above
            };

            self.circumcircle_contain((self.tr(edge).points.0, self.tr(edge).points.1, third_point), next_candidat_pid)
        }
    }

    fn remove_edge_r(&mut self, edge: EdgeIndex){
        let o = self.tr(edge).neighbors[0];
        debug_assert!(self.tr(o).points.2.is_some());
        let p = self.tr(edge).neighbors[2];

        let o_id = self.find_neighbor_id(o, edge);

        let on = self.tr(o).neighbors[(o_id+1)%3];
        let op = self.tr(o).neighbors[(o_id+2)%3];

        let on_id = self.find_neighbor_id(on, o);

        self.tr_mut(edge).neighbors[0] = on;
        self.tr_mut(on).neighbors[on_id] = edge.into();

        self.tr_mut(edge).neighbors[2] = o;
        self.tr_mut(o).neighbors[1] = edge.into();

        self.tr_mut(o).neighbors[2] = p;
        self.tr_mut(p).neighbors[1] = o;

        self.tr_mut(o).neighbors[0] = op;

        let (o_points, edge_points) = match o_id{ // since opposite triangle should contain all three points unwrap() is safe
            0 => (
                (self.tr(o).points.1, self.tr(o).points.0, None),
                (self.tr(o).points.0, self.tr(o).points.2.unwrap(), None)
            ),
            1 => (
                (self.tr(o).points.2.unwrap(), self.tr(o).points.1, None),
                (self.tr(o).points.1, self.tr(o).points.0, None)
            ),
            _ => ( // cound be only 0, 1 or 2
                (self.tr(o).points.0, self.tr(o).points.2.unwrap(), None),
                (self.tr(o).points.2.unwrap(), self.tr(o).points.1, None)
            )
        };

        self.tr_mut(o).points = o_points;
        self.tr_mut(edge).points = edge_points;
    }

    fn remove_edge_l(&mut self, edge: EdgeIndex){
        let o = self.tr(edge).neighbors[0];
        debug_assert!(self.tr(o).points.2.is_some());
        let n = self.tr(edge).neighbors[1];

        let o_id = self.find_neighbor_id(o, edge);

        let on = self.tr(o).neighbors[(o_id+1)%3];
        let op = self.tr(o).neighbors[(o_id+2)%3];

        let op_id = self.find_neighbor_id(op, o);

        self.tr_mut(edge).neighbors[0] = op;
        self.tr_mut(op).neighbors[op_id] = edge.into();

        self.tr_mut(edge).neighbors[1] = o;
        self.tr_mut(o).neighbors[2] = edge.into();

        self.tr_mut(o).neighbors[1] = n;
        self.tr_mut(n).neighbors[2] = o;

        self.tr_mut(o).neighbors[0] = on;

        let (o_points, edge_points) = match o_id{ // since opposite triangle should contain all three points unwrap() is safe
            0 => (
                (self.tr(o).points.0, self.tr(o).points.2.unwrap(), None),
                (self.tr(o).points.1, self.tr(o).points.0, None)
            ),
            1 => (
                (self.tr(o).points.1, self.tr(o).points.0, None),
                (self.tr(o).points.2.unwrap(), self.tr(o).points.1, None)
            ),
            _ => ( // cound be only 0, 1 or 2
                (self.tr(o).points.2.unwrap(), self.tr(o).points.1, None),
                (self.tr(o).points.0, self.tr(o).points.2.unwrap(), None)
            )
        };

        self.tr_mut(o).points = o_points;
        self.tr_mut(edge).points = edge_points;
    }


    #[inline]
    fn find_neighbor_id<T1: TrIndex+Copy, T2: TrIndex+Copy>(&self, triangle: T1, neighbor: T2)->usize{
        for i in 0..3{
            if self.tr(triangle).neighbors[i].id() == neighbor.id(){
                return i;
            }
        }
        panic!("Neighbor is not neighbor of triangle");
    }

    fn circumcircle_contain(&self, (a, b, c): (PointIndex, PointIndex, PointIndex), d: PointIndex)->bool{
        circumcircle_contain((self.p(a), self.p(b), self.p(c)), self.p(d))
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

/// check if 3 points clockwise
fn is_clockwise(a: &Point, b: &Point, c: &Point)->bool{
    cross_product(&(b-a), &(c-a)) < 0.0
}

/// calculate z-component of cross product of vectors extended with z=0
fn cross_product(a: &Vector, b: &Vector) -> f64{
    a.x * b.y - a.y * b.x
}

fn circumcircle_contain((a, b, c): (&Point, &Point, &Point), d: &Point)->bool{
    let mat = cgmath::Matrix3::<f64>::new(
        a.x - d.x, a.y - d.y, a.x.powi(2) - d.x.powi(2) + a.y.powi(2) - d.y.powi(2),
        b.x - d.x, b.y - d.y, b.x.powi(2) - d.x.powi(2) + b.y.powi(2) - d.y.powi(2),
        c.x - d.x, c.y - d.y, c.x.powi(2) - d.x.powi(2) + c.y.powi(2) - d.y.powi(2)
    );
    let det = mat.determinant();

    if is_counterclockwise(a, b, c){
        det > 0.0
    }else{
        det < 0.0
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    extern crate rand;
    use self::rand::distributions::{IndependentSample, Range};

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

        assert_eq!(d.triangles[2].neighbors[0], TriangleIndex(1));
        assert_eq!(d.triangles[2].neighbors[1], TriangleIndex(1));

        assert_eq!(d.triangles[6].neighbors[0], TriangleIndex(5));
        assert_eq!(d.triangles[6].neighbors[1], TriangleIndex(5));

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

        assert_eq!(d.triangles[4].neighbors[0], TriangleIndex(2));
        assert_eq!(d.triangles[4].neighbors[1], TriangleIndex(1));

        assert_eq!(d.triangles[10].neighbors[0], TriangleIndex(8));
        assert_eq!(d.triangles[10].neighbors[1], TriangleIndex(8));

        // test for build3 line
        let points = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(2.0, 0.0),

            Point::new(3.0, -1.0),
            Point::new(3.0, 0.0),
            Point::new(3.0, 1.0)
        ];
        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(12, TriangleLike::default());
        d.build_3points(0);
        d.build_3points(3);

        assert_eq!(d.triangles[4].neighbors[0], TriangleIndex(3));
        assert_eq!(d.triangles[4].neighbors[1], TriangleIndex(2));

        assert_eq!(d.triangles[10].neighbors[0], TriangleIndex(9));
        assert_eq!(d.triangles[10].neighbors[1], TriangleIndex(8));


        // test following invariant for merge
        let mut d = prepare_diagram2();

        d.merge(0, 3, 6);
        assert_eq!(d.triangles[10].neighbors[0], TriangleIndex(5));
        assert_eq!(d.triangles[10].neighbors[1], TriangleIndex(7));
    }

    #[test]
    fn test_find_lower_tangent(){
        let points = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(2.0, 0.0),

            Point::new(3.0, -1.0),
            Point::new(3.0, 0.0),
            Point::new(3.0, 1.0)
        ];
        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(12, TriangleLike::default());
        d.build_3points(0);
        d.build_3points(3);

        let l = d.find_rightmost_edge(0..3);
        let r = d.find_leftmost_edge(3..6);
        
        let (l, r) = d.find_lower_tangent(l, r);
        assert_eq!(l, EdgeIndex(3));
        assert_eq!(r, EdgeIndex(6));
    }


    #[test]
    fn test_circumcircle_contains(){
        let p1 = Point::new(2.0, 0.0);
        let p2 = Point::new(-2.0, 0.0);
        let p3 = Point::new(0.0, 4.0);

        let tp0 = Point::new(0.0, 0.0);
        assert!(circumcircle_contain((&p1, &p2, &p3), &tp0));
        let tp1 = Point::new(0.0, -0.999999);
        assert!(circumcircle_contain((&p1, &p2, &p3), &tp1));
        let tp2 = Point::new(0.0, -1.0);
        assert!(!circumcircle_contain((&p1, &p2, &p3), &tp2));
        let tp3 = Point::new(1.0, -1.0);
        assert!(!circumcircle_contain((&p1, &p2, &p3), &tp3));
    }

    #[test]
    fn test_circumcircle_contains_next_candidate(){
        let points = vec![
            Point::new(2.0, 2.0),
            Point::new(2.0, 3.0),
            Point::new(3.0, 3.0),

            Point::new(4.0, -2.0),
            Point::new(4.0, 4.0),
            Point::new(6.0, 1.0)
        ];

        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(12, TriangleLike::default());
        d.build_3points(0);
        d.build_3points(3);

        assert!(d.circumcircle_contain_next_candidate(EdgeIndex(8), PointIndex(0)));
        assert!(!d.circumcircle_contain_next_candidate(EdgeIndex(3), PointIndex(3)));
    }


    #[test]
    fn test_remove_edge_r(){
        let points = vec![
            Point::new(2.0, 2.0),
            Point::new(2.0, 3.0),
            Point::new(3.0, 3.0),

            Point::new(4.0, -2.0),
            Point::new(4.0, 4.0),
            Point::new(6.0, 1.0)
        ];

        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(12, TriangleLike::default());
        d.build_3points(0);
        d.build_3points(3);

        d.remove_edge_r(EdgeIndex(8));

        assert_eq!(d.tr(TriangleIndex(8)).neighbors, [TriangleIndex(9), TriangleIndex(9), TriangleIndex(6)]);
        assert_eq!(d.tr(TriangleIndex(6)).neighbors, [TriangleIndex(7), TriangleIndex(8), TriangleIndex(7)]);
        assert_eq!(d.tr(TriangleIndex(7)).neighbors, [TriangleIndex(6), TriangleIndex(6), TriangleIndex(9)]);
        assert_eq!(d.tr(TriangleIndex(9)).neighbors, [TriangleIndex(8), TriangleIndex(7), TriangleIndex(8)]);

        assert_eq!(d.tr(TriangleIndex(6)).points, (PointIndex(4), PointIndex(5), None));
        assert_eq!(d.tr(TriangleIndex(8)).points, (PointIndex(5), PointIndex(3), None));
    }


    #[test]
    fn test_remove_edge_l(){
        let points = vec![
            Point::new(2.0, 2.0),
            Point::new(3.0, 2.0),
            Point::new(3.0, 3.0),

            Point::new(4.0, -2.0),
            Point::new(4.0, 4.0),
            Point::new(6.0, 1.0)
        ];
        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(12, TriangleLike::default());
        d.build_3points(0);
        d.build_3points(3);

        d.remove_edge_l(EdgeIndex(1));

        assert_eq!(d.tr(TriangleIndex(1)).neighbors, [TriangleIndex(3), TriangleIndex(0), TriangleIndex(3)]);
        assert_eq!(d.tr(TriangleIndex(0)).neighbors, [TriangleIndex(2), TriangleIndex(2), TriangleIndex(1)]);
        assert_eq!(d.tr(TriangleIndex(2)).neighbors, [TriangleIndex(0), TriangleIndex(3), TriangleIndex(0)]);
        assert_eq!(d.tr(TriangleIndex(3)).neighbors, [TriangleIndex(1), TriangleIndex(1), TriangleIndex(2)]);

        assert_eq!(d.tr(TriangleIndex(0)).points, (PointIndex(0), PointIndex(2), None));
        assert_eq!(d.tr(TriangleIndex(1)).points, (PointIndex(1), PointIndex(0), None));
    }

    #[test]
    fn test_find_next_candidate(){
        // TODO: make it clear
        let mut d = prepare_diagram1();
        let lower_tangent_id = EdgeIndex(4);
        let merge_edge_id = EdgeIndex(5);

        let lower_tangent = TriangleLike{
            points: (PointIndex(1), PointIndex(3), None),
            neighbors: [merge_edge_id.into(), TriangleIndex(9), TriangleIndex(3)]
        };
        let merge_edge = TriangleLike{
            points: (PointIndex(3), PointIndex(1), None),
            neighbors: [lower_tangent_id.into(), TriangleIndex(1), TriangleIndex(8)]
        };

        *d.tr_mut(lower_tangent_id) = lower_tangent;
        *d.tr_mut(merge_edge_id) = merge_edge;

        d.tr_mut(TriangleIndex(8)).neighbors[1] = merge_edge_id.into();
        d.tr_mut(TriangleIndex(1)).neighbors[2] = merge_edge_id.into();

        // and finally test
        let l = d.next_left_candidate(merge_edge_id);
        let r = d.next_right_candidate(merge_edge_id);

        assert_eq!(l, Some(EdgeIndex(1)));
        assert_eq!(r, Some(EdgeIndex(8)));

        let l = l.unwrap();
        let r = r.unwrap();

        assert_eq!(d.tr(l).neighbors, [TriangleIndex(0), TriangleIndex(2), TriangleIndex(5)]);
        assert_eq!(d.tr(l).points, (PointIndex(1), PointIndex(2), None));
        assert_eq!(d.tr(r).neighbors, [TriangleIndex(9), TriangleIndex(5), TriangleIndex(6)]);
        assert_eq!(d.tr(r).points, (PointIndex(5), PointIndex(3), None));

        // another one

        let mut d = prepare_diagram1();

        let lower_tangent_id = EdgeIndex(4);
        let merge_edge_id = EdgeIndex(5);

        let lower_tangent = TriangleLike{
            points: (PointIndex(1), PointIndex(3), None),
            neighbors: [merge_edge_id.into(), TriangleIndex(9), TriangleIndex(3)]
        };
        let merge_edge = TriangleLike{
            points: (PointIndex(4), PointIndex(1), None),
            neighbors: [lower_tangent_id.into(), TriangleIndex(1), TriangleIndex(7)]
        };
        *d.tr_mut(lower_tangent_id) = lower_tangent;
        *d.tr_mut(merge_edge_id) = merge_edge;

        d.tr_mut(TriangleIndex(7)).neighbors[1] = merge_edge_id.into();
        d.tr_mut(TriangleIndex(1)).neighbors[2] = merge_edge_id.into();

        let l = d.next_left_candidate(merge_edge_id);
        let r = d.next_right_candidate(merge_edge_id);

        assert_eq!(l, Some(EdgeIndex(1)));
        assert_eq!(r, None);
        
    }

    fn prepare_diagram1() -> Delaunay{
        let points = vec![
            Point::new(2.0, 2.0),
            Point::new(3.0, 2.0),
            Point::new(3.0, 3.0),

            Point::new(4.0, -2.0),
            Point::new(4.0, 4.0),
            Point::new(6.0, 1.0)
        ];
        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(12, TriangleLike::default());
        d.build_3points(0);
        d.build_3points(3);

        d
    }

    fn prepare_diagram2() -> Delaunay{
        let points = vec![
            Point::new(2.0, 2.0),
            Point::new(3.0, 2.0),
            Point::new(3.0, 3.0),

            Point::new(4.0, -2.0),
            Point::new(4.0, 5.0),
            Point::new(6.0, 1.0)
        ];
        let mut d = Delaunay{points: points, triangles: vec![]};
        d.triangles.resize(12, TriangleLike::default());
        d.build_3points(0);
        d.build_3points(3);

        d
    }

    #[test]
    fn test_merge(){
        let mut d = prepare_diagram2();

        d.merge(0, 3, 6);

        // triangles
        assert_eq!(d.tr(TriangleIndex(0)), &TriangleLike{
            points: (PointIndex(0), PointIndex(1), Some(PointIndex(2))),
            neighbors: [TriangleIndex(1), TriangleIndex(2), TriangleIndex(3)]
        });
        assert_eq!(d.tr(TriangleIndex(1)), &TriangleLike{
            points: (PointIndex(5), PointIndex(2), Some(PointIndex(1))),
            neighbors: [TriangleIndex(0), TriangleIndex(8), TriangleIndex(6)]
        });
        assert_eq!(d.tr(TriangleIndex(2)), &TriangleLike{
            points: (PointIndex(4), PointIndex(0), Some(PointIndex(2))),
            neighbors: [TriangleIndex(0), TriangleIndex(6), TriangleIndex(5)]
        });
        assert_eq!(d.tr(TriangleIndex(3)), &TriangleLike{
            points: (PointIndex(3), PointIndex(1), Some(PointIndex(0))),
            neighbors: [TriangleIndex(0), TriangleIndex(4), TriangleIndex(8)]
        });
        assert_eq!(d.tr(TriangleIndex(6)), &TriangleLike{
            points: (PointIndex(2), PointIndex(5), Some(PointIndex(4))),
            neighbors: [TriangleIndex(7), TriangleIndex(2), TriangleIndex(1)]
        });
        assert_eq!(d.tr(TriangleIndex(8)), &TriangleLike{
            points: (PointIndex(1), PointIndex(3), Some(PointIndex(5))),
            neighbors: [TriangleIndex(9), TriangleIndex(1), TriangleIndex(3)]
        });

        // edges
        assert_eq!(d.tr(TriangleIndex(4)), &TriangleLike{
            points: (PointIndex(0), PointIndex(3), None),
            neighbors: [TriangleIndex(3), TriangleIndex(9), TriangleIndex(5)]
        });
        assert_eq!(d.tr(TriangleIndex(5)), &TriangleLike{
            points: (PointIndex(4), PointIndex(0), None),
            neighbors: [TriangleIndex(2), TriangleIndex(4), TriangleIndex(7)]
        });
        assert_eq!(d.tr(TriangleIndex(7)), &TriangleLike{
            points: (PointIndex(5), PointIndex(4), None),
            neighbors: [TriangleIndex(6), TriangleIndex(5), TriangleIndex(9)]
        });
        assert_eq!(d.tr(TriangleIndex(9)), &TriangleLike{
            points: (PointIndex(3), PointIndex(5), None),
            neighbors: [TriangleIndex(8), TriangleIndex(7), TriangleIndex(4)]
        });

        test_for_delaunay_triangulation(&d);
    }

    fn test_for_delaunay_triangulation(d: &Delaunay){
        let mut edges_count = 0;
        let mut some_edge = None;
        for i in 0..(d.points.len()*2-2){
            let tr_id = TriangleIndex(i);
            if d.tr(tr_id).points.2.is_some(){
                // Triangle should contain no point
                for i in 0..d.points.len(){
                    let triangle = d.tr(tr_id).points;
                    let triangle = (triangle.0, triangle.1, triangle.2.unwrap());
                    assert!(!d.circumcircle_contain(triangle, PointIndex(i)), "'Triangle should contain not point' constrait violated");
                }
            }else{
                some_edge = Some(tr_id);
                edges_count += 1;
            }
        }

        let some_edge = EdgeIndex::new(d, some_edge.unwrap());
        let mut current = some_edge;
        loop{
            current = d.edge_clockwise(current);
            edges_count -= 1;
            if current == some_edge{
                break;
            }
        }
        assert_eq!(edges_count, 0);
    }

    fn test_random_delaunay_of_size(size: usize){
        let iter_count = 100;
        for _ in 0..iter_count{
            let points = random_point_set((0, 1000, 0, 1000), size);
            let d = Delaunay::new(points);

            test_for_delaunay_triangulation(&d);
        }
    }

    #[test]
    fn test_random_delaunay(){
        test_random_delaunay_of_size(5);
        test_random_delaunay_of_size(10);
        test_random_delaunay_of_size(20);
        test_random_delaunay_of_size(50);
    }


    fn random_point_set(rect: (i64, i64, i64, i64), size: usize)->Vec<Point>{
        let mut v = Vec::with_capacity(size);

        let x_between = Range::new(rect.0, rect.1);
        let y_between = Range::new(rect.2, rect.3);
        let mut rng = rand::thread_rng();
    
        for _ in 0..size{
            let x = x_between.ind_sample(&mut rng);
            let y = y_between.ind_sample(&mut rng);
            v.push(Point::new(x as f64, y as f64)); 
        }
        v
    }
    
    #[test]
    fn test_new(){
        let points = vec![
            Point::new(2.0, 1.0),
            Point::new(3.0, 8.0),

            Point::new(4.0, 8.0),
            Point::new(5.0, 1.0),
            Point::new(6.0, 5.0)
        ];

        //let d = Delaunay::new(points);
        //test_for_delaunay_triangulation(&d);


        let points = vec![
            Point::new(2.41, 3.67),
            Point::new(2.68, 4.30),
            Point::new(3.11, 1.46),
            Point::new(5.15, 1.47),
            Point::new(5.63, 0.80),
            Point::new(6.32, 7.64),
            Point::new(7.02, 0.05),
            Point::new(7.26, 0.59),
            Point::new(8.34, 1.76),
            Point::new(9.31, 5.10)
        ];
        let d = Delaunay::new(points);
        test_for_delaunay_triangulation(&d);
    }

}
