use std::ops::Range;
use std::cmp::Ordering;


extern crate cgmath;
use cgmath::prelude::*;

type Point = cgmath::Point2<f64>;


/// Represent point id in delaunay data
#[derive(Debug, Clone, Copy)]
struct PointIndex(usize);
/// Represent triangle id in delaunay data
trait TrIndex{
    fn id(&self)->usize;
}

/// Represent any triangle index (further check required)
#[derive(Debug, Clone, Copy)]
struct TriangleIndex(usize);
impl TrIndex for TriangleIndex{
    fn id(&self)->usize{
        self.0
    }
}

/// Represent ghost triangle index
#[derive(Debug, Clone, Copy)]
struct EdgeIndex(usize);
impl TrIndex for EdgeIndex{
    fn id(&self)->usize{
        self.0
    }
}


/// Represent triangle in delaunay data
///
/// Contain 2 cases:
/// Regular triangle
///     points.3 is Some(PointIndex)
///     points kept from leftmost in counterclockwise order
///     neighbors[i] - triangle having same edge as opposite to points[i]
/// Ghost triangle
///     points.3 is None
///     represent hull edge
///     neighbors[0] - triangle on the opposite of edge
///     neighbors[1] - next counterclockwise edge on hull
///     neighbors[2] - previous counterclockwise edge on hull (next clockwise)
///     points kept in counterclockwise order on hull
#[derive(Debug, Clone)]
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
        let mut d = Delaunay{points: points, triangles: vec![]};
        d.sort_points();
        // since sort_points not only sorting points, but also remove duplicates
        // len calculation should be exactly here
        let len = d.points.len();
        d.triangles.resize(len*2, TriangleLike::default());
        d.build(0..len);
        d
    }


    /// sort points in from left to right order
    ///
    /// if several points has similar x-coordinate, then sorting done by y-coordinate
    /// removes duplicate points
    fn sort_points(&mut self){
        self.points.sort_by(|a, b|{
            debug_assert!(
                a.x.is_finite()
                && a.y.is_finite()
                && b.x.is_finite()
                && b.y.is_finite()
                , "Delaunay do not support infinite point coordinates."
            );
            // since none coordinate is NaN unwrap is correct
            let x_cmp = a.x.partial_cmp(&b.x).unwrap();
            if x_cmp == Ordering::Equal{
                a.y.partial_cmp(&b.y).unwrap()
            }else{
                x_cmp
            }
        });
        self.points.dedup();
    }

    /// Construct triangulation on selected slice of points
    ///
    /// space for triangles must be preallocated
    /// triangles kept in range.start*2..range..(end*2-2)
    /// triangles (end*2-2)..end*2 reserved
    fn build(&mut self, range: Range<usize>){
        match range.len(){
            2 => self.build_2points(range),
            3 => self.build_3points(range),
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
    fn build_2points(&mut self, range: Range<usize>){
        debug_assert!(range.len() == 2);

        let i1 = TriangleIndex(range.start * 2); // index to t1
        let i2 = TriangleIndex(range.start * 2 + 1); // index to t2

        let t1 = TriangleLike{
            neighbors: [i2, i2, i2],
            points: (PointIndex(range.start), PointIndex(range.start+1), None)
        };
        let t2 = TriangleLike{
            neighbors: [i1, i1, i1],
            points: (PointIndex(range.start+1), PointIndex(range.start), None)
        };

        *self.tr_mut(i1) = t1;
        *self.tr_mut(i2) = t2;
    }

    /// Construct triangulation on 3 points
    ///
    /// consist of triangle and
    /// hull through 'ghost' triangles
    fn build_3points(&mut self, range: Range<usize>){
        // TODO
    }

    fn tr_mut<T: TrIndex>(&mut self, id: T)->&mut TriangleLike{
        &mut self.triangles[id.id()]
    }

    /// merge to partial triangulation
    ///
    /// left triangulation consist from points in from..sep
    /// right from points in sep..to
    fn merge(&mut self, from: usize, sep: usize, to: usize){
        // TODO
    }

}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sort_points(){
        let points = vec![
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

        let mut d = Delaunay{points: points, triangles: vec![]};
        d.sort_points();
        assert_eq!(d.points, expected);
    }
}
