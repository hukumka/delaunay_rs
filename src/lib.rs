extern crate cgmath;
use cgmath::prelude::*;

type Point = cgmath::Point2<f64>;


/// Represent point id in delaunay data
struct PointIndex(usize);
/// Represent triangle id in delaunay data
trait TrIndex{
    fn id(&self)->usize;
}

/// Represent any triangle index (further check required)
struct TriangleIndex(usize);
impl TrIndex for TriangleIndex{
    fn id(&self)->usize{
        self.0
    }
}

/// Represent ghost triangle index
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
struct TriangleLike{
    neighbors: [TriangleIndex; 3],
    points: (PointIndex, PointIndex, Option<PointIndex>)
}


/// Represent Delaunay triangulation
struct Delaunay{
    points: Vec<Point>,
    triangles: Vec<TriangleLike>
}


impl Delaunay{
}


#[cfg(test)]
mod tests {
}
