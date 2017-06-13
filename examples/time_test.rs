extern crate delaunay2;
extern crate rand;
extern crate cgmath;
extern crate time;

use time::{PreciseTime, Duration};
use delaunay2::{Delaunay, Point, test_for_delaunay_triangulation};
use rand::distributions::{IndependentSample, Range};

use std::fs::File;
use std::io::{BufWriter, Write};


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


fn test_uniform_size(size: usize)->(Duration, Duration){
    let mut points = random_point_set((0, 100000, 0, 100000), size);

    let sort_start = PreciseTime::now();
    Delaunay::sort_points(&mut points);
    let start = PreciseTime::now();
    let d = Delaunay::from_sorted(points);
    (sort_start.to(start), start.to(PreciseTime::now()))
}


fn main(){
    let mut f = File::create("time.csv").unwrap();
    for i in 10..22{
        let size = 1<<i;
        let time = test_uniform_size(size);
        writeln!(f, "{}, {:?}", size, time);
        println!("{}, {:?}", size, time)
    }
}