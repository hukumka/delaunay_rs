#[macro_use]
extern crate quick_error;

extern crate delaunay2;
extern crate cgmath;

use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Write;

use delaunay2::{Delaunay, Point};


fn main(){
    let args: Vec<_> = env::args().collect();
    if args.len() < 3 || args[1] == "-h" || args[1] == "--help"{
        println!("triangulate <input> <output>");
        println!("  <input> - name of input file with following format: new line for ");
        println!("  each point where coordinates separated with coma. for example:");
        println!("      -21.356, 13.0");
        println!("      -14.356, 15.0");
        println!("  <output> - name of output file containing edges by indexes");
        println!("  where index is line number for point in input file");
    }else{
        let input = match File::open(&args[1]){
            Ok(file) => fil,
            Err(..) => panic!("Unable to open input file {}", args[1])
        };

        let mut output = match File::create(&args[2]){
            Ok(file) => file,
            Err(..) => panic!("Unable to use output file {}", args[2])
        };

        let points = read_data(&input);
        let shuffle = Delaunay::sort_point_indexes(&points);
        let d = Delaunay::new(points);

        let triangles = d.data();
        let mut edges = vec![];

        for (a, b, c) in triangles{
            let (a, b, c) = (shuffle[a], shuffle[b], shuffle[c]);
            let pair = |a, b|{
                if a < b{
                    (a, b)
                }else{
                    (b, a)
                }
            };
            edges.push(pair(a, b));
            edges.push(pair(b, c));
            edges.push(pair(c, a));
        }
        edges.sort();
        edges.dedup();

        for (a, b) in edges{
            writeln!(output, "{}, {}", a, b).unwrap();
        }
    }
}


fn read_data(input: &File)->Vec<Point>{
    let mut points = vec![];

    let mut reader = BufReader::new(input);
    let mut buf = String::new();
    let mut counter = 1;
    while let Ok(len) = reader.read_line(&mut buf){
        if len == 0{
            break;
        }
        let point = match read_point(&buf){
            Ok(point) => point,
            Err(err) => {panic!("Error in line {}: {}", counter, err)}
        };
        points.push(point);
        buf.clear();
        counter += 1;
    }
    points
}


quick_error!{
    #[derive(Debug)]
    enum ParseError{
        ParseFloat(err: std::num::ParseFloatError){
            from()
        }
        Other(err: String){
            description(err)
        }
    }
}


fn read_point(buf: &str)->Result<Point, ParseError>{
    let strings: Vec<_> = buf.split(",").collect();
    if strings.len() != 2{
        Err(ParseError::Other(format!("{}: 2 coordinates must be separated with coma", buf)))
    }else{
        let x: f64 = strings[0].trim().parse() ?;
        let y: f64 = strings[1].trim().parse() ?;
        Ok(Point::new(x, y))
    }
}
