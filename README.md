# delaunay_rs

This library could be used to build [delaunay triangulation](https://github.com/hukumka/delaunay_rs/) of given set of points.

## warning!
This library been developed for studing purpose, might contain several issues, 
probaby won't be maintained and strongly not recommended to use in serious projects.

## requiments

1. [rust 1.18+](https://www.rust-lang.org)

## build details

### as rust library
to use this in your rust project add following line to your [dependencies]

delaunay2 = {git = "https://github.com/hukumka/delaunay_rs/", rev = "91263b1"}

usage might be seen in examples

### as external executable
build example triangulate

1. `git clone https://github.com/hukumka/delaunay_rs/`
2. `cd delaunay_rs`
3. `cargo build --release --example triangulate`

Executable will appear as ./target/release/examples/triangulate

executable usage:

`triangulate <in_file> <out_file>`

there <in_file> contains set of points, each represented by line with two floating point values separated by coma. For example:
`10.32, 1.09`

After execution file <out_file> will containt list of triangulation edges as pairs of points id separated with coma.
Point id correspond to point position in input file (numeration starts with 0)
