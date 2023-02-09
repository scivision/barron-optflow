[![Build status](https://ci.appveyor.com/api/projects/status/y45eymwuq53pgnsa?svg=true)](https://ci.appveyor.com/project/scivision/barron-optflow)

# Barron Optical Flow Code

[Original code](http://www.csd.uwo.ca/faculty/barron/FTP/)
modified to run with modern compilers.

[Examples](https://scivision.dev/barron1994opticalflow/)

## Build

A C compiler is required.
Some of the programs require a Fortran compiler and/or BLAS library.

-   Linux: `apt install gcc gfortran libopenblas-dev`
-   Mac: `brew install gcc openblas`

```sh
cmake -Bbuild
cmake --build build
```

## Horn-Schunck

1. Run Horn-Schunck optical flow on Tree image set
   ```sh
   ./horn new2binarytreet. 0.5 1.5 20 100 ../TESTDATA/TREE_DATA/TRANS out/ -B 150 150 -MH -T 5.0
   ```
2. convert output to Postscript and view (can use other viewer beside Evince)
   ```sh
   ./flow2ps out/horn.modified.new2binarytreet.F-5.00 out/hornOF.ps

   evince out/hornOF.ps
   ```

## Lucas-Kanade

1. Run Horn-Schunck optical flow on Tree image set
   ```sh
   ./lucas new2binarytreet. 1.5 20 1.0 ../TESTDATA/TREE_DATA/TRANS out/ -B 150 150
   ```
2. convert output to Postscript and view (can use other viewer beside Evince)
   ```sh
   ./flow2ps out/lucas.new2binarytreet.20F-1.00-1.5 out/lucas20F.ps

   evince out/lucas20F.ps
   ```
