#!/bin/csh

xyz2pipe -in d2o_hcch_tocsy.ft1 -x -verb \
| nmrPipe  -fn TP                                             \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5            \
| nmrPipe  -fn FT                                             \
| nmrPipe  -fn PS -p0 0.0 -p1  0.0                           \
| nmrPipe  -fn CS  -ls 2ppm   -sw                           \
| nmrPipe  -fn EXT -x1 -1ppm -xn 6.5ppm -di -sw     \
| nmrPipe  -fn TP                                             \
| nmrPipe  -fn POLY -auto -ord 0                                    \
| nmrPipe  -fn TP                                    \
| nmrPipe  -fn POLY -auto -ord 0                        \
| pipe2xyz -out data/test%03d.ft2 -x -ov

xyz2pipe -in data/test%03d.ft2 -z -verb                        \
| nmrPipe  -fn SP -off 0.50 -end 0.99 -pow 1 -c 0.5            \
| nmrPipe  -fn FT                                             \
| nmrPipe  -fn PS -p0  0.0 -p1 0.0 -di                       \
| nmrPipe  -fn POLY -ord 0                        \
| pipe2xyz -out ft/test%03d.ft3 -z -ov

/bin/rm -rf ft2

proj3D.tcl -in ft/test%03d.ft3 -axis -outDir ftproj

