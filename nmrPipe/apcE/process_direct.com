#!/bin/csh

xyz2pipe -in data/test%03d.fid -x -verb              \
| nmrPipe  -fn SP -size 512 -off 0.5 -end 0.99 -pow 2 -c 0.5            \
| nmrPipe  -fn FT                                             \
| nmrPipe  -fn PS -p0 187.0 -p1 0.0                         \
| nmrPipe  -fn EXT -x1 -2ppm -xn 7.5ppm -di -sw     \
| pipe2xyz -out data/test%03d.ft1 -y -ov -verb      \

xyz2pipe  -in data/test%03d.ft1 -y -verb            \
  > ./d2o_hcch_tocsy.ft1


exit
    
