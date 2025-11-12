#!/bin/csh

xyz2pipe -in data/test%03d.fid -x -verb              \
| nmrPipe  -fn POLY -time                            \
| nmrPipe  -fn SP -size 512 -off 0.5 -end 0.98 -pow 2 -c 0.5   \
| nmrPipe  -fn FT                                    \
| nmrPipe  -fn PS -p0  174.0 -p1 0.0                   \
| nmrPipe  -fn EXT -x1 12.0ppm -xn 5.0ppm -di -sw      \
| pipe2xyz -out data/test%03d.ft1 -y -ov -verb      \

xyz2pipe  -in data/test%03d.ft1 -y -verb            \
  > ./csr4_n15_noesy1.ft


exit
    
