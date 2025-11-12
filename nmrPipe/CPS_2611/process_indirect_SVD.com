#!/bin/csh

xyz2pipe -in Nnoesy_recon_SVD.ft3 -x -verb \
| nmrPipe  -fn TP                                    \
| nmrPipe  -fn SP -off 0.50 -end 0.98 -pow 2 -c 0.5   \
| nmrPipe  -fn FT                                    \
| nmrPipe  -fn PS -p0 43.0 -p1 0.0 -di                \
| nmrPipe  -fn POLY -ord 0                           \
| nmrPipe  -fn TP                                    \
| nmrPipe  -fn POLY -ord 0                           \
| pipe2xyz -out ft2/test%03d.ft2 -x -ov -verb      \

xyz2pipe  -in ft2/test%03d.ft2 -z -verb            \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5   \
| nmrPipe  -fn FT                                    \
| nmrPipe  -fn PS -p0 0.0 -p1 0.0 -di                \
| nmrPipe  -fn POLY -ord 0                           \
| pipe2xyz -out ft_SVD/test%03d.ft3 -z -ov

/bin/rm -rf ft2

proj3D.tcl -in ft_SVD/test%03d.ft3 -axis -outDir ftproj_SVD

