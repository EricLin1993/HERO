#!/bin/csh

xyz2pipe -in Hnoesy_4D.ft1 -x -verb \
| nmrPipe -fn TP \
| nmrPipe -fn SP -off 0.50 -end 0.95 -pow 1 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn ZF -size 64 -auto \
| nmrPipe -fn FT \
| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \
| nmrPipe -fn TP \
| nmrPipe -fn ZTP \
| nmrPipe -fn SP -off 0.50 -end 0.95 -pow 1 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn ZF -size 32 -auto \
| nmrPipe -fn FT \
| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \
| nmrPipe -fn ZTP \
| pipe2xyz -out ft3/test%05d.ft3 -x -ov

xyz2pipe -in ft3/test%05d.ft3 -a -verb \
| nmrPipe -fn SP -off 0.50 -end 0.95 -pow 1 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn ZF -size 16 -auto \
| nmrPipe -fn FT \
| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \
| pipe2xyz -out ft/test%03d.ft4 -a -ov

/bin/rm -rf ft3

proj4D.tcl -in ft/test%03d.ft4 -axis -outDir ftproj

