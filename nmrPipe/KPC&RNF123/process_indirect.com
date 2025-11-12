#!/bin/csh

xyz2pipe -in D2O_Cnoesy.ft1 -x -verb \
| nmrPipe -fn TP \
| nmrPipe -fn SP -off 0.50 -end 0.95 -pow 1 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn FT \
| nmrPipe -fn PS -p0 -94.0 -p1 0.0 -di \
| nmrPipe -fn TP \
| pipe2xyz -out ft2/test%05d.ft2 -x -ov

xyz2pipe -in ft2/test%05d.ft2 -z -verb \
| nmrPipe -fn SP -off 0.50 -end 0.95 -pow 1 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn ZF -zf 1 -auto \
| nmrPipe -fn FT \
| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \
| pipe2xyz -out ft/test%03d.ft3 -z -ov

/bin/rm -rf ft2

proj3D.tcl -in ft/test%03d.ft3 -axis -outDir ftproj

