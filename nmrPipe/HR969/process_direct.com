#!/bin/csh -f

echo
echo Processing YZ dimensions
xyz2pipe -in data/test%04d.fid -x -verb \
| nmrPipe -fn SP -off 0.5 -end 0.95 -pow 2 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn FT \
| nmrPipe -fn PS -p0 22.0 -p1 0.0 -di \
| nmrPipe -fn POLY -auto -ord 0                    \
| nmrPipe -fn TP \
| nmrPipe -fn SP -off 0.50 -end 0.95 -pow 1 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn ZF -size 64 -auto \
| nmrPipe -fn TP \
| nmrPipe -fn ZTP \
| nmrPipe -fn SP -off 0.50 -end 0.95 -pow 1 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn ZF -size 32 -auto \
| nmrPipe -fn ZTP \
| nmrPipe -fn ATP \
| nmrPipe -fn SP -off 0.50 -end 0.95 -pow 1 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn ZF -size 16 -auto \
| nmrPipe -fn ATP \
| pipe2xyz -out data1/test%04d.ft1 -y -ov -verb      \

xyz2pipe  -in data1/test%04d.ft1 -y -verb            \
  > ./Hnoesy_4D.ft1


exit
    
