#!/bin/csh -f

echo
echo Processing YZ dimensions
xyz2pipe -in fid/test%03d.fid -x -verb \
| nmrPipe -fn SP -off 0.5 -end 0.95 -pow 2 -elb 0.0 -glb 0.0 -c 0.5 \
| nmrPipe -fn FT \
| nmrPipe -fn PS -p0 -195.0 -p1 0.0 -di \
| nmrPipe  -fn POLY -auto -ord 0                    \
| pipe2xyz -out data1/test%03d.ft1 -y -ov -verb      \

xyz2pipe  -in data1/test%03d.ft1 -y -verb            \
  > ./D2O_Cnoesy.ft1


exit
    
