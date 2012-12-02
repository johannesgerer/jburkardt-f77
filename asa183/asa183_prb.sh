#!/bin/bash
#
gfortran -c -g asa183_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa183_prb.f"
  exit
fi
rm compiler.txt
#
gfortran asa183_prb.o -L$HOME/libf77/$ARCH -lasa183
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa183_prb.o"
  exit
fi
rm asa183_prb.o
#
mv a.out asa183_prb
./asa183_prb > asa183_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa183_prb"
  exit
fi
rm asa183_prb
#
echo "Test results written to asa183_prb_output.txt."
