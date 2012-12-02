#!/bin/bash
#
gfortran -c -g asa005_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa005_prb.f"
  exit
fi
rm compiler.txt
#
gfortran asa005_prb.o -L$HOME/libf77/$ARCH -lasa005
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa005_prb.o"
  exit
fi
rm asa005_prb.o
#
mv a.out asa005_prb
./asa005_prb > asa005_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa005_prb"
  exit
fi
rm asa005_prb
#
echo "Test results written to asa005_prb_output.txt."
