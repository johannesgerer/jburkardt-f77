#!/bin/bash
#
gfortran -c -g asa103_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa103_prb.f"
  exit
fi
rm compiler.txt
#
gfortran asa103_prb.o -L$HOME/libf77/$ARCH -lasa103
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa103_prb.o"
  exit
fi
rm asa103_prb.o
#
mv a.out asa103_prb
./asa103_prb > asa103_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa103_prb"
  exit
fi
rm asa103_prb
#
echo "Test results written to asa103_prb_output.txt."
