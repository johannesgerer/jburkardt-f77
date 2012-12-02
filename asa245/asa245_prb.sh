#!/bin/bash
#
gfortran -c -g asa245_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa245_prb.f"
  exit
fi
rm compiler.txt
#
gfortran asa245_prb.o -L$HOME/libf77/$ARCH -lasa245
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa245_prb.o"
  exit
fi
rm asa245_prb.o
#
mv a.out asa245_prb
./asa245_prb > asa245_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa245_prb"
  exit
fi
rm asa245_prb
#
echo "Test results written to asa245_prb_output.txt."
