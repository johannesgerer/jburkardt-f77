#!/bin/bash
#
gfortran -c -g asa006_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa006_prb.f"
  exit
fi
rm compiler.txt
#
gfortran asa006_prb.o -L$HOME/libf77/$ARCH -lasa006
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa006_prb.o"
  exit
fi
rm asa006_prb.o
#
mv a.out asa006_prb
./asa006_prb > asa006_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa006_prb"
  exit
fi
rm asa006_prb
#
echo "Test results written to asa006_prb_output.txt."
