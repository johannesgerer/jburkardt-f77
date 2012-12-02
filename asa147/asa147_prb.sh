#!/bin/bash
#
gfortran -c -g asa147_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa147_prb.f"
  exit
fi
rm compiler.txt
#
gfortran asa147_prb.o -L$HOME/libf77/$ARCH -lasa147
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa147_prb.o"
  exit
fi
rm asa147_prb.o
#
mv a.out asa147_prb
./asa147_prb > asa147_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa147_prb"
  exit
fi
rm asa147_prb
#
echo "Test results written to asa147_prb_output.txt."
