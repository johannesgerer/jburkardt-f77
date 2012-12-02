#!/bin/bash
#
gfortran -c -g asa299_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa299_prb.f"
  exit
fi
rm compiler.txt
#
gfortran asa299_prb.o -L$HOME/libf77/$ARCH -lasa299
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa299_prb.o"
  exit
fi
rm asa299_prb.o
#
mv a.out asa299_prb
./asa299_prb > asa299_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa299_prb"
  exit
fi
rm asa299_prb
#
echo "Test results written to asa299_prb_output.txt."
