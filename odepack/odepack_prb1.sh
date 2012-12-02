#!/bin/bash
#
gfortran -c -g odepack_prb1.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling odepack_prb1.f"
  exit
fi
rm compiler.txt
#
gfortran odepack_prb1.o -L$HOME/libf77/$ARCH -lodepack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading odepack_prb1.o"
  exit
fi
rm odepack_prb1.o
#
mv a.out odepack_prb1
./odepack_prb1 > odepack_prb1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running odepack_prb1"
  exit
fi
rm odepack_prb1
#
echo "Test results written to odepack_prb1_output.txt."
