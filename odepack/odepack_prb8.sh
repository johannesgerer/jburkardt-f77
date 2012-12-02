#!/bin/bash
#
gfortran -c -g odepack_prb8.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling odepack_prb8.f"
  exit
fi
rm compiler.txt
#
gfortran odepack_prb8.o -L$HOME/libf77/$ARCH -lodepack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading odepack_prb8.o"
  exit
fi
rm odepack_prb8.o
#
mv a.out odepack_prb8
./odepack_prb8 > odepack_prb8_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running odepack_prb8"
  exit
fi
rm odepack_prb8
#
echo "Test results written to odepack_prb8_output.txt."
