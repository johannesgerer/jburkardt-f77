#!/bin/bash
#
gfortran -c -g odepack_prb4.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling odepack_prb4.f"
  exit
fi
rm compiler.txt
#
gfortran odepack_prb4.o -L$HOME/libf77/$ARCH -lodepack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading odepack_prb4.o"
  exit
fi
rm odepack_prb4.o
#
mv a.out odepack_prb4
./odepack_prb4 > odepack_prb4_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running odepack_prb4"
  exit
fi
rm odepack_prb4
#
echo "Test results written to odepack_prb4_output.txt."
