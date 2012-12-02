#!/bin/bash
#
gfortran -c -g odepack_prb9.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling odepack_prb9.f"
  exit
fi
rm compiler.txt
#
gfortran odepack_prb9.o -L$HOME/libf77/$ARCH -lodepack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading odepack_prb9.o"
  exit
fi
rm odepack_prb9.o
#
mv a.out odepack_prb9
./odepack_prb9 > odepack_prb9_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running odepack_prb9"
  exit
fi
rm odepack_prb9
#
echo "Test results written to odepack_prb9_output.txt."
