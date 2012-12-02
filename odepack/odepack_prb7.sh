#!/bin/bash
#
gfortran -c -g odepack_prb7.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling odepack_prb7.f"
  exit
fi
rm compiler.txt
#
gfortran odepack_prb7.o -L$HOME/libf77/$ARCH -lodepack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading odepack_prb7.o"
  exit
fi
rm odepack_prb7.o
#
mv a.out odepack_prb7
./odepack_prb7 > odepack_prb7_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running odepack_prb7"
  exit
fi
rm odepack_prb7
#
echo "Test results written to odepack_prb7_output.txt."
