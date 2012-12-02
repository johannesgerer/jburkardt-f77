#!/bin/bash
#
gfortran -c -g odepack_prb3.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling odepack_prb3.f"
  exit
fi
rm compiler.txt
#
gfortran odepack_prb3.o -L$HOME/libf77/$ARCH -lodepack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading odepack_prb3.o"
  exit
fi
rm odepack_prb3.o
#
mv a.out odepack_prb3
./odepack_prb3 > odepack_prb3_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running odepack_prb3"
  exit
fi
rm odepack_prb3
#
echo "Test results written to odepack_prb3_output.txt."
