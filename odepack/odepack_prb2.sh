#!/bin/bash
#
gfortran -c -g odepack_prb2.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling odepack_prb2.f"
  exit
fi
rm compiler.txt
#
gfortran odepack_prb2.o -L$HOME/libf77/$ARCH -lodepack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading odepack_prb2.o"
  exit
fi
rm odepack_prb2.o
#
mv a.out odepack_prb2
./odepack_prb2 > odepack_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running odepack_prb2"
  exit
fi
rm odepack_prb2
#
echo "Test results written to odepack_prb2_output.txt."
