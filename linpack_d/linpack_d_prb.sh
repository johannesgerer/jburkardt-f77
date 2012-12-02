#!/bin/bash
#
gfortran -c -g linpack_d_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_d_prb.f"
  exit
fi
rm compiler.txt
#
gfortran linpack_d_prb.o -L$HOME/libf77/$ARCH -llinpack_d -lblas1_d
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linpack_d_prb.o"
  exit
fi
rm linpack_d_prb.o
#
mv a.out linpack_d_prb
./linpack_d_prb > linpack_d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running linpack_d_prb"
  exit
fi
rm linpack_d_prb
#
echo "Test results written to linpack_d_prb_output.txt."
