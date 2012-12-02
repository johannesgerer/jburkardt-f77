#!/bin/bash
#
gfortran -c -g linpack_c_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_c_prb.f"
  exit
fi
rm compiler.txt
#
gfortran linpack_c_prb.o -L$HOME/libf77/$ARCH -llinpack_c -lblas1_c -lblas1_s
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linpack_c_prb.o"
  exit
fi
rm linpack_c_prb.o
#
mv a.out linpack_c_prb
./linpack_c_prb > linpack_c_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running linpack_c_prb"
  exit
fi
rm linpack_c_prb
#
echo "Test results written to linpack_c_prb_output.txt."
