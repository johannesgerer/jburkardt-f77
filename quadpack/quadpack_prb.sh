#!/bin/bash
#
gfortran -c -g quadpack_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadpack_prb.f"
  exit
fi
rm compiler.txt
#
gfortran quadpack_prb.o -L$HOME/libf77/$ARCH -lquadpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadpack_prb.o"
  exit
fi
rm quadpack_prb.o
#
mv a.out quadpack_prb
./quadpack_prb > quadpack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quadpack_prb"
  exit
fi
rm quadpack_prb
#
echo "Test results written to quadpack_prb_output.txt."
