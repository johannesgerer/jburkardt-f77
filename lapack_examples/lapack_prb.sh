#!/bin/bash
#
gfortran -c lapack_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lapack_prb.f"
  exit
fi
rm compiler.txt
#
gfortran lapack_prb.o -L$HOME/libf77/$ARCH -llapack -lblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lapack_prb.o"
  exit
fi
rm lapack_prb.o
#
mv a.out lapack_prb
./lapack_prb > lapack_prb_output.txt
rm lapack_prb
#
echo "Program output written to lapack_prb_output.txt"
