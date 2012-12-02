#!/bin/bash
#
gfortran -c -g sparse_interp_nd_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_interp_nd_prb.f"
  exit
fi
rm compiler.txt
#
gfortran sparse_interp_nd_prb.o -L$HOME/libf77/$ARCH -lsparse_interp_nd -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_interp_nd_prb.o"
  exit
fi
rm sparse_interp_nd_prb.o
#
mv a.out sparse_interp_nd_prb
./sparse_interp_nd_prb > sparse_interp_nd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparse_interp_nd_prb"
  exit
fi
rm sparse_interp_nd_prb
#
echo "Test program output written to sparse_interp_nd_prb_output.txt."
