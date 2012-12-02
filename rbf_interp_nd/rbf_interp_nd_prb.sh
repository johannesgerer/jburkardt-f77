#!/bin/bash
#
gfortran -c -g rbf_interp_nd_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rbf_interp_nd_prb.f"
  exit
fi
rm compiler.txt
#
gfortran rbf_interp_nd_prb.o -L$HOME/libf77/$ARCH -lrbf_interp_nd -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rbf_interp_nd_prb.o"
  exit
fi
rm rbf_interp_nd_prb.o
#
mv a.out rbf_interp_nd_prb
./rbf_interp_nd_prb > rbf_interp_nd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rbf_interp_nd_prb"
  exit
fi
rm rbf_interp_nd_prb
#
echo "Test results written to rbf_interp_nd_prb_output.txt."
