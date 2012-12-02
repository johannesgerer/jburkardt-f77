#!/bin/bash
#
gfortran -c -g cg_plus_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cg_plus_prb.f"
  exit
fi
rm compiler.txt
#
gfortran cg_plus_prb.o -L$HOME/libf77/$ARCH -lcg_plus
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cg_plus_prb.o"
  exit
fi
rm cg_plus_prb.o
#
mv a.out cg_plus_prb
./cg_plus_prb > cg_plus_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cg_plus_prb"
  exit
fi
rm cg_plus_prb
#
echo "Test results written to cg_plus_prb_output.txt."
