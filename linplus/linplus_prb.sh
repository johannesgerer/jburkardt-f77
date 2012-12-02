#!/bin/bash
#
gfortran -c -g linplus_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linplus_prb.f"
  exit
fi
rm compiler.txt
#
gfortran linplus_prb.o -L$HOME/libf77/$ARCH -llinplus
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linplus_prb.o"
  exit
fi
rm linplus_prb.o
#
mv a.out linplus_prb
./linplus_prb > linplus_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running linplus_prb"
  exit
fi
rm linplus_prb
#
echo "Test results written to linplus_prb_output.txt."
