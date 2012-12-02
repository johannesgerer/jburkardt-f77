#!/bin/bash
#
gfortran -c -g toms661_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms661_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms661_prb.o -L$HOME/libf77/$ARCH -ltoms661
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms661_prb.o"
  exit
fi
rm toms661_prb.o
#
mv a.out toms661_prb
./toms661_prb > toms661_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms661_prb"
  exit
fi
rm toms661_prb
#
echo "Program output written to toms661_prb_output.txt"
