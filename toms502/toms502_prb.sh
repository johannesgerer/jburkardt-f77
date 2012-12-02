#!/bin/bash
#
gfortran -c -g toms502_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms502_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms502_prb.o -L$HOME/libf77/$ARCH -ltoms502
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms502_prb.o"
  exit
fi
rm toms502_prb.o
#
mv a.out toms502_prb
./toms502_prb > toms502_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms502_prb"
  exit
fi
rm toms502_prb
#
echo "Test results written to toms502_prb_output.txt."
