#!/bin/bash
#
gfortran -c -g toms772_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms772_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms772_prb.o -L$HOME/libf77/$ARCH -ltoms772
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms772_prb.o"
  exit
fi
rm toms772_prb.o
#
mv a.out toms772_prb
./toms772_prb > toms772_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms772_prb"
  exit
fi
rm toms772_prb
#
echo "Test results written to toms772_prb_output.txt."
