#!/bin/bash
#
gfortran -c -g toms461_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms461_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms461_prb.o -L$HOME/libf77/$ARCH -ltoms461
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms461_prb.o"
  exit
fi
rm toms461_prb.o
#
mv a.out toms461_prb
./toms461_prb > toms461_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms461_prb"
  exit
fi
rm toms461_prb
#
echo "Test results written to toms461_prb_output.txt."
