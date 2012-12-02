#!/bin/bash
#
gfortran -c -g toms429_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms429_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms429_prb.o -L$HOME/libf77/$ARCH -ltoms429
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms429_prb.o"
  exit
fi
rm toms429_prb.o
#
mv a.out toms429_prb
./toms429_prb > toms429_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms429_prb"
  exit
fi
rm toms429_prb
#
echo "Test results written to toms429_prb_output.txt."
