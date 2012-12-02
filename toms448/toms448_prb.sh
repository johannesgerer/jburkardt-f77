#!/bin/bash
#
gfortran -c -g toms448_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms448_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms448_prb.o -L$HOME/libf77/$ARCH -ltoms448
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms448_prb.o"
  exit
fi
rm toms448_prb.o
#
mv a.out toms448_prb
./toms448_prb > toms448_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms448_prb"
  exit
fi
rm toms448_prb
#
echo "Test results written to toms448_prb_output.txt."
