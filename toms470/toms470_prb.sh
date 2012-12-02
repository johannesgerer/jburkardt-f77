#!/bin/bash
#
gfortran -c -g toms470_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms470_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms470_prb.o -L$HOME/libf77/$ARCH -ltoms470
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms470_prb.o"
  exit
fi
rm toms470_prb.o
#
mv a.out toms470_prb
./toms470_prb > toms470_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms470_prb"
  exit
fi
rm toms470_prb
#
echo "Test results written to toms470_prb_output.txt."
