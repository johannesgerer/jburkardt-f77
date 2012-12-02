#!/bin/bash
#
gfortran -c -g toms423_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms423_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms423_prb.o -L$HOME/libf77/$ARCH -ltoms423
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms423_prb.o"
  exit
fi
rm toms423_prb.o
#
mv a.out toms423_prb
./toms423_prb > toms423_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms423_prb"
  exit
fi
rm toms423_prb
#
echo "Test results written to toms423_prb_output.txt."
