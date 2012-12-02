#!/bin/bash
#
gfortran -c -g toms467_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms467_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms467_prb.o -L$HOME/libf77/$ARCH -ltoms467
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms467_prb.o"
  exit
fi
rm toms467_prb.o
#
mv a.out toms467_prb
./toms467_prb > toms467_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms467_prb"
  exit
fi
rm toms467_prb
#
echo "Test results written to toms467_prb_output.txt."
