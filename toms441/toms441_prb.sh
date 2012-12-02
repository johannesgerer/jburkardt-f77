#!/bin/bash
#
gfortran -c -g toms441_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms441_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms441_prb.o -L$HOME/libf77/$ARCH -ltoms441
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms441_prb.o"
  exit
fi
rm toms441_prb.o
#
mv a.out toms441_prb
./toms441_prb > toms441_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms441_prb"
  exit
fi
rm toms441_prb
#
echo "Test results written to toms441_prb_output.txt."
