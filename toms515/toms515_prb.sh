#!/bin/bash
#
gfortran -c -g toms515_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms515_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms515_prb.o -L$HOME/libf77/$ARCH -ltoms515
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms515_prb.o"
  exit
fi
rm toms515_prb.o
#
mv a.out toms515_prb
./toms515_prb > toms515_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms515_prb"
  exit
fi
rm toms515_prb
#
echo "Test results written to toms515_prb_output.txt."
