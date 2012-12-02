#!/bin/bash
#
gfortran -c -g toms322_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms322_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms322_prb.o -L$HOME/libf77/$ARCH -ltoms322
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms322_prb.o"
  exit
fi
rm toms322_prb.o
#
mv a.out toms322_prb
./toms322_prb > toms322_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms322_prb"
  exit
fi
rm toms322_prb
#
echo "Test results written to toms322_prb_output.txt."
