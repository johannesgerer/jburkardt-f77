#!/bin/bash
#
gfortran -c -g toms571_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms571_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms571_prb.o -L$HOME/libf77/$ARCH -ltoms571
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms571_prb.o"
  exit
fi
rm toms571_prb.o
#
mv a.out toms571_prb
./toms571_prb > toms571_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms571_prb"
  exit
fi
rm toms571_prb
#
echo "Test results written to toms571_prb_output.txt."
