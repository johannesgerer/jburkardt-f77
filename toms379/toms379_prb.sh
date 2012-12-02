#!/bin/bash
#
gfortran -c -g toms379_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms379_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms379_prb.o -L$HOME/libf77/$ARCH -ltoms379
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms379_prb.o"
  exit
fi
rm toms379_prb.o
#
mv a.out toms379_prb
./toms379_prb > toms379_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms379_prb"
  exit
fi
rm toms379_prb
#
echo "Test results written to toms379_prb_output.txt."
