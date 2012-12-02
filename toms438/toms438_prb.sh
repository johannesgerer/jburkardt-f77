#!/bin/bash
#
gfortran -c -g toms438_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms438_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms438_prb.o -L$HOME/libf77/$ARCH -ltoms438
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms438_prb.o"
  exit
fi
rm toms438_prb.o
#
mv a.out toms438_prb
./toms438_prb > toms438_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms438_prb"
  exit
fi
rm toms438_prb
#
echo "Test results written to toms438_prb_output.txt."
