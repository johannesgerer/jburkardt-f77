#!/bin/bash
#
gfortran -c -g toms347_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms347_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms347_prb.o -L$HOME/libf77/$ARCH -ltoms347
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms347_prb.o"
  exit
fi
rm toms347_prb.o
#
mv a.out toms347_prb
./toms347_prb > toms347_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms347_prb"
  exit
fi
rm toms347_prb
#
echo "Test results written to toms347_prb_output.txt."
