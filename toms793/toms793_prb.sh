#!/bin/bash
#
gfortran -c -g toms793_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms793_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms793_prb.o -L$HOME/libf77/$ARCH -ltoms793
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms793_prb.o"
  exit
fi
rm toms793_prb.o
#
mv a.out toms793_prb
./toms793_prb > toms793_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms793_prb"
  exit
fi
rm toms793_prb
#
echo "Test results written to toms793_prb_output.txt."
