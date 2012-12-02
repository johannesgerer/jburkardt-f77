#!/bin/bash
#
gfortran -c -g toms352_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms352_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms352_prb.o -L$HOME/libf77/$ARCH -ltoms352
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms352_prb.o"
  exit
fi
rm toms352_prb.o
#
mv a.out toms352_prb
./toms352_prb > toms352_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms352_prb"
  exit
fi
rm toms352_prb
#
echo "Test results written to toms352_prb_output.txt."
