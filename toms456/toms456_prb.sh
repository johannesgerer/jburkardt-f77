#!/bin/bash
#
gfortran -c -g toms456_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms456_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms456_prb.o -L$HOME/libf77/$ARCH -ltoms456
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms456_prb.o"
  exit
fi
rm toms456_prb.o
#
mv a.out toms456_prb
./toms456_prb > toms456_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms456_prb"
  exit
fi
rm toms456_prb
#
echo "Test results written to toms456_prb_output.txt."
