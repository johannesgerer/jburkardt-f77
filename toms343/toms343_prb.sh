#!/bin/bash
#
gfortran -c -g toms343_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms343_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms343_prb.o -L$HOME/libf77/$ARCH -ltoms343
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms343_prb.o"
  exit
fi
rm toms343_prb.o
#
mv a.out toms343_prb
./toms343_prb > toms343_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms343_prb"
  exit
fi
rm toms343_prb
#
echo "Test results written to toms343_prb_output.txt."
