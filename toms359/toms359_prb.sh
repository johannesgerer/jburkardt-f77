#!/bin/bash
#
gfortran -c -g toms359_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms359_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms359_prb.o -L$HOME/libf77/$ARCH -ltoms359
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms359_prb.o"
  exit
fi
rm toms359_prb.o
#
mv a.out toms359_prb
./toms359_prb > toms359_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms359_prb"
  exit
fi
rm toms359_prb
#
echo "Test results written to toms359_prb_output.txt."
