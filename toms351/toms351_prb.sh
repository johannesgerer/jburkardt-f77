#!/bin/bash
#
gfortran -c -g toms351_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms351_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms351_prb.o -L$HOME/libf77/$ARCH -ltoms351
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms351_prb.o"
  exit
fi
rm toms351_prb.o
#
mv a.out toms351_prb
./toms351_prb > toms351_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms351_prb"
  exit
fi
rm toms351_prb
#
echo "Test results written to toms351_prb_output.txt."
