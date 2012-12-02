#!/bin/bash
#
gfortran -c -g toms370_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms370_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms370_prb.o -L$HOME/libf77/$ARCH -ltoms370
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms370_prb.o"
  exit
fi
rm toms370_prb.o
#
mv a.out toms370_prb
./toms370_prb > toms370_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms370_prb"
  exit
fi
rm toms370_prb
#
echo "Test results written to toms370_prb_output.txt."
