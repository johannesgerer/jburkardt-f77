#!/bin/bash
#
gfortran -c -g toms435_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms435_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms435_prb.o -L$HOME/libf77/$ARCH -ltoms435
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms435_prb.o"
  exit
fi
rm toms435_prb.o
#
mv a.out toms435_prb
./toms435_prb > toms435_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms435_prb"
  exit
fi
rm toms435_prb
#
echo "Test results written to toms435_prb_output.txt."
