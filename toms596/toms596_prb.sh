#!/bin/bash
#
gfortran -c -g toms596_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms596_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms596_prb.o -L$HOME/libf77/$ARCH -ltoms596
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms596_prb.o"
  exit
fi
rm toms596_prb.o
#
mv a.out toms596_prb
./toms596_prb > toms596_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms596_prb"
  exit
fi
rm toms596_prb
#
echo "Test results written to toms596_prb_output.txt."
