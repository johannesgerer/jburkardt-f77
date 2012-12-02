#!/bin/bash
#
gfortran -c -g toms332_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms332_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms332_prb.o -L$HOME/libf77/$ARCH -ltoms332
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms332_prb.o"
  exit
fi
rm toms332_prb.o
#
mv a.out toms332_prb
./toms332_prb > toms332_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms332_prb"
  exit
fi
rm toms332_prb
#
echo "Test results written to toms332_prb_output.txt."
