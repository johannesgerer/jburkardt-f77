#!/bin/bash
#
gfortran -c -g toms581_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms581_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms581_prb.o -L$HOME/libf77/$ARCH -ltoms581
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms581_prb.o"
  exit
fi
rm toms581_prb.o
#
mv a.out toms581_prb
./toms581_prb > toms581_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms581_prb"
  exit
fi
rm toms581_prb
#
echo "Test results written to toms581_prb_output.txt."
