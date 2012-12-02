#!/bin/bash
#
gfortran -c -g toms612_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms612_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms612_prb.o -L$HOME/libf77/$ARCH -ltoms612
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms612_prb.o"
  exit
fi
rm toms612_prb.o
#
mv a.out toms612_prb
./toms612_prb > toms612_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms612_prb"
  exit
fi
rm toms612_prb
#
echo "Test results written to toms612_prb_output.txt."
