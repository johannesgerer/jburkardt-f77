#!/bin/bash
#
gfortran -c -g toms494_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms494_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms494_prb.o -L$HOME/libf77/$ARCH -ltoms494
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms494_prb.o"
  exit
fi
rm toms494_prb.o
#
mv a.out toms494_prb
./toms494_prb > toms494_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms494_prb"
  exit
fi
rm toms494_prb
#
echo "Test results written to toms494_prb_output.txt."
