#!/bin/bash
#
gfortran -c -g toms385_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms385_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms385_prb.o -L$HOME/libf77/$ARCH -ltoms385
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms385_prb.o"
  exit
fi
rm toms385_prb.o
#
mv a.out toms385_prb
./toms385_prb > toms385_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms385_prb"
  exit
fi
rm toms385_prb
#
echo "Test results written to toms385_prb_output.txt."
