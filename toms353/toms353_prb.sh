#!/bin/bash
#
gfortran -c -g toms353_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms353_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms353_prb.o -L$HOME/libf77/$ARCH -ltoms353
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms353_prb.o"
  exit
fi
rm toms353_prb.o
#
mv a.out toms353_prb
./toms353_prb > toms353_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms353_prb"
  exit
fi
rm toms353_prb
#
echo "Test results written to toms353_prb_output.txt."
