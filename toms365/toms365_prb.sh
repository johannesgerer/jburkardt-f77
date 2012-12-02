#!/bin/bash
#
gfortran -c -g toms365_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms365_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms365_prb.o -L$HOME/libf77/$ARCH -ltoms365
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms365_prb.o"
  exit
fi
rm toms365_prb.o
#
mv a.out toms365_prb
./toms365_prb > toms365_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms365_prb"
  exit
fi
rm toms365_prb
#
echo "Test results written to toms365_prb_output.txt."
