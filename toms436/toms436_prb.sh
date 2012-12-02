#!/bin/bash
#
gfortran -c -g toms436_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms436_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms436_prb.o -L$HOME/libf77/$ARCH -ltoms436
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms436_prb.o"
  exit
fi
rm toms436_prb.o
#
mv a.out toms436_prb
./toms436_prb > toms436_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms436_prb"
  exit
fi
rm toms436_prb
#
echo "Test results written to toms436_prb_output.txt."
