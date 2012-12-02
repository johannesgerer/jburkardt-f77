#!/bin/bash
#
gfortran -c -g toms443_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms443_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms443_prb.o -L$HOME/libf77/$ARCH -ltoms443
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms443_prb.o"
  exit
fi
rm toms443_prb.o
#
mv a.out toms443_prb
./toms443_prb > toms443_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms443_prb"
  exit
fi
rm toms443_prb
#
echo "Test results written to toms443_prb_output.txt."
