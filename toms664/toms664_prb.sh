#!/bin/bash
#
gfortran -c -g toms664_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms664_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms664_prb.o -L$HOME/libf77/$ARCH -ltoms664
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms664_prb.o"
  exit
fi
rm toms664_prb.o
#
mv a.out toms664_prb
./toms664_prb > toms664_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms664_prb"
  exit
fi
rm toms664_prb
#
echo "Test results written to toms664_prb_output.txt."
