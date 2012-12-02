#!/bin/bash
#
gfortran -c -g toms626_prb1.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms626_prb1.f"
  exit
fi
rm compiler.txt
#
gfortran toms626_prb1.o -L$HOME/libf77/$ARCH -ltoms626 -lcalcomp
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms626_prb1.o"
  exit
fi
rm toms626_prb1.o
#
mv a.out toms626_prb1
./toms626_prb1 > toms626_prb1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms626_prb1"
  exit
fi
rm toms626_prb1
#
echo "Test results written to toms626_prb1_output.txt."
