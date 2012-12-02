#!/bin/bash
#
gfortran -c -g toms706_prb2.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms706_prb2.f"
  exit
fi
rm compiler.txt
#
gfortran toms706_prb2.o -L$HOME/libf77/$ARCH -ltoms706
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms706_prb2.o"
  exit
fi
rm toms706_prb2.o
#
mv a.out toms706_prb2
./toms706_prb2 > toms706_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms706_prb2"
  exit
fi
rm toms706_prb2
#
echo "Test results written to toms706_prb2_output.txt."
