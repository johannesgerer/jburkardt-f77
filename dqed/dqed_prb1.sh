#!/bin/bash
#
gfortran -c -g dqed_prb1.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dqed_prb1.f"
  exit
fi
rm compiler.txt
#
gfortran dqed_prb1.o -L$HOME/libf77/$ARCH -ldqed
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dqed_prb1.o"
  exit
fi
rm dqed_prb1.o
#
mv a.out dqed_prb1
./dqed_prb1 > dqed_prb1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dqed_prb1"
  exit
fi
rm dqed_prb1
#
echo "Test results written to dqed_prb1_output.txt."
