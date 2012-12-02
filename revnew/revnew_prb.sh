#!/bin/bash
#
gfortran -c -g revnew_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling revnew_prb.f"
  exit
fi
rm compiler.txt
#
gfortran revnew_prb.o -L$HOME/libf77/$ARCH -lrevnew
if [ $? -ne 0 ]; then
  echo "Errors linking and loading revnew_prb.o"
  exit
fi
rm revnew_prb.o
#
mv a.out revnew_prb
./revnew_prb > revnew_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running revnew_prb"
  exit
fi
rm revnew_prb
#
echo "Test results written to revnew_prb_output.txt."
