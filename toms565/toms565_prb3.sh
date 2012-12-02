#!/bin/bash
#
gfortran -c -g toms565_prb3.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms565_prb3.f"
  exit
fi
rm compiler.txt
#
gfortran toms565_prb3.o -L$HOME/libf77/$ARCH -ltoms565
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms565_prb3.o"
  exit
fi
rm toms565_prb3.o
#
mv a.out toms565_prb3
./toms565_prb3 > toms565_prb3_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms565_prb3"
  exit
fi
rm toms565_prb3
#
echo "Test results written to toms565_prb3_output.txt."
