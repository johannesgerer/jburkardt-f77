#!/bin/bash
#
gfortran -c -g toms708_prb1.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms708_prb1.f"
  exit
fi
rm compiler.txt
#
gfortran toms708_prb1.o -L$HOME/libf77/$ARCH -ltoms708
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms708_prb1.o"
  exit
fi
rm toms708_prb1.o
#
mv a.out toms708_prb1
./toms708_prb1 > toms708_prb1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms708_prb1"
  exit
fi
rm toms708_prb1
#
echo "Test results written to toms708_prb1_output.txt."
