#!/bin/bash
#
gfortran -c -g lamp_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lamp_prb.f"
  exit
fi
rm compiler.txt
#
gfortran lamp_prb.o -L$HOME/libf77/$ARCH -llamp
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lamp_prb.o"
  exit
fi
rm lamp_prb.o
#
mv a.out lamp_prb
./lamp_prb > lamp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lamp_prb"
  exit
fi
rm lamp_prb
#
echo "Test results written to lamp_prb_output.txt."
