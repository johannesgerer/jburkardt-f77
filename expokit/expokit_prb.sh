#!/bin/bash
#
gfortran -c -g expokit_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling expokit_prb.f"
  exit
fi
rm compiler.txt
#
gfortran expokit_prb.o -L$HOME/libf77/$ARCH -lexpokit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading expokit_prb.o"
  exit
fi
rm expokit_prb.o
#
mv a.out expokit_prb
./expokit_prb > expokit_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running expokit_prb"
  exit
fi
rm expokit_prb
#
echo "Test results written to expokit_prb_output.txt."
