#!/bin/bash
#
gfortran -c -g praxis_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling praxis_prb.f"
  exit
fi
rm compiler.txt
#
gfortran praxis_prb.o -L$HOME/libf77/$ARCH -lpraxis
if [ $? -ne 0 ]; then
  echo "Errors linking and loading praxis_prb.o"
  exit
fi
rm praxis_prb.o
#
mv a.out praxis_prb
./praxis_prb > praxis_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running praxis_prb"
  exit
fi
rm praxis_prb
#
echo "Test results written to praxis_prb_output.txt."
