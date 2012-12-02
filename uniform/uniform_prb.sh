#!/bin/bash
#
gfortran -c -g uniform_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling uniform_prb.f"
  exit
fi
rm compiler.txt
#
gfortran uniform_prb.o -L$HOME/libf77/$ARCH -luniform
if [ $? -ne 0 ]; then
  echo "Errors linking and loading uniform_prb.o"
  exit
fi
rm uniform_prb.o
#
mv a.out uniform_prb
./uniform_prb > uniform_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running uniform_prb"
  exit
fi
rm uniform_prb
#
echo "Test results written to uniform_prb_output.txt."
