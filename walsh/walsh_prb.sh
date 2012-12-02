#!/bin/bash
#
gfortran -c -g walsh_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling walsh_prb.f"
  exit
fi
rm compiler.txt
#
gfortran walsh_prb.o -L$HOME/libf77/$ARCH -lwalsh
if [ $? -ne 0 ]; then
  echo "Errors linking and loading walsh_prb.o"
  exit
fi
rm walsh_prb.o
#
mv a.out walsh_prb
./walsh_prb > walsh_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running walsh_prb"
  exit
fi
rm walsh_prb
#
echo "Program output written to walsh_prb_output.txt"
