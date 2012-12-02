#!/bin/bash
#
gfortran -c -g unicycle_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling unicycle_prb.f"
  exit
fi
rm compiler.txt
#
gfortran unicycle_prb.o -L$HOME/libf77/$ARCH -lunicycle
if [ $? -ne 0 ]; then
  echo "Errors linking and loading unicycle_prb.o"
  exit
fi
rm unicycle_prb.o
#
mv a.out unicycle_prb
./unicycle_prb > unicycle_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running unicycle_prb"
  exit
fi
rm unicycle_prb
#
echo "Test program output written to unicycle_prb_output.txt."
