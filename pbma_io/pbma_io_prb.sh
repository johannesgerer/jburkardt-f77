#!/bin/bash
#
gfortran -c -g pbma_io_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbma_io_prb.f"
  exit
fi
rm compiler.txt
#
gfortran pbma_io_prb.o -L$HOME/libf77/$ARCH -lpbma_io
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pbma_io_prb.o"
  exit
fi
rm pbma_io_prb.o
#
mv a.out pbma_io_prb
./pbma_io_prb > pbma_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pbma_io_prb"
  exit
fi
rm pbma_io_prb
#
echo "Test results written to pbma_io_prb_output.txt."
