#!/bin/bash
#
gfortran -c -g owens_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling owens_prb.f"
  exit
fi
rm compiler.txt
#
gfortran owens_prb.o -L$HOME/libf77/$ARCH -lowens
if [ $? -ne 0 ]; then
  echo "Errors linking and loading owens_prb.o"
  exit
fi
rm owens_prb.o
#
mv a.out owens_prb
./owens_prb > owens_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running owens_prb"
  exit
fi
rm owens_prb
#
echo "Test results written to owens_prb_output.txt."
