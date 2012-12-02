#!/bin/bash
#
gfortran -c -g lawson_prb6.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lawson_prb6.f"
  exit
fi
rm compiler.txt
#
gfortran lawson_prb6.o -L$HOME/libf77/$ARCH -llawson
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lawson_prb6.o"
  exit
fi
rm lawson_prb6.o
#
mv a.out lawson_prb6
./lawson_prb6 > lawson_prb6_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lawson_prb6"
  exit
fi
rm lawson_prb6
#
echo "Test results written to lawson_prb6_output.txt."
