#!/bin/bash
#
gfortran -c -g lawson_prb5.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lawson_prb5.f"
  exit
fi
rm compiler.txt
#
gfortran lawson_prb5.o -L$HOME/libf77/$ARCH -llawson
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lawson_prb5.o"
  exit
fi
rm lawson_prb5.o
#
mv a.out lawson_prb5
./lawson_prb5 > lawson_prb5_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lawson_prb5"
  exit
fi
rm lawson_prb5
#
echo "Test results written to lawson_prb5_output.txt."
