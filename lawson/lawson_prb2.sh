#!/bin/bash
#
gfortran -c -g lawson_prb2.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lawson_prb2.f"
  exit
fi
rm compiler.txt
#
gfortran lawson_prb2.o -L$HOME/libf77/$ARCH -llawson
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lawson_prb2.o"
  exit
fi
rm lawson_prb2.o
#
mv a.out lawson_prb2
./lawson_prb2 > lawson_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lawson_prb2"
  exit
fi
rm lawson_prb2
#
echo "Test results written to lawson_prb2_output.txt."
