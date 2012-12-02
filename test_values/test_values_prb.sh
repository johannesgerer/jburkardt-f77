#!/bin/bash
#
gfortran -c -g test_values_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_values_prb.f"
  exit
fi
rm compiler.txt
#
gfortran test_values_prb.o -L$HOME/libf77/$ARCH -ltest_values
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_values_prb.o"
  exit
fi
rm test_values_prb.o
#
mv a.out test_values_prb
./test_values_prb > test_values_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_values_prb"
  exit
fi
rm test_values_prb
#
echo "Test results written to test_values_prb_output.txt."
