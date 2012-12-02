#!/bin/bash
#
gfortran -c -g test_optimization_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_optimization_prb.f"
  exit
fi
rm compiler.txt
#
gfortran test_optimization_prb.o -L$HOME/libf77/$ARCH -ltest_optimization
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_optimization_prb.o"
  exit
fi
rm test_optimization_prb.o
#
mv a.out test_optimization_prb
./test_optimization_prb > test_optimization_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_optimization_prb"
  exit
fi
rm test_optimization_prb
#
echo "Test program output written to test_optimization_prb_output.txt."
