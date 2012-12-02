#!/bin/bash
#
gfortran -c -g stochastic_rk_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stochastic_rk_prb.f"
  exit
fi
rm compiler.txt
#
gfortran stochastic_rk_prb.o -L$HOME/libf77/$ARCH -lstochastic_rk
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stochastic_rk_prb.o"
  exit
fi
rm stochastic_rk_prb.o
#
mv a.out stochastic_rk_prb
./stochastic_rk_prb > stochastic_rk_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running stochastic_rk_prb"
  exit
fi
rm stochastic_rk_prb
#
echo "Test results written to stochastic_rk_prb_output.txt."
