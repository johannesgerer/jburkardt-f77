#!/bin/bash
#
gfortran -c -g quadrule_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrule_prb.f"
  exit
fi
rm compiler.txt
#
gfortran quadrule_prb.o -L$HOME/libf77/$ARCH -lquadrule
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadrule_prb.o"
  exit
fi
rm quadrule_prb.o
#
mv a.out quadrule_prb
./quadrule_prb > quadrule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quadrule_prb"
  exit
fi
rm quadrule_prb
#
echo "Test results written to quadrule_prb_output.txt."
