#!/bin/bash
#
gfortran -c anyplt_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling anyplt_prb.f"
  exit
fi
rm compiler.txt
#
gfortran anyplt_prb.o -L$HOME/libf77/$ARCH -lanytty
if [ $? -ne 0 ]; then
  echo "Errors linking and loading anyplt_prb.o"
  exit
fi
rm anyplt_prb.o
#
mv a.out anytty_prb
./anytty_prb > anytty_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running anytty_prb"
fi
rm anytty_prb
#
echo "Program output written to anytty_prb_output.txt"
