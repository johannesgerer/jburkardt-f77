#!/bin/bash
#
gfortran -c -g geompack_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geompack_prb.f"
  exit
fi
rm compiler.txt
#
gfortran geompack_prb.o -L$HOME/libf77/$ARCH -lgeompack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading geompack_prb.o"
  exit
fi
rm geompack_prb.o
#
mv a.out geompack_prb
./geompack_prb > geompack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running geompack_prb"
  exit
fi
rm geompack_prb
#
echo "Test program output written to geompack_prb_output.txt."
