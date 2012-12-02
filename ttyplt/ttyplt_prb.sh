#!/bin/bash
#
gfortran -c -g ttyplt_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ttyplt_prb.f"
  exit
fi
rm compiler.txt
#
gfortran ttyplt_prb.o -L$HOME/libf77/$ARCH -lttyplt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ttyplt_prb.o"
  exit
fi
rm ttyplt_prb.o
#
mv a.out ttyplt_prb
./ttyplt_prb > ttyplt_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ttyplt_prb"
  exit
fi
rm ttyplt_prb
#
echo "Program output written to ttyplt_prb_output.txt"
