#!/bin/bash
#
gfortran -c -g special_functions_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling special_functions_prb.f"
  exit
fi
rm compiler.txt
#
gfortran special_functions_prb.o -L$HOME/libf77/$ARCH -lspecial_functions
if [ $? -ne 0 ]; then
  echo "Errors linking and loading special_functions_prb.o"
  exit
fi
rm special_functions_prb.o
#
mv a.out special_functions_prb
./special_functions_prb > special_functions_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running special_functions_prb"
  exit
fi
rm special_functions_prb
#
echo "Program output written to special_functions_prb_output.txt"
