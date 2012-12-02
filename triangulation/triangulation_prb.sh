#!/bin/bash
#
gfortran -c -g triangulation_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_prb.f"
  exit
fi
rm compiler.txt
#
gfortran triangulation_prb.o -L$HOME/libf77/$ARCH -ltriangulation
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_prb.o"
  exit
fi
rm triangulation_prb.o
#
mv a.out triangulation_prb
./triangulation_prb > triangulation_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangulation_prb"
  exit
fi
rm triangulation_prb
#
echo "Test results written to triangulation_prb_output.txt."
