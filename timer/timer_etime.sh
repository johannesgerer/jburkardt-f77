#!/bin/bash
#
gfortran -c timer_etime.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timer_etime.f"
  exit
fi
rm compiler.txt
#
gfortran timer_etime.o
if [ $? -ne 0 ]; then
  echo "Errors loading timer_etime.o"
  exit
fi
rm timer_etime.o
#
mv a.out timer_etime
./timer_etime > timer_etime_output.txt
rm timer_etime
#
echo "Program output written to timer_etime_output.txt"
