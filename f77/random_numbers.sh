#!/bin/bash
#
gfortran -c random_numbers.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling random_numbers.f"
  exit
fi
rm compiler.txt
#
gfortran random_numbers.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading random_numbers.o"
  exit
fi
rm random_numbers.o
#
mv a.out random_numbers
./random_numbers > random_numbers_output.txt
rm random_numbers
#
echo "Program output written to random_numbers_output.txt"
