#!/bin/bash
#
gfortran -c -g orbital.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling orbital.f"
  exit
fi
rm compiler.txt
#
gfortran orbital.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading orbital.o."
  exit
fi
#
rm orbital.o
#
mv a.out orbital
./orbital > orbital_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running orbital."
  exit
fi
rm orbital
#
echo "Program output written to orbital_output.txt"
