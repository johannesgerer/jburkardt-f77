#!/bin/bash
#
gfortran -c fd1d_predator_prey.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_predator_prey.f"
  exit
fi
rm compiler.txt
#
gfortran fd1d_predator_prey.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fd1d_predator_prey.o"
  exit
fi
rm fd1d_predator_prey.o
#
mv a.out fd1d_predator_prey
./fd1d_predator_prey > fd1d_predator_prey_output.txt
rm fd1d_predator_prey
#
echo "Program output written to fd1d_predator_prey_output.txt"
