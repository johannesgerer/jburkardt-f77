#!/bin/bash
#
gfortran -c -g ising_2d_simulation.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ising_2d_simulation.f"
  exit
fi
rm compiler.txt
#
gfortran ising_2d_simulation.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ising_2d_simulation.o"
  exit
fi
rm ising_2d_simulation.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/ising_2d_simulation
#
echo "The program has been installed as ~/binf77/$ARCH/ising_2d_simulation."
