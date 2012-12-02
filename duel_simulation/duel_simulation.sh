#!/bin/bash
#
gfortran -c -g duel_simulation.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling duel_simulation.f"
  exit
fi
rm compiler.txt
#
gfortran duel_simulation.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading duel_simulation.o"
  exit
fi
rm duel_simulation.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/duel_simulation
#
echo "Executable installed as ~/binf77/$ARCH/duel_simulation"
