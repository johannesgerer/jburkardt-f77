#!/bin/bash
#
gfortran -c tcell_mass.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tcell_mass.f"
  exit
fi
rm compiler.txt
#
gfortran tcell_mass.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tcell_mass.o"
  exit
fi
rm tcell_mass.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/tcell_mass
#
echo "Executable installed as ~/binf77/$ARCH/tcell_mass"
