#!/bin/bash
#
gfortran -c -g tcell.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tcell.f"
  exit
fi
rm compiler.txt
#
gfortran tcell.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tcell.o"
  exit
fi
rm tcell.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/tcell
#
echo "Program installed as ~/binf77/$ARCH/tcell"
