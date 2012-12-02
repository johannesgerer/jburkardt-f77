#!/bin/bash
#
gfortran -c -g netode.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling netode.f"
  exit
fi
rm compiler.txt
#
gfortran netode.o -L$HOME/libf77/$ARCH -llsodi -lmachine
if [ $? -ne 0 ]; then
  echo "Errors loading args.o"
  exit
fi
rm netode.o
#
mv a.out ~/binf77/$ARCH/netode
#
echo "Executable installed as ~/binf77/$ARCH/netode."
