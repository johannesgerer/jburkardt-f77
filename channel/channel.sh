#!/bin/bash
#
gfortran -c -g channel.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling channel.f"
  exit
fi
rm compiler.txt
#
gfortran channel.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading channel.o"
  exit
fi
rm channel.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/channel
#
echo "Executable installed as ~/binf77/$ARCH/channel"
