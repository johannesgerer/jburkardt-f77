#!/bin/bash
#
gfortran -c -g cvt_main.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_main.f"
  exit
fi
rm compiler.txt
#
gfortran cvt_main.o -L$HOME/libf77/$ARCH -lsandia_cvt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_main.o"
  exit
fi
rm cvt_main.o
#
mv a.out ~/binf77/$ARCH/cvt_main
#
echo "Executable installed as ~/binf77/$ARCH/cvt_main"
