#!/bin/bash
#
gfortran -c -g subpak_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subpak_prb.f"
  exit
fi
rm compiler.txt
#
gfortran subpak_prb.o -L$HOME/libf77/$ARCH -lsubpak
if [ $? -ne 0 ]; then
  echo "Errors linking and loading subpak_prb.o"
  exit
fi
rm subpak_prb.o
#
mv a.out subpak_prb
./subpak_prb > subpak_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running subpak_prb"
  exit
fi
rm subpak_prb
#
echo "Test results written to subpak_prb_output.txt."
