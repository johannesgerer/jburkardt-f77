#!/bin/bash
#
gfortran -c -g hregion_05.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hregion_05.f"
  exit
fi
rm compiler.txt
#
gfortran hregion_05.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hregion_05.o"
  exit
fi
rm hregion_05.o
#
mv a.out ~/binf77/$ARCH/hregion_05
#
echo "A new version of hregion_05 has been created."
