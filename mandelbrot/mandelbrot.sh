#!/bin/bash
#
gfortran -c -g mandelbrot.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mandelbrot.f"
  exit
fi
rm compiler.txt
#
gfortran mandelbrot.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mandelbrot.o"
  exit
fi
rm mandelbrot.o
#
mv a.out ~/binf77/$ARCH/mandelbrot
#
echo "Executable installed as ~/binf77/$ARCH/mandelbrot"
