#!/bin/bash
#
cp npbparams_C.h npbparams.h
#
gfortran -c -g mg_serial.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mg_serial.f"
  exit
fi
rm compiler.txt
#
gfortran -c print_results.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling print_results.f"
  exit
fi
rm compiler.txt
#
gfortran -c randdp.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling randdp.f"
  exit
fi
rm compiler.txt
#
gfortran -c timers.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timers.f"
  exit
fi
rm compiler.txt
#
gcc -c wtime.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wtime.c"
  exit
fi
rm compiler.txt
#
gfortran mg_serial.o print_results.o randdp.o timers.o wtime.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mg_serial.o + print_results.o + randdp.o + timers.o + wtime.o"
  exit
fi
rm mg_serial.o
rm print_results.o
rm randdp.o
rm timers.o
rm wtime.o
rm npbparams.h
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/mg_serial_C
#
echo "Executable installed as ~/binf77/$ARCH/mg_serial_C"
