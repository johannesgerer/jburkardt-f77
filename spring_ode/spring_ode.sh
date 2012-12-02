#!/bin/bash
#
gfortran -c spring_ode.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spring_ode.f"
  exit
fi
rm compiler.txt
#
gfortran spring_ode.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spring_ode.o"
  exit
fi
rm spring_ode.o
#
mv a.out ~/binf77/$ARCH/spring_ode
#
echo "Executable installed as ~/binf77/$ARCH/spring_ode"
