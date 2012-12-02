#!/bin/bash
#
gfortran -c -g fem1d_project.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_project.f"
  exit
fi
rm compiler.txt
#
gfortran fem1d_project.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_project.o"
  exit
fi
#
rm fem1d_project.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/fem1d_project
#
echo "Program installed as ~/binf77/$ARCH/fem1d_project"
