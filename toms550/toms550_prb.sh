#!/bin/bash
#
gfortran -c -g toms550_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms550_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms550_prb.o -L$HOME/libf77/$ARCH -ltoms550
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms550_prb.o"
  exit
fi
rm toms550_prb.o
#
mv a.out toms550_prb
./toms550_prb < cutout_cube.txt > cutout_cube_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms550_prb"
  exit
fi
rm toms550_prb
#
echo "Test results written to cutout_cube_output.txt."
