#!/bin/bash
#
gfortran -c -g driv_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling driv_prb.f"
  exit
fi
rm compiler.txt
#
gfortran driv_prb.o -L$HOME/libf77/$ARCH -ldriv
if [ $? -ne 0 ]; then
  echo "Errors linking and loading driv_prb.o"
  exit
fi
rm driv_prb.o
#
mv a.out driv_prb
./driv_prb > driv_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running driv_prb"
  exit
fi
rm driv_prb
#
echo "Test results written to driv_prb_output.txt."
