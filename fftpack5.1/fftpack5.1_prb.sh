#!/bin/bash
#
gfortran -c -O3 fftpack5.1_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fftpack5.1_prb.f"
  exit
fi
rm compiler.txt
#
gfortran fftpack5.1_prb.o -L$HOME/libf77/$ARCH -lfftpack5.1
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fftpack5.1_prb.o"
  exit
fi
rm fftpack5.1_prb.o
#
mv a.out fftpack5.1_prb
./fftpack5.1_prb > fftpack5.1_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fftpack5.1_prb"
  exit
fi
rm fftpack5.1_prb
#
echo "Test results written to fftpack5.1_prb_output.txt."
