#!/bin/bash
#
gfortran -c -g band_qr_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling band_qr_prb.f"
  exit
fi
rm compiler.txt
#
gfortran band_qr_prb.o -L$HOME/libf77/$ARCH -lband_qr -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading band_qr_prb.o"
  exit
fi
rm band_qr_prb.o
#
mv a.out band_qr_prb
./band_qr_prb > band_qr_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running band_qr_prb"
  exit
fi
rm band_qr_prb
#
echo "Test results written to band_qr_prb_output.txt."
