#!/bin/bash
#
gfortran -c -g toms648_nsd_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms648_nsd_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toms648_nsd_prb.o -L$HOME/libf77/$ARCH -ltoms648_nsd
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms648_nsd_prb.o"
  exit
fi
rm toms648_nsd_prb.o
#
mv a.out toms648_nsd_prb
./toms648_nsd_prb > toms648_nsd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms648_nsd_prb"
  exit
fi
rm toms648_nsd_prb
#
echo "Test results written to toms648_nsd_prb_output.txt."
