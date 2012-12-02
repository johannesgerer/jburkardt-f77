#!/bin/bash
#
gfortran -c -g beta_nc_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling beta_nc_prb.f"
  exit
fi
rm compiler.txt
#
gfortran beta_nc_prb.o -L$HOME/libf77/$ARCH -lbeta_nc
if [ $? -ne 0 ]; then
  echo "Errors linking and loading beta_nc_prb.o"
  exit
fi
rm beta_nc_prb.o
#
mv a.out beta_nc_prb
./beta_nc_prb > beta_nc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running beta_nc_prb"
  exit
fi
rm beta_nc_prb
#
echo "Test results written to beta_nc_prb_output.txt."
