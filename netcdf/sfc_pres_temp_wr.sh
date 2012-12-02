#!/bin/bash
#
cp ~/include/netcdf.inc .
gfortran -c -g sfc_pres_temp_wr.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sfc_pres_temp_wr.f"
  exit
fi
rm compiler.txt
rm netcdf.inc
#
gfortran sfc_pres_temp_wr.o -L$HOME/lib/$ARCH -lnetcdf
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sfc_pres_temp_wr.o"
  exit
fi
rm sfc_pres_temp_wr.o
#
mv a.out sfc_pres_temp_wr
./sfc_pres_temp_wr > sfc_pres_temp_wr_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sfc_pres_temp_wr"
  exit
fi
rm sfc_pres_temp_wr
#
echo "Test results written to sfc_pres_temp_wr_output.txt."
