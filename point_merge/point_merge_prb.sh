#!/bin/bash
#
gfortran -c -g point_merge_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling point_merge_prb.f"
  exit
fi
rm compiler.txt
#
gfortran point_merge_prb.o -L$HOME/libf77/$ARCH -lpoint_merge
if [ $? -ne 0 ]; then
  echo "Errors linking and loading point_merge_prb.o"
  exit
fi
rm point_merge_prb.o
#
mv a.out point_merge_prb
./point_merge_prb > point_merge_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running point_merge_prb"
  exit
fi
rm point_merge_prb
#
echo "Test results written to point_merge_prb_output.txt."
