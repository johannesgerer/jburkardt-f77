#!/bin/bash
#
gfortran -c -g fem2d_pack_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_pack_prb.f"
  exit
fi
rm compiler.txt
#
gfortran fem2d_pack_prb.o -L$HOME/libf77/$ARCH -lfem2d_pack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_pack_prb.o"
  exit
fi
rm fem2d_pack_prb.o
#
mv a.out fem2d_pack_prb
./fem2d_pack_prb > fem2d_pack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem2d_pack_prb"
  exit
fi
rm fem2d_pack_prb
#
echo "Test results written to fem2d_pack_prb_output.txt."
