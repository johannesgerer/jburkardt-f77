#!/bin/csh
#
F77 -c -g starpac_prb.f >& compiler.txt
if ( $status != 0 ) then
  echo "Errors compiling starpac_prb.f"
  exit
endif
rm compiler.txt
#
F77 starpac_prb.o -L$HOME/libf77/$ARCH -lstarpac -lblas1 -lmachine -lnl2sol
if ( $status != 0 ) then
  echo "Errors linking and loading starpac_prb.o"
  exit
endif
rm starpac_prb.o
#
mv a.out starpac_prb
./starpac_prb > starpac_prb_output.txt
if ( $status != 0 ) then
  echo "Errors running starpac_prb"
  exit
endif
rm starpac_prb
#
echo "Program output written to starpac_prb_output.txt"
