#!/bin/bash
#
gfortran -c -g clebsch_gordan_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling clebsch_gordan_prb.f"
  exit
fi
rm compiler.txt
#
gfortran clebsch_gordan_prb.o -L$HOME/libf77/$ARCH -lclebsch_gordan
if [ $? -ne 0 ]; then
  echo "Errors linking and loading clebsch_gordan_prb.o"
  exit
fi
rm clebsch_gordan_prb.o
#
mv a.out clebsch_gordan_prb
./clebsch_gordan_prb > clebsch_gordan_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running clebsch_gordan_prb"
  exit
fi
rm clebsch_gordan_prb
#
echo "Test results written to clebsch_gordan_prb_output.txt."
